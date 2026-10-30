using System;
using System.Numerics;
using System.Runtime.CompilerServices;
using Graphs.Models.Potts;
using Graphs.Primitives;
using Maths.Rng;

namespace Clustering.Graphical.SPC.Runtime.Core.Sampler;

/// <summary>
/// Single-temperature Swendsen-Wang sampler for the q-state Potts model.
/// Generic over an <see cref="ISwConfig"/> struct so the JIT can specialize
/// the hot path and dead-code-eliminate the disabled accumulator branch.
/// </summary>
/// <remarks>
/// <para><b>Hot-path discipline.</b> No per-cycle allocation. All scratch
/// arrays are fields, sized at construction. Bond probabilities are precomputed
/// once; the inner loop is a table lookup. The bond-formation and energy
/// accumulation share a single edge walk because both decisions key off the
/// same spin-equality test. The <c>j &lt;= i</c> skip handles symmetric CSR
/// (each undirected edge stored as both (i,j) and (j,i)). The FK accumulator
/// widens to <see langword="long"/> on multiply to survive large N.</para>
///
/// <para><b>Color semantics.</b> <c>_spins[i]</c> is a Potts color in
/// <c>[0, q)</c>, not a cluster ID. Bond formation reads the OLD spins (Pass 1)
/// and the energy is evaluated against that same configuration. The cluster
/// flip (Pass 2) overwrites <c>_spins</c> with NEW colors, one per UF root.
/// Two distinct bond clusters may independently draw the same color; this is
/// the "color collision" that distinguishes Blatt (which counts color buckets)
/// from FK (which counts UF root sizes directly).</para>
///
/// <para><b>Always-on accumulators.</b> FK χ, specific-heat moments, Blatt
/// magnetization, the cluster-size histogram and the per-draw record (E, m,
/// Σ|c|², c_max) are all populated every cycle regardless of
/// <typeparamref name="TConfig"/>. They are free
/// byproducts of the SW union-find pass — gating them would double the
/// specialization count to save a handful of float ops per cycle, which is
/// the wrong trade. Only the per-edge currency gates
/// (<see cref="ISwConfig.Affinities"/> / <see cref="ISwConfig.Alignments"/>) earn
/// JIT-specialization gates — they own the only per-bond scatter writes, independently.</para>
/// </remarks>
internal sealed class SwendsenWang<TConfig> : ISwEngine
    where TConfig : struct, ISwConfig
{
    // ── Inputs (read-only after construction) ─────────────────────────────
    private readonly CsrGraph _graph;
    private readonly double _temperature;
    private readonly int _q;
    private readonly int _n;

    // ── State ─────────────────────────────────────────────────────────────
    private readonly int[] _spins;             // [N] Potts color of each node, in [0, q)
    private readonly UnionFind _uf;            // bond-cluster structure, rebuilt every cycle
    private readonly double[] _bondTable;      // [E] p_e = 1 - exp(-J_e / T), precomputed once
    private readonly Xoshiro256PlusPlus _rng;

    // ── Scratch buffers (reused across cycles; never reallocated) ─────────
    private readonly int[] _rootSizeScratch;   // [N] compacted UF root sizes after WriteRootSizesTo
    private readonly int[] _clusterMap;        // [N] root index -> new color, sentinel -1
    private readonly int[] _colorBucketScratch; // [Q] color counts (Blatt is always on)

    // ── Per-edge accumulators (currency specializations only) ─────────────
    // Indexed by CSR entry. Symmetric CSR stores each undirected edge twice;
    // only the j > i half is ever incremented (the SW loop's j ≤ i skip).
    private readonly int[]? _bondFormedCount;     // [E] cycles where bond formed at edge e
    private readonly int[]? _spinAgreementCount;  // [E] cycles where _spins[i] == _spins[j] at edge e

    // ── Per-node thermodynamic landscapes (runtime-gated; null when off) ──
    // Un-reduced 0-forms: the same per-node quantities the sweep collapses to
    // global χ (Σ_c |c|²) and the M order parameter, kept un-collapsed so a
    // downstream resolution step can ascend them. Gated by a plain runtime flag,
    // NOT a TConfig specialization — these are O(N) post-cycle folds, so a branch
    // is negligible and they don't earn a JIT gate (that's the per-edge path only).
    private readonly double[]? _sumClusterSizePerNode;      // [N] Σ_draws |cluster(i)|       (un-reduced χ)
    private readonly double[]? _sumInGiantClusterPerNode;   // [N] Σ_draws 1{i ∈ giant cluster} (un-reduced M)

    // ── Co-membership accumulator (runtime-gated; null when off) ─────────
    // O(E) post-pass per draw: for each undirected edge (i,j), j>i, increment when
    // Find(i)==Find(j). Lower variance than BondFormedCount because it captures
    // transitive co-clustering via multi-hop paths (⟨n_ij⟩, eq. 4 / Niedermayer).
    private readonly int[]? _coMembershipCount;

    // ── Observables ───────────────────────────────────────────────────────
    private readonly long[] _clusterSizeHistogram; // log-binned; bin k = sizes [2^k, 2^(k+1)); 64-bit, as a bin gains up to N a cycle
    private readonly long[]? _exactClusterSizeHistogram; // [N+1] entry s = clusters of size s (runtime-gated; null when off)

    private double _runningSumSqClusterSizes;       // Σ_cycles Σ_c |c|², the FK accumulator
    private double _runningSumSqClusterSizesExcl;   // Σ_cycles (Σ_c |c|² − max_c|c|²), percolation diagnostic
    private double _runningSumEnergy;               // Σ_cycles ⟨H⟩
    private double _runningSumEnergySq;             // Σ_cycles ⟨H²⟩, for specific heat
    private double _runningSumMag;                  // Σ_cycles m  (Blatt; always populated)
    private double _runningSumMagSq;                // Σ_cycles m² (Blatt; always populated)

    // ── Per-draw record (always on) ───────────────────────────────────────
    // Entry t is draw t since the last reset, so the record's length is _drawCount and the
    // running sums above are its sums in order. Run and BurnIn reserve their draw counts up
    // front, so a budgeted run never copies; a bare Draw past the capacity doubles it.
    private double[] _energyPerDraw = [];             // E of the spins entering the draw
    private double[] _magPerDraw = [];                // m of the spins leaving the draw
    private long[]   _sumSqClusterSizesPerDraw = [];  // Σ|c|² of the entering spins' FK clusters
    private int[]    _maxClusterSizePerDraw = [];     // c_max of the same clusters

    // ── Burn-in trace (runtime-gated; null when off) ─────────────────────
    // The per-draw record of every draw a BurnIn discards, appended just before the reset.
    private double[]? _burnInEnergyPerDraw;
    private double[]? _burnInMagPerDraw;
    private long[]?   _burnInSumSqClusterSizesPerDraw;
    private int[]?    _burnInMaxClusterSizePerDraw;

    private int _drawCount;

    public SwendsenWang(
        CsrGraph graph,
        double temperature,
        int q,
        AccumulationSpec accumulation,
        int? seed)
    {
        if (graph.NodeCount <= 0) throw new ArgumentException("Graph must have at least one node.", nameof(graph));
        if (graph.Targets is null || graph.Weights is null || graph.RowPointers is null)
            throw new ArgumentException("Graph CSR arrays are not initialized.", nameof(graph));
        if (temperature <= 0.0) throw new ArgumentOutOfRangeException(nameof(temperature), "Temperature must be strictly positive.");
        if (q < 2) throw new ArgumentOutOfRangeException(nameof(q), "Q must be at least 2.");

        _graph = graph;
        _temperature = temperature;
        _q = q;
        _n = graph.NodeCount;

        _spins = new int[_n];
        _uf = new UnionFind(_n);
        _bondTable = new double[graph.Targets.Length];
        _rng = new Xoshiro256PlusPlus(seed);

        // Precompute p_e = 1 - exp(-J_e / T). For symmetric CSR, (i,j) and (j,i)
        // entries hold the same p; only one is touched per cycle (j <= i skip).
        // The redundancy is a memory-vs-cache tradeoff: matching the CSR access
        // pattern keeps the table walk linear and prefetcher-friendly.
        for (int e = 0; e < _bondTable.Length; e++)
        {
            _bondTable[e] = FkKernel.BondProbability(graph.Weights[e], _temperature);
        }

        _rootSizeScratch = new int[_n];
        _clusterMap = new int[_n];

        // Blatt color-bucket scratch is always allocated (Blatt is unconditional).
        // O(Q) memory, Q typically ≤ 20.
        _colorBucketScratch = new int[_q];

        // Per-edge currency allocation; each array folds away independently when
        // its currency gate is false (Affinities ⟂ Alignments).
        _bondFormedCount    = TConfig.Affinities ? new int[_bondTable.Length] : null;
        _spinAgreementCount = TConfig.Alignments ? new int[_bondTable.Length] : null;

        // Per-node landscapes: allocated only when the run asks for them. The two
        // flags are independent capability bits (the sinks toggle separately).
        _sumClusterSizePerNode    = accumulation.ClusterSizeLandscape ? new double[_n] : null;
        _sumInGiantClusterPerNode = accumulation.OrderLandscape       ? new double[_n] : null;

        // Co-membership accumulator: allocated when the run requests it.
        _coMembershipCount = accumulation.CoMembership ? new int[_bondTable.Length] : null;

        // Burn-in trace: present (empty until a BurnIn) when the run asks for it.
        if (accumulation.BurnInTrace)
        {
            _burnInEnergyPerDraw            = [];
            _burnInMagPerDraw               = [];
            _burnInSumSqClusterSizesPerDraw = [];
            _burnInMaxClusterSizePerDraw    = [];
        }

        int maxBin = BitOperations.Log2((uint)Math.Max(1, _n)) + 1;
        _clusterSizeHistogram = new long[maxBin];
        _exactClusterSizeHistogram = accumulation.ExactClusterSizeHistogram ? new long[_n + 1] : null;

        // Random initial configuration (high-T disorder).
        for (int i = 0; i < _n; i++)
            _spins[i] = _rng.NextInt(_q);
    }

    // ── ISwEngine read surface ─────────────────────────────────────────

    public int DrawCount => _drawCount;
    public double Temperature => _temperature;
    public int Q => _q;
    public int N => _n;

    public ReadOnlySpan<int> Spins => _spins;
    public ReadOnlySpan<long> ClusterSizeHistogram => _clusterSizeHistogram;

    public double RunningSumSqClusterSizes     => _runningSumSqClusterSizes;
    public double RunningSumSqClusterSizesExcl => _runningSumSqClusterSizesExcl;
    public double RunningSumEnergy             => _runningSumEnergy;
    public double RunningSumEnergySq           => _runningSumEnergySq;
    public double RunningSumMag                => _runningSumMag;
    public double RunningSumMagSq              => _runningSumMagSq;

    public ReadOnlySpan<double> EnergyPerDraw            => _energyPerDraw.AsSpan(0, _drawCount);
    public ReadOnlySpan<double> MagPerDraw               => _magPerDraw.AsSpan(0, _drawCount);
    public ReadOnlySpan<long>   SumSqClusterSizesPerDraw => _sumSqClusterSizesPerDraw.AsSpan(0, _drawCount);
    public ReadOnlySpan<int>    MaxClusterSizePerDraw    => _maxClusterSizePerDraw.AsSpan(0, _drawCount);

    // ── Drivers ───────────────────────────────────────────────────────────

    public void Run(int drawCount)
    {
        if (drawCount < 0) throw new ArgumentOutOfRangeException(nameof(drawCount));
        ReserveRecord(checked(_drawCount + drawCount));
        for (int i = 0; i < drawCount; i++)
            Draw();
    }

    public void BurnIn(int cycles)
    {
        if (cycles < 0) throw new ArgumentOutOfRangeException(nameof(cycles));
        ReserveRecord(checked(_drawCount + cycles));
        for (int i = 0; i < cycles; i++)
            Draw();
        if (_burnInEnergyPerDraw is not null)
            AppendBurnInTrace();
        ResetAccumulators();
    }

    /// <summary>
    /// Append the per-draw record so far — every draw the coming reset discards — to the burn-in trace.
    /// </summary>
    private void AppendBurnInTrace()
    {
        int kept = _burnInEnergyPerDraw!.Length;
        int total = checked(kept + _drawCount);
        Array.Resize(ref _burnInEnergyPerDraw, total);
        Array.Resize(ref _burnInMagPerDraw, total);
        Array.Resize(ref _burnInSumSqClusterSizesPerDraw, total);
        Array.Resize(ref _burnInMaxClusterSizePerDraw, total);
        Array.Copy(_energyPerDraw, 0, _burnInEnergyPerDraw, kept, _drawCount);
        Array.Copy(_magPerDraw, 0, _burnInMagPerDraw, kept, _drawCount);
        Array.Copy(_sumSqClusterSizesPerDraw, 0, _burnInSumSqClusterSizesPerDraw, kept, _drawCount);
        Array.Copy(_maxClusterSizePerDraw, 0, _burnInMaxClusterSizePerDraw, kept, _drawCount);
    }

    /// <summary>
    /// One Swendsen-Wang cycle: bond formation (with energy as a free byproduct),
    /// cluster flip, and observable accumulation.
    /// </summary>
    public void Draw()
    {
        _uf.Reset();

        // ── Pass 1: single edge walk. Bond formation + energy in the same loop.
        // For each undirected edge (i,j) where the CURRENT spins match:
        //   * subtract J_ij from this cycle's energy (Potts H = -Σ J_ij δ),
        //   * with probability bondTable[e], union i and j in the UF.
        // Mismatched edges contribute neither energy nor a bond opportunity.
        double cycleEnergy = 0.0;
        for (int i = 0; i < _n; i++)
        {
            int rowEnd = _graph.RowPointers[i + 1];
            for (int e = _graph.RowPointers[i]; e < rowEnd; e++)
            {
                int j = _graph.Targets[e];
                if (j <= i) continue;  // undirected: visit each edge once

                if (_spins[i] == _spins[j])
                {
                    if (TConfig.Alignments) _spinAgreementCount![e]++;
                    cycleEnergy -= _graph.Weights[e];
                    if (_rng.NextDouble() < _bondTable[e])
                    {
                        if (TConfig.Affinities) _bondFormedCount![e]++;
                        _uf.Union(i, j);
                    }
                }
            }
        }

        // ── Pass 2: cluster flip. Each UF root draws ONE color uniformly from
        // [0, q); every member of that bond cluster inherits it. _clusterMap is
        // reset with a -1 sentinel, then used as a root -> color cache so the
        // RNG fires exactly once per cluster, not once per node.
        _clusterMap.AsSpan().Fill(-1);
        for (int i = 0; i < _n; i++)
        {
            int root = _uf.Find(i);
            int newColor = _clusterMap[root];
            if (newColor < 0)
            {
                newColor = _rng.NextInt(_q);
                _clusterMap[root] = newColor;
            }
            _spins[i] = newColor;
        }

        // ── Pass 3a: FK susceptibility accumulator (always on).
        // Σ_c |c|² over BOND CLUSTERS — read from the UF, not from spin buckets.
        // This is the §10.1 fix: immune to coincidental color collisions because
        // it never references _spins. long widening on the multiply lets us
        // safely accumulate sums for N well past int range.
        int rootCount = _uf.WriteRootSizesTo(_rootSizeScratch);
        ReadOnlySpan<int> rootSizes = _rootSizeScratch.AsSpan(0, rootCount);
        long sumSq = 0L;
        long maxSize = 0L;
        for (int r = 0; r < rootCount; r++)
        {
            long sz = rootSizes[r];
            sumSq += sz * sz;
            if (sz > maxSize) maxSize = sz;
        }
        _runningSumSqClusterSizes     += sumSq;
        _runningSumSqClusterSizesExcl += sumSq - maxSize * maxSize;

        // ── Pass 3b: Blatt 1996 magnetization estimator (always on).
        // m = (N_max · q / N - 1) / (q - 1), where N_max is the largest COLOR
        // BUCKET in the current spin configuration. O(N + Q) per cycle; with Q
        // typically ≤ 20, this is essentially free relative to Pass 1/2.
        int nMax = MaxColorBucket(_spins, _q, _colorBucketScratch);
        double m = ((double)nMax * _q / _n - 1.0) / (_q - 1.0);
        _runningSumMag   += m;
        _runningSumMagSq += m * m;

        // ── Pass 3c: log-binned cluster size histogram scatter. Scalar writes
        // (scatter doesn't vectorize), but the loop is short — rootCount is
        // typically << N. BitOperations.Log2 is a single hardware instruction.
        for (int r = 0; r < rootCount; r++)
        {
            int bin = BitOperations.Log2((uint)rootSizes[r]);
            _clusterSizeHistogram[bin]++;
        }
        // The exact histogram beside it (runtime-gated): one count at each root's own size.
        if (_exactClusterSizeHistogram is not null)
        {
            for (int r = 0; r < rootCount; r++)
                _exactClusterSizeHistogram[rootSizes[r]]++;
        }

        // ── Pass 3d: per-node thermodynamic landscapes (runtime-gated). The
        // un-collapsed 0-forms of Pass 3a's χ sum and giant-cluster size: for
        // each node, fold its FK-cluster size and its giant-cluster membership.
        // Find is O(1) here — Pass 2 already path-compressed every node — and
        // maxSize is reused from Pass 3a. A null array skips its branch.
        if (_sumClusterSizePerNode != null || _sumInGiantClusterPerNode != null)
        {
            // maxSize is a value, not a root; resolve a single canonical giant
            // root (the lowest-index root attaining maxSize) so a size tie can't
            // count every max-size cluster as "the giant". One giant per draw.
            int giantRoot = -1;
            if (_sumInGiantClusterPerNode != null)
            {
                for (int i = 0; i < _n; i++)
                {
                    int root = _uf.Find(i);
                    if (_uf.Size(root) == maxSize && (giantRoot < 0 || root < giantRoot))
                        giantRoot = root;
                }
            }
            for (int i = 0; i < _n; i++)
            {
                int root = _uf.Find(i);
                if (_sumClusterSizePerNode != null)
                    _sumClusterSizePerNode[i] += _uf.Size(root);
                if (_sumInGiantClusterPerNode != null && root == giantRoot)
                    _sumInGiantClusterPerNode[i] += 1.0;
            }
        }

        // ── Pass 3e: co-membership post-pass (runtime-gated). O(E) walk: for each
        // undirected edge (i,j), j>i, count draws where i and j are in the same bond
        // cluster. Path-compression from Pass 2 makes Find O(1) here.
        if (_coMembershipCount != null)
        {
            for (int i = 0; i < _n; i++)
            {
                int rowEnd = _graph.RowPointers[i + 1];
                for (int e = _graph.RowPointers[i]; e < rowEnd; e++)
                {
                    int j = _graph.Targets[e];
                    if (j <= i) continue;
                    if (_uf.Find(i) == _uf.Find(j))
                        _coMembershipCount[e]++;
                }
            }
        }

        _runningSumEnergy   += cycleEnergy;
        _runningSumEnergySq += cycleEnergy * cycleEnergy;

        // ── Record: the draw's entry in the per-draw series — the same values the sums
        // above just folded, so each sum stays its series summed in order.
        if (_drawCount == _energyPerDraw.Length)
            ReserveRecord((int)Math.Min(Array.MaxLength, Math.Max(16L, 2L * _drawCount)));
        _energyPerDraw[_drawCount]            = cycleEnergy;
        _magPerDraw[_drawCount]               = m;
        _sumSqClusterSizesPerDraw[_drawCount] = sumSq;
        _maxClusterSizePerDraw[_drawCount]    = (int)maxSize;
        _drawCount++;
    }

    /// <summary>
    /// Grow the per-draw record to hold <paramref name="capacity"/> draws, keeping the entries so far.
    /// </summary>
    private void ReserveRecord(int capacity)
    {
        if (capacity <= _energyPerDraw.Length) return;
        _energyPerDraw            = Grown(_energyPerDraw, capacity);
        _magPerDraw               = Grown(_magPerDraw, capacity);
        _sumSqClusterSizesPerDraw = Grown(_sumSqClusterSizesPerDraw, capacity);
        _maxClusterSizePerDraw    = Grown(_maxClusterSizePerDraw, capacity);
    }

    /// <summary>A copy of the live entries [0, draw count) in a series of the new capacity — nothing past them.</summary>
    private T[] Grown<T>(T[] series, int capacity)
    {
        var grown = new T[capacity];
        Array.Copy(series, grown, _drawCount);
        return grown;
    }

    // ── Checkpoint surface ────────────────────────────────────────────────

    public Accumulator GetCheckpoint()
    {
        var (s0, s1, s2, s3) = _rng.SaveState();

        return new Accumulator
        {
            Temperature              = _temperature,
            Q                        = _q,
            DrawCount                = _drawCount,
            Spins                    = (int[])_spins.Clone(),
            ClusterSizeHistogram     = (long[])_clusterSizeHistogram.Clone(),
            ExactClusterSizeHistogram = (long[]?)_exactClusterSizeHistogram?.Clone(),
            RngState0                = s0,
            RngState1                = s1,
            RngState2                = s2,
            RngState3                = s3,
            RunningSumSqClusterSizes     = _runningSumSqClusterSizes,
            RunningSumSqClusterSizesExcl = _runningSumSqClusterSizesExcl,
            RunningSumEnergy             = _runningSumEnergy,
            RunningSumEnergySq           = _runningSumEnergySq,
            RunningSumMag                = _runningSumMag,
            RunningSumMagSq              = _runningSumMagSq,
            EnergyPerDraw            = EnergyPerDraw.ToArray(),
            MagPerDraw               = MagPerDraw.ToArray(),
            SumSqClusterSizesPerDraw = SumSqClusterSizesPerDraw.ToArray(),
            MaxClusterSizePerDraw    = MaxClusterSizePerDraw.ToArray(),
            // The burn-in trace travels the same way — runtime-gated, null when off.
            BurnInEnergyPerDraw            = (double[]?)_burnInEnergyPerDraw?.Clone(),
            BurnInMagPerDraw               = (double[]?)_burnInMagPerDraw?.Clone(),
            BurnInSumSqClusterSizesPerDraw = (long[]?)_burnInSumSqClusterSizesPerDraw?.Clone(),
            BurnInMaxClusterSizePerDraw    = (int[]?)_burnInMaxClusterSizePerDraw?.Clone(),
            // Per-edge currency arrays travel in the same accumulator — each null when
            // its currency gate is off (Affinities ⟂ Alignments ⟂ CoMembership).
            BondFormedCount    = TConfig.Affinities ? (int[])_bondFormedCount!.Clone() : null,
            SpinAgreementCount = TConfig.Alignments ? (int[])_spinAgreementCount!.Clone() : null,
            CoMembershipCount  = _coMembershipCount is null ? null : (int[])_coMembershipCount.Clone(),
            // Per-node landscapes travel the same way — runtime-gated, null when off.
            SumClusterSizePerNode    = _sumClusterSizePerNode    is null ? null : (double[])_sumClusterSizePerNode.Clone(),
            SumInGiantClusterPerNode = _sumInGiantClusterPerNode is null ? null : (double[])_sumInGiantClusterPerNode.Clone(),
        };
    }

    public void Restore(Accumulator result)
    {
        // Capability validation. Each per-edge currency array's presence IS its gate
        // now: a snapshot taken with a currency can only restore into a specialization
        // that tracks it, and vice versa — validated independently.
        if (TConfig.Affinities != (result.BondFormedCount is not null))
            throw new InvalidOperationException(
                $"Snapshot Affinities presence ({result.BondFormedCount is not null}) does not " +
                $"match model capability ({TConfig.Affinities}).");
        if (TConfig.Alignments != (result.SpinAgreementCount is not null))
            throw new InvalidOperationException(
                $"Snapshot Alignments presence ({result.SpinAgreementCount is not null}) does not " +
                $"match model capability ({TConfig.Alignments}).");
        // Per-node landscapes gate on presence too: a snapshot taken with a
        // landscape on can only restore into a model tracking that landscape.
        if ((_sumClusterSizePerNode is not null) != (result.SumClusterSizePerNode is not null))
            throw new InvalidOperationException(
                $"Snapshot per-node cluster-size presence ({result.SumClusterSizePerNode is not null}) does not " +
                $"match model capability ({_sumClusterSizePerNode is not null}).");
        if ((_sumInGiantClusterPerNode is not null) != (result.SumInGiantClusterPerNode is not null))
            throw new InvalidOperationException(
                $"Snapshot per-node giant-participation presence ({result.SumInGiantClusterPerNode is not null}) does not " +
                $"match model capability ({_sumInGiantClusterPerNode is not null}).");
        // The burn-in trace gates on presence too, its four series present together with one length.
        int? traceLength = result.BurnInEnergyPerDraw?.Length;
        if (result.BurnInMagPerDraw?.Length != traceLength
            || result.BurnInSumSqClusterSizesPerDraw?.Length != traceLength
            || result.BurnInMaxClusterSizePerDraw?.Length != traceLength)
            throw new InvalidOperationException(
                $"Snapshot burn-in trace series disagree (E {result.BurnInEnergyPerDraw?.Length.ToString() ?? "null"}, " +
                $"m {result.BurnInMagPerDraw?.Length.ToString() ?? "null"}, Σ|c|² {result.BurnInSumSqClusterSizesPerDraw?.Length.ToString() ?? "null"}, " +
                $"c_max {result.BurnInMaxClusterSizePerDraw?.Length.ToString() ?? "null"}).");
        if ((_burnInEnergyPerDraw is not null) != traceLength.HasValue)
            throw new InvalidOperationException(
                $"Snapshot burn-in trace presence ({traceLength.HasValue}) does not " +
                $"match model capability ({_burnInEnergyPerDraw is not null}).");
        // The exact histogram gates on presence, with one count per size 0..N.
        if ((_exactClusterSizeHistogram is not null) != (result.ExactClusterSizeHistogram is not null))
            throw new InvalidOperationException(
                $"Snapshot exact cluster-size histogram presence ({result.ExactClusterSizeHistogram is not null}) does not " +
                $"match model capability ({_exactClusterSizeHistogram is not null}).");
        if (_exactClusterSizeHistogram is not null && result.ExactClusterSizeHistogram!.Length != _exactClusterSizeHistogram.Length)
            throw new InvalidOperationException(
                $"Snapshot exact cluster-size histogram length ({result.ExactClusterSizeHistogram.Length}) does not " +
                $"match model length N + 1 ({_exactClusterSizeHistogram.Length}).");
        if (result.Temperature != _temperature)
            throw new InvalidOperationException(
                $"Snapshot temperature {result.Temperature} does not match model temperature {_temperature}.");
        if (result.Q != _q)
            throw new InvalidOperationException(
                $"Snapshot Q {result.Q} does not match model Q {_q}.");
        if (result.Spins is null || result.Spins.Length != _n)
            throw new InvalidOperationException(
                $"Snapshot spin count {result.Spins?.Length ?? 0} does not match model node count {_n}.");
        if (result.ClusterSizeHistogram is null || result.ClusterSizeHistogram.Length != _clusterSizeHistogram.Length)
            throw new InvalidOperationException(
                $"Snapshot histogram bin count {result.ClusterSizeHistogram?.Length ?? 0} does not match model bin count {_clusterSizeHistogram.Length}.");
        // The per-draw record holds one entry per draw, so a restored run appends at DrawCount.
        if (result.EnergyPerDraw?.Length != result.DrawCount
            || result.MagPerDraw?.Length != result.DrawCount
            || result.SumSqClusterSizesPerDraw?.Length != result.DrawCount
            || result.MaxClusterSizePerDraw?.Length != result.DrawCount)
            throw new InvalidOperationException(
                $"Snapshot per-draw record lengths (E {result.EnergyPerDraw?.Length ?? 0}, m {result.MagPerDraw?.Length ?? 0}, " +
                $"Σ|c|² {result.SumSqClusterSizesPerDraw?.Length ?? 0}, c_max {result.MaxClusterSizePerDraw?.Length ?? 0}) " +
                $"do not match its draw count {result.DrawCount}.");

        // Copy mutable state. Array.Copy because we want elementwise copy,
        // not reference assignment — the DTO keeps its arrays, the runner
        // keeps its arrays, neither aliases the other.
        Array.Copy(result.Spins, _spins, _n);
        Array.Copy(result.ClusterSizeHistogram, _clusterSizeHistogram, _clusterSizeHistogram.Length);
        if (_exactClusterSizeHistogram is not null)
            Array.Copy(result.ExactClusterSizeHistogram!, _exactClusterSizeHistogram, _exactClusterSizeHistogram.Length);

        _rng.LoadState(result.RngState0, result.RngState1, result.RngState2, result.RngState3);

        _drawCount                   = result.DrawCount;
        _runningSumSqClusterSizes     = result.RunningSumSqClusterSizes;
        _runningSumSqClusterSizesExcl = result.RunningSumSqClusterSizesExcl;
        _runningSumEnergy             = result.RunningSumEnergy;
        _runningSumEnergySq           = result.RunningSumEnergySq;
        _runningSumMag                = result.RunningSumMag;
        _runningSumMagSq              = result.RunningSumMagSq;

        // The record is copied, never aliased; the next Run reserves past it.
        _energyPerDraw            = (double[])result.EnergyPerDraw.Clone();
        _magPerDraw               = (double[])result.MagPerDraw.Clone();
        _sumSqClusterSizesPerDraw = (long[])result.SumSqClusterSizesPerDraw.Clone();
        _maxClusterSizePerDraw    = (int[])result.MaxClusterSizePerDraw.Clone();
        if (_burnInEnergyPerDraw is not null)
        {
            _burnInEnergyPerDraw            = (double[])result.BurnInEnergyPerDraw!.Clone();
            _burnInMagPerDraw               = (double[])result.BurnInMagPerDraw!.Clone();
            _burnInSumSqClusterSizesPerDraw = (long[])result.BurnInSumSqClusterSizesPerDraw!.Clone();
            _burnInMaxClusterSizePerDraw    = (int[])result.BurnInMaxClusterSizePerDraw!.Clone();
        }

        // Restore the per-edge counts elementwise (presence already validated
        // above). Without this copy a resumed run keeps the snapshot's DrawCount
        // but zeroed bond/agreement counts — silently corrupted frequencies.
        if (TConfig.Affinities)
        {
            if (result.BondFormedCount!.Length != _bondFormedCount!.Length)
                throw new InvalidOperationException(
                    $"Snapshot Affinities array length ({result.BondFormedCount.Length}) does not " +
                    $"match model edge count ({_bondFormedCount.Length}).");
            Array.Copy(result.BondFormedCount, _bondFormedCount, _bondFormedCount.Length);
        }
        if (TConfig.Alignments)
        {
            if (result.SpinAgreementCount!.Length != _spinAgreementCount!.Length)
                throw new InvalidOperationException(
                    $"Snapshot Alignments array length ({result.SpinAgreementCount.Length}) does not " +
                    $"match model edge count ({_spinAgreementCount.Length}).");
            Array.Copy(result.SpinAgreementCount, _spinAgreementCount, _spinAgreementCount.Length);
        }
        if (_coMembershipCount is not null)
        {
            if (result.CoMembershipCount is null || result.CoMembershipCount.Length != _coMembershipCount.Length)
                throw new InvalidOperationException(
                    $"Snapshot CoMembership presence/length mismatch (snapshot={result.CoMembershipCount?.Length.ToString() ?? "null"}, " +
                    $"model={_coMembershipCount.Length}).");
            Array.Copy(result.CoMembershipCount, _coMembershipCount, _coMembershipCount.Length);
        }

        // Restore the per-node landscape sums elementwise (presence already
        // validated above). Faithful resume: the running fields continue folding
        // from where the snapshot left off, matching the scalar-accumulator treatment.
        if (_sumClusterSizePerNode is not null)
            Array.Copy(result.SumClusterSizePerNode!, _sumClusterSizePerNode, _n);
        if (_sumInGiantClusterPerNode is not null)
            Array.Copy(result.SumInGiantClusterPerNode!, _sumInGiantClusterPerNode, _n);
    }

    /// <summary>
    /// Zero every running accumulator and reset the cycle counter. Called
    /// at the tail of <see cref="BurnIn"/>; the model continues to hold
    /// post-burn-in spins and RNG state, but starts collecting statistics
    /// fresh from this point. The per-draw record restarts with the counter
    /// (its length is the draw count), so the measured series begins here.
    /// </summary>
    private void ResetAccumulators()
    {
        _runningSumSqClusterSizes     = 0.0;
        _runningSumSqClusterSizesExcl = 0.0;
        _runningSumEnergy             = 0.0;
        _runningSumEnergySq           = 0.0;
        _runningSumMag                = 0.0;
        _runningSumMagSq              = 0.0;
        _drawCount                   = 0;
        _clusterSizeHistogram.AsSpan().Clear();
        _exactClusterSizeHistogram?.AsSpan().Clear();
        if (TConfig.Affinities) _bondFormedCount!.AsSpan().Clear();
        if (TConfig.Alignments) _spinAgreementCount!.AsSpan().Clear();
        _coMembershipCount?.AsSpan().Clear();
        _sumClusterSizePerNode?.AsSpan().Clear();
        _sumInGiantClusterPerNode?.AsSpan().Clear();
    }

    /// <summary>
    /// Returns the size of the largest color bucket in <paramref name="spins"/>.
    /// O(N + Q) with Q small (typically ≤ 20).
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int MaxColorBucket(ReadOnlySpan<int> spins, int q, Span<int> scratch)
    {
        scratch.Clear();
        for (int i = 0; i < spins.Length; i++)
            scratch[spins[i]]++;

        int max = 0;
        for (int c = 0; c < q; c++)
            if (scratch[c] > max) max = scratch[c];
        return max;
    }
}
