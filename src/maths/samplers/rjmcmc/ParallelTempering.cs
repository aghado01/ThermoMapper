using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Maths.Samplers.Ensemble;

namespace Maths.Samplers.Rjmcmc;

/// <summary>
/// A ladder's restart state: every replica's <see cref="ChainSnapshot{TState}"/> (its rung is its temperature, so a
/// resized ladder's rungs ride here), the swap stream, the scan count, which state each rung holds, the round-trip
/// bookkeeping and the swap statistics with the scans they cover, including those gathered since the last adaptation,
/// the growth stream new replicas draw their streams from, for a ladder with an exact reference rung, that rung's
/// stream, for a ladder tracing its rungs, copies of the traces, and, for a ladder with per-replica blocks, each
/// replica's block stream and block counters.
/// </summary>
public sealed record LadderSnapshot<TState>(
    IReadOnlyList<ChainSnapshot<TState>> Replicas, TemperingPath Path, SwapSchedule Schedule,
    ulong SwapRng0, ulong SwapRng1, ulong SwapRng2, ulong SwapRng3,
    long Scans, long SwapAttempts, long SwapAccepts, long RoundTrips, long AnnealedRestarts,
    int[] LabelAt, sbyte[] LastEnd, long[] PairAttempts, long[] PairAccepts, double[] RejectionSum, long PairScans,
    double[] RoundRejectionSum, long RoundScans, (ulong, ulong, ulong, ulong) GrowthRng, (ulong, ulong, ulong, ulong)? ReferenceRng = null,
    IReadOnlyList<ChainTrace>? RungTraces = null, bool TracesSlopes = false,
    IReadOnlyList<(ulong, ulong, ulong, ulong)>? BlockRngs = null, long[][]? BlockAttempts = null,
    long[][]? BlockAccepts = null);

/// <summary>
/// The functionals a ladder records at every rung once it traces them (<see cref="ParallelTempering{TState}.TraceRungs"/>),
/// the per-rung draws a ladder evidence estimate reads (<see cref="Evidence.LadderEvidence"/>).
/// </summary>
public static class RungFunctional
{
    /// <summary>∂_b log q_b at the rung's own b (<see cref="Temperature.LogDensitySlope"/>), path sampling's integrand;
    /// NaN when the target supplies no slope.</summary>
    public const int Slope = 0;

    /// <summary><c>log q_c − log q_b</c> with c the next colder rung's b, stepping-stone's log importance ratio; NaN at
    /// the cold rung.</summary>
    public const int LogRatio = 1;

    internal const int Count = 2;
}

/// <summary>
/// Parallel tempering (replica exchange) over a ladder of reversible-jump chains at inverse temperatures
/// 1 = b₀ &gt; b₁ &gt; … &gt; b_{L-1}, all on one declared <see cref="TemperingPath"/>: the powered posterior
/// <c>(prior · likelihood)^b</c> or the likelihood-tempered path <c>prior · likelihood^b</c>. After every scan, adjacent
/// replicas attempt Metropolis state swaps on the pairs the <see cref="SwapSchedule"/> chooses, each replica's density
/// evaluated at the other's state from the cached evaluations — so the hot replicas, which cross modal barriers freely,
/// ferry the cold replica (b = 1, the true posterior) between modes. The cold chain is the sampler; the rest are
/// scaffolding. For multimodal targets (e.g. a multi-peak T_c posterior) where a single chain would lock into one mode.
/// <para>The ladder reports what Syed et al. (2022) measure it by: round trips of states from the hot end to the cold
/// end and back, and each pair's rejection rate, from which <see cref="AdaptLadder"/> equalizes the rungs.</para>
/// <para>With an <see cref="IReferenceSampler{TState}"/> the hot end is b = 0 and its rung takes an exact independent
/// draw from the reference every scan instead of a local step, so each state reaching the hot end is regenerated
/// (OQ2). On the powered path this is the only way to place a rung at b = 0, the normalized reference measure.</para>
/// <para>After burn-in a ladder can trace its rungs (<see cref="TraceRungs"/>): each rung's slope and log ratio to the
/// next rung, every scan, from the evaluations it already holds, for an evidence over the ladder.</para>
/// <para>A replica can also run blocks between its moves and the swaps: <see cref="IBlock{TState}"/>s over the replica
/// chain, each leaving the replica's tempered target at its <see cref="ReversibleJumpChain{TState}.Temperature"/>
/// invariant and handing back the new state with its evaluation (<see cref="ReversibleJumpChain{TState}.SetState"/>),
/// such as a slice step on a continuous coordinate that no reversible-jump move updates. Each replica's blocks draw from
/// a stream of their own, after the ladder's others in the <see cref="SeedTree"/>, so a ladder without blocks keeps every
/// stream. The same blocks serve a lone chain through a <see cref="BlockScan{TState}"/> beside its
/// <see cref="MoveStep{TState}"/>.</para>
/// <para>During burn-in a ladder can be resized (<see cref="Resize"/>) to rungs its caller places, by the barrier
/// (<see cref="TemperingLadder.Place"/>) or for a readout's resolution: kept replicas keep their states and streams, and
/// new ones start from forks of their sources' states on streams drawn from one growth stream, the last in the
/// <see cref="SeedTree"/>, so a ladder that never resizes keeps every stream.</para>
/// As an <see cref="IChain{TState}"/> the ladder presents as a single chain whose draw is the cold state and whose
/// acceptance is the cold replica's, so the ensemble shell can drive a tempered replica like any other.
/// </summary>
public sealed class ParallelTempering<TState> : IChain<TState>
{
    private ReversibleJumpChain<TState>[] _chains;
    private readonly ISamplingTarget<TState> _target;   // new replicas' (Resize)
    private readonly Xoshiro256PlusPlus _swapRng;
    private readonly Xoshiro256PlusPlus _growthRng;     // new replicas' move and block streams are drawn from it
    private readonly IReferenceSampler<TState>? _reference;
    private readonly Xoshiro256PlusPlus? _referenceRng;
    private double[] _logA;                        // this scan's log swap ratio per pair
    private long[] _pairAttempts, _pairAccepts;
    private double[] _rejectionSum;                // Σ (1 − α) per pair, every scan since the rungs' count was set
    private long _pairScans;                       // the scans it covers
    private double[] _roundRejectionSum;           // the same since the last adaptation or resize
    private long _roundScans;
    private int[] _labelAt;                        // which state (by its starting rung) each rung holds
    private sbyte[] _lastEnd;                      // per state: +1 last at the hot end, −1 at the cold end, 0 neither yet
    private ChainTrace[]? _rungTraces;             // per rung, cold first, once traced
    private bool _tracesSlopes;
    private IReadOnlyList<IBlock<ReversibleJumpChain<TState>>>[]? _blocks;   // per replica, when given
    private Xoshiro256PlusPlus[]? _blockRngs;
    private long[][]? _blockAttempts, _blockAccepts;

    /// <summary>Adjacent-pair swap attempts so far.</summary>
    public long SwapAttempts { get; private set; }

    /// <summary>Adjacent-pair swaps accepted so far; <c>SwapAccepts / SwapAttempts</c> is the exchange rate.</summary>
    public long SwapAccepts { get; private set; }

    /// <param name="movesForReplica">Builds replica i's own move set (i = 0 is the cold replica). It is called once
    /// per replica, so any mutable move or kernel state is private to that replica.</param>
    /// <param name="target">The target; replica i samples it at <c>Temperature(b_i, path)</c>.</param>
    /// <param name="start">Common initial state for all replicas.</param>
    /// <param name="inverseTemperatures">The ladder: finite, strictly descending, cold value exactly 1 at index 0. With a
    /// <paramref name="reference"/> the hot end is exactly 0. Without one the hot end is above 0 on the powered path,
    /// where b = 0 is the reference measure itself; on the likelihood path b = 0 is the prior, sampled by local moves,
    /// which a target with an improper nuisance prior refuses. The array is not retained.</param>
    /// <param name="path">The tempering path every replica takes.</param>
    /// <param name="masterSeed">Seeds the per-replica streams, the swap stream, the reference rung's stream, the
    /// per-replica block streams and the growth stream via the <see cref="SeedTree"/>; each prefix is that of a ladder
    /// without what follows.</param>
    /// <param name="schedule">Which pairs each scan attempts; deterministic even-odd by default.</param>
    /// <param name="reference">Exact draws from the target's b = 0 member on <paramref name="path"/>, which the hot rung
    /// takes every scan in place of its local moves.</param>
    /// <param name="blocksForReplica">Builds replica i's blocks, run in order after its moves on every scan their cadence
    /// (<see cref="IBlock{TState}.Every"/>) falls due, on a stream of their own. It is called once per replica, as
    /// <paramref name="movesForReplica"/> is. The reference rung's exact draw runs none.</param>
    public ParallelTempering(Func<int, IReadOnlyList<IRjMove<TState>>> movesForReplica, ISamplingTarget<TState> target,
                             TState start, double[] inverseTemperatures, TemperingPath path, int masterSeed,
                             SwapSchedule schedule = SwapSchedule.DeterministicEvenOdd, IReferenceSampler<TState>? reference = null,
                             Func<int, IReadOnlyList<IBlock<ReversibleJumpChain<TState>>>>? blocksForReplica = null)
    {
        ArgumentNullException.ThrowIfNull(movesForReplica);
        ArgumentNullException.ThrowIfNull(target);
        ValidateLadder(inverseTemperatures, path, reference);
        if (schedule is not (SwapSchedule.DeterministicEvenOdd or SwapSchedule.StochasticEvenOdd or SwapSchedule.Sequential))
            throw new ArgumentOutOfRangeException(nameof(schedule));

        Path = path;
        Schedule = schedule;
        _target = target;
        int levels = inverseTemperatures.Length;
        // Prefix-stable: the replica streams, the swap stream and the reference stream come first, as without blocks, and
        // the growth stream last, which a ladder that never resizes never draws.
        int[] seeds = SeedTree.Derive(masterSeed, levels + 2 + (blocksForReplica is null ? 0 : levels) + 1);
        _chains = new ReversibleJumpChain<TState>[levels];
        for (int i = 0; i < levels; i++)
            _chains[i] = new ReversibleJumpChain<TState>(
                movesForReplica(i), target, start, new Xoshiro256PlusPlus(seeds[i]), new Temperature(inverseTemperatures[i], path));
        _swapRng = new Xoshiro256PlusPlus(seeds[levels]);
        _growthRng = new Xoshiro256PlusPlus(seeds[^1]);
        if (blocksForReplica is not null)
        {
            (_blocks, _blockAttempts, _blockAccepts) = BuildBlocks(blocksForReplica, levels);
            _blockRngs = new Xoshiro256PlusPlus[levels];
            for (int i = 0; i < levels; i++) _blockRngs[i] = new Xoshiro256PlusPlus(seeds[levels + 2 + i]);
        }
        if (reference is not null)
        {
            _reference = reference;
            _referenceRng = new Xoshiro256PlusPlus(seeds[levels + 1]);
            (TState x, ITargetEvaluation e) = reference.Draw(_referenceRng);   // the hot rung starts on an exact draw
            _chains[levels - 1].SetState(x, e);
        }
        Restart(levels);
    }

    // A ladder rebuilt from a snapshot: the replicas and the swap and growth streams are restored, the bookkeeping copied
    // after.
    private ParallelTempering(ReversibleJumpChain<TState>[] chains, ISamplingTarget<TState> target, Xoshiro256PlusPlus swapRng,
                              Xoshiro256PlusPlus growthRng, TemperingPath path, SwapSchedule schedule, IReferenceSampler<TState>? reference,
                              Xoshiro256PlusPlus? referenceRng, Func<int, IReadOnlyList<IBlock<ReversibleJumpChain<TState>>>>? blocksForReplica = null,
                              Xoshiro256PlusPlus[]? blockRngs = null)
    {
        _chains = chains;
        _target = target;
        _swapRng = swapRng;
        _growthRng = growthRng;
        _reference = reference;
        _referenceRng = referenceRng;
        Path = path;
        Schedule = schedule;
        int levels = chains.Length;
        if (blocksForReplica is not null)
        {
            (_blocks, _blockAttempts, _blockAccepts) = BuildBlocks(blocksForReplica, levels);
            _blockRngs = blockRngs;
        }
        Restart(levels);
    }

    /// <summary>This ladder's restart state; the ladder itself is unchanged.</summary>
    public LadderSnapshot<TState> Snapshot()
    {
        var replicas = new ChainSnapshot<TState>[_chains.Length];
        for (int i = 0; i < replicas.Length; i++) replicas[i] = _chains[i].Snapshot();
        (ulong s0, ulong s1, ulong s2, ulong s3) = _swapRng.SaveState();
        return new LadderSnapshot<TState>(replicas, Path, Schedule, s0, s1, s2, s3, Scans, SwapAttempts, SwapAccepts,
            RoundTrips, AnnealedRestarts, (int[])_labelAt.Clone(), (sbyte[])_lastEnd.Clone(), (long[])_pairAttempts.Clone(),
            (long[])_pairAccepts.Clone(), (double[])_rejectionSum.Clone(), _pairScans, (double[])_roundRejectionSum.Clone(), _roundScans,
            _growthRng.SaveState(), _referenceRng?.SaveState(), _rungTraces is null ? null : Array.ConvertAll(_rungTraces, t => t.Copy()), _tracesSlopes,
            _blockRngs is null ? null : Array.ConvertAll(_blockRngs, r => r.SaveState()),
            _blockAttempts is null ? null : Array.ConvertAll(_blockAttempts, a => (long[])a.Clone()),
            _blockAccepts is null ? null : Array.ConvertAll(_blockAccepts, a => (long[])a.Clone()));
    }

    /// <summary>
    /// A ladder that continues exactly where <paramref name="snapshot"/>'s left off, over the move sets, target,
    /// reference and blocks it ran with, on the rungs it held, after any resize. Each replica's current state is evaluated
    /// once, against <paramref name="target"/>.
    /// </summary>
    public static ParallelTempering<TState> Restore(Func<int, IReadOnlyList<IRjMove<TState>>> movesForReplica,
                                                    ISamplingTarget<TState> target, LadderSnapshot<TState> snapshot,
                                                    IReferenceSampler<TState>? reference = null,
                                                    Func<int, IReadOnlyList<IBlock<ReversibleJumpChain<TState>>>>? blocksForReplica = null)
    {
        ArgumentNullException.ThrowIfNull(movesForReplica);
        ArgumentNullException.ThrowIfNull(snapshot);
        if ((reference is null) != (snapshot.ReferenceRng is null))
            throw new ArgumentException("The ladder ran " + (reference is null ? "with" : "without") + " a reference rung.", nameof(reference));
        if ((blocksForReplica is null) != (snapshot.BlockRngs is null))
            throw new ArgumentException("The ladder ran " + (blocksForReplica is null ? "with" : "without") + " per-replica blocks.",
                                        nameof(blocksForReplica));
        int levels = snapshot.Replicas.Count;
        var chains = new ReversibleJumpChain<TState>[levels];
        for (int i = 0; i < levels; i++)
            chains[i] = ReversibleJumpChain<TState>.Restore(movesForReplica(i), target, snapshot.Replicas[i]);
        Xoshiro256PlusPlus? referenceRng = snapshot.ReferenceRng is (ulong, ulong, ulong, ulong) r
            ? Xoshiro256PlusPlus.Restore(r.Item1, r.Item2, r.Item3, r.Item4)
            : null;
        Xoshiro256PlusPlus[]? blockRngs = null;
        if (snapshot.BlockRngs is IReadOnlyList<(ulong, ulong, ulong, ulong)> words)
        {
            if (words.Count != levels || snapshot.BlockAttempts?.Length != levels || snapshot.BlockAccepts?.Length != levels)
                throw new ArgumentException("One block stream and one set of block counters per replica are required.", nameof(snapshot));
            blockRngs = new Xoshiro256PlusPlus[levels];
            for (int i = 0; i < levels; i++) blockRngs[i] = Xoshiro256PlusPlus.Restore(words[i].Item1, words[i].Item2, words[i].Item3, words[i].Item4);
        }
        (ulong g0, ulong g1, ulong g2, ulong g3) = snapshot.GrowthRng;
        var ladder = new ParallelTempering<TState>(chains, target,
            Xoshiro256PlusPlus.Restore(snapshot.SwapRng0, snapshot.SwapRng1, snapshot.SwapRng2, snapshot.SwapRng3),
            Xoshiro256PlusPlus.Restore(g0, g1, g2, g3), snapshot.Path, snapshot.Schedule, reference, referenceRng, blocksForReplica, blockRngs)
        {
            Scans = snapshot.Scans,
            SwapAttempts = snapshot.SwapAttempts,
            SwapAccepts = snapshot.SwapAccepts,
            RoundTrips = snapshot.RoundTrips,
            AnnealedRestarts = snapshot.AnnealedRestarts,
        };
        snapshot.LabelAt.CopyTo(ladder._labelAt, 0);
        snapshot.LastEnd.CopyTo(ladder._lastEnd, 0);
        snapshot.PairAttempts.CopyTo(ladder._pairAttempts, 0);
        snapshot.PairAccepts.CopyTo(ladder._pairAccepts, 0);
        snapshot.RejectionSum.CopyTo(ladder._rejectionSum, 0);
        ladder._pairScans = snapshot.PairScans;
        snapshot.RoundRejectionSum.CopyTo(ladder._roundRejectionSum, 0);
        ladder._roundScans = snapshot.RoundScans;
        if (snapshot.RungTraces is IReadOnlyList<ChainTrace> traces)
        {
            if (traces.Count != levels) throw new ArgumentException("One rung trace per replica is required.", nameof(snapshot));
            ladder._rungTraces = new ChainTrace[levels];
            for (int i = 0; i < levels; i++) ladder._rungTraces[i] = traces[i].Copy();
            ladder._tracesSlopes = snapshot.TracesSlopes;
        }
        if (ladder._blocks is not null)
            for (int i = 0; i < levels; i++)
            {
                if (snapshot.BlockAttempts![i].Length != ladder._blocks[i].Count || snapshot.BlockAccepts![i].Length != ladder._blocks[i].Count)
                    throw new ArgumentException($"Replica {i}'s snapshot was taken over a different block set.", nameof(snapshot));
                snapshot.BlockAttempts[i].CopyTo(ladder._blockAttempts![i], 0);
                snapshot.BlockAccepts[i].CopyTo(ladder._blockAccepts![i], 0);
            }
        return ladder;
    }

    /// <summary>The path every replica takes.</summary>
    public TemperingPath Path { get; }

    /// <summary>Whether the hot rung draws exactly from a reference at b = 0.</summary>
    public bool HasReferenceRung => _reference is not null;

    /// <summary>Which pairs each scan attempts.</summary>
    public SwapSchedule Schedule { get; }

    /// <summary>The current cold-chain (b = 1) state — the posterior sample.</summary>
    public TState Cold => _chains[0].Current;

    /// <summary>The cold replica (b = 1) itself — the chain whose draws are the posterior, exposed so a consumer
    /// can read its acceptance / per-move diagnostics. Its state advances both by its own moves and by swaps.</summary>
    public ReversibleJumpChain<TState> ColdChain => _chains[0];

    /// <summary>Every replica, cold first, for per-rung diagnostics. Each one's <see cref="ReversibleJumpChain{TState}.Temperature"/>
    /// is its rung. Stepping a replica directly takes it out of the ladder's schedule.</summary>
    public IReadOnlyList<ReversibleJumpChain<TState>> Replicas => _chains;

    /// <summary>The current ladder, cold first.</summary>
    public double[] InverseTemperatures
    {
        get
        {
            var b = new double[_chains.Length];
            for (int i = 0; i < b.Length; i++) b[i] = _chains[i].Temperature.InverseTemperature;
            return b;
        }
    }

    /// <summary>Scans so far: one <see cref="Step"/> each.</summary>
    public long Scans { get; private set; }

    /// <summary>Each rung's trace of the <see cref="RungFunctional"/>s, cold first, one draw per scan since
    /// <see cref="TraceRungs"/>; null until then.</summary>
    public IReadOnlyList<ChainTrace>? RungTraces => _rungTraces;

    /// <summary>Whether the traces hold slopes: when asked for (<see cref="TraceRungs"/>), always on the powered path,
    /// and on the likelihood path when the target supplies ∂ℓ_b/∂b.</summary>
    public bool TracesSlopes => _rungTraces is not null && _tracesSlopes;

    /// <summary>
    /// Opt in to per-rung traces, after burn-in: from the next scan on, every <see cref="Step"/> records at each rung, on
    /// the state it holds after the scan's swaps, the <see cref="RungFunctional"/>s, read from the evaluations the rungs
    /// already hold. With <paramref name="traceSlopes"/> false the slopes are not read and trace as NaN, for an estimator
    /// that reads only the log ratios: a slope costs the target a query per rung per scan. A second call starts the traces
    /// afresh. The ladder is then fixed, so <see cref="AdaptLadder"/> and <see cref="Resize"/> refuse. Whether the target
    /// supplies the slope is probed once, on the cold rung's state.
    /// </summary>
    public void TraceRungs(bool traceSlopes = true)
    {
        _tracesSlopes = traceSlopes && (Path == TemperingPath.Powered || SuppliesSlope(_chains[0]));
        _rungTraces = new ChainTrace[_chains.Length];
        for (int i = 0; i < _rungTraces.Length; i++) _rungTraces[i] = new ChainTrace(RungFunctional.Count);
    }

    /// <summary>
    /// Round trips so far, summed over states: a state that leaves the hot end, reaches the cold end and returns
    /// completes one (Syed et al.'s definition, with the hot end as the reference). <c>RoundTrips / Scans</c> is the
    /// round-trip rate, the measure of how well the ladder moves.
    /// </summary>
    public long RoundTrips { get; private set; }

    /// <summary>Arrivals at the cold end of states last seen at the hot end ("annealed restarts"): each carries what
    /// the hot rungs explored down to the posterior.</summary>
    public long AnnealedRestarts { get; private set; }

    /// <summary>Per adjacent pair, cold first: attempts, accepts and the mean rejection probability over every scan since
    /// the ladder took its count of rungs (<see cref="PairScans"/>): its start, or its last <see cref="Resize"/>.</summary>
    public IReadOnlyList<SwapPairStat> PairStats
    {
        get
        {
            var stats = new SwapPairStat[_chains.Length - 1];
            for (int i = 0; i < stats.Length; i++)
                stats[i] = new SwapPairStat(_chains[i].Temperature.InverseTemperature, _chains[i + 1].Temperature.InverseTemperature,
                                            _pairAttempts[i], _pairAccepts[i], _pairScans > 0 ? _rejectionSum[i] / _pairScans : double.NaN);
            return stats;
        }
    }

    /// <summary>The scans <see cref="PairStats"/> covers: every scan, or those since the last <see cref="Resize"/>.</summary>
    public long PairScans => _pairScans;

    /// <summary>
    /// Per adjacent pair, cold first, the mean rejection probability over the scans since the last adaptation or resize:
    /// the round <see cref="AdaptLadder"/> equalizes from, and a caller places a resize's rungs from
    /// (<see cref="TemperingLadder.Place"/>). NaN before a scan.
    /// </summary>
    public double[] RoundRejection
    {
        get
        {
            var rejection = new double[_chains.Length - 1];
            for (int i = 0; i < rejection.Length; i++) rejection[i] = _roundScans > 0 ? _roundRejectionSum[i] / _roundScans : double.NaN;
            return rejection;
        }
    }

    /// <summary>Whether each replica runs blocks between its moves and the swaps.</summary>
    public bool HasBlocks => _blocks is not null;

    /// <summary>Replica <paramref name="replica"/>'s blocks' cumulative attempts and acceptances, in order; empty for a
    /// ladder without blocks.</summary>
    public IReadOnlyList<MoveStat> BlockStats(int replica)
    {
        if ((uint)replica >= (uint)_chains.Length) throw new ArgumentOutOfRangeException(nameof(replica));
        if (_blocks is null) return [];
        IReadOnlyList<IBlock<ReversibleJumpChain<TState>>> blocks = _blocks[replica];
        var stats = new MoveStat[blocks.Count];
        for (int j = 0; j < stats.Length; j++) stats[j] = new MoveStat(blocks[j].Name, _blockAttempts![replica][j], _blockAccepts![replica][j]);
        return stats;
    }

    /// <summary>
    /// Advance every replica one step (local exploration) and run its blocks, then attempt the scan's swaps. A swap of
    /// x_i and x_j is accepted with <c>log α = d_i(x_j) + d_j(x_i) − d_i(x_i) − d_j(x_j)</c>, d_i replica i's log
    /// density, which on the powered path reduces to <c>(b_i − b_j)(L_j − L_i)</c>. Every pair's rejection probability
    /// 1 − α is recorded every scan, attempted or not. Returns the cold state.
    /// </summary>
    public TState Step()
    {
        int last = _chains.Length - 1;
        for (int i = 0; i < last; i++)
        {
            _chains[i].Step();
            RunBlocks(i);
        }
        if (_reference is null)
        {
            _chains[last].Step();
            RunBlocks(last);
        }
        else
        {
            (TState x, ITargetEvaluation e) = _reference.Draw(_referenceRng!);
            _chains[last].SetState(x, e);
        }
        Scans++;
        _pairScans++;
        _roundScans++;
        int pairs = _chains.Length - 1;

        if (Schedule == SwapSchedule.Sequential)
        {
            // Every pair, cold to hot; each sees the swaps before it.
            for (int i = 0; i < pairs; i++)
            {
                double logA = Record(i);
                Attempt(i, logA);
            }
        }
        else
        {
            // The even and the odd pairs are each disjoint, so this scan's ratios are all taken before any swap.
            for (int i = 0; i < pairs; i++) _logA[i] = Record(i);
            int parity = Schedule == SwapSchedule.DeterministicEvenOdd
                ? (int)(Scans % 2 == 0 ? 0 : 1)
                : _swapRng.NextDouble() < 0.5 ? 0 : 1;
            for (int i = parity; i < pairs; i += 2) Attempt(i, _logA[i]);
        }

        TrackEnds();
        if (_rungTraces is not null) RecordRungs();
        return _chains[0].Current;
    }

    /// <summary>
    /// Re-space the interior rungs so that adjacent rejection rates are equal (<see cref="TemperingLadder.Equalize"/>),
    /// from the rates observed since the last adaptation or resize, and return those rates' barrier estimate Λ̂. The cold and hot
    /// ends stay where they are. Retempering is a query on each replica's cached evaluation: no state is re-evaluated.
    /// Adapt during burn-in only; a ladder that keeps changing does not leave the posterior invariant.
    /// </summary>
    public double AdaptLadder()
    {
        if (_rungTraces is not null) throw new InvalidOperationException("The ladder is fixed while its rungs are traced.");
        if (_roundScans == 0) throw new InvalidOperationException("No scans since the last adaptation.");
        double[] rejection = RoundRejection;
        double[] next = TemperingLadder.Equalize(InverseTemperatures, rejection);
        ValidateLadder(next, Path, _reference);
        for (int i = 0; i < _chains.Length; i++) _chains[i].SetTemperature(new Temperature(next[i], Path));
        Array.Clear(_roundRejectionSum);
        _roundScans = 0;
        return TemperingLadder.Barrier(rejection);
    }

    /// <summary>
    /// Resize the ladder during burn-in to <paramref name="inverseTemperatures"/>, which its caller places: by the barrier
    /// (<see cref="TemperingLadder.Place"/> from <see cref="RoundRejection"/>, with <see cref="TemperingLadder.Sources"/>)
    /// or for a readout's resolution.
    /// <list type="bullet">
    /// <item><b>Kept replicas.</b> A replica that a new rung names in <paramref name="sources"/> is kept, with its state,
    /// streams, moves, blocks and counters, at the naming rung nearest its own b (the colder on a tie), and moves there as
    /// <see cref="AdaptLadder"/> moves rungs, by a query on its cached evaluation. A replica no rung names is dropped.</item>
    /// <item><b>New replicas.</b> Each further rung that names a replica starts a new one from a fork of its state, with
    /// moves and blocks built for it, and its move stream, then its block stream, seeded by the next draws of the growth
    /// stream, cold first. It samples the ladder's target, evaluated once at the fork.</item>
    /// <item><b>Statistics.</b> The pairs' swap statistics, the round-trip labels and the adaptation round restart, since
    /// the old pairs do not map onto the new ones; the totals (scans, swaps, round trips, annealed restarts) carry on.</item>
    /// </list>
    /// Refused while the rungs are traced, as <see cref="AdaptLadder"/> is; the burn-in's freeze is the caller's to keep.
    /// A refused resize leaves the ladder as it was. Returns, per new rung, cold first, the replica (by its old rung) kept
    /// there, or −1 where a new one started.
    /// </summary>
    /// <param name="inverseTemperatures">The new ladder, valid as the constructor's: strictly descending from exactly 1, and with a
    /// reference rung ending at exactly 0. The array is not retained.</param>
    /// <param name="sources">Per new rung, cold first, the current rung whose state it starts from: the cold rung for the
    /// cold rung and, with a reference rung, the reference rung for the last.</param>
    /// <param name="fork">An independent copy of a state, which can be mutated without touching the original, as
    /// resampling's offspring are (<see cref="Smc.ISmcParticles{TState}.Fork"/>).</param>
    /// <param name="movesForReplica">Builds a new replica's move set, called with its new rung, as the constructor's is.</param>
    /// <param name="blocksForReplica">Builds a new replica's blocks, as the constructor's does: given exactly when the
    /// ladder runs blocks.</param>
    public int[] Resize(double[] inverseTemperatures, int[] sources, Func<TState, TState> fork, Func<int, IReadOnlyList<IRjMove<TState>>> movesForReplica,
                        Func<int, IReadOnlyList<IBlock<ReversibleJumpChain<TState>>>>? blocksForReplica = null)
    {
        if (_rungTraces is not null) throw new InvalidOperationException("The ladder is fixed while its rungs are traced.");
        ArgumentNullException.ThrowIfNull(sources);
        ArgumentNullException.ThrowIfNull(fork);
        ArgumentNullException.ThrowIfNull(movesForReplica);
        ValidateLadder(inverseTemperatures, Path, _reference);
        if ((blocksForReplica is null) != (_blocks is null))
            throw new ArgumentException("The ladder runs " + (_blocks is null ? "no" : "per-replica") + " blocks.", nameof(blocksForReplica));
        int levels = inverseTemperatures.Length, current = _chains.Length;
        if (sources.Length != levels) throw new ArgumentException("One source per new rung is required.", nameof(sources));
        foreach (int s in sources)
            if ((uint)s >= (uint)current) throw new ArgumentOutOfRangeException(nameof(sources), $"The ladder has {current} rungs, not a rung {s}.");
        if (sources[0] != 0) throw new ArgumentException("The cold rung keeps the cold replica.", nameof(sources));
        if (_reference is not null && sources[^1] != current - 1)
            throw new ArgumentException("The reference rung keeps the reference rung's replica.", nameof(sources));

        // Each named replica is kept at the naming rung nearest its own b, the colder on a tie.
        var keptAt = new int[current];
        Array.Fill(keptAt, -1);
        for (int k = 0; k < levels; k++)
        {
            int s = sources[k];
            double b = _chains[s].Temperature.InverseTemperature;
            if (keptAt[s] < 0 || Math.Abs(inverseTemperatures[k] - b) < Math.Abs(inverseTemperatures[keptAt[s]] - b)) keptAt[s] = k;
        }
        var kept = new int[levels];
        Array.Fill(kept, -1);
        for (int s = 0; s < current; s++)
            if (keptAt[s] >= 0) kept[keptAt[s]] = s;

        var chains = new ReversibleJumpChain<TState>[levels];
        IReadOnlyList<IBlock<ReversibleJumpChain<TState>>>[]? blocks = _blocks is null ? null : new IReadOnlyList<IBlock<ReversibleJumpChain<TState>>>[levels];
        Xoshiro256PlusPlus[]? blockRngs = _blocks is null ? null : new Xoshiro256PlusPlus[levels];
        long[][]? blockAttempts = _blocks is null ? null : new long[levels][], blockAccepts = _blocks is null ? null : new long[levels][];
        (ulong, ulong, ulong, ulong) growth = _growthRng.SaveState();
        try
        {
            for (int k = 0; k < levels; k++)
            {
                if (kept[k] is int s and >= 0)
                {
                    chains[k] = _chains[s];
                    if (blocks is not null)
                        (blocks[k], blockRngs![k], blockAttempts![k], blockAccepts![k]) = (_blocks![s], _blockRngs![s], _blockAttempts![s], _blockAccepts![s]);
                    continue;
                }
                TState start = fork(_chains[sources[k]].Current);
                var moveRng = new Xoshiro256PlusPlus(NextGrowthSeed());
                chains[k] = new ReversibleJumpChain<TState>(movesForReplica(k), _target, start, moveRng, new Temperature(inverseTemperatures[k], Path));
                if (blocks is null) continue;
                blockRngs![k] = new Xoshiro256PlusPlus(NextGrowthSeed());
                blocks[k] = BlockSet(blocksForReplica!, k);
                blockAttempts![k] = new long[blocks[k].Count];
                blockAccepts![k] = new long[blocks[k].Count];
            }
        }
        catch
        {
            _growthRng.LoadState(growth.Item1, growth.Item2, growth.Item3, growth.Item4);
            throw;
        }

        for (int k = 0; k < levels; k++)
            if (kept[k] >= 0) chains[k].SetTemperature(new Temperature(inverseTemperatures[k], Path));
        _chains = chains;
        (_blocks, _blockRngs, _blockAttempts, _blockAccepts) = (blocks, blockRngs, blockAttempts, blockAccepts);
        Restart(levels);
        return kept;
    }

    // IChain view: the ladder advances by one swept Step; its acceptance is the cold replica's own move acceptance
    // (swaps move state via SetState and are tracked separately as SwapAccepts/SwapAttempts).
    long IChain<TState>.Accepted => _chains[0].Accepted;
    long IChain<TState>.Attempts => _chains[0].Attempts;

    /// <summary>Geometric inverse-temperature ladder from 1 (cold) down to <paramref name="hotEnd"/>.</summary>
    public static double[] GeometricLadder(int levels, double hotEnd)
    {
        if (levels < 2) throw new ArgumentOutOfRangeException(nameof(levels), "Need at least two levels.");
        if (!(hotEnd > 0.0 && hotEnd < 1.0)) throw new ArgumentOutOfRangeException(nameof(hotEnd), "The hot end must lie in (0, 1).");
        var b = new double[levels];
        for (int i = 0; i < levels; i++)
            b[i] = Math.Pow(hotEnd, i / (double)(levels - 1));   // 1 … the hot end
        return b;
    }

    /// <summary>
    /// A ladder ending at an exact reference rung: <paramref name="levels"/> − 1 geometric rungs from 1 down to
    /// <paramref name="hotEnd"/>, then b = 0.
    /// </summary>
    public static double[] ReferenceLadder(int levels, double hotEnd)
    {
        if (levels < 3) throw new ArgumentOutOfRangeException(nameof(levels), "A reference ladder needs at least three levels.");
        double[] geometric = GeometricLadder(levels - 1, hotEnd);
        var b = new double[levels];
        geometric.CopyTo(b, 0);
        return b;
    }

    // The per-pair and per-state bookkeeping of a ladder of `levels` rungs, from scratch: no swap statistics, no round yet,
    // and every state at its starting rung, the hottest's first round trip starting at its first visit to the hot end.
    [MemberNotNull(nameof(_logA), nameof(_pairAttempts), nameof(_pairAccepts), nameof(_rejectionSum), nameof(_roundRejectionSum),
                   nameof(_labelAt), nameof(_lastEnd))]
    private void Restart(int levels)
    {
        _logA = new double[levels - 1];
        _pairAttempts = new long[levels - 1];
        _pairAccepts = new long[levels - 1];
        _rejectionSum = new double[levels - 1];
        _pairScans = 0;
        _roundRejectionSum = new double[levels - 1];
        _roundScans = 0;
        _labelAt = new int[levels];
        _lastEnd = new sbyte[levels];
        for (int i = 0; i < levels; i++) _labelAt[i] = i;
        _lastEnd[levels - 1] = +1;
    }

    // A new replica's stream seed, the growth stream's next draw (its high 32 bits).
    private int NextGrowthSeed() => unchecked((int)(_growthRng.NextUInt64() >> 32));

    // Replica i's blocks whose cadence falls due on the scan being taken, in order, on its block stream.
    private void RunBlocks(int i)
    {
        if (_blocks is null) return;
        IReadOnlyList<IBlock<ReversibleJumpChain<TState>>> blocks = _blocks[i];
        long scan = Scans + 1;
        for (int j = 0; j < blocks.Count; j++)
        {
            if (scan % blocks[j].Every != 0) continue;
            _blockAttempts![i][j]++;
            if (blocks[j].Update(_chains[i], _blockRngs![i])) _blockAccepts![i][j]++;
        }
    }

    private static (IReadOnlyList<IBlock<ReversibleJumpChain<TState>>>[] Blocks, long[][] Attempts, long[][] Accepts) BuildBlocks(
        Func<int, IReadOnlyList<IBlock<ReversibleJumpChain<TState>>>> blocksForReplica, int levels)
    {
        var blocks = new IReadOnlyList<IBlock<ReversibleJumpChain<TState>>>[levels];
        var attempts = new long[levels][];
        var accepts = new long[levels][];
        for (int i = 0; i < levels; i++)
        {
            blocks[i] = BlockSet(blocksForReplica, i);
            attempts[i] = new long[blocks[i].Count];
            accepts[i] = new long[blocks[i].Count];
        }
        return (blocks, attempts, accepts);
    }

    private static IReadOnlyList<IBlock<ReversibleJumpChain<TState>>> BlockSet(
        Func<int, IReadOnlyList<IBlock<ReversibleJumpChain<TState>>>> blocksForReplica, int i)
    {
        IReadOnlyList<IBlock<ReversibleJumpChain<TState>>> set = blocksForReplica(i)
            ?? throw new ArgumentException($"Replica {i} has no block set.", nameof(blocksForReplica));
        foreach (IBlock<ReversibleJumpChain<TState>> block in set)
        {
            ArgumentNullException.ThrowIfNull(block, nameof(blocksForReplica));
            if (block.Every < 1) throw new ArgumentException($"Block '{block.Name}' has a cadence below 1.", nameof(blocksForReplica));
        }
        return set;
    }

    // Pair i's log swap ratio on the current states, with its rejection probability recorded.
    private double Record(int i)
    {
        ReversibleJumpChain<TState> a = _chains[i], b = _chains[i + 1];
        double logA = a.Temperature.LogDensity(b.CurrentEvaluation) + b.Temperature.LogDensity(a.CurrentEvaluation)
                    - a.CurrentLogDensity - b.CurrentLogDensity;
        double accept = logA >= 0.0 ? 1.0 : logA < 0.0 ? Math.Exp(logA) : 0.0;   // NaN: no swap
        _rejectionSum[i] += 1.0 - accept;
        _roundRejectionSum[i] += 1.0 - accept;
        return logA;
    }

    private void Attempt(int i, double logA)
    {
        SwapAttempts++;
        _pairAttempts[i]++;
        // Drawn on every attempt, so the swap stream does not advance by the decision: between states whose targets
        // agree to rounding, rounding decides a tie's sign, and it must not decide whether a uniform is consumed.
        double u = _swapRng.NextDouble();
        if (!(logA >= 0.0 || u < Math.Exp(logA))) return;

        SwapAccepts++;
        _pairAccepts[i]++;
        ReversibleJumpChain<TState> a = _chains[i], b = _chains[i + 1];
        TState xa = a.Current, xb = b.Current;
        ITargetEvaluation ea = a.CurrentEvaluation, eb = b.CurrentEvaluation;
        a.SetState(xb, eb);
        b.SetState(xa, ea);
        (_labelAt[i], _labelAt[i + 1]) = (_labelAt[i + 1], _labelAt[i]);
    }

    // Each rung's slope at its own b and its log ratio to the next colder rung, on the state it now holds.
    private void RecordRungs()
    {
        Span<double> values = stackalloc double[RungFunctional.Count];
        for (int i = 0; i < _chains.Length; i++)
        {
            ReversibleJumpChain<TState> rung = _chains[i];
            values[RungFunctional.Slope] = _tracesSlopes ? rung.Temperature.LogDensitySlope(rung.CurrentEvaluation) : double.NaN;
            values[RungFunctional.LogRatio] = i == 0
                ? double.NaN
                : _chains[i - 1].Temperature.LogDensity(rung.CurrentEvaluation) - rung.CurrentLogDensity;
            _rungTraces![i].Add(values);
        }
    }

    private static bool SuppliesSlope(ReversibleJumpChain<TState> chain)
    {
        try
        {
            _ = chain.Temperature.LogDensitySlope(chain.CurrentEvaluation);
            return true;
        }
        catch (NotSupportedException)
        {
            return false;
        }
    }

    // A state at the hot end that was last at the cold end completes a round trip; one at the cold end that was last
    // at the hot end is an annealed restart. A state that has not yet been to the hot end is at neither.
    private void TrackEnds()
    {
        int hot = _chains.Length - 1;
        int atHot = _labelAt[hot], atCold = _labelAt[0];
        if (_lastEnd[atHot] == -1) RoundTrips++;
        _lastEnd[atHot] = +1;
        if (_lastEnd[atCold] == +1)
        {
            AnnealedRestarts++;
            _lastEnd[atCold] = -1;
        }
    }

    private static void ValidateLadder(double[] inverseTemperatures, TemperingPath path, IReferenceSampler<TState>? reference)
    {
        ArgumentNullException.ThrowIfNull(inverseTemperatures);
        if (inverseTemperatures.Length < 2) throw new ArgumentException("Tempering needs at least two levels.", nameof(inverseTemperatures));
        if (inverseTemperatures[0] != 1.0) throw new ArgumentException("The cold rung must be exactly b = 1.", nameof(inverseTemperatures));
        for (int i = 1; i < inverseTemperatures.Length; i++)
            if (!(inverseTemperatures[i] < inverseTemperatures[i - 1]))
                throw new ArgumentException("The ladder must be strictly descending.", nameof(inverseTemperatures));
        double hot = inverseTemperatures[^1];
        if (reference is not null)
        {
            if (reference.Path != path) throw new ArgumentException("The reference samples another path's b = 0.", nameof(reference));
            if (hot != 0.0) throw new ArgumentException("A ladder with a reference rung ends at exactly b = 0.", nameof(inverseTemperatures));
            return;
        }
        if (!(path == TemperingPath.Powered ? hot > 0.0 : hot >= 0.0))
            throw new ArgumentException(path == TemperingPath.Powered
                ? "On the powered path every b must be positive, unless the hot rung draws from a reference."
                : "Every b must be non-negative.", nameof(inverseTemperatures));
    }
}
