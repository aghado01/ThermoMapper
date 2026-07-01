// src/clustering/hdbscan/HdbscanRunner.cs
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using Clustering.Dendrograms;
using Graphs.Distance;
using Graphs.Primitives;
using Graphs.Primitives.Mst;

namespace Clustering.Graphical.HdbScan;

/// <summary>
/// Full HDBSCAN pipeline. Pre-allocates all scratch at construction time and
/// reuses it across <see cref="Run{TMetric}"/> calls on data of the same N.
///
/// Pipeline phases:
///   1. Core-distance computation  — k-NN distance per point (minPts neighbour).
///   2. Implicit Prim's MST        — <see cref="Prim.ComputeMutualReachabilityMst{TMetric}"/>
///      over the mutual reachability graph; zero edge-list allocation.
///   3. Edge materialisation       — N−1 <see cref="MstEdge"/> structs, sorted
///      ascending by weight (Kruskal order).
///   4. Dendrogram construction    — single-linkage merge tree via UnionFind;
///      emits <see cref="DendrogramNode"/>[N−1].
///   5. Condensation + extraction  — λ-stability / excess-of-mass pass;
///      produces <see cref="HdbscanResult"/>.
/// </summary>
public sealed class HdbscanRunner
{
    private readonly int      _n;
    private readonly UnionFind _uf;         // capacity 2N−1; Reset() between runs
    private readonly double[]  _coreDist;
    private readonly bool[]    _visited;
    private readonly double[]  _minWeight;
    private readonly int[]     _parent;
    private readonly MstEdge[] _mstEdges;

    public HdbscanRunner(int n)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), "Must be >= 2.");
        _n = n;
        _uf = new UnionFind(2 * n - 1);
        _coreDist = new double[n];
        _visited = new bool[n];
        _minWeight = new double[n];
        _parent = new int[n];
        _mstEdges = new MstEdge[n - 1];
    }

    /// <param name="data">Flat row-major buffer, length N × dim.</param>
    /// <param name="dim">Dimensionality of each row in <paramref name="data"/>.</param>
    /// <param name="minPts">Smoothing parameter for core-distance computation (≥ 2).</param>
    /// <param name="metric">Distance metric; struct-generic so the JIT inlines the call.</param>
    /// <param name="minClusterSize">Minimum size for a subtree to be treated as
    /// a "real" cluster during condensation. Smaller subtrees fall out of their
    /// parent. Defaults to <paramref name="minPts"/> when null. Controls cluster
    /// granularity, not the cluster count: larger values → fewer/bigger clusters.</param>
    /// <param name="allowSingleCluster">If true, the root cluster can be selected
    /// by EOM — useful when the input is one dense blob with outliers (mapper-style
    /// cover patches). If false, datasets with no real splits return all-noise
    /// (sklearn default behaviour).</param>
    /// <param name="selectionMethod">Condensed-tree selection rule: excess-of-mass
    /// (default) or leaf (finest stable clusters; recovers structure EOM
    /// under-segments).</param>
    /// <param name="clusterSelectionEpsilon">Mutual-reachability distance floor;
    /// clusters born below it are merged upward (0 = off). See
    /// <see cref="HdbscanSettings.ClusterSelectionEpsilon"/>.</param>
    public HdbscanResult Run<TMetric>(
        ReadOnlySpan<double>    data,
        int                     dim,
        int                     minPts,
        TMetric                 metric,
        int?                    minClusterSize          = null,
        bool                    allowSingleCluster      = true,
        ClusterSelectionMethod  selectionMethod         = ClusterSelectionMethod.Eom,
        double                  clusterSelectionEpsilon = 0.0)
        where TMetric : struct, IDistanceMetric
    {
        if (minPts < 2)
            throw new ArgumentOutOfRangeException(nameof(minPts), "Must be >= 2.");

        int effMinClusterSize = minClusterSize ?? minPts;
        if (effMinClusterSize < 2)
            throw new ArgumentOutOfRangeException(nameof(minClusterSize), "Must be >= 2.");

        int n = _n;

        // ── Phase 1: core distances ───────────────────────────────────────────
        CoreDistances.Compute(data, n, dim, minPts, metric, _coreDist.AsSpan());

        // ── Phase 2: implicit Prim's MST + edge materialisation ─────────────
        for (int i = 0; i < n; i++)
        {
            _visited[i] = false;
            _minWeight[i] = double.PositiveInfinity;
            _parent[i] = -1;
        }

        Prim.ComputeMutualReachabilityMst(
            data, n, dim,
            _coreDist.AsSpan(),
            _visited.AsSpan(),
            _minWeight.AsSpan(),
            _parent.AsSpan(),
            metric);

        for (int v = 1; v < n; v++)
            _mstEdges[v - 1] = new MstEdge(_parent[v], v, _minWeight[v]);

        Array.Sort(_mstEdges, 0, n - 1);

        // ── Phase 4: build dendrogram ─────────────────────────────────────────
        DendrogramNode[] tree = DendrogramBuilder.BuildSingleLinkageDendrogram(
            _mstEdges.AsSpan(0, n - 1), n, _uf);

        // Wrap the raw merge sequence in the shared Dendrogram DTO so the
        // persistence without reproducing the build pass. CostAxis names
        // the y-axis units (mutual-reachability distance, the same value
        // HDBSCAN's condensation pass inverts to λ = 1/d).
        var dendrogram = new Dendrogram(
            Merges:    tree,
            LeafCount: n,
            CostAxis:  "mutual_reachability_distance");

        // ── Phase 5: condense + extract clusters ─────────────────────────────
        return ExtractClusters(dendrogram, effMinClusterSize, allowSingleCluster, selectionMethod, clusterSelectionEpsilon);
    }

    /// <summary>
    /// Approximate HDBSCAN over a k-nearest-neighbour graph — the
    /// <see cref="MstAlgorithm.SparseKnn"/> path. Core distances stay <b>exact</b>
    /// (the minPts-th neighbour from an exact kNN); only the MST edge set is
    /// approximate: mutual-reachability edges are restricted to the kNN adjacency
    /// (<see cref="Kruskal"/>), and <see cref="Boruvka.AddMinimalBridges"/> bridges
    /// any leftover components into one spanning tree. Phases 4–5 (dendrogram,
    /// condensation, selection) are shared with the dense path unchanged.
    /// <para>Takes <c>double[]</c> (not a span) so the per-pair distance closure
    /// the kNN / bridge passes need can capture it. With <paramref name="graphK"/>
    /// = n−1 every pair is a candidate ⇒ the result equals the dense MST exactly.</para>
    /// </summary>
    /// <param name="graphK">kNN neighbour count for the candidate-edge graph;
    /// must be in <c>[minPts, n-1]</c>.</param>
    public HdbscanResult RunSparse<TMetric>(
        double[]                data,
        int                     dim,
        int                     minPts,
        TMetric                 metric,
        int                     graphK,
        int?                    minClusterSize          = null,
        bool                    allowSingleCluster      = true,
        ClusterSelectionMethod  selectionMethod         = ClusterSelectionMethod.Eom,
        double                  clusterSelectionEpsilon = 0.0)
        where TMetric : struct, IDistanceMetric
    {
        ArgumentNullException.ThrowIfNull(data);
        if (minPts < 2)
            throw new ArgumentOutOfRangeException(nameof(minPts), "Must be >= 2.");
        int effMinClusterSize = minClusterSize ?? minPts;
        if (effMinClusterSize < 2)
            throw new ArgumentOutOfRangeException(nameof(minClusterSize), "Must be >= 2.");

        int n = _n;
        if (graphK < minPts)
            throw new ArgumentOutOfRangeException(nameof(graphK), "graphK must be >= minPts.");
        if (graphK > n - 1) graphK = n - 1;

        // ── Phases 1+2: one exact kNN pass — per-point ascending (dist,index)
        // buffer of size graphK. Core distance = the minPts-th entry (exact, same
        // as the dense CoreDistances). Neighbour rows stored flat for the edge
        // pass below (mutual-reachability needs both endpoints' core, so edges
        // can only be weighted after every core is known). O(n²·dim), like dense.
        double[] core = _coreDist;
        var nbDst = new double[n * graphK];
        var nbIdx = new int[n * graphK];
        var nbLen = new int[n];
        for (int i = 0; i < n; i++)
        {
            ReadOnlySpan<double> rowI = data.AsSpan(i * dim, dim);
            int @base = i * graphK;
            for (int p = 0; p < graphK; p++) nbDst[@base + p] = double.PositiveInfinity;
            int filled = 0;
            double worst = double.PositiveInfinity;
            for (int j = 0; j < n; j++)
            {
                if (j == i) continue;
                double d = metric.Distance(rowI, data.AsSpan(j * dim, dim));
                if (filled == graphK && d >= worst) continue;
                int pos = filled < graphK ? filled : graphK - 1;
                nbDst[@base + pos] = d;
                nbIdx[@base + pos] = j;
                while (pos > 0 && nbDst[@base + pos] < nbDst[@base + pos - 1])
                {
                    (nbDst[@base + pos], nbDst[@base + pos - 1]) = (nbDst[@base + pos - 1], nbDst[@base + pos]);
                    (nbIdx[@base + pos], nbIdx[@base + pos - 1]) = (nbIdx[@base + pos - 1], nbIdx[@base + pos]);
                    pos--;
                }
                if (filled < graphK) filled++;
                worst = nbDst[@base + filled - 1];
            }
            nbLen[i] = filled;
            core[i]  = nbDst[@base + minPts - 1];   // minPts-th NN = core distance
        }

        // ── Phase 3: mutual-reachability edges over the kNN adjacency (dedup). ──
        var seen  = new HashSet<long>(n * graphK);
        var edges = new List<MstEdge>(n * graphK);
        for (int i = 0; i < n; i++)
        {
            int @base = i * graphK;
            for (int t = 0; t < nbLen[i]; t++)
            {
                int j  = nbIdx[@base + t];
                int lo = Math.Min(i, j), hi = Math.Max(i, j);
                if (!seen.Add(((long)lo << 32) | (uint)hi)) continue;
                double w = Math.Max(core[i], Math.Max(core[j], nbDst[@base + t]));
                edges.Add(new MstEdge(lo, hi, w));
            }
        }
        edges.Sort();   // ascending by weight (MstEdge : IComparable)

        // ── MST: Kruskal over the sparse candidates, bridge leftover components. ──
        int written = Kruskal.BuildMinimumSpanningTree(
            CollectionsMarshal.AsSpan(edges), n, _mstEdges.AsSpan(0, n - 1));

        if (written < n - 1)
        {
            var components = new UnionFind(n);
            for (int e = 0; e < written; e++)
                components.Union(_mstEdges[e].U, _mstEdges[e].V);

            List<Boruvka.BridgeEdge> bridges = Boruvka.AddMinimalBridges(n,
                (i, j) =>
                {
                    double d = metric.Distance(data.AsSpan(i * dim, dim), data.AsSpan(j * dim, dim));
                    return Math.Max(core[i], Math.Max(core[j], d));
                },
                components);

            for (int b = 0; b < bridges.Count; b++)
                _mstEdges[written++] = new MstEdge(bridges[b].LoIndex, bridges[b].HiIndex, bridges[b].Weight);
        }

        if (written != n - 1)
            throw new InvalidOperationException(
                $"Sparse MST produced {written} edges, expected {n - 1} — graph not connectable.");

        Array.Sort(_mstEdges, 0, n - 1);

        DendrogramNode[] tree = DendrogramBuilder.BuildSingleLinkageDendrogram(
            _mstEdges.AsSpan(0, n - 1), n, _uf);
        var dendrogram = new Dendrogram(
            Merges:    tree,
            LeafCount: n,
            CostAxis:  "mutual_reachability_distance");

        return ExtractClusters(dendrogram, effMinClusterSize, allowSingleCluster, selectionMethod, clusterSelectionEpsilon);
    }

    // ── Phase 5 ───────────────────────────────────────────────────────────────

    /// <summary>
    /// HDBSCAN cluster extraction, delegated to the shared resolution spine:
    /// <see cref="Condensation.Condense"/> collapses the single-linkage tree into
    /// a <see cref="CondensedTree"/> of "real" clusters (both sides clear
    /// <paramref name="minClusterSize"/> at a split), then the selector picks
    /// clusters off it — excess-of-mass (few large persistent clusters) or leaf
    /// (the finest stable clusters). <see cref="CondensedTree.ResolveLabeled"/>
    /// walks each leaf to its nearest selected ancestor and emits the λ-ratio
    /// membership probability.
    /// <para>The condensation, EOM, and λ-ratio membership formerly inlined here
    /// now live once on the spine (<c>Clustering.Dendrograms</c>) so SPC / GMM
    /// share the same selector axis; the runner is a consumer of it.</para>
    /// </summary>
    private static HdbscanResult ExtractClusters(
        Dendrogram             dendrogram,
        int                    minClusterSize,
        bool                   allowSingleCluster,
        ClusterSelectionMethod selectionMethod,
        double                 clusterSelectionEpsilon)
    {
        CondensedTree condensed = Condensation.Condense(dendrogram, minClusterSize);

        bool[] selected = selectionMethod == ClusterSelectionMethod.Leaf
            ? condensed.SelectByLeaf(allowSingleCluster)
            : condensed.SelectByExcessOfMass(allowSingleCluster);

        if (clusterSelectionEpsilon > 0.0)
            selected = condensed.SelectByEpsilon(selected, clusterSelectionEpsilon, allowSingleCluster);

        var (labels, memberProb, clusterCount) = condensed.ResolveLabeled(selected);
        return new HdbscanResult(labels, memberProb, clusterCount, dendrogram);
    }
}
/// <summary>
/// Output of a completed <see cref="HdbscanRunner.Run{TMetric}"/> call.
///
/// <see cref="Labels"/>: cluster index per point in [0, ClusterCount), or
///   −1 for noise points.
/// <see cref="MembershipProbabilities"/>: soft assignment score ∈ [0, 1] for
///   each point. For a point x labelled to cluster L:
///   λ_fallout(x) / λ_death(L) when x's deepest condensed cluster is L itself
///   (the point may have left L early via a small-side falls-out); 1.0 when L
///   is a strict ancestor of x's deepest cluster (x stayed in L until L died).
///   0.0 for noise points.
/// </summary>
public sealed class HdbscanResult(
    int[]      labels,
    double[]   membershipProbabilities,
    int        clusterCount,
    Dendrogram dendrogram)
{
    public int[]      Labels                   { get; } = labels;
    public double[]   MembershipProbabilities  { get; } = membershipProbabilities;
    public int        ClusterCount             { get; } = clusterCount;

    /// <summary>
    /// Raw single-linkage dendrogram produced in Phase 4. Cost axis is
    /// mutual-reachability distance; λ = 1/cost is the persistence
    /// scalar HDBSCAN's condensation pass consumes. Preserved on the
    /// result so downstream plotting / re-analysis (and the
    /// <c>userrepl hdbscan</c> persistence layer) can render the merge
    /// tree without re-running the pipeline.
    /// </summary>
    public Dendrogram Dendrogram               { get; } = dendrogram;

    public bool IsNoise(int pointIndex) => Labels[pointIndex] < 0;
}

// ── Internal helpers ──────────────────────────────────────────────────────────
//
// The single-linkage dendrogram record (<c>DendrogramNode</c>) now lives in
// <see cref="Clustering.Dendrograms"/> so other agglomerative algorithms (GMM
// entropy-merge, future SPC threshold-sweep hierarchies) can emit the same
// shape. The HDBSCAN convention still applies: nodes 0..N-1 are leaves; nodes
// N..2N-2 are internal merge nodes assigned by the UnionFind during Kruskal's
// pass; Distance is the mutual-reachability weight, and λ = 1/Distance is the
// persistence value the excess-of-mass selection pass consumes.
