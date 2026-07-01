namespace Clustering.Graphical.HdbScan;

/// <summary>
/// Condensed-tree cluster-selection rule — the shared selector axis over the
/// merge tree (the downstream half of the resolution spine). Maps to
/// <see cref="Clustering.Dendrograms.CondensedTree"/>'s two selectors.
/// </summary>
public enum ClusterSelectionMethod
{
    /// <summary>Excess-of-mass: few large, persistent clusters
    /// (<see cref="Clustering.Dendrograms.CondensedTree.SelectByExcessOfMass"/>).
    /// sklearn's <c>'eom'</c> — the HDBSCAN default.</summary>
    Eom,

    /// <summary>Leaf: the finest stable clusters (condensed-tree leaves;
    /// <see cref="Clustering.Dendrograms.CondensedTree.SelectByLeaf"/>). sklearn's
    /// <c>'leaf'</c> — recovers structure EOM under-segments (multi-class data).</summary>
    Leaf,
}

/// <summary>
/// Mutual-reachability MST construction — the intrinsic geometry substrate (the
/// metric/density model, NOT preprocessing). Both build the SAME single-linkage
/// tree the condensation spine consumes; they differ only in which candidate
/// edges the MST is drawn from.
/// </summary>
public enum MstAlgorithm
{
    /// <summary>Exact dense MST over all O(n²) pairs (implicit Prim). The faithful
    /// default — every pair is a candidate, so the MST is the true minimum.</summary>
    Dense,

    /// <summary>Approximate MST over a k-nearest-neighbour graph: exact core
    /// distances + mutual-reachability edges restricted to the kNN adjacency
    /// (Kruskal), with <see cref="Graphs.Primitives.Mst.Boruvka"/> bridging any
    /// disconnected components into one spanning tree. Reuses the shared kNN /
    /// MST primitives; the substrate intrinsic geodesics would also ride on
    /// (see the intrinsic-reachability note). Approximate — the kNN restriction
    /// can miss a globally-cheaper edge; raise <see cref="HdbscanSettings.GraphNeighbors"/>
    /// to tighten toward the dense result (k = n−1 reproduces it exactly).</summary>
    SparseKnn,
}

/// <summary>
/// All knobs that define an HDBSCAN run. Pure declarative data — no fluent
/// surface, mirroring the strict-core / fluent-shell split used by the graph
/// engine (the fluent/CLI layer builds one of these and hands it to
/// <see cref="HdbscanSession"/>.Run). Algorithm <i>intrinsics</i> only — feature
/// scaling / DR are model-agnostic Transforms the shell applies upstream, never a
/// field here (see <c>.discussion/issues/hdbscan/preprocessing-intrinsics-boundary.md</c>).
/// </summary>
public sealed record HdbscanSettings
{
    /// <summary>Core-distance neighbour count (≥ 2). The session clamps the
    /// effective value to <c>[2, n-1]</c> so small inputs still resolve a kth
    /// neighbour.</summary>
    public int MinPts { get; init; } = 5;

    /// <summary>Minimum subtree size to be treated as a real cluster during
    /// condensation. <c>null</c> ⇒ defaults to <see cref="MinPts"/> (the runner's
    /// convention).</summary>
    public int? MinClusterSize { get; init; }

    /// <summary>When true, the root cluster can be selected by excess-of-mass —
    /// the mapper-friendly "one dense blob + outliers ⇒ one cluster" behaviour.
    /// When false, datasets with no real splits return all-noise (sklearn
    /// default).</summary>
    public bool AllowSingleCluster { get; init; } = true;

    /// <summary>Condensed-tree selection rule: <see cref="ClusterSelectionMethod.Eom"/>
    /// (default; few large persistent clusters) or
    /// <see cref="ClusterSelectionMethod.Leaf"/> (finest stable clusters — recovers
    /// the structure EOM under-segments on multi-class data).</summary>
    public ClusterSelectionMethod ClusterSelectionMethod { get; init; } = ClusterSelectionMethod.Eom;

    /// <summary>Mutual-reachability distance floor (sklearn's
    /// <c>cluster_selection_epsilon</c>): clusters split off below this distance
    /// are merged upward into the nearest coarser cluster — a DBSCAN-like floor
    /// that tames over-segmentation (notably leaf's). <c>0</c> (default) ⇒ off;
    /// composes with both <see cref="ClusterSelectionMethod"/> values.</summary>
    public double ClusterSelectionEpsilon { get; init; }

    /// <summary>Mutual-reachability MST construction: <see cref="MstAlgorithm.Dense"/>
    /// (default; exact, O(n²)) or <see cref="MstAlgorithm.SparseKnn"/> (kNN-graph
    /// approximate MST — cheaper edge set, the geometry substrate).</summary>
    public MstAlgorithm MstAlgorithm { get; init; } = MstAlgorithm.Dense;

    /// <summary>kNN neighbour count for <see cref="MstAlgorithm.SparseKnn"/>'s
    /// candidate-edge graph (ignored when dense). Clamped to <c>[MinPts, n-1]</c>.
    /// <c>null</c> ⇒ <c>max(MinPts, 10)</c>. Larger ⇒ closer to the dense MST
    /// (n−1 reproduces it exactly), at more edges.</summary>
    public int? GraphNeighbors { get; init; }

    /// <summary>Distance-metric spec:
    /// <c>euclidean</c> | <c>manhattan</c>/<c>l1</c> | <c>minkowski:p=N</c> |
    /// <c>hamming</c> | <c>poincare</c> | <c>cosine</c>. The Minkowski exponent
    /// folds into the spec string (no separate field) — same canonical
    /// <see cref="Graphs.Distance.MetricRegistry"/> vocabulary SPC's
    /// <c>DistanceMetricFactory</c> uses, which is what kills the A3 drift.</summary>
    public string Metric { get; init; } = "euclidean";
}
