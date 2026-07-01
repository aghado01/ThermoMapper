namespace Clustering.Dendrograms;

/// <summary>
/// One merge event in an agglomerative single-linkage-style dendrogram.
/// Leaves are points <c>0..N-1</c>; internal merge nodes are assigned
/// dense ids <c>N..2N-2</c> in build order. <see cref="LeftChild"/> and
/// <see cref="RightChild"/> may reference either a leaf or an internal id.
/// </summary>
/// <remarks>
/// <para><b>Cost axis.</b> <see cref="Distance"/> is the scalar at which
/// the two subtrees merge. The build pass is expected to emit merges in
/// monotone non-decreasing <see cref="Distance"/> order — true for both
/// Kruskal single-linkage (HDBSCAN's mutual-reachability MST) and greedy
/// entropy merging (GMM agglomerative). The
/// <see cref="Dendrogram.CostAxis"/> field on the enclosing wrapper
/// documents the physical interpretation (distance, ΔH, dissimilarity).</para>
///
/// <para><b>λ convention.</b> HDBSCAN's excess-of-mass pass consumes
/// <c>λ = 1/Distance</c>; that conversion lives in the consumer (e.g.
/// the condensation pass), not in this DTO — keeping the record neutral
/// across cost-axis interpretations.</para>
/// </remarks>
public readonly record struct DendrogramNode(
    int    LeftChild,
    int    RightChild,
    double Distance,
    int    Size);
