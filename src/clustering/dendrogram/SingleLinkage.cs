using System;
using Graphs.Primitives;
using Graphs.Primitives.Mst;

namespace Clustering.Dendrograms;

/// <summary>
/// Single-linkage clustering of a <i>coupling</i> graph — CSR weights are
/// couplings <c>J</c> (higher = closer), the convention SPC and the proximity
/// builders use. The dual of the distance-native builders: weights are negated
/// so ascending MST weight = descending coupling, i.e. the strongest couplings
/// merge first.
/// </summary>
/// <remarks>
/// Thin composition over the canonical pipeline (<see cref="Kruskal"/> →
/// <see cref="DendrogramBuilder.BuildSingleLinkageDendrogram"/>); it exists so
/// callers — and the PKWang sampler's validation harness — have one reusable
/// entry point rather than re-deriving the MST→dendrogram chain. This is the
/// ground-truth method PKWang's MeanField rung reduces to (Lemma B).
/// </remarks>
public static class SingleLinkage
{
    /// <summary>
    /// Build the single-linkage dendrogram of a connected coupling graph. Cut it
    /// with <see cref="Dendrogram.CutToK"/> / <see cref="Dendrogram.CutAt"/>.
    /// </summary>
    /// <exception cref="InvalidOperationException">The graph is disconnected; a
    /// single-linkage dendrogram requires a spanning tree.</exception>
    public static Dendrogram FromCouplingGraph(CsrGraph graph)
    {
        if (graph.NodeCount < 2)
            throw new ArgumentOutOfRangeException(nameof(graph), "Single-linkage needs at least 2 nodes.");

        int[] rowPtr = graph.RowPointers;
        int[] targets = graph.Targets;
        double[] weights = graph.Weights;
        int n = graph.NodeCount;

        int edgeCount = 0;
        for (int i = 0; i < n; i++)
            for (int e = rowPtr[i]; e < rowPtr[i + 1]; e++)
                if (targets[e] > i) edgeCount++;

        var edges = new MstEdge[edgeCount];
        int k = 0;
        for (int i = 0; i < n; i++)
            for (int e = rowPtr[i]; e < rowPtr[i + 1]; e++)
                if (targets[e] > i)
                    edges[k++] = new MstEdge(i, targets[e], -weights[e]); // negate: strongest first

        Array.Sort(edges); // ascending weight = descending coupling

        var mst = new MstEdge[n - 1];
        int treeEdges = Kruskal.BuildMinimumSpanningTree(edges, n, mst);
        if (treeEdges != n - 1)
            throw new InvalidOperationException(
                $"Coupling graph is disconnected ({n - treeEdges} components); " +
                "single-linkage requires a connected graph.");

        var uf = new UnionFind(2 * n - 1);
        DendrogramNode[] merges = DendrogramBuilder.BuildSingleLinkageDendrogram(mst.AsSpan(0, treeEdges), n, uf);
        return new Dendrogram(merges, n, CostAxis: "neg_coupling");
    }
}
