using System;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Strategies;

/// <summary>
/// The shared affinity→partition close: connected components of the active-bond
/// subgraph (<c>g[e] &gt; theta</c>) over a per-edge affinity field, densified
/// into a <see cref="Assignment"/>. Every inference method's affinity — SW's
/// bond-frequency reduction, PKWang's closed-form field — is partitioned
/// identically here, decoupling the partition close from how the affinity was
/// produced. The single home of the threshold-and-connect step the per-T sweep
/// loop closes on.
/// </summary>
internal static class AffinityThreshold
{
    /// <summary>
    /// Union every undirected edge whose affinity <paramref name="g"/> exceeds
    /// <paramref name="theta"/> (walking the <c>j &gt; i</c> CSR half), then
    /// densify the components. <paramref name="g"/> is indexed by CSR slot,
    /// parallel to <see cref="CsrGraph.Targets"/>.
    /// </summary>
    /// <param name="peripheralCapture">
    /// When <see langword="true"/>, after the threshold pass each node is also
    /// unioned with its single highest-affinity neighbor regardless of θ
    /// (Domany1999 step 2 — captures cluster periphery whose density decreases
    /// toward the perimeter). O(E) extra pass; default off.
    /// </param>
    internal static Assignment Connect(CsrGraph graph, double[] g, double theta, bool peripheralCapture = false)
    {
        if (theta < 0.0 || theta > 1.0)
            throw new InvalidOperationException($"theta ({theta}) must lie in [0, 1].");

        int n = graph.NodeCount;
        var uf = new UnionFind(n);
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
            if (g[edge.Slot] > theta) uf.Union(edge.Source, edge.Target);

        if (peripheralCapture)
        {
            // Affinities.G is only meaningful at upper-triangle (j > i) CSR slots — lower-triangle
            // slots carry 0. Walk upper-triangle edges and let each edge update both endpoints'
            // max-G tracker so that lower-index nodes see their edge values correctly.
            var bestNeighbor = new int[n];
            var bestGPerNode = new double[n];
            Array.Fill(bestNeighbor, -1);

            foreach (UndirectedEdge edge in graph.UndirectedEdges())
            {
                double gij = g[edge.Slot];
                if (gij > bestGPerNode[edge.Source]) { bestGPerNode[edge.Source] = gij; bestNeighbor[edge.Source] = edge.Target; }
                if (gij > bestGPerNode[edge.Target]) { bestGPerNode[edge.Target] = gij; bestNeighbor[edge.Target] = edge.Source; }
            }

            for (int i = 0; i < n; i++)
                if (bestNeighbor[i] >= 0)
                    uf.Union(i, bestNeighbor[i]);
        }

        return UnionFindLabeler.Densify(uf, n);
    }
}
