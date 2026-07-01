using System;
using Graphs.Primitives;

namespace Graphs.Primitives.Mst;

/// <summary>
/// Kruskal's MST algorithm over an explicit edge list.
/// </summary>
public static class Kruskal
{
    /// <summary>
    /// Builds a minimum spanning tree from an edge list sorted by ascending
    /// weight. Returns the number of tree edges written to <paramref name="output"/>.
    /// </summary>
    /// <param name="sortedEdges">Input edges sorted by non-decreasing weight.</param>
    /// <param name="nodeCount">Number of nodes in the graph.</param>
    /// <param name="output">Caller-provided destination for the selected MST edges.
    /// Must have capacity at least <c>nodeCount - 1</c>.</param>
    /// <returns>The number of edges written to <paramref name="output"/>.</returns>
    public static int BuildMinimumSpanningTree(
        ReadOnlySpan<MstEdge> sortedEdges,
        int                   nodeCount,
        Span<MstEdge>         output)
    {
        if (nodeCount < 1)
            throw new ArgumentOutOfRangeException(nameof(nodeCount), "nodeCount must be positive.");
        if (output.Length < Math.Max(0, nodeCount - 1))
            throw new ArgumentException("Output span must have capacity at least nodeCount - 1.", nameof(output));

        var uf = new UnionFind(nodeCount);
        int written = 0;

        for (int i = 0; i < sortedEdges.Length && written < nodeCount - 1; i++)
        {
            var edge = sortedEdges[i];
            int ra = uf.Find(edge.U);
            int rb = uf.Find(edge.V);
            if (ra == rb)
                continue;

            uf.Union(ra, rb);
            output[written++] = edge;
        }

        return written;
    }
}
