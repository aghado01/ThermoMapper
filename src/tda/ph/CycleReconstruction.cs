#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Reconstructs full H1 cycle representatives for infinite bars (birth edge + shortest u→v path).
/// </summary>
public static class CycleReconstruction
{
    /// <summary>
    /// Birth edge plus shortest-path edges among the edges before it in the filtration. An edge tied with it and after it
    /// stays out: a triangle tied with both can make the two a boundary, which is no cycle of the bar. Any path among the
    /// earlier edges gives a cycle of the bar's class, so "shortest" only picks a small one: edge values read as lengths when
    /// they are all non-negative, as a Rips filtration's distances are. Otherwise a value is not a length and an undirected
    /// negative edge has no shortest path, so a path's length is its number of edges.
    /// </summary>
    public static int[] ReconstructH1Cycle(SimplicialFiltration filtration, int birthEdgeIndex)
    {
        ArgumentNullException.ThrowIfNull(filtration);

        Simplex birth = filtration.Simplices[birthEdgeIndex];
        if (birth.Dimension != 1 || birth.Vertices.Length != 2)
            throw new ArgumentException(
                $"Simplex at index {birthEdgeIndex} is not a 1-simplex edge.",
                nameof(birthEdgeIndex));

        int u = birth.Vertices[0];
        int v = birth.Vertices[1];

        bool valuesAreLengths = true;
        for (int i = 0; i < birthEdgeIndex; i++)
        {
            Simplex s = filtration.Simplices[i];
            if (s.Dimension == 1)
                valuesAreLengths &= s.FiltrationValue >= 0.0;
        }

        var adjacency = new Dictionary<int, List<(int Neighbor, int EdgeIndex, double Weight)>>();
        for (int i = 0; i < birthEdgeIndex; i++)
        {
            Simplex s = filtration.Simplices[i];
            if (s.Dimension != 1)
                continue;

            int a = s.Vertices[0];
            int b = s.Vertices[1];
            double w = valuesAreLengths ? s.FiltrationValue : 1.0;
            AddAdjacency(adjacency, a, b, i, w);
            AddAdjacency(adjacency, b, a, i, w);
        }

        List<int> pathEdges = ShortestPathEdgeIndices(adjacency, filtration, u, v);
        if (pathEdges.Count == 0 && u != v)
            throw new InvalidOperationException(
                $"No u→v path among the edges before birth edge index {birthEdgeIndex}, so it closes no cycle.");

        var cycle = new SortedSet<int> { birthEdgeIndex };
        foreach (int edgeIdx in pathEdges)
            cycle.Add(edgeIdx);

        var result = new int[cycle.Count];
        cycle.CopyTo(result);
        return result;
    }

    static void AddAdjacency(
        Dictionary<int, List<(int Neighbor, int EdgeIndex, double Weight)>> adjacency,
        int from,
        int to,
        int edgeIndex,
        double weight)
    {
        if (!adjacency.TryGetValue(from, out List<(int, int, double)>? list))
        {
            list = new List<(int, int, double)>();
            adjacency[from] = list;
        }

        list.Add((to, edgeIndex, weight));
    }

    static List<int> ShortestPathEdgeIndices(
        Dictionary<int, List<(int Neighbor, int EdgeIndex, double Weight)>> adjacency,
        SimplicialFiltration filtration,
        int start,
        int end)
    {
        if (start == end)
            return new List<int>();

        var dist = new Dictionary<int, double> { [start] = 0.0 };
        var parentEdge = new Dictionary<int, int>();
        var pq = new PriorityQueue<int, double>();
        pq.Enqueue(start, 0.0);

        while (pq.Count > 0)
        {
            pq.TryDequeue(out int cur, out double d);
            if (d > dist.GetValueOrDefault(cur, double.PositiveInfinity))
                continue;

            if (cur == end)
                break;

            if (!adjacency.TryGetValue(cur, out List<(int Neighbor, int EdgeIndex, double Weight)>? neighbors))
                continue;

            foreach (var (neighbor, edgeIndex, weight) in neighbors)
            {
                // Huge lengths can sum to infinity, which still reaches a vertex not reached before.
                double nd = d + weight;
                if (dist.TryGetValue(neighbor, out double known) && nd >= known)
                    continue;

                dist[neighbor] = nd;
                parentEdge[neighbor] = edgeIndex;
                pq.Enqueue(neighbor, nd);
            }
        }

        if (!parentEdge.ContainsKey(end))
            return new List<int>();

        var edges = new List<int>();
        int at = end;
        while (at != start)
        {
            int edgeIdx = parentEdge[at];
            edges.Add(edgeIdx);
            at = OtherVertex(filtration.Simplices[edgeIdx], at);
        }

        edges.Reverse();
        return edges;
    }

    static int OtherVertex(Simplex edge, int vertex)
    {
        int a = edge.Vertices[0];
        int b = edge.Vertices[1];
        return a == vertex ? b : a;
    }
}
