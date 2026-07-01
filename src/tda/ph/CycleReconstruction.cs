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
    /// Birth edge plus shortest-path edges in the 1-skeleton at the birth filtration value.
    /// </summary>
    public static int[] ReconstructH1Cycle(SimplicialFiltration filtration, int birthEdgeIndex)
    {
        ArgumentNullException.ThrowIfNull(filtration);

        Simplex birth = filtration.Simplices[birthEdgeIndex];
        if (birth.Dimension != 1 || birth.Vertices.Length != 2)
            throw new ArgumentException(
                $"Simplex at index {birthEdgeIndex} is not a 1-simplex edge.",
                nameof(birthEdgeIndex));

        double t = birth.FiltrationValue;
        int u = birth.Vertices[0];
        int v = birth.Vertices[1];

        var adjacency = new Dictionary<int, List<(int Neighbor, int EdgeIndex, double Weight)>>();
        for (int i = 0; i < filtration.Simplices.Count; i++)
        {
            if (i == birthEdgeIndex)
                continue;

            Simplex s = filtration.Simplices[i];
            if (s.Dimension != 1 || s.FiltrationValue > t)
                continue;

            int a = s.Vertices[0];
            int b = s.Vertices[1];
            double w = s.FiltrationValue;
            AddAdjacency(adjacency, a, b, i, w);
            AddAdjacency(adjacency, b, a, i, w);
        }

        List<int> pathEdges = ShortestPathEdgeIndices(adjacency, filtration, u, v);
        if (pathEdges.Count == 0 && u != v)
            throw new InvalidOperationException(
                $"No u→v path in 1-skeleton at filtration ≤ {t} (birth edge index {birthEdgeIndex}).");

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
                double nd = d + weight;
                if (nd >= dist.GetValueOrDefault(neighbor, double.PositiveInfinity))
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
