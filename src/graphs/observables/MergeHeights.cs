using System;
using System.Collections.Generic;
using Graphs.Primitives;

namespace Graphs.Observables;

/// <summary>
/// One join of two components as the coupling threshold falls: the weakest edge so far, in log space, and the sizes
/// of the components it merged.
/// </summary>
/// <param name="LogCoupling">ln J of the joining edge: the height of the merge.</param>
/// <param name="SizeA">Nodes in the first component before the join.</param>
/// <param name="SizeB">Nodes in the second component before the join.</param>
/// <param name="Source">The joining edge's lower endpoint.</param>
/// <param name="Target">The joining edge's higher endpoint.</param>
public readonly record struct LinkageMerge(double LogCoupling, int SizeA, int SizeB, int Source, int Target);

/// <summary>
/// The single-linkage merge heights of a coupling graph: Kruskal's pass over the edges strongest first, so each merge's
/// height is the bottleneck coupling at which two components first connect. Read in
/// <see cref="CsrGraph.LogWeights"/> when the graph has them, so a far gap keeps its height past the double range (a
/// Gaussian gap at d/a = 60 merges at −1800, not at a tie of zeros); else in the log of the weights, with zero weight no
/// coupling. An edge at −∞ joins nothing, so a graph held together only by such edges returns a forest.
/// </summary>
/// <remarks>
/// The observable of a graph, not of a temperature: it is the maximum-spanning-tree profile on a log axis, the
/// stage structure a coupling graph offers a thermal single-linkage reading, and the cumulative-field solver's cut for
/// an edge sits at T = Hcum/ln 2, where Hcum is the ladder's running sum and not the bare coupling.
/// </remarks>
public static class MergeHeights
{
    /// <summary>The merges in order of falling height (strongest first); ties by edge endpoints, so the order is deterministic.</summary>
    public static IReadOnlyList<LinkageMerge> Compute(CsrGraph graph)
    {
        int n = graph.NodeCount;
        var edges = new List<(double Log, int U, int V)>();
        for (int u = 0; u < n; u++)
            for (int s = graph.RowPointers[u]; s < graph.RowPointers[u + 1]; s++)
            {
                int v = graph.Targets[s];
                if (v <= u) continue;
                double log = graph.LogWeights is not null
                    ? graph.LogWeights[s]
                    : graph.Weights[s] > 0.0 ? Math.Log(graph.Weights[s]) : double.NegativeInfinity;
                if (double.IsNaN(log) || double.IsNegativeInfinity(log)) continue;
                edges.Add((log, u, v));
            }

        edges.Sort(static (a, b) =>
        {
            int byHeight = b.Log.CompareTo(a.Log);
            if (byHeight != 0) return byHeight;
            int byU = a.U.CompareTo(b.U);
            return byU != 0 ? byU : a.V.CompareTo(b.V);
        });

        var uf = new UnionFind(n);
        var merges = new List<LinkageMerge>(Math.Max(0, n - 1));
        foreach (var (log, u, v) in edges)
        {
            int ru = uf.Find(u), rv = uf.Find(v);
            if (ru == rv) continue;
            merges.Add(new LinkageMerge(log, uf.Size(ru), uf.Size(rv), u, v));
            uf.Union(u, v);
        }
        return merges;
    }
}
