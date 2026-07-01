using System;
using Graphs.Neighbors;
using Graphs.Primitives;

namespace Graphs.Pipeline.Filters;

/// <summary>
/// Stage 2 — the OR-rule symmetrizer. Edge (i, j) exists if either
/// <c>i</c> is in <c>j</c>'s top-K *or* <c>j</c> is in <c>i</c>'s top-K.
/// </summary>
public sealed class PassThroughFilter : ITopologyFilter
{
    public NeighborSelection Filter(NeighborSelection directed, int n, Func<int, int, double> pairDistance)
    {
        if (IsAlreadySymmetric(directed.AllNeighbors))
            return directed;

        return Symmetrization.OrUnion(directed, n);
    }

    /// <summary>
    /// Fast symmetry heuristic: samples only the first edge from the first
    /// non-empty row and checks for its reverse. Full symmetry is not verified.
    /// </summary>
    private static bool IsAlreadySymmetric(Neighbor[][] neighbors)
    {
        for (int i = 0; i < neighbors.Length; i++)
        {
            if (neighbors[i].Length == 0) continue;
            int j = neighbors[i][0].Index;
            foreach (var nb in neighbors[j])
            {
                if (nb.Index == i) return true;
            }
            return false;
        }
        return true;
    }
}
