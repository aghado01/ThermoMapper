using System;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Core.Solver;

/// <summary>
/// Per-site field: each node sorts its own incident couplings and cumulative-
/// sums them into a site-local energy ladder, so an edge's <c>Hcum</c> is
/// calibrated to its source node's neighborhood ranking rather than the global
/// pool. Directed — every CSR slot carries the value from its source row's
/// perspective, and the two directions of each undirected edge are reconciled by
/// the chosen <see cref="SymmetrizationRule"/> before clustering.
/// </summary>
/// <remarks>
/// O(N·K log K) vs MeanField's O(N·K log N·K). Density robustness is realized in
/// concert with adaptive-bandwidth graph construction upstream; per-site ranking
/// on a fixed-bandwidth graph reorders cuts but does not by itself rescue
/// small-coupling regions. Not robust to frustration (a topological, not
/// density, property — that needs the full Potts spin consistency SW keeps).
/// </remarks>
internal readonly struct LocalField : IField
{
    public static bool DirectedSymmetrize => true;

    public static double[] BuildHcum(CsrGraph graph)
    {
        int[] rowPtr = graph.RowPointers;
        double[] weights = graph.Weights;
        int n = graph.NodeCount;

        var hcum = new double[graph.Targets.Length];

        int maxDeg = 0;
        for (int i = 0; i < n; i++)
            maxDeg = Math.Max(maxDeg, rowPtr[i + 1] - rowPtr[i]);
        if (maxDeg == 0) return hcum;

        // Scratch reused across rows; only the [0, deg) prefix is touched per row.
        var localSlot = new int[maxDeg];
        var localCoupling = new double[maxDeg];

        for (int i = 0; i < n; i++)
        {
            int start = rowPtr[i];
            int deg = rowPtr[i + 1] - start;
            if (deg == 0) continue;

            for (int t = 0; t < deg; t++)
            {
                localSlot[t] = start + t;
                localCoupling[t] = weights[start + t];
            }

            // Ascending sort of this row's couplings, then a rank-inclusive
            // cumulative sum — the same ladder as MeanField but scoped to row i.
            Array.Sort(localCoupling, localSlot, 0, deg);
            double cum = 0.0;
            for (int t = 0; t < deg; t++)
            {
                cum += localCoupling[t];
                hcum[localSlot[t]] = cum;
            }
        }
        return hcum;
    }
}
