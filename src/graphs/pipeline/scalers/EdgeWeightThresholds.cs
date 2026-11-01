using Graphs.Pipeline;
using Graphs.Primitives;

namespace Graphs.Pipeline.Scalers;

/// <summary>
/// The near-zero level the delta-collapse detectors count against.
/// </summary>
internal static class EdgeWeightThresholds
{
    /// <summary>A coupling under this counts as collapsed in the delta-collapse detectors.</summary>
    public const double NearZero = 1e-8;

    /// <summary>
    /// Fraction of the edges the kernel produced whose coupling is under
    /// <see cref="NearZero"/>. The <paramref name="cutEdges"/> (undirected) the floor already
    /// dropped count as collapsed, so the ratio reads the edge set before the cut. The graph
    /// stores each edge in two CSR slots, so the cut edges count twice, as the kept ones do.
    /// Edges in <paramref name="excluded"/> (the repair bridges, weak by construction) are left
    /// out of both counts.
    /// </summary>
    public static double NearZeroRatio(CsrGraph graph, int cutEdges, RepairEdgeSet? excluded = null)
    {
        long collapsed = 2L * cutEdges;
        long total = collapsed;
        for (int u = 0; u < graph.NodeCount; u++)
            for (int slot = graph.RowPointers[u]; slot < graph.RowPointers[u + 1]; slot++)
            {
                if (excluded is not null && excluded.Contains(u, graph.Targets[slot])) continue;
                total++;
                if (graph.Weights[slot] < NearZero) collapsed++;
            }

        return total == 0 ? 0.0 : (double)collapsed / total;
    }
}
