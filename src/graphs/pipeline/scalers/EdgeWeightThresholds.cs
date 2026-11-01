using Graphs.Primitives;

namespace Graphs.Pipeline.Scalers;

/// <summary>
/// The two coupling thresholds the scalers and the collapse detectors share: the
/// floor under which a kernel coupling is dropped, and the near-zero level the
/// delta-collapse detectors count against.
/// </summary>
internal static class EdgeWeightThresholds
{
    /// <summary>A kernel coupling at or under this is dropped from the graph.</summary>
    public const double CouplingFloor = 1e-10;

    /// <summary>A coupling under this counts as collapsed in the delta-collapse detectors.</summary>
    public const double NearZero = 1e-8;

    /// <summary>
    /// Fraction of the edges the kernel produced whose coupling is under
    /// <see cref="NearZero"/>. The <paramref name="cutEdges"/> (undirected) the floor already
    /// dropped count as collapsed, so the ratio reads the edge set before the cut. The graph
    /// stores each edge in two CSR slots, so the cut edges count twice, as the kept ones do.
    /// </summary>
    public static double NearZeroRatio(CsrGraph graph, int cutEdges)
    {
        long cutSlots = 2L * cutEdges;
        long collapsed = cutSlots;
        double[] weights = graph.Weights;
        for (int i = 0; i < weights.Length; i++)
            if (weights[i] < NearZero) collapsed++;

        long total = weights.Length + cutSlots;
        return total == 0 ? 0.0 : (double)collapsed / total;
    }
}
