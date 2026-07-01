namespace Graphs.Models.Potts.Observables;

/// <summary>
/// Magnetic susceptibility χ of the q-state Potts model — the normalized second
/// moment of FK cluster sizes, χ = ⟨Σ sᵢ²⟩ / N. The model's observable vocabulary:
/// what χ <i>means</i> as a Gibbs-measure quantity, independent of the inference
/// strategy that produced the accumulator and of the temperature sweep that
/// assembles per-draw values into a χ(T) curve (that assembly lives in profiling).
/// </summary>
public static class Susceptibility
{
    /// <summary>
    /// Per-draw-set reduction: χ = <paramref name="runningSumSqClusterSizes"/> /
    /// (<paramref name="draws"/> · <paramref name="siteCount"/>). Returns 0 when
    /// there are no draws or no sites.
    /// </summary>
    /// <param name="runningSumSqClusterSizes">Σ over retained draws of the summed squared FK cluster sizes.</param>
    /// <param name="draws">Number of retained MC draws averaged (0 ⇒ 0).</param>
    /// <param name="siteCount">N, the number of sites.</param>
    public static double Reduce(double runningSumSqClusterSizes, long draws, int siteCount)
        => draws > 0 && siteCount > 0
            ? runningSumSqClusterSizes / (draws * (double)siteCount)
            : 0.0;
}
