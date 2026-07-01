using System;

namespace Clustering.Graphical.SPC.Runtime.Core.Sampler;

/// <summary>
/// Reduces a Swendsen–Wang <see cref="Accumulator"/>'s per-node sums to per-node <b>landscapes</b> —
/// the 0-form siblings of <see cref="SwCurrencies"/>' per-edge 1-form reductions. Each landscape is a
/// length-N <c>double[]</c>: an existing sweep scalar kept un-collapsed as a per-node field, then
/// averaged over draws (<c>landscape[i] = sum[i] / DrawCount</c>). <c>SumClusterSizePerNode →</c>
/// mean FK cluster size (the un-reduced χ; high on a dense core); <c>SumInGiantClusterPerNode →</c>
/// giant-cluster participation in <c>[0,1]</c> (the un-reduced M order parameter, per node). A
/// downstream resolution step ascends these surfaces.
/// </summary>
/// <remarks>
/// These exist only where the SW pass materialized the per-node arrays
/// (<see cref="AccumulationSpec.ClusterSizeLandscape"/> / <see cref="AccumulationSpec.OrderLandscape"/>);
/// on an accumulator that did not track them the source array is null and the reduction throws.
/// </remarks>
public static class SwLandscapes
{
    /// <summary>
    /// Reduce the per-node FK-cluster-size sums to the mean-cluster-size landscape:
    /// <c>landscape[i] = SumClusterSizePerNode[i] / DrawCount</c>.
    /// </summary>
    public static double[] MeanClusterSize(Accumulator accumulator)
    {
        ArgumentNullException.ThrowIfNull(accumulator);
        double[] sums = accumulator.SumClusterSizePerNode
            ?? throw new InvalidOperationException(
                "Accumulator carries no SumClusterSizePerNode — MeanClusterSize requires AccumulationSpec.ClusterSizeLandscape = true.");

        return PerNodeRate(sums, accumulator.DrawCount);
    }

    /// <summary>
    /// Reduce the per-node giant-participation sums to the order-parameter landscape:
    /// <c>landscape[i] = SumInGiantClusterPerNode[i] / DrawCount</c>, in <c>[0,1]</c>.
    /// </summary>
    public static double[] GiantParticipation(Accumulator accumulator)
    {
        ArgumentNullException.ThrowIfNull(accumulator);
        double[] sums = accumulator.SumInGiantClusterPerNode
            ?? throw new InvalidOperationException(
                "Accumulator carries no SumInGiantClusterPerNode — GiantParticipation requires AccumulationSpec.OrderLandscape = true.");

        return PerNodeRate(sums, accumulator.DrawCount);
    }

    private static double[] PerNodeRate(double[] sums, int draws)
    {
        if (draws <= 0)
            throw new InvalidOperationException(
                $"DrawCount must be positive to reduce a landscape; was {draws}.");

        double inv = 1.0 / draws;
        var landscape = new double[sums.Length];
        for (int i = 0; i < sums.Length; i++)
            landscape[i] = sums[i] * inv;
        return landscape;
    }
}
