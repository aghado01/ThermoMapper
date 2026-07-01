using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Primitives;

namespace Clustering.Graphical.SPC.Profiling;

/// <summary>Which per-node sink supplies the height field.</summary>
public enum SwLandscapeSink
{
    /// <summary>Mean FK cluster size per node (the un-reduced χ) — the canonical climbing surface.</summary>
    MeanClusterSize,

    /// <summary>Giant-cluster participation per node in [0,1] (the un-reduced M) — a detector field
    /// (plateau-heavy; poorly conditioned for ascent).</summary>
    GiantParticipation,
}

/// <summary>
/// Assembles the thermal <see cref="Landscape"/> carrier from a rich sweep's
/// accumulators: per temperature, the per-node sums and draw counts are
/// POOLED across replicas and divided once (the sinks are linear per node, so
/// pooling is exact — no Jensen hazard), then the per-T columns stack
/// grid-major along the ascending temperature axis. The sweep-tier
/// pool-then-reduce twin of <see cref="SweepProfile.From"/>, at degree 0 per
/// node instead of fully collapsed.
/// </summary>
public static class SweepLandscapes
{
    /// <summary>
    /// Builds the thermal landscape from sweep frames. Every frame must carry
    /// the selected sink's per-node array (rich sweeps post-purge do:
    /// <see cref="AccumulationSpec.ClusterSizeLandscape"/> /
    /// <see cref="AccumulationSpec.OrderLandscape"/>); throws otherwise.
    /// </summary>
    public static Landscape FromFrames(
        IReadOnlyList<Accumulator> frames,
        SwLandscapeSink sink,
        string graphId = "unspecified")
    {
        ArgumentNullException.ThrowIfNull(frames);
        if (frames.Count == 0)
            throw new ArgumentException("At least one frame is required.", nameof(frames));

        var groups = frames
            .GroupBy(f => f.Temperature)
            .OrderBy(g => g.Key)
            .ToArray();

        var grid    = new double[groups.Length];
        var columns = new double[groups.Length][];

        for (int t = 0; t < groups.Length; t++)
        {
            grid[t] = groups[t].Key;
            double[]? pooled = null;
            long draws = 0;

            foreach (Accumulator frame in groups[t])
            {
                double[] sums = sink switch
                {
                    SwLandscapeSink.MeanClusterSize => frame.SumClusterSizePerNode
                        ?? throw new InvalidOperationException(
                            $"Frame at T={frame.Temperature:G4} carries no SumClusterSizePerNode — " +
                            "MeanClusterSize requires AccumulationSpec.ClusterSizeLandscape = true."),
                    SwLandscapeSink.GiantParticipation => frame.SumInGiantClusterPerNode
                        ?? throw new InvalidOperationException(
                            $"Frame at T={frame.Temperature:G4} carries no SumInGiantClusterPerNode — " +
                            "GiantParticipation requires AccumulationSpec.OrderLandscape = true."),
                    _ => throw new ArgumentOutOfRangeException(nameof(sink)),
                };

                pooled ??= new double[sums.Length];
                if (pooled.Length != sums.Length)
                    throw new InvalidOperationException(
                        $"Inconsistent node counts across frames at T={frame.Temperature:G4}.");
                for (int i = 0; i < sums.Length; i++) pooled[i] += sums[i];
                draws += frame.DrawCount;
            }

            if (draws <= 0)
                throw new InvalidOperationException(
                    $"Pooled draw count at T={grid[t]:G4} must be positive.");

            double inv = 1.0 / draws;
            for (int i = 0; i < pooled!.Length; i++) pooled[i] *= inv;
            columns[t] = pooled;
        }

        return Landscape.Create(
            axis: "temperature",
            grid: grid,
            valuesByGridPoint: columns,
            provenance: new LandscapeProvenance(sink.ToString(), graphId));
    }
}
