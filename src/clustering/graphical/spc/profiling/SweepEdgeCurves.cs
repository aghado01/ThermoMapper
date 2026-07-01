using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;

namespace Clustering.Graphical.SPC.Profiling;

/// <summary>
/// Stacks a rich sweep's per-edge co-membership counts into grid-major eq-4
/// discriminant columns over the ascending temperature grid — the
/// structure-side twin of <see cref="SweepLandscapes"/>. Counts and draws are
/// POOLED across replicas per temperature and divided once; the eq-4
/// transform <c>δ̄ = ((q−1)·⟨n_ij⟩ + 1)/q</c> (affine, so pooling commutes)
/// puts the columns on the same θ scale as the chosen-T
/// <c>ThresholdCoMembership</c> cut.
/// </summary>
public static class SweepEdgeCurves
{
    public static (double[] Temperatures, double[][] DeltaByGridPoint) CoMembershipDelta(
        IReadOnlyList<Accumulator> frames)
    {
        ArgumentNullException.ThrowIfNull(frames);
        if (frames.Count == 0)
            throw new ArgumentException("At least one frame is required.", nameof(frames));

        int q = frames[0].Q;
        var groups = frames
            .GroupBy(f => f.Temperature)
            .OrderBy(g => g.Key)
            .ToArray();

        var temps   = new double[groups.Length];
        var columns = new double[groups.Length][];

        for (int t = 0; t < groups.Length; t++)
        {
            temps[t] = groups[t].Key;
            long[]? pooled = null;
            long draws = 0;

            foreach (Accumulator frame in groups[t])
            {
                int[] counts = frame.CoMembershipCount ?? throw new InvalidOperationException(
                    $"Frame at T={frame.Temperature:G4} carries no CoMembershipCount — the thermal " +
                    "composition needs AccumulationSpec.CoMembership (CLI: --accumulation comembership).");
                if (frame.Q != q)
                    throw new InvalidOperationException($"Frames disagree on q ({frame.Q} vs {q}).");
                pooled ??= new long[counts.Length];
                if (pooled.Length != counts.Length)
                    throw new InvalidOperationException(
                        $"Inconsistent slot counts across frames at T={frame.Temperature:G4}.");
                for (int e = 0; e < counts.Length; e++) pooled[e] += counts[e];
                draws += frame.DrawCount;
            }

            if (draws <= 0)
                throw new InvalidOperationException($"Pooled draw count at T={temps[t]:G4} must be positive.");

            var delta = new double[pooled!.Length];
            double inv = 1.0 / draws;
            for (int e = 0; e < delta.Length; e++)
            {
                double rate = pooled[e] * inv;
                delta[e] = ((q - 1.0) * rate + 1.0) / q;
            }
            columns[t] = delta;
        }

        return (temps, columns);
    }
}
