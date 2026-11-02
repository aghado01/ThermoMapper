using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Dendrograms;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Profiling;

/// <summary>The closure gap at each temperature of a sweep, ascending.</summary>
/// <param name="Theta">The cut on the co-membership scale, (qθ − 1)/(q − 1) for θ on the δ̄ scale.</param>
/// <param name="Temperatures">The ascending grid.</param>
/// <param name="Results">The gap from everything pooled at each temperature, its costs read off the estimate that
/// decided, so sampling noise inflates it.</param>
/// <param name="CrossFitGap">The gap with one half of the replicas deciding and the other evaluating, averaged over the
/// two ways round, which removes that inflation. NaN at a temperature with one replica.</param>
/// <param name="Integrated">The pooled gap integrated over θ, which noise also inflates, by about a pair's variance each.</param>
public sealed record ClosureGapCurve(
    double Theta,
    double[] Temperatures,
    ClosureGapResult[] Results,
    double[] CrossFitGap,
    double[] Integrated);

/// <summary>
/// Reduces a sweep's co-membership to <see cref="ClosureGap"/> curves over temperature, on the same cut
/// <see cref="Partitions.Hierarchical.SweepDenseTStack"/> makes: how much each level of the T-stack pays for closing
/// its θ-graph.
/// <list type="bullet">
/// <item><see cref="OverEdges"/> reads every draw, through <c>AccumulationSpec.CoMembership</c>.</item>
/// <item><see cref="OverPairs"/> reads the kept draws, through <c>AccumulationSpec.PartitionEvery</c>, so under thinning
/// its closure can differ from the T-stack's near θ.</item>
/// </list>
/// </summary>
public static class SweepClosureGap
{
    /// <summary>The cut on the co-membership scale that matches <paramref name="theta"/> on the δ̄ scale.</summary>
    public static double CoMembershipTheta(double theta, int q)
    {
        if (q < 2)
            throw new ArgumentOutOfRangeException(nameof(q), q, "The δ̄ scale needs q ≥ 2.");
        if (theta < 1.0 / q || theta > 1.0)
            throw new ArgumentOutOfRangeException(nameof(theta), theta,
                $"Every δ̄ is at least 1/q = {1.0 / q:G4}, so a θ below it passes every edge and reads as no cost ratio.");
        return (q * theta - 1.0) / (q - 1.0);
    }

    /// <summary>The curve over the graph's edges, from frames carrying co-membership counts.</summary>
    /// <param name="graph">The CSR graph the frames were sampled on.</param>
    /// <param name="frames">A rich sweep's accumulators, one or more per temperature.</param>
    /// <param name="theta">The cut on the δ̄ scale, as <see cref="Partitions.Hierarchical.SweepDenseTStack.Build"/> takes it.</param>
    public static ClosureGapCurve OverEdges(CsrGraph graph, IReadOnlyList<Accumulator> frames, double theta = 0.5)
        => Reduce<double[]>(
            frames, theta,
            field: Column,
            pooled: (group, _, _) => Column(group),
            gap: (field, cut, evaluate) => ClosureGap.OverEdges(graph, field, cut, evaluate),
            integrated: field => ClosureGap.IntegratedOverEdges(graph, field));

    /// <summary>The curve over all pairs, from frames carrying kept partitions.</summary>
    /// <param name="graph">The CSR graph the frames were sampled on.</param>
    /// <param name="frames">A sweep's accumulators with kept partitions, one or more per temperature.</param>
    /// <param name="theta">The cut on the δ̄ scale, as <see cref="Partitions.Hierarchical.SweepDenseTStack.Build"/> takes it.</param>
    public static ClosureGapCurve OverPairs(CsrGraph graph, IReadOnlyList<Accumulator> frames, double theta = 0.5)
        => Reduce<DenseCoMembership>(
            frames, theta,
            field: group => Dense(graph, group),
            pooled: (group, first, second) => first is null ? Dense(graph, group) : DenseCoMembership.Pool(first, second!),
            gap: (field, cut, evaluate) => ClosureGap.OverPairs(graph, field, cut, evaluate),
            integrated: field => ClosureGap.IntegratedOverPairs(graph, field));

    private static ClosureGapCurve Reduce<TField>(
        IReadOnlyList<Accumulator> frames,
        double theta,
        Func<Accumulator[], TField> field,
        Func<Accumulator[], TField?, TField?, TField> pooled,
        Func<TField, double, TField?, ClosureGapResult> gap,
        Func<TField, double> integrated)
        where TField : class
    {
        ArgumentNullException.ThrowIfNull(frames);
        if (frames.Count == 0)
            throw new ArgumentException("At least one frame is required.", nameof(frames));

        int q = frames[0].Q;
        foreach (Accumulator frame in frames)
            if (frame.Q != q)
                throw new InvalidOperationException($"Frames disagree on q ({frame.Q} vs {q}).");
        double cut = CoMembershipTheta(theta, q);

        Accumulator[][] groups = frames
            .GroupBy(f => f.Temperature)
            .OrderBy(g => g.Key)
            .Select(g => g.OrderBy(f => f.ReplicaIndex).ToArray())
            .ToArray();

        var temperatures = new double[groups.Length];
        var results = new ClosureGapResult[groups.Length];
        var crossFit = new double[groups.Length];
        var integrals = new double[groups.Length];
        for (int t = 0; t < groups.Length; t++)
        {
            Accumulator[] group = groups[t];
            temperatures[t] = group[0].Temperature;

            // The replicas split alternately, so each half keeps its own draws.
            TField? first = null, second = null;
            if (group.Length > 1)
            {
                first = field(group.Where((_, r) => r % 2 == 0).ToArray());
                second = field(group.Where((_, r) => r % 2 == 1).ToArray());
            }

            TField all = pooled(group, first, second);
            results[t] = gap(all, cut, null);
            crossFit[t] = first is null ? double.NaN : 0.5 * (gap(first, cut, second).Gap + gap(second!, cut, first).Gap);
            integrals[t] = integrated(all);
        }
        return new ClosureGapCurve(cut, temperatures, results, crossFit, integrals);
    }

    private static double[] Column(Accumulator[] frames)
    {
        long[]? pooled = null;
        long draws = 0;
        foreach (Accumulator frame in frames)
        {
            int[] counts = frame.CoMembershipCount ?? throw new InvalidOperationException(
                $"Frame at T={frame.Temperature:G4} carries no CoMembershipCount. The edge gap needs " +
                "AccumulationSpec.CoMembership (CLI: --accumulation comembership).");
            pooled ??= new long[counts.Length];
            if (pooled.Length != counts.Length)
                throw new InvalidOperationException($"Inconsistent slot counts across frames at T={frame.Temperature:G4}.");
            for (int e = 0; e < counts.Length; e++) pooled[e] += counts[e];
            draws += frame.DrawCount;
        }
        if (draws <= 0)
            throw new InvalidOperationException($"Pooled draw count at T={frames[0].Temperature:G4} must be positive.");

        var column = new double[pooled!.Length];
        for (int e = 0; e < column.Length; e++) column[e] = (double)pooled[e] / draws;
        return column;
    }

    private static DenseCoMembership Dense(CsrGraph graph, Accumulator[] frames)
    {
        var runs = new int[frames.Length][];
        for (int r = 0; r < frames.Length; r++)
            runs[r] = frames[r].KeptPartitions ?? throw new InvalidOperationException(
                $"Frame at T={frames[r].Temperature:G4} keeps no partitions. The pair gap needs " +
                "AccumulationSpec.PartitionEvery (CLI: --accumulation partitions=K).");
        return DenseCoMembership.FromPartitions(graph.NodeCount, runs);
    }
}
