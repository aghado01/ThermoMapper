using System;
using System.Collections.Generic;
using Clustering.Primitives;
using Graphs;
using Graphs.Primitives;

namespace Clustering.Graphical.CumulativeField;

/// <summary>
/// Focused cumulative-field clustering backend: prepare the energy ladder once, then
/// emit a partition at each temperature in a grid. Draw-free and embarrassingly parallel —
/// the closed-form kernel makes each temperature an O(E) evaluation after the
/// one-time O(E log E) prepare, so a whole sweep is cheap and the prepared context
/// is immutable (safe to fan out per temperature). The thermomapper drives
/// clustering through this; phase selection over the sweep is a separate (Phase C)
/// concern.
/// </summary>
public static class CumulativeFieldSweep
{
    /// <summary>
    /// Prepare once for <paramref name="field"/>/<paramref name="symmetrization"/>,
    /// then cluster at every temperature in <paramref name="temperatures"/>
    /// (caller order preserved). Each temperature is an independent closed-form
    /// evaluation; non-positive temperatures are rejected by the kernel.
    /// </summary>
    public static CumulativeFieldSweepResult Run(
        CsrGraph graph,
        EdgeWeightKind weightKind,
        IReadOnlyList<double> temperatures,
        Field field,
        SymmetrizationRule rule = SymmetrizationRule.Mutual,
        double theta = 0.5)
    {
        ArgumentNullException.ThrowIfNull(temperatures);
        if (temperatures.Count == 0)
            throw new ArgumentException("Temperature grid must be non-empty.", nameof(temperatures));

        CumulativeFieldContext context = CumulativeFieldSolver.Prepare(graph, weightKind, field, rule);

        var temps = new double[temperatures.Count];
        var partitions = new Assignment[temperatures.Count];
        for (int i = 0; i < temperatures.Count; i++)
        {
            temps[i] = temperatures[i];
            partitions[i] = CumulativeFieldSolver.Cluster(context, temperatures[i], theta); // validates T > 0
        }

        return new CumulativeFieldSweepResult
        {
            Field = field,
            Symmetrization = rule,
            Theta = theta,
            Temperatures = temps,
            Partitions = partitions,
        };
    }
}
