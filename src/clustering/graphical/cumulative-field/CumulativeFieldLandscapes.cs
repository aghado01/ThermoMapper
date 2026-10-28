using System;
using System.Collections.Generic;
using Clustering.Primitives;
using Graphs.Observables;
using Graphs.Primitives;

namespace Clustering.Graphical.CumulativeField;

/// <summary>Which per-node marginal supplies the height field.</summary>
public enum CumulativeFieldLandscapeSink
{
    /// <summary>Σ_j G_ij(T) — the attachment field; ascend.</summary>
    BondMass,

    /// <summary>Σ_j J_ij·(1−G_ij(T)) — the frustration field; descend.</summary>
    LocalEnergy,
}

/// <summary>
/// The solver-side <see cref="Landscape"/> producer: evaluates the cumulative-field
/// closed-form affinity at each requested temperature and reduces it to a
/// per-node marginal column — no draws, no accumulator, exact. The second
/// independent mint of the landscape currency (Swendsen–Wang's
/// <c>SweepLandscapes.FromFrames</c> is the sampled one), proving the
/// two-producer seam in code.
/// </summary>
/// <remarks>
/// Deliberately written against the stable <see cref="CumulativeFieldSolver.Prepare"/> /
/// <see cref="CumulativeFieldSolver.Solve"/> lip only — the solver's internals can be
/// rebuilt under it.
/// </remarks>
public static class CumulativeFieldLandscapes
{
    /// <summary>
    /// Mints the thermal landscape over an ascending temperature grid.
    /// </summary>
    public static Landscape OverGrid(
        CumulativeFieldContext context,
        IReadOnlyList<double> temperatures,
        CumulativeFieldLandscapeSink sink,
        string graphId = "unspecified")
    {
        ArgumentNullException.ThrowIfNull(context);
        ArgumentNullException.ThrowIfNull(temperatures);
        if (temperatures.Count == 0)
            throw new ArgumentException("At least one temperature is required.", nameof(temperatures));

        var grid = new double[temperatures.Count];
        for (int t = 0; t < grid.Length; t++) grid[t] = temperatures[t];

        var columns = new double[grid.Length][];
        for (int t = 0; t < grid.Length; t++)
        {
            Affinities affinities = CumulativeFieldSolver.Solve(context, grid[t]);
            columns[t] = sink switch
            {
                CumulativeFieldLandscapeSink.BondMass    => AffinityNodeMarginals.BondMass(context.Graph, affinities),
                CumulativeFieldLandscapeSink.LocalEnergy => AffinityNodeMarginals.LocalEnergy(context.Graph, affinities),
                _ => throw new ArgumentOutOfRangeException(nameof(sink)),
            };
        }

        return Landscape.Create(
            axis: "temperature",
            grid: grid,
            valuesByGridPoint: columns,
            provenance: new LandscapeProvenance(sink.ToString(), graphId, GaugeNote: "cumulative-field:closed-form"));
    }
}
