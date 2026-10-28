using System;
using System.Collections.Generic;
using Clustering.Primitives;
using Graphs.Observables;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Core.Solver;

/// <summary>Which per-node marginal supplies the height field.</summary>
public enum PKWangLandscapeSink
{
    /// <summary>Σ_j G_ij(T) — the attachment field; ascend.</summary>
    BondMass,

    /// <summary>Σ_j J_ij·(1−G_ij(T)) — the frustration field; descend.</summary>
    LocalEnergy,
}

/// <summary>
/// The solver-side <see cref="Landscape"/> producer: evaluates PKWang's
/// closed-form affinity at each requested temperature and reduces it to a
/// per-node marginal column — no draws, no accumulator, exact. The second
/// independent mint of the landscape currency (SW's
/// <c>SweepLandscapes.FromFrames</c> is the sampled one), proving the
/// two-producer seam in code.
/// </summary>
/// <remarks>
/// Deliberately written against the stable <see cref="PKWang.Prepare"/> /
/// <see cref="PKWang.Solve"/> lip only — the solver's internals can be
/// rebuilt under it (see
/// <c>.discussion/issues/spc-samplers/pkwang-solver-review-followup.md</c>).
/// </remarks>
public static class PKWangLandscapes
{
    /// <summary>
    /// Mints the thermal landscape over an ascending temperature grid.
    /// </summary>
    public static Landscape OverGrid(
        PKWangContext context,
        IReadOnlyList<double> temperatures,
        PKWangLandscapeSink sink,
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
            Affinities affinities = PKWang.Solve(context, grid[t]);
            columns[t] = sink switch
            {
                PKWangLandscapeSink.BondMass    => AffinityNodeMarginals.BondMass(context.Graph, affinities),
                PKWangLandscapeSink.LocalEnergy => AffinityNodeMarginals.LocalEnergy(context.Graph, affinities),
                _ => throw new ArgumentOutOfRangeException(nameof(sink)),
            };
        }

        return Landscape.Create(
            axis: "temperature",
            grid: grid,
            valuesByGridPoint: columns,
            provenance: new LandscapeProvenance(sink.ToString(), graphId, GaugeNote: "pkwang:closed-form"));
    }
}
