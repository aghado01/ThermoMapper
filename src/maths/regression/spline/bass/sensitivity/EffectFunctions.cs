using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>Which effect functions to compute (O8), and on what grid.</summary>
public sealed record EffectSpec
{
    /// <summary>The effects, each one input or two by name; when null, every input's main effect and the two-way effect of
    /// every pair some posterior term holds.</summary>
    public IReadOnlyList<string[]>? Subsets { get; init; }

    /// <summary>A continuous input's grid points, evenly spaced over its measure's support, ends included (R's 100).</summary>
    public int GridPoints { get; init; } = 100;

    /// <summary>The most values the effects may hold, draws × cells summed over the effects (and coordinates): more is
    /// refused before anything is computed.</summary>
    public long CellBudget { get; init; } = 10_000_000;

    /// <summary>On the sampled route, the most term values the pinned blocks may evaluate, summed over the batches of draws: each
    /// cell's block is the base rows with its inputs pinned, read for every term its batch holds. More is refused before anything is
    /// computed.</summary>
    public long EvaluationBudget { get; init; } = 2_000_000_000;
}

/// <summary>An input's effect grid: points in the input's own units for a continuous input, level indices with their names
/// for a categorical one.</summary>
public sealed record EffectGrid(double[] Values, IReadOnlyList<string>? Levels);

/// <summary>
/// One effect function per draw: for one input, the main effect f_v(x) = E[f | x_v = x] − f₀; for two, the two-way
/// effect f_vw = E[f | x_v, x_w] − f_v − f_w − f₀, the surface's ANOVA component on the pair. Values are draws × cells,
/// the cells running over the grids with the last input fastest.
/// </summary>
public sealed record Effect(int[] Inputs, IReadOnlyList<EffectGrid> Grids, double[,] Values)
{
    /// <summary>The sampled route's Monte Carlo standard error of each value, draws × cells; null on the exact route.</summary>
    public double[,]? Errors { get; init; }
}

/// <summary>A surface's effect functions under the declared measure and degree, with each draw's mean f₀ = E[f].</summary>
public sealed record EffectFunctions(IReadOnlyList<string> Inputs, IReadOnlyList<Effect> Effects, double[] Mean)
{
    /// <summary>The sampled route's standard error of each draw's mean; null on the exact route.</summary>
    public double[]? MeanError { get; init; }
}

/// <summary>Effect functions of a BASS fit (<see cref="BassSensitivity.Effects"/>): one set for a scalar fit or the integrated
/// functional form, one per coordinate for the conditioned form, the draws in <see cref="SobolDecomposition"/>'s order.</summary>
public sealed record EffectAnalysis(FunctionalForm? Form, ProductMeasure Measure, IReadOnlyList<EffectFunctions> At,
                                    IReadOnlyList<double[]>? Coordinates, IReadOnlyList<string>? CoordinateVariables)
{
    /// <summary>The route that ran (the sampled-sensitivity brief, decision 1).</summary>
    public SensitivityRoute Route { get; init; } = SensitivityRoute.Exact;

    /// <summary>The sampled route's rows N in each of the two base sets the Sobol route draws from the same seed, the effects averaging
    /// over both, 2N rows; 0 on the exact route.</summary>
    public int Points { get; init; }
}
