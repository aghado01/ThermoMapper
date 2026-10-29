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
}

/// <summary>An input's effect grid: points in the input's own units for a continuous input, level indices with their names
/// for a categorical one.</summary>
public sealed record EffectGrid(double[] Values, IReadOnlyList<string>? Levels);

/// <summary>
/// One effect function per draw: for one input, the main effect f_v(x) = E[f | x_v = x] − f₀; for two, the two-way
/// effect f_vw = E[f | x_v, x_w] − f_v − f_w − f₀, the surface's ANOVA component on the pair. Values are draws × cells,
/// the cells running over the grids with the last input fastest.
/// </summary>
public sealed record Effect(int[] Inputs, IReadOnlyList<EffectGrid> Grids, double[,] Values);

/// <summary>A surface's effect functions under the declared measure and degree, with each draw's mean f₀ = E[f].</summary>
public sealed record EffectFunctions(IReadOnlyList<string> Inputs, IReadOnlyList<Effect> Effects, double[] Mean);

/// <summary>Effect functions of a BASS fit (<see cref="BassSensitivity.Effects"/>): one set for a scalar fit or the integrated
/// functional form, one per coordinate for the conditioned form, the draws in <see cref="SobolDecomposition"/>'s order.</summary>
public sealed record EffectAnalysis(FunctionalForm? Form, ProductMeasure Measure, IReadOnlyList<EffectFunctions> At,
                                    IReadOnlyList<double[]>? Coordinates);
