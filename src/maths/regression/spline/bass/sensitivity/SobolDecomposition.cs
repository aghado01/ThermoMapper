using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>What sensitivity decomposes (O-X1): the surface's functional of the response.</summary>
public enum DecomposedSurface
{
    /// <summary>E[y | x]: for Gaussian BASS the latent surface Ba itself (D08).</summary>
    Mean,
}

/// <summary>How a direct functional fit's output variables enter sensitivity (D08, JSS §3.1).</summary>
public enum FunctionalForm
{
    /// <summary>The output variables are more inputs, integrated under their declared measures: one decomposition.</summary>
    Integrate,

    /// <summary>The surface is held at each declared output coordinate: one decomposition per coordinate, indices as
    /// functions of t.</summary>
    Condition,
}

/// <summary>A draw's decomposition status.</summary>
[Flags]
public enum DecompositionStatus
{
    /// <summary>The indices are defined, every contribution nonnegative.</summary>
    Defined = 0,

    /// <summary>The total variance is zero within rounding: the indices are undefined (NaN), never 0 or 1.</summary>
    ZeroVariance = 1,

    /// <summary>A contribution was negative within the rounding tolerance relative to the total, and is clamped to 0.</summary>
    Clamped = 2,

    /// <summary>A contribution is negative beyond the tolerance; it is kept as computed.</summary>
    NegativeBeyondTolerance = 4,
}

/// <summary>
/// A declared sensitivity analysis (D08): the product measure, the order and budget of the subsets, the decomposed surface,
/// and a functional fit's form. Every field defaults, so <c>new SensitivitySpec()</c> is R's analysis with D08's
/// corrections.
/// </summary>
public sealed record SensitivitySpec
{
    /// <summary>The predictors' measures by name, in their own units; a predictor left out takes the uniform measure on its
    /// training range, or uniform level probabilities.</summary>
    public IReadOnlyDictionary<string, InputMeasure>? Measures { get; init; }

    /// <summary>The output variables' measures by name, for <see cref="FunctionalForm.Integrate"/>; uniform over the training
    /// range when left out, <see cref="EmpiricalMeasure"/> for the training output grid.</summary>
    public IReadOnlyDictionary<string, InputMeasure>? OutputMeasures { get; init; }

    /// <summary>The largest subset whose main index is reported; the posterior's largest term (its variable count) when null,
    /// so every order is reported and the main indices sum to 1.</summary>
    public int? MaxOrder { get; init; }

    /// <summary>The most subsets the analysis enumerates: the posterior's terms spanning more, up to the order, are refused
    /// before anything is computed.</summary>
    public int SubsetBudget { get; init; } = 10_000;

    /// <summary>The decomposed surface.</summary>
    public DecomposedSurface Surface { get; init; } = DecomposedSurface.Mean;

    /// <summary>A direct functional fit's form; ignored for a scalar response.</summary>
    public FunctionalForm Functional { get; init; } = FunctionalForm.Integrate;

    /// <summary>The output coordinates of <see cref="FunctionalForm.Condition"/>, each one value per output variable in its own
    /// units; the training output locations when null.</summary>
    public IReadOnlyList<double[]>? Coordinates { get; init; }

    /// <summary>The rounding tolerance: a total variance at most this times its terms' moment magnitudes is zero, and a
    /// contribution above −this times the total is clamped to 0.</summary>
    public double RoundingTolerance { get; init; } = 1e-10;
}

/// <summary>
/// The variance decomposition of one surface under the product measure, per posterior draw (D08, O6): each subset's main
/// (Sobol) variance, the variance of the surface's ANOVA component on it, and each input's total-effect variance, both in the
/// response's units, with the total variance and the draw's status. The indices are their ratios to the total, views
/// computed on demand; the draws, with <see cref="Weights"/>, are the indices' posterior distribution.
/// </summary>
public sealed class SobolDecomposition
{
    internal SobolDecomposition(IReadOnlyList<string> inputs, int[][] subsets, double[,] main, double[,] total, double[] variance,
                                DecompositionStatus[] status, int[] chains, int[] draws, double[]? weights)
    {
        Inputs = inputs;
        Subsets = Array.AsReadOnly(subsets);
        MainVariance = main;
        TotalVariance = total;
        Variance = variance;
        Status = status;
        Chains = chains;
        Draws = draws;
        Weights = weights;
    }

    /// <summary>The inputs' names: the predictors in the data's order, then the integrated output variables.</summary>
    public IReadOnlyList<string> Inputs { get; }

    /// <summary>The reported subsets, as ascending input indices, smaller first: those inside some posterior term, up to the
    /// order. A subset inside no term has main variance 0.</summary>
    public IReadOnlyList<int[]> Subsets { get; }

    /// <summary>Draws × subsets: each subset's main variance.</summary>
    public double[,] MainVariance { get; }

    /// <summary>Draws × inputs: each input's total-effect variance.</summary>
    public double[,] TotalVariance { get; }

    /// <summary>Each draw's total variance of the surface.</summary>
    public double[] Variance { get; }

    /// <summary>Each draw's status.</summary>
    public DecompositionStatus[] Status { get; }

    /// <summary>Each draw's chain, the draws in chain order and each chain's in draw order.</summary>
    public int[] Chains { get; }

    /// <summary>Each draw's index within its chain.</summary>
    public int[] Draws { get; }

    /// <summary>The draws' weights; null for equal weights (MCMC).</summary>
    public double[]? Weights { get; }

    /// <summary>The draws.</summary>
    public int Count => Variance.Length;

    /// <summary>Whether draw <paramref name="draw"/>'s indices are defined.</summary>
    public bool IsDefined(int draw) => (Status[draw] & DecompositionStatus.ZeroVariance) == 0;

    /// <summary>Draw <paramref name="draw"/>'s main index of subset <paramref name="subset"/>; NaN where undefined.</summary>
    public double MainIndex(int draw, int subset) => IsDefined(draw) ? MainVariance[draw, subset] / Variance[draw] : double.NaN;

    /// <summary>Draw <paramref name="draw"/>'s total index of input <paramref name="input"/>; NaN where undefined.</summary>
    public double TotalIndex(int draw, int input) => IsDefined(draw) ? TotalVariance[draw, input] / Variance[draw] : double.NaN;

    /// <summary>Each subset's posterior mean main index over the defined draws, under the weights; NaN with none defined.</summary>
    public double[] MeanMainIndices() => Means(Subsets.Count, MainIndex);

    /// <summary>Each input's posterior mean total index over the defined draws, under the weights; NaN with none defined.</summary>
    public double[] MeanTotalIndices() => Means(Inputs.Count, TotalIndex);

    private double[] Means(int count, Func<int, int, double> index)
    {
        var mean = new double[count];
        double mass = 0.0;
        for (int d = 0; d < Count; d++)
        {
            if (!IsDefined(d)) continue;
            double w = Weights?[d] ?? 1.0;
            mass += w;
            for (int k = 0; k < count; k++) mean[k] += w * index(d, k);
        }
        for (int k = 0; k < count; k++) mean[k] = mass > 0.0 ? mean[k] / mass : double.NaN;
        return mean;
    }
}

/// <summary>
/// A sensitivity analysis of a BASS fit (<see cref="BassSensitivity.Of"/>): one <see cref="SobolDecomposition"/> for a scalar
/// fit or the integrated functional form, one per coordinate for the conditioned form.
/// </summary>
/// <param name="Form">The functional form; null for a scalar response.</param>
/// <param name="Measure">The product measure, as declared with defaults filled in and resolved.</param>
/// <param name="MaxOrder">The order the main indices were reported to.</param>
/// <param name="Decompositions">The decomposition, or one per coordinate.</param>
/// <param name="Coordinates">The conditioned form's output coordinates in their own units; null otherwise.</param>
/// <param name="CoordinatesOutside">Whether each coordinate lies outside the training range, where the output hinges
/// extrapolate; null otherwise.</param>
/// <param name="MomentError">The largest quadrature error estimate among the one-variable moments (0 unless a Student
/// component was integrated).</param>
public sealed record SobolAnalysis(
    FunctionalForm? Form,
    ProductMeasure Measure,
    int MaxOrder,
    IReadOnlyList<SobolDecomposition> Decompositions,
    IReadOnlyList<double[]>? Coordinates,
    IReadOnlyList<bool>? CoordinatesOutside,
    double MomentError);
