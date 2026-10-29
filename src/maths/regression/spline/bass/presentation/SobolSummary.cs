using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass.Sensitivity;

namespace Maths.Regression.Spline.Bass.Presentation;

/// <summary>
/// A decomposition's posterior summary (R's <c>plot.bassSob</c>; P11), reducing its draws to values under its weights:
/// <list type="bullet">
/// <item>each subset's main index and each input's total index over the defined draws, the zero-variance draws left out
/// (their indices are undefined), whose means are the decomposition's own (<see cref="SobolDecomposition.MeanMainIndices"/>,
/// <see cref="SobolDecomposition.MeanTotalIndices"/>), the same sums in the same order;</item>
/// <item>each subset's main variance, each input's total-effect variance and the total variance, in the response's units,
/// over every draw, since a zero variance is a value: R's partitioned variance;</item>
/// <item>the draws, the defined ones, and those clamped or negative beyond the tolerance.</item>
/// </list>
/// The subsets are input indices, smaller first, into <see cref="Inputs"/>.
/// </summary>
public sealed record SobolSummary(
    IReadOnlyList<string> Inputs, IReadOnlyList<int[]> Subsets,
    IReadOnlyList<PosteriorSummary> MainIndex, IReadOnlyList<PosteriorSummary> TotalIndex,
    IReadOnlyList<PosteriorSummary> MainVariance, IReadOnlyList<PosteriorSummary> TotalVariance, PosteriorSummary Variance,
    int Draws, int Defined, int Clamped, int NegativeBeyondTolerance)
{
    /// <summary>The summary of <paramref name="decomposition"/> at <paramref name="probabilities"/>
    /// (<see cref="PosteriorSummary.DefaultProbabilities"/> when null).</summary>
    public static SobolSummary Of(SobolDecomposition decomposition, IReadOnlyList<double>? probabilities = null)
        => Over(decomposition, PosteriorSummary.Checked(probabilities));

    internal static SobolSummary Over(SobolDecomposition decomposition, IReadOnlyList<double> p)
    {
        ArgumentNullException.ThrowIfNull(decomposition);
        int n = decomposition.Count;
        int[] defined = [.. Enumerable.Range(0, n).Where(decomposition.IsDefined)];
        double[]? weights = decomposition.Weights;
        double[]? definedWeights = weights is null ? null : [.. defined.Select(d => weights[d])];

        PosteriorSummary[] Indices(int count, Func<int, int, double> index)
        {
            var summaries = new PosteriorSummary[count];
            for (int k = 0; k < count; k++) summaries[k] = PosteriorSummary.Over([.. defined.Select(d => index(d, k))], definedWeights, p);
            return summaries;
        }
        PosteriorSummary[] Variances(int count, double[,] table)
        {
            var summaries = new PosteriorSummary[count];
            for (int k = 0; k < count; k++) summaries[k] = PosteriorSummary.Over([.. Enumerable.Range(0, n).Select(d => table[d, k])], weights, p);
            return summaries;
        }

        int subsets = decomposition.Subsets.Count, inputs = decomposition.Inputs.Count;
        int clamped = decomposition.Status.Count(s => (s & DecompositionStatus.Clamped) != 0);
        int negative = decomposition.Status.Count(s => (s & DecompositionStatus.NegativeBeyondTolerance) != 0);
        return new SobolSummary(decomposition.Inputs, decomposition.Subsets,
            Indices(subsets, decomposition.MainIndex), Indices(inputs, decomposition.TotalIndex),
            Variances(subsets, decomposition.MainVariance), Variances(inputs, decomposition.TotalVariance),
            PosteriorSummary.Over(decomposition.Variance, weights, p), n, defined.Length, clamped, negative);
    }
}

/// <summary>
/// A direct fit's sensitivity summary (R's <c>plot.bassSob</c>; P11): one <see cref="SobolSummary"/> for a scalar fit or the
/// integrated functional form, one per coordinate for the conditioned form, with the analysis's form, measure, order and
/// moment error, and the conditioned form's coordinates, one value per output variable in <see cref="OutputVariables"/>'
/// order.
/// </summary>
public sealed record SobolAnalysisSummary(
    FunctionalForm? Form, ProductMeasure Measure, int MaxOrder, IReadOnlyList<SobolSummary> Decompositions,
    IReadOnlyList<double[]>? Coordinates, IReadOnlyList<bool>? CoordinatesOutside, IReadOnlyList<string> OutputVariables, double MomentError)
{
    internal static SobolAnalysisSummary Of(BassFit fit, SobolAnalysis analysis, IReadOnlyList<double>? probabilities)
    {
        IReadOnlyList<double> p = PosteriorSummary.Checked(probabilities);
        return new SobolAnalysisSummary(analysis.Form, analysis.Measure, analysis.MaxOrder,
            Array.AsReadOnly(analysis.Decompositions.Select(d => SobolSummary.Over(d, p)).ToArray()), analysis.Coordinates,
            analysis.CoordinatesOutside, Array.AsReadOnly(fit.Transform.OutputVariables.Select(v => v.Name).ToArray()), analysis.MomentError);
    }
}

/// <summary>
/// A basis fit's sensitivity summary (R's <c>plot.bassSob</c> on <c>sobolBasis</c>; P11), in output units: one
/// <see cref="SobolSummary"/> per output location (<see cref="BasisSobolAnalysis.AtLocation"/>, R's functional panels) and
/// the integrated form's (<see cref="BasisSobolAnalysis.Integrated"/>), with the analysis's measure, orders and moment
/// error.
/// </summary>
public sealed record BasisSobolSummary(
    ProductMeasure Measure, int MaxOrder, int IntegratedMaxOrder, IReadOnlyList<SobolSummary> PerLocation, SobolSummary Integrated,
    double MomentError)
{
    /// <summary>The summary of <paramref name="analysis"/> at <paramref name="probabilities"/>
    /// (<see cref="PosteriorSummary.DefaultProbabilities"/> when null), each location's decomposition built in turn.</summary>
    public static BasisSobolSummary Of(BasisSobolAnalysis analysis, IReadOnlyList<double>? probabilities = null)
    {
        ArgumentNullException.ThrowIfNull(analysis);
        IReadOnlyList<double> p = PosteriorSummary.Checked(probabilities);
        var perLocation = new SobolSummary[analysis.Locations];
        for (int l = 0; l < perLocation.Length; l++) perLocation[l] = SobolSummary.Over(analysis.AtLocation(l), p);
        return new BasisSobolSummary(analysis.Measure, analysis.MaxOrder, analysis.IntegratedMaxOrder, Array.AsReadOnly(perLocation),
                                     SobolSummary.Over(analysis.Integrated, p), analysis.MomentError);
    }
}
