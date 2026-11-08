using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// The posterior's partial derivatives at given cases (KHAOS slice 4, E8's input; <see cref="BassFit.Gradient"/>): for each
/// continuous predictor v, ∂f/∂z_v of each selected draw's surface at each case and, for a direct functional fit, each training output location, in
/// the predictor's own units. The terms' derivatives come through the design kernels (O-X4), and the chain rule through the range
/// scale divides by the box's width. At a hinge's knot the derivative is the mean of its one-sided limits.
/// <para>A pointwise readout, under no measure: a derivative moment under any measure, a dependent one included, is its average over
/// the measure's points, which is the route E8 and E25 take where the moment service's product-measure closed forms do not hold.
/// Under the product uniform law on a Legendre fit's boxes, the Jacobi-weighted measure has a closed form
/// (<see cref="Sensitivity.PolynomialChaosSensitivity.Jacobi"/>).</para>
/// </summary>
/// <param name="Inputs">The continuous predictors, in the data's order.</param>
/// <param name="Partials">One draw × case × location array per continuous predictor, in <paramref name="Inputs"/>' order.</param>
/// <param name="Weights">Each selected draw's weight, normalized to sum to 1.</param>
/// <param name="CaseOutsideRange">Whether each case lies outside the training range of some continuous predictor.</param>
public sealed record BassGradient(IReadOnlyList<string> Inputs, IReadOnlyList<double[,,]> Partials, IReadOnlyList<double> Weights,
                                  IReadOnlyList<bool> CaseOutsideRange)
{
    /// <summary>The partials of the predictor named <paramref name="name"/>.</summary>
    public double[,,] Of(string name)
    {
        ArgumentNullException.ThrowIfNull(name);
        for (int v = 0; v < Inputs.Count; v++)
            if (Inputs[v] == name) return Partials[v];
        throw new KeyNotFoundException($"No continuous predictor is named \"{name}\".");
    }
}

/// <summary>The gradient readout's computation (<see cref="BassGradient"/>).</summary>
internal static class BassGradients
{
    public static BassGradient Of(BassModel model, BassTermTable terms, IReadOnlyList<IReadOnlyList<BassDraw>> draws, IReadOnlyList<double>? drawWeights,
                                  PredictionCases cases, IReadOnlyList<int>? selection)
    {
        ArgumentNullException.ThrowIfNull(cases);
        (List<BassDraw> pooled, _) = BassDrawFit.Pool(draws, terms, false, nameof(draws));
        (int[] selected, double[] weights) = Draws3.Select(selection, pooled.Count, drawWeights, nameof(selection));
        (CaseInputs caseInputs, bool[] caseOutside) = cases.Prepare(model.Transform, model.Inputs, UnseenLevel.Refuse, model.Data.Grouping?.Name);   // u_g is flat in x
        OutputInputs? locations = null;
        if (model.Data.Response is FunctionalResponse functional)
        {
            ContinuousColumn[] given = Array.ConvertAll([.. functional.OutputVariables], v => new ContinuousColumn(v.Name, v.Values.ToArray()));
            (locations, _) = PredictionCases.PrepareLocations(model.Transform, given);
        }
        var surface = new BassSurface(model.Basis, terms, caseInputs, locations);

        // Each continuous predictor's position among the continuous predictors, its name and its box's width.
        var names = new List<string>();
        var widths = new List<double>();
        for (int k = 0; k < model.Data.Predictors.Count; k++)
        {
            if (model.Data.Predictors[k] is not ContinuousPredictor c) continue;
            var scale = (RangeScale)model.Transform.Predictors[k];
            names.Add(c.Name);
            widths.Add(scale.Max - scale.Min);
        }

        int s = selected.Length, n = surface.Cases, m = surface.Locations;
        var partials = new double[names.Count][,,];
        var f = new double[n * m];
        for (int v = 0; v < names.Count; v++)
        {
            var values = new double[s, n, m];
            double chain = 1.0 / widths[v];
            for (int d = 0; d < s; d++)
            {
                BassDraw draw = pooled[selected[d]];
                surface.EvaluateDerivative(draw.Structure.TermIds, draw.Coefficients, v, f);
                for (int i = 0; i < n; i++)
                    for (int l = 0; l < m; l++) values[d, i, l] = chain * f[i * m + l];
            }
            partials[v] = values;
        }
        return new BassGradient(names, partials, weights, caseOutside);
    }
}
