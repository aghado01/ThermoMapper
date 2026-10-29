using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass.Presentation;

/// <summary>Where a training fit's intervals come from (D07).</summary>
public enum IntervalSource
{
    /// <summary>The predictive mixture Σₛ wₛ·N(fₛ, σₛ²), its quantile solved exactly (<see cref="BassPrediction.Quantile"/>).</summary>
    PredictiveMixture,

    /// <summary>The weighted quantile of the predictive draws, the mixture's up to Monte Carlo error
    /// (<see cref="BassBasisPrediction.Quantile"/>).</summary>
    PredictiveDraws,
}

/// <summary>
/// Fitted against observed at the training inputs (R's <c>plot.bass</c> and <c>plot.bassBasis</c> fit and residual panels;
/// P10, P11): the observed responses, the posterior mean of the latent surface, the predictive intervals' quantiles and the
/// residuals, each case × location, a scalar response with one location. Every value is the prediction view's at the
/// training cases, so nothing is recomputed differently: the mean is <c>Mean()</c>, the quantiles <c>Quantile(p)</c>, of
/// the predictive mixture, never of per-draw interval endpoints as R's plot takes them (D07). A basis fit's values are in
/// physical units.
/// </summary>
/// <param name="Observed">The training responses, case × location; a basis fit's are those its decomposition represents,
/// y_m + y_s ∘ (K w + e) with the in-sample residual (or the residual supplied with the basis).</param>
/// <param name="Mean">The posterior mean of the latent surface, case × location.</param>
/// <param name="Probabilities">The quantiles' probabilities, in the order given.</param>
/// <param name="Quantiles">The predictive quantiles, probability × case × location.</param>
/// <param name="Intervals">How the quantiles were computed.</param>
/// <param name="Residuals">Observed minus the posterior mean, case × location: R's posterior mean residuals.</param>
/// <param name="NoiseSd">A direct fit's posterior mean of σ under the draws' weights, the sd of R's reference normal density
/// for the residuals; null for a basis fit, whose noise varies by location.</param>
/// <param name="LocationCoordinates">A direct functional fit's output variables at the training locations, by name; null
/// for a scalar fit and a basis fit, whose locations are indices.</param>
public sealed record BassTrainingFit(
    double[,] Observed, double[,] Mean, IReadOnlyList<double> Probabilities, double[,,] Quantiles, IntervalSource Intervals,
    double[,] Residuals, double? NoiseSd, IReadOnlyList<ContinuousColumn>? LocationCoordinates)
{
    /// <summary>The cases.</summary>
    public int Cases => Observed.GetLength(0);

    /// <summary>The locations: 1 for a scalar response.</summary>
    public int Locations => Observed.GetLength(1);

    // A direct fit's: its prediction at the training cases and locations, latent draws only, the mixture solved exactly.
    internal static BassTrainingFit Of(BassFit fit, BassData data, IReadOnlyList<double>? probabilities)
    {
        IReadOnlyList<double> p = PosteriorSummary.Checked(probabilities, open: true);
        BassPrediction prediction = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Predictive = false });
        int n = prediction.Cases, m = prediction.Locations;
        ReadOnlySpan<double> y = data.Response.Values;
        var observed = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++) observed[i, l] = y[i * m + l];
        double sd = 0.0;
        for (int d = 0; d < prediction.DrawCount; d++) sd += prediction.Weights[d] * Math.Sqrt(prediction.NoiseVariances[d]);
        IReadOnlyList<ContinuousColumn>? coordinates = data.Response is FunctionalResponse functional
            ? Array.AsReadOnly(Array.ConvertAll([.. functional.OutputVariables], v => new ContinuousColumn(v.Name, v.Values.ToArray())))
            : null;
        return Assemble(observed, prediction.Mean(), p, prediction.Quantile, IntervalSource.PredictiveMixture, sd, coordinates);
    }

    // A basis fit's: its prediction over the joint draws at the training cases, predictive draws from the seed.
    internal static BassTrainingFit Of(BassBasisFit fit, IReadOnlyList<double>? probabilities, int seed)
    {
        IReadOnlyList<double> p = PosteriorSummary.Checked(probabilities, open: true);
        BassBasisPrediction prediction = fit.Predict(PredictionCases.Of(fit.Predictors), new BassPredictOptions { Seed = seed });
        BassDecomposition decomposition = fit.Decomposition;
        int n = decomposition.Cases, m = decomposition.Locations, r = decomposition.Rank;
        double[,] residuals = decomposition.Truncation.InSample ?? decomposition.Truncation.Values;
        var observed = new double[n, m];
        var w = new double[r];
        var e = new double[m];
        var row = new double[m];
        for (int i = 0; i < n; i++)
        {
            for (int k = 0; k < r; k++) w[k] = decomposition.Coefficient(k, i);
            for (int l = 0; l < m; l++) e[l] = residuals[i, l];
            decomposition.Map.Reconstruct(w, e, row);
            for (int l = 0; l < m; l++) observed[i, l] = row[l];
        }
        return Assemble(observed, prediction.Mean(), p, prediction.Quantile, IntervalSource.PredictiveDraws, null, null);
    }

    private static BassTrainingFit Assemble(double[,] observed, double[,] mean, IReadOnlyList<double> p, Func<double, double[,]> quantile,
                                            IntervalSource source, double? noiseSd, IReadOnlyList<ContinuousColumn>? coordinates)
    {
        int n = observed.GetLength(0), m = observed.GetLength(1);
        var quantiles = new double[p.Count, n, m];
        for (int q = 0; q < p.Count; q++)
        {
            double[,] at = quantile(p[q]);
            for (int i = 0; i < n; i++)
                for (int l = 0; l < m; l++) quantiles[q, i, l] = at[i, l];
        }
        var residuals = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++) residuals[i, l] = observed[i, l] - mean[i, l];
        return new BassTrainingFit(observed, mean, p, quantiles, source, residuals, noiseSd, coordinates);
    }
}
