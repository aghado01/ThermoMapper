using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Rng;

namespace Maths.Regression.Spline.Bars;

public sealed partial class BarsEnsemble
{
    /// <summary>
    /// BARS's candidate density (E16's seam, the conformal layer's slice 5): log p(y | x_t, draw) for every retained draw and every value y
    /// on each new case's grid, beside <see cref="PointwiseLogDensity"/>'s retained-parameters matrix over the same draws in the same rows,
    /// with its chain lengths and latent report, so the conformal layer (<see cref="Conformal.ConformalBayes"/>) reads both. Each draw's
    /// curve at the case is its knots' design at x_t times its coefficients, and the density is the training terms' own, each observation's
    /// latent integrated as there: the Gaussian, a robust weight's marginal, a noise's law at the draw's σ² and skew, the Poisson mass or the
    /// Whittle ordinate's exponential. A value outside the response's support, a count that is not a whole number or a negative ordinate,
    /// takes −∞. The Gibbs-loss model is refused, as its pointwise density is: its likelihood is a pseudo-likelihood. On
    /// <paramref name="cdf"/> each draw's log CDF rides beside the density, from the same law; a count's CDF steps, so the Poisson model's
    /// is refused, its PIT not uniform without a randomized tie, which the predictive distribution does not take.
    /// </summary>
    /// <param name="result">A run of this ensemble, whose retained draws are read.</param>
    /// <param name="x">The new cases' inputs, each in [0, 1] as the training inputs are.</param>
    /// <param name="grids">Each case's candidate values, ascending.</param>
    /// <param name="strata">Each case's stratum in the training density's partition; null puts every case in stratum 0.</param>
    /// <param name="quadrature">Takes a noise whose marginal is an integral over v, draws × candidates integrals; refused without it, as
    /// <see cref="PointwiseLogDensity"/> refuses it.</param>
    /// <param name="cdf">Adds each draw's log CDF at each candidate, for the conformal predictive distribution.</param>
    /// <param name="observationLatents">CB-4's reading of each observation's latent, the robust weight or the noise's scale: integrated, the
    /// default; or <see cref="LatentTreatment.Conditioned"/>, the latent-posterior reading (the conformal layer's slice 12b), beside the
    /// training density's conditioned content. Under it each candidate's latent is drawn per draw from its exact conditional given the
    /// candidate's value, p(v′ | θ_s, y), on a stream from <paramref name="seed"/>, and ln N(y | η + m·v′, σ²/w′) is
    /// <see cref="Maths.Regression.CandidateLogDensity.ConditionedValues"/>; the values stay the integrated terms, the add-one-in weights. A
    /// noise without an exact conditional draw (the generalized beta prime's slice step) is refused, and so is the CDF.</param>
    /// <param name="seed">The candidates' latent draws' seed under the conditioned reading.</param>
    public CandidateLogDensity CandidateLogDensity(BarsResult result, double[] x, IReadOnlyList<double[]> grids, int[]? strata = null,
                                                   bool quadrature = false, bool cdf = false,
                                                   LatentTreatment observationLatents = LatentTreatment.Integrated, int seed = 0)
    {
        ArgumentNullException.ThrowIfNull(result);
        ArgumentNullException.ThrowIfNull(x);
        bool conditioned = CheckObservationLatents(observationLatents, cdf);
        CheckQuadrature(quadrature);
        Func<double, double, double, double> term = RetainedTerm();
        if (cdf) CheckLogCdf();
        int cases = x.Length;
        int[] offsets = CandidateSeam.Offsets(grids, cases);
        Func<double, bool> inSupport = _model switch
        {
            PoissonModel => yi => yi >= 0.0 && yi == Math.Floor(yi),
            WhittleModel => yi => yi >= 0.0,
            _ => _ => true,
        };

        int draws = 0;
        foreach (IReadOnlyList<BarsDraw> chain in result.Draws) draws += chain.Count;
        var values = new double[draws, offsets[cases]];
        double[,]? logCdfs = cdf ? new double[draws, offsets[cases]] : null;
        Func<BarsDraw, double[], Xoshiro256PlusPlus, double[]>? drawLatents = conditioned ? LatentDraws() : null;
        double[,]? conditionedValues = drawLatents is null ? null : new double[draws, offsets[cases]];
        Xoshiro256PlusPlus? rng = drawLatents is null ? null : new Xoshiro256PlusPlus(seed);
        int row = 0;
        foreach (IReadOnlyList<BarsDraw> chain in result.Draws)
            foreach (BarsDraw d in chain)
            {
                double[,] z = _basis.Design(d.Knots, x);
                Func<double, double, double>? logCdf = cdf ? LogCdfAt(d) : null;
                for (int t = 0; t < cases; t++)
                {
                    double eta = Linear(z, t, d.Coefficients);
                    double[] grid = grids[t];
                    for (int g = 0; g < grid.Length; g++)
                    {
                        values[row, offsets[t] + g] = inSupport(grid[g]) ? Term(term, grid[g], eta, d) : double.NegativeInfinity;
                        if (logCdf is not null) logCdfs![row, offsets[t] + g] = logCdf(grid[g], eta);
                    }
                    if (drawLatents is null) continue;
                    var residuals = new double[grid.Length];
                    for (int g = 0; g < grid.Length; g++) residuals[g] = grid[g] - eta;
                    double[] weights = drawLatents(d, residuals, rng!);
                    for (int g = 0; g < grid.Length; g++) conditionedValues![row, offsets[t] + g] = ConditionedTerm(d, residuals[g], weights[g]);
                }
                row++;
            }
        return new CandidateLogDensity(values, grids, RetainedLatents(observationLatents), strata, ChainLengths.Of(result.Draws))
        {
            LogCdf = logCdfs,
            ConditionedValues = conditionedValues,
        };
    }

    // Whether the model holds an observation latent: a robust weight, or a noise's scale.
    private bool HasObservationLatents => _model is WeightedNormalModel && (_resampler is not null || _noise is not null);

    // Integrated or conditioned, the latter refusing the CDF, which the predictive distribution reads with the latents integrated; true when
    // conditioned and the model holds a latent to condition on.
    private bool CheckObservationLatents(LatentTreatment observationLatents, bool cdf)
    {
        if (observationLatents is not (LatentTreatment.Integrated or LatentTreatment.Conditioned))
            throw new ArgumentOutOfRangeException(nameof(observationLatents), "An observation's latent is read integrated or conditioned.");
        if (observationLatents == LatentTreatment.Integrated || !HasObservationLatents) return false;
        if (cdf) throw new ArgumentException("The CDF is read with each observation's latent integrated, as the predictive distribution reads it.", nameof(cdf));
        return true;
    }

    // ln f(e | θ, the observation's latent) at a draw: N(e; m·v, σ²/w), with w the latent's precision weight, the robust weight or 1/(c·v)
    // under a noise, and m·v the noise's skew shift, m = √σ²·β. The training terms and the candidates' conditioned terms share it.
    internal double ConditionedTerm(BarsDraw d, double residual, double weight)
    {
        double shift = _noise is { } noise ? Math.Sqrt(d.Sigma2) * d.Skew / (noise.VarianceScale * weight) : 0.0;
        double variance = d.Sigma2 / weight, e = residual - shift;
        return -0.5 * Math.Log(2.0 * Math.PI * variance) - 0.5 * e * e / variance;
    }

    // Each candidate's latent drawn from its exact conditional given its residual at a draw, as a precision weight: the robust weights'
    // resampler on the grid's residuals at once, or the noise's mixing law, v′ given r with s² = σ²c and m = √σ²·β, read as w′ = 1/(c·v′). A
    // law whose conditional is a slice step from the current value has no exact draw, and is refused.
    internal Func<BarsDraw, double[], Xoshiro256PlusPlus, double[]> LatentDraws()
    {
        if (_noise is not { } noise) return (d, residuals, rng) => _resampler!.Resample(residuals, d.Sigma2, rng);
        Func<double, double, double, Xoshiro256PlusPlus, double> draw = noise.Mixing switch
        {
            StudentTMixing t => (r, s2, m, rng) => t.SampleConditional(r, s2, m, 1.0, rng),
            ExponentialMixing e => (r, s2, m, rng) => e.SampleConditional(r, s2, m, 1.0, rng),
            GigMixing g => (r, s2, m, rng) => g.SampleConditional(r, s2, m, 1.0, rng),
            _ => throw new NotSupportedException($"{noise.Name}'s latent has no exact draw from its conditional, only a step from its current value, "
                                                 + "so the latent-posterior reading cannot draw a candidate's latent."),
        };
        double c = noise.VarianceScale;
        return (d, residuals, rng) =>
        {
            double s2 = d.Sigma2 * c, m = Math.Sqrt(d.Sigma2) * d.Skew;
            var weights = new double[residuals.Length];
            for (int g = 0; g < residuals.Length; g++) weights[g] = 1.0 / (c * draw(residuals[g], s2, m, rng));
            return weights;
        };
    }

    // ηᵢ = zᵢβ, summed in column order: the training terms' and the candidates' one arithmetic.
    private static double Linear(double[,] z, int i, double[] beta)
    {
        double eta = 0.0;
        for (int j = 0; j < beta.Length; j++) eta += z[i, j] * beta[j];
        return eta;
    }

    // Under a noise the law's density with vᵢ integrated, at the draw's σ² and skew; otherwise the model's term.
    private double Term(Func<double, double, double, double> term, double yi, double eta, BarsDraw d) => _noise is { } noise
        ? noise.Mixing.MarginalLogDensity(yi - eta, d.Sigma2 * noise.VarianceScale, Math.Sqrt(d.Sigma2) * d.Skew)
        : term(yi, eta, d.Sigma2);

    // (yᵢ, ηᵢ = zᵢβ, σ²) → the retained parameters' log term, a robust weight integrated under its prior.
    private Func<double, double, double, double> RetainedTerm() => _model switch
    {
        WeightedNormalModel when _resampler is IWeightResampler resampler
            => (yi, eta, s2) => resampler.MarginalLogDensity(yi - eta, s2),
        WeightedNormalModel => (yi, eta, s2) => -0.5 * Math.Log(2.0 * Math.PI * s2) - 0.5 * (yi - eta) * (yi - eta) / s2,
        PoissonModel => (yi, eta, _) => yi * eta - Math.Exp(eta) - SpecialFunctions.LogGamma(yi + 1.0),
        WhittleModel => (yi, eta, _) => -eta - yi * Math.Exp(-eta),
        _ => throw new NotSupportedException(
            $"{_model.GetType().Name} has no pointwise density: its likelihood is a pseudo-likelihood, not a density of a new observation."),
    };

    // The retained parameters' latents: knots, coefficients and σ² global and conditioned on, each observation's weight or scale integrated, or
    // conditioned under the latent-posterior reading.
    private List<LatentReport> RetainedLatents(LatentTreatment observation = LatentTreatment.Integrated)
    {
        var latents = new List<LatentReport>
        {
            new("knot configuration", LatentScope.Global, LatentTreatment.Conditioned),
            new("coefficients β", LatentScope.Global, LatentTreatment.Conditioned),
        };
        if (_model is WeightedNormalModel)
        {
            latents.Add(new LatentReport("noise variance σ²", LatentScope.Global, LatentTreatment.Conditioned));
            if (_resampler is not null) latents.Add(new LatentReport("robust weights wᵢ", LatentScope.Observation, observation));
            if (_noise is { } noise)
            {
                latents.Add(new LatentReport("latent scales vᵢ", LatentScope.Observation, observation));
                if (noise.SkewPrior is not null) latents.Add(new LatentReport("skew β", LatentScope.Global, LatentTreatment.Conditioned));
            }
        }
        return latents;
    }

    private void CheckQuadrature(bool quadrature)
    {
        if (_noise is { } law && law.Mixing.MarginalByQuadrature && !quadrature)
            throw new NotSupportedException($"{law.Name}'s marginal is an integral over v per observation and draw; ask for it with quadrature: true.");
    }

    private void CheckLogCdf()
    {
        if (_model is PoissonModel)
            throw new NotSupportedException("A count's CDF steps, so its PIT is not uniform without a randomized tie, which the conformal predictive distribution does not take.");
        if (_model is not (WeightedNormalModel or WhittleModel))
            throw new NotSupportedException($"{_model.GetType().Name} has no density of a new observation, so no CDF.");
    }

    // (yᵢ, ηᵢ) → ln F(yᵢ | ηᵢ, draw), each tail kept, from the term's own law: a noise's or a robust weight's marginal, the Gaussian, or the
    // Whittle ordinate's exponential with mean e^η, ln(1 − exp(−y·e^(−η))).
    private Func<double, double, double> LogCdfAt(BarsDraw d)
    {
        if (_model is WhittleModel)
            return (yi, eta) => yi > 0.0 ? LogFunctions.Log1mExp(-yi * Math.Exp(-eta)) : double.NegativeInfinity;
        if (_noise is { } noise)
        {
            Func<double, double> residual = noise.LogCdfAt(d.Sigma2, d.Skew);
            return (yi, eta) => residual(yi - eta);
        }
        double s2 = d.Sigma2;
        if (_resampler is { } resampler) return (yi, eta) => resampler.MarginalLogCdf(yi - eta, s2);
        double sigma = Math.Sqrt(s2);
        return (yi, eta) => Normal.StandardLogCdf((yi - eta) / sigma);
    }
}
