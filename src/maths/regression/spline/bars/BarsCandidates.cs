using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Regression.GPrior;

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
    public CandidateLogDensity CandidateLogDensity(BarsResult result, double[] x, IReadOnlyList<double[]> grids, int[]? strata = null,
                                                   bool quadrature = false, bool cdf = false)
    {
        ArgumentNullException.ThrowIfNull(result);
        ArgumentNullException.ThrowIfNull(x);
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
                }
                row++;
            }
        return new CandidateLogDensity(values, grids, RetainedLatents(), strata, ChainLengths.Of(result.Draws)) { LogCdf = logCdfs };
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

    // The retained parameters' latents: knots, coefficients and σ² global and conditioned on, each observation's weight or scale integrated.
    private List<LatentReport> RetainedLatents()
    {
        var latents = new List<LatentReport>
        {
            new("knot configuration", LatentScope.Global, LatentTreatment.Conditioned),
            new("coefficients β", LatentScope.Global, LatentTreatment.Conditioned),
        };
        if (_model is WeightedNormalModel)
        {
            latents.Add(new LatentReport("noise variance σ²", LatentScope.Global, LatentTreatment.Conditioned));
            if (_resampler is not null) latents.Add(new LatentReport("robust weights wᵢ", LatentScope.Observation, LatentTreatment.Integrated));
            if (_noise is { } noise)
            {
                latents.Add(new LatentReport("latent scales vᵢ", LatentScope.Observation, LatentTreatment.Integrated));
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
