using System;
using Maths.LinAlg;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// Whittle spectral model (MRA2015 ch4): the spline is the log-spectral-density <c>η = Bβ</c>, and the response
/// is the periodogram <c>I</c>, which under Whittle's approximation is exponential with mean <c>f = exp(η)</c>.
/// That is a Gamma(shape 1) / log-link GLM whose IRLS weights are constant (<c>W = 1</c>), so <c>BᵀB</c> factors
/// once and only the working response iterates. Marginal = <c>ℓ(β̂) − (ν/2) log m</c> with
/// <c>ℓ = Σ(−η − I·exp(−η))</c>. Feed (Fourier frequency scaled to [0,1], periodogram) pairs; the recovered
/// spline is the log-spectrum, its peak the dominant frequency. Up to the (k,ξ)-independent constant.
/// </summary>
public sealed class WhittleModel : IObservationModel
{
    private const int MaxIterations = 100;
    private const double Tolerance = 1e-9;
    private const double EtaClamp = 30.0;

    /// <inheritdoc/>
    public int EffectiveDimension(double[,] design, double[]? weights)
    {
        ArgumentNullException.ThrowIfNull(design);
        return design.GetLength(1);
    }

    /// <inheritdoc/>
    public double LogMarginalLikelihood(double[,] design, double[] y, double[]? weights)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        int m = design.GetLength(0);
        int nu = design.GetLength(1);
        if (y.Length != m)
            throw new ArgumentException("Periodogram length must match the design row count.", nameof(y));
        if (m <= nu)
            return double.NegativeInfinity;

        double[] beta = FitGlm(design, y);
        double logLik = 0.0;
        for (int i = 0; i < m; i++)
        {
            double eta = Clamp(LinearPredictor(design, beta, i));
            logLik += -eta - y[i] * Math.Exp(-eta);
        }
        return logLik - 0.5 * nu * Math.Log(m);
    }

    /// <inheritdoc/>
    public double[] PosteriorMeanCoefficients(double[,] design, double[] y, double[]? weights)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        return FitGlm(design, y);
    }

    // Gamma(shape 1) / log-link IRLS: weights are 1, so A = BᵀB is constant — factor once, iterate z only.
    private static double[] FitGlm(double[,] design, double[] periodogram)
    {
        int nu = design.GetLength(1);
        var bd = new BandedDesign(design);
        int m = bd.Rows;

        var band = new double[bd.Bandwidth + 1, nu];
        bd.Accumulate(null, band, null, null);   // A = BᵀB only (W = 1)
        var chol = new BandCholesky(nu, bd.Bandwidth, BandFactorization.Ldlt);
        chol.DecomposeBanded(band);               // factored once

        var beta = new double[nu];
        var z = new double[m];
        var b = new double[nu];
        for (int iter = 0; iter < MaxIterations; iter++)
        {
            for (int i = 0; i < m; i++)
            {
                double eta = Clamp(LinearPredictor(design, beta, i));
                z[i] = eta + periodogram[i] * Math.Exp(-eta) - 1.0;   // working response, W = 1
            }
            Array.Clear(b);
            bd.AccumulateRhs(null, z, b);

            double[] sol = chol.Solve(b);
            double shift = 0.0;
            for (int p = 0; p < nu; p++)
            {
                double d = sol[p] - beta[p];
                shift += d * d;
                beta[p] = sol[p];
            }
            if (shift < Tolerance) break;
        }
        return beta;
    }

    private static double LinearPredictor(double[,] design, double[] beta, int row)
    {
        double eta = 0.0;
        for (int j = 0; j < beta.Length; j++) eta += design[row, j] * beta[j];
        return eta;
    }

    private static double Clamp(double eta) => eta > EtaClamp ? EtaClamp : eta < -EtaClamp ? -EtaClamp : eta;
}
