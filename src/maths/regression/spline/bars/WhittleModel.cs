using System;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

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
    private const double Tolerance = 1e-9;
    private const double EtaClamp = 30.0;
    private readonly int _maxIterations;

    public WhittleModel() : this(maxIterations: 100) { }

    // The IRLS iteration cap; below the default only to exercise the non-convergence rejection.
    internal WhittleModel(int maxIterations) => _maxIterations = maxIterations;

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

        (double[] beta, bool usable, _) = FitGlm(design, y);
        if (!usable)
            return double.NegativeInfinity;   // inadmissible BᵀB, or IRLS did not converge
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
        return FitGlm(design, y).Beta;
    }

    /// <inheritdoc/>
    public CoefficientUncertainty Uncertainty => CoefficientUncertainty.Approximate;

    /// <summary>
    /// A draw from the Laplace approximation <c>N(β̂, (BᵀB)⁻¹)</c>: the Whittle log-likelihood's Fisher information is
    /// BᵀB (unit IRLS weights). Approximate; falls back to β̂ when the fit is unusable.
    /// </summary>
    public double[] SampleCoefficients(double[,] design, double[] y, double[]? weights, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        (double[] beta, bool usable, BandCholesky factor) = FitGlm(design, y);
        if (!usable) return beta;
        int nu = beta.Length;
        var z = new double[nu];
        for (int j = 0; j < nu; j++) z[j] = Normal.SampleStandard(rng);
        var innov = new double[nu];
        factor.SampleInnovation(z, innov);
        for (int j = 0; j < nu; j++) beta[j] += innov[j];
        return beta;
    }

    /// <summary>The Laplace variance <c>z(t)ᵀ(BᵀB)⁻¹z(t)</c> of the log-spectrum at each readout row; NaN when the fit
    /// is unusable.</summary>
    public double[]? ConditionalVariance(double[,] design, double[] y, double[]? weights, double[,] readoutDesign)
    {
        ArgumentNullException.ThrowIfNull(readoutDesign);
        (_, bool usable, BandCholesky factor) = FitGlm(design, y);
        return usable ? ReadoutQuadratics.Of(factor, readoutDesign, 1.0) : ReadoutQuadratics.Unavailable(readoutDesign);
    }

    // Gamma(shape 1) / log-link IRLS: weights are 1, so A = BᵀB is constant — factor once, iterate z only. Usable is
    // false when BᵀB fails DesignAdmission or the iteration cap is reached before the step settles.
    private (double[] Beta, bool Usable, BandCholesky Factor) FitGlm(double[,] design, double[] periodogram)
    {
        int nu = design.GetLength(1);
        var bd = new BandedDesign(design);
        int m = bd.Rows;

        var band = new double[bd.Bandwidth + 1, nu];
        bd.Accumulate(null, band, null, null);   // A = BᵀB only (W = 1)
        var chol = new BandCholesky(nu, bd.Bandwidth, BandFactorization.Ldlt);
        chol.DecomposeBanded(band);               // factored once

        var beta = new double[nu];
        if (!DesignAdmission.Admits(chol.Status))
            return (beta, false, chol);
        var z = new double[m];
        var b = new double[nu];
        for (int iter = 0; iter < _maxIterations; iter++)
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
            if (shift < Tolerance) return (beta, true, chol);
        }
        return (beta, false, chol);
    }

    private static double LinearPredictor(double[,] design, double[] beta, int row)
    {
        double eta = 0.0;
        for (int j = 0; j < beta.Length; j++) eta += design[row, j] * beta[j];
        return eta;
    }

    private static double Clamp(double eta) => eta > EtaClamp ? EtaClamp : eta < -EtaClamp ? -EtaClamp : eta;
}
