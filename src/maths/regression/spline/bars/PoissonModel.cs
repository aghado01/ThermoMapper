using System;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// Poisson observation model for count / intensity data (DMGK §6; WLS2008 — peri-stimulus spike histograms):
/// the spline is the log-intensity <c>η = Bβ</c>, <c>μ = exp(η)</c>, and the marginal is the BIC/Laplace
/// approximation <c>log p(y|Z) ≈ ℓ(β̂) − (ν/2) log m</c> with β̂ the Poisson-GLM MLE (Fisher-scoring IRLS). Same
/// reversible-jump engine, different codomain — the spike-train / density case (Poisson-process intensity is a
/// density up to <c>∫λ</c>). Returned up to the (k,ξ)-independent <c>Σ log yᵢ!</c> constant, which cancels in
/// every marginal ratio. (v1: non-robust — <c>weights</c> are ignored.)
/// </summary>
public sealed class PoissonModel : IObservationModel
{
    private const double Tolerance = 1e-9;
    private const double EtaClamp = 30.0;
    private readonly int _maxIterations;

    public PoissonModel() : this(maxIterations: 100) { }

    // The IRLS iteration cap; below the default only to exercise the non-convergence rejection.
    internal PoissonModel(int maxIterations) => _maxIterations = maxIterations;

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
            throw new ArgumentException("Response length must match the design row count.", nameof(y));
        if (m <= nu)
            return double.NegativeInfinity;

        (double[] beta, bool usable) = FitGlm(design, y);
        if (!usable)
            return double.NegativeInfinity;   // inadmissible weighted Gram, or IRLS did not converge

        double logLik = 0.0;
        for (int i = 0; i < m; i++)
        {
            double eta = Clamp(LinearPredictor(design, beta, i));
            logLik += y[i] * eta - Math.Exp(eta);
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

    /// <inheritdoc/>
    public double[] SampleCoefficients(double[,] design, double[] y, double[]? weights, Xoshiro256PlusPlus rng)
        => SampleCoefficients(design, y, rng);

    /// <summary>
    /// An approximate draw from the conditional posterior <c>π(β | ξ, y) ∝ Poisson likelihood × N(0, 100·I)</c> by
    /// independence Metropolis–Hastings, proposing from the Laplace approximation <c>N(β̂, J⁻¹)</c> with
    /// <c>J = BᵀŴB</c> at the MLE (WLS2008 §5.4). The chain starts from a proposal draw and runs
    /// <paramref name="mhIterations"/> accept/reject steps. It approaches the posterior as the steps grow (geometrically
    /// when π/q is bounded), so any finite count is an approximation; the steps matter most for skewed, low-count
    /// posteriors. The proposal exponent is <c>−½‖z‖²</c> by construction (<c>innov = L⁻ᵀD^{−½}z</c>), so no quadratic
    /// form is needed.
    /// </summary>
    public double[] SampleCoefficients(double[,] design, double[] y, Xoshiro256PlusPlus rng, int mhIterations = 5)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(rng);
        if (mhIterations < 1) throw new ArgumentOutOfRangeException(nameof(mhIterations));

        int nu = design.GetLength(1);
        double[] betaHat = FitGlm(design, y).Beta;
        BandCholesky chol = LaplaceFactor(design, betaHat);

        var z = new double[nu];
        var innov = new double[nu];

        // One proposal draw, with its target and proposal log-densities.
        (double[] Beta, double LogTarget, double LogProposal) Propose()
        {
            double ss = 0.0;
            for (int j = 0; j < nu; j++) { z[j] = Normal.SampleStandard(rng); ss += z[j] * z[j]; }
            chol.SampleInnovation(z, innov);             // innov ~ N(0, J⁻¹)
            var beta = new double[nu];
            for (int j = 0; j < nu; j++) beta[j] = betaHat[j] + innov[j];
            return (beta, LogPosterior(design, y, beta), -0.5 * ss);
        }

        (double[] cur, double logTargetCur, double logPropCur) = Propose();
        for (int it = 0; it < mhIterations; it++)
        {
            (double[] cand, double logTargetCand, double logPropCand) = Propose();
            double r = (logTargetCand - logTargetCur) + (logPropCur - logPropCand);   // independence-MH ratio
            if (r >= 0.0 || Math.Log(1.0 - rng.NextDouble()) < r)
            {
                cur = cand;
                logTargetCur = logTargetCand;
                logPropCur = logPropCand;
            }
        }
        return cur;
    }

    /// <summary>The Laplace variance <c>z(t)ᵀJ⁻¹z(t)</c> of the log-intensity at each readout row, with
    /// <c>J = BᵀŴB</c> at the MLE; NaN when the fit is unusable.</summary>
    public double[]? ConditionalVariance(double[,] design, double[] y, double[]? weights, double[,] readoutDesign)
    {
        ArgumentNullException.ThrowIfNull(readoutDesign);
        (double[] betaHat, bool usable) = FitGlm(design, y);
        if (!usable) return ReadoutQuadratics.Unavailable(readoutDesign);
        BandCholesky j = LaplaceFactor(design, betaHat);
        return DesignAdmission.Admits(j.Status)
            ? ReadoutQuadratics.Of(j, readoutDesign, 1.0)
            : ReadoutQuadratics.Unavailable(readoutDesign);
    }

    // J = BᵀŴB at β̂, the Laplace precision of the coefficients.
    private static BandCholesky LaplaceFactor(double[,] design, double[] betaHat)
    {
        int nu = design.GetLength(1);
        var bd = new BandedDesign(design);
        int m = bd.Rows;
        var w = new double[m];
        for (int i = 0; i < m; i++)
        {
            double mu = Math.Exp(Clamp(LinearPredictor(design, betaHat, i)));
            w[i] = mu < 1e-10 ? 1e-10 : mu;
        }
        var jBand = new double[bd.Bandwidth + 1, nu];
        bd.Accumulate(w, jBand, null, null);
        var chol = new BandCholesky(nu, bd.Bandwidth, BandFactorization.Ldlt);
        chol.DecomposeBanded(jBand);
        return chol;
    }

    // Poisson log-likelihood + a weak Gaussian ridge (propriety; the skew the MH corrects is the likelihood's).
    private static double LogPosterior(double[,] design, double[] y, double[] beta)
    {
        int m = design.GetLength(0);
        double ll = 0.0;
        for (int i = 0; i < m; i++)
        {
            double eta = Clamp(LinearPredictor(design, beta, i));
            ll += y[i] * eta - Math.Exp(eta);
        }
        double prior = 0.0;
        for (int j = 0; j < beta.Length; j++) prior -= 0.5 * beta[j] * beta[j] / 100.0;
        return ll + prior;
    }

    // Fisher-scoring IRLS for the canonical (log-link) Poisson GLM: working response z = η + (y−μ)/μ,
    // weight w = μ, solving the banded weighted normal equations (BᵀWB) β⁺ = BᵀW z each iteration. Usable is false
    // when a weighted Gram fails DesignAdmission or the iteration cap is reached before the step settles.
    private (double[] Beta, bool Usable) FitGlm(double[,] design, double[] y)
    {
        int nu = design.GetLength(1);
        var bd = new BandedDesign(design);
        int m = bd.Rows;
        var beta = new double[nu];
        var w = new double[m];
        var z = new double[m];
        var band = new double[bd.Bandwidth + 1, nu];
        var b = new double[nu];
        var chol = new BandCholesky(nu, bd.Bandwidth, BandFactorization.Ldlt);

        for (int iter = 0; iter < _maxIterations; iter++)
        {
            for (int i = 0; i < m; i++)
            {
                double eta = Clamp(LinearPredictor(design, beta, i));
                double mu = Math.Exp(eta);
                if (mu < 1e-10) mu = 1e-10;
                w[i] = mu;
                z[i] = eta + (y[i] - mu) / mu;   // working response
            }
            Array.Clear(band);
            Array.Clear(b);
            bd.Accumulate(w, band, z, b);

            chol.DecomposeBanded(band);
            if (!DesignAdmission.Admits(chol.Status))
                return (beta, false);
            double[] sol = chol.Solve(b);

            double shift = 0.0;
            for (int p = 0; p < nu; p++)
            {
                double d = sol[p] - beta[p];
                shift += d * d;
                beta[p] = sol[p];
            }
            if (shift < Tolerance) return (beta, true);
        }
        return (beta, false);
    }

    private static double LinearPredictor(double[,] design, double[] beta, int row)
    {
        double eta = 0.0;
        for (int j = 0; j < beta.Length; j++) eta += design[row, j] * beta[j];
        return eta;
    }

    private static double Clamp(double eta) => eta > EtaClamp ? EtaClamp : eta < -EtaClamp ? -EtaClamp : eta;
}
