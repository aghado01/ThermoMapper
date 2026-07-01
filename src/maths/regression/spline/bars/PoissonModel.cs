using System;
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
            throw new ArgumentException("Response length must match the design row count.", nameof(y));
        if (m <= nu)
            return double.NegativeInfinity;

        double[] beta = FitGlm(design, y);

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
        return FitGlm(design, y);
    }

    /// <inheritdoc/>
    public double[] SampleCoefficients(double[,] design, double[] y, double[]? weights, Xoshiro256PlusPlus rng)
        => SampleCoefficients(design, y, rng);

    /// <summary>
    /// Draws a coefficient vector from the (possibly skewed) conditional posterior <c>π(β | ξ, y) ∝ Poisson-lik ×
    /// weak prior</c> by independence Metropolis–Hastings, proposing from the Laplace/normal approximation
    /// <c>N(β̂, (XᵀŴX)⁻¹)</c> at the MLE (WLS2008 §5.4). For large Poisson means the normal approximation is good
    /// and the first draw is kept (governed by <paramref name="suspectThreshold"/> — accept early when the log MH
    /// ratio exceeds it); for small means the posterior is skewed and the MH steps correct toward the true
    /// likelihood — the refinement that matters for low-count neuronal data. (<see cref="PosteriorMeanCoefficients"/>
    /// is the cheaper Laplace-mean readout.) The proposal exponent is <c>−½‖z‖²</c> by construction
    /// (<c>innov = L⁻ᵀD^{−½}z</c>), so no quadratic form is needed.
    /// </summary>
    public double[] SampleCoefficients(double[,] design, double[] y, Xoshiro256PlusPlus rng,
                                       int mhIterations = 5, double suspectThreshold = -10.0)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(rng);
        if (mhIterations < 1) throw new ArgumentOutOfRangeException(nameof(mhIterations));

        int nu = design.GetLength(1);
        double[] betaHat = FitGlm(design, y);

        var bd = new BandedDesign(design);
        int m = bd.Rows;
        var w = new double[m];
        for (int i = 0; i < m; i++)
        {
            double mu = Math.Exp(Clamp(LinearPredictor(design, betaHat, i)));
            w[i] = mu < 1e-10 ? 1e-10 : mu;
        }
        var jBand = new double[bd.Bandwidth + 1, nu];
        bd.Accumulate(w, jBand, null, null);              // J = XᵀŴX at the MLE — the normal-approx precision
        var chol = new BandCholesky(nu, bd.Bandwidth, BandFactorization.Ldlt);
        chol.DecomposeBanded(jBand);

        var z = new double[nu];
        var innov = new double[nu];
        var cand = new double[nu];
        var cur = (double[])betaHat.Clone();
        double logTargetCur = LogPosterior(design, y, cur);
        double logPropCur = 0.0;                          // cur = β̂ ⇒ proposal exponent 0

        for (int it = 0; it < mhIterations; it++)
        {
            double ss = 0.0;
            for (int j = 0; j < nu; j++) { z[j] = StandardNormal(rng); ss += z[j] * z[j]; }
            chol.SampleInnovation(z, innov);             // innov ~ N(0, J⁻¹)
            for (int j = 0; j < nu; j++) cand[j] = betaHat[j] + innov[j];

            double logTargetCand = LogPosterior(design, y, cand);
            double logPropCand = -0.5 * ss;
            double r = (logTargetCand - logTargetCur) + (logPropCur - logPropCand);   // independence-MH ratio

            if (it == 0 && r > suspectThreshold)
                return cand;                             // normal approx adequate — keep the first draw

            if (r >= 0.0 || Math.Log(1.0 - rng.NextDouble()) < r)
            {
                Array.Copy(cand, cur, nu);
                logTargetCur = logTargetCand;
                logPropCur = logPropCand;
            }
        }
        return cur;
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

    private static double StandardNormal(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    // Fisher-scoring IRLS for the canonical (log-link) Poisson GLM: working response z = η + (y−μ)/μ,
    // weight w = μ, solving the banded weighted normal equations (BᵀWB) β⁺ = BᵀW z each iteration.
    private static double[] FitGlm(double[,] design, double[] y)
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

        for (int iter = 0; iter < MaxIterations; iter++)
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
