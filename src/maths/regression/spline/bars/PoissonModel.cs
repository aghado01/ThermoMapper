using System;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>The knot marginal a <see cref="PoissonModel"/> scores configurations by (step 3.10, OQ12).</summary>
public enum PoissonMarginal
{
    /// <summary>
    /// The default. The Laplace marginal under the unit-information g-prior β ~ N(β₀, g·(BᵀW₀B)⁻¹), taken at the
    /// posterior mode: g = m, W₀ = ȳ·I is the Fisher weight of the null fit, and β₀ projects the null fit's constant
    /// log ȳ onto the basis. It is the GLM analogue of the Normal model's Zellner g = m (Sabanés Bové and Held 2011;
    /// Li and Clyde 2018), centred on the null fit where theirs leaves the intercept flat. The prior is proper, so the
    /// mode stays finite where every count under a basis function is zero.
    /// </summary>
    GPrior,

    /// <summary>
    /// DiMatteo, Genovese and Kass's original BIC marginal ℓ(β̂) − (ν/2)·log m at the MLE, with the ridge N(0, 100·I)
    /// in the coefficient draw. It is the base of He, Yang and Kang's EBIC for non-normal responses. Over a run of
    /// zero counts its MLE lies at −∞, and IRLS stops at the η clamp (OQ12).
    /// </summary>
    Bic,
}

/// <summary>
/// Poisson observation model for count / intensity data (DMGK §6; WLS2008 — peri-stimulus spike histograms):
/// the spline is the log-intensity <c>η = Bβ</c>, <c>μ = exp(η)</c>, and β is fitted by Fisher-scoring IRLS. The knot
/// marginal is a setting (<see cref="PoissonMarginal"/>): by default the Laplace marginal under the unit-information
/// g-prior at its posterior mode, or DMGK's BIC at the MLE. Same reversible-jump engine, different codomain — the
/// spike-train / density case (Poisson-process intensity is a density up to <c>∫λ</c>). Returned up to the
/// (k,ξ)-independent <c>Σ log yᵢ!</c> constant, which cancels in every marginal ratio. (v1: non-robust —
/// <c>weights</c> are ignored.)
/// </summary>
public sealed class PoissonModel : IObservationModel
{
    private const double Tolerance = 1e-9;
    private const double EtaClamp = 30.0;
    private readonly int _maxIterations;

    /// <param name="marginal">The knot marginal: the g-prior's Laplace marginal (default) or DMGK's BIC.</param>
    public PoissonModel(PoissonMarginal marginal = PoissonMarginal.GPrior) : this(maxIterations: 100, marginal) { }

    // The IRLS iteration cap; below the default only to exercise the non-convergence rejection.
    internal PoissonModel(int maxIterations, PoissonMarginal marginal = PoissonMarginal.GPrior)
    {
        if (marginal is not (PoissonMarginal.GPrior or PoissonMarginal.Bic)) throw new ArgumentOutOfRangeException(nameof(marginal));
        _maxIterations = maxIterations;
        Marginal = marginal;
    }

    /// <summary>The knot marginal this model scores by.</summary>
    public PoissonMarginal Marginal { get; }

    // The g-prior N(β₀, Λ⁻¹): its precision Λ = Scale·BᵀB, with BᵀB in band storage and its log-determinant; the mean
    // β₀ and its fit Bβ₀; and Λβ₀, the prior's share of the IRLS right-hand side.
    private sealed record GPriorTerms(double Scale, double[,] GramBand, double LogDetGram, double[] Mean, double[] NullFit,
                                      double[] PrecisionMean);

    /// <inheritdoc/>
    public int EffectiveDimension(double[,] design, double[]? weights)
    {
        ArgumentNullException.ThrowIfNull(design);
        return design.GetLength(1);
    }

    /// <inheritdoc/>
    public double LogMarginalLikelihood(double[,] design, double[] y, double[]? weights) => Score(design, y).Value;

    // The knot marginal: BIC at the MLE, or the g-prior's Laplace marginal at the posterior mode,
    // ℓ(β̂) − ½·(β̂−β₀)ᵀΛ(β̂−β₀) + ½·log|Λ| − ½·log|H| with H = BᵀŴB + Λ (the 2π terms cancel).
    private (double Value, DesignRejection Rejection) Score(double[,] design, double[] y)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        int m = design.GetLength(0);
        int nu = design.GetLength(1);
        if (y.Length != m)
            throw new ArgumentException("Response length must match the design row count.", nameof(y));
        if (m <= nu)
            return (double.NegativeInfinity, DesignRejection.TooFewObservations);

        (double[] beta, DesignRejection rejection, GPriorTerms? prior) = FitGlm(design, y);
        if (rejection != DesignRejection.None)
            return (double.NegativeInfinity, rejection);   // inadmissible Gram, IRLS did not converge, or no counts

        if (prior is null)
        {
            double logLik = 0.0;
            for (int i = 0; i < m; i++)
            {
                double eta = Clamp(LinearPredictor(design, beta, i));
                logLik += y[i] * eta - Math.Exp(eta);
            }
            return (logLik - 0.5 * nu * Math.Log(m), DesignRejection.None);
        }

        BandCholesky h = LaplaceFactor(design, beta, prior);
        if (!DesignAdmission.Admits(h)) return (double.NegativeInfinity, DesignRejection.Inadmissible);
        double logNormalizer = 0.5 * (nu * Math.Log(prior.Scale) + prior.LogDetGram);
        return (LogPosterior(design, y, beta, prior) + logNormalizer - 0.5 * h.LogDet, DesignRejection.None);
    }

    /// <summary>
    /// The linear tempered marginal, <c>b·ℓ₁</c>: the Laplace marginal raised to b, with the knot prior untempered. Integrating
    /// β against the tempered likelihood would need a fit at every rung, and at every swap.
    /// </summary>
    public ITemperedMarginal TemperedMarginal(double[,] design, double[] y, double[]? weights)
    {
        (double value, DesignRejection rejection) = Score(design, y);
        return LinearMarginal.Of(value, rejection);
    }

    /// <inheritdoc/>
    public bool CanTemperLikelihood => true;

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
    /// <remarks>The Poisson likelihood has no noise scale, so σ² is NaN.</remarks>
    public ConditionalDraw SampleConditional(double[,] design, double[] y, double[]? weights, Xoshiro256PlusPlus rng)
        => new(SampleCoefficients(design, y, rng), double.NaN);

    /// <summary>
    /// An approximate draw from the conditional posterior <c>π(β | ξ, y) ∝ Poisson likelihood × prior</c> by
    /// independence Metropolis–Hastings (WLS2008 §5.4). Under the g-prior it proposes from <c>N(β̂, H⁻¹)</c>, with β̂ the
    /// posterior mode and <c>H = BᵀŴB + Λ</c>; under BIC the prior is the ridge N(0, 100·I) and the proposal
    /// <c>N(β̂, J⁻¹)</c>, with β̂ the MLE and <c>J = BᵀŴB</c>. The chain starts from a proposal draw and runs
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
        (double[] betaHat, _, GPriorTerms? prior) = FitGlm(design, y);
        if (prior is null && Marginal == PoissonMarginal.GPrior) return betaHat;   // no g-prior: ȳ = 0 or BᵀB refused
        BandCholesky chol = LaplaceFactor(design, betaHat, prior);
        if (!DesignAdmission.Admits(chol)) return betaHat;   // no Laplace proposal: fall back to β̂, as Whittle does

        var z = new double[nu];
        var innov = new double[nu];

        // One proposal draw, with its target and proposal log-densities.
        (double[] Beta, double LogTarget, double LogProposal) Propose()
        {
            double ss = 0.0;
            for (int j = 0; j < nu; j++) { z[j] = Normal.SampleStandard(rng); ss += z[j] * z[j]; }
            chol.SampleInnovation(z, innov);             // innov ~ N(0, H⁻¹)
            var beta = new double[nu];
            for (int j = 0; j < nu; j++) beta[j] = betaHat[j] + innov[j];
            return (beta, LogPosterior(design, y, beta, prior), -0.5 * ss);
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

    /// <summary>The Laplace variance <c>z(t)ᵀH⁻¹z(t)</c> of the log-intensity at each readout row: under the g-prior
    /// <c>H = BᵀŴB + Λ</c> at the posterior mode, under BIC <c>H = J = BᵀŴB</c> at the MLE; NaN when the fit is
    /// unusable.</summary>
    public double[]? ConditionalVariance(double[,] design, double[] y, double[]? weights, double[,] readoutDesign)
    {
        ArgumentNullException.ThrowIfNull(readoutDesign);
        (double[] betaHat, DesignRejection rejection, GPriorTerms? prior) = FitGlm(design, y);
        if (rejection != DesignRejection.None) return ReadoutQuadratics.Unavailable(readoutDesign);
        BandCholesky j = LaplaceFactor(design, betaHat, prior);
        return DesignAdmission.Admits(j)
            ? ReadoutQuadratics.Of(j, readoutDesign, 1.0)
            : ReadoutQuadratics.Unavailable(readoutDesign);
    }

    // The Laplace precision of the coefficients at β̂: BᵀŴB, plus the g-prior's Λ when there is one.
    private static BandCholesky LaplaceFactor(double[,] design, double[] betaHat, GPriorTerms? prior)
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
        if (prior is not null) AddScaled(jBand, prior.GramBand, prior.Scale);
        var chol = new BandCholesky(nu, bd.Bandwidth, BandFactorization.Ldlt, strict: true);
        chol.DecomposeBanded(jBand);
        return chol;
    }

    // Poisson log-likelihood + the coefficient prior's exponent: the g-prior's −½·(β−β₀)ᵀΛ(β−β₀) =
    // −½·Scale·‖Bβ − Bβ₀‖², or under BIC a weak Gaussian ridge (propriety; the skew the MH corrects is the likelihood's).
    private static double LogPosterior(double[,] design, double[] y, double[] beta, GPriorTerms? prior)
    {
        int m = design.GetLength(0);
        double ll = 0.0, squares = 0.0;
        for (int i = 0; i < m; i++)
        {
            double f = LinearPredictor(design, beta, i);
            double eta = Clamp(f);
            ll += y[i] * eta - Math.Exp(eta);
            if (prior is not null) squares += (f - prior.NullFit[i]) * (f - prior.NullFit[i]);
        }
        if (prior is not null) return ll - 0.5 * prior.Scale * squares;
        double ridge = 0.0;
        for (int j = 0; j < beta.Length; j++) ridge -= 0.5 * beta[j] * beta[j] / 100.0;
        return ll + ridge;
    }

    private static void AddScaled(double[,] band, double[,] addend, double scale)
    {
        for (int d = 0; d < band.GetLength(0); d++)
            for (int j = 0; j < band.GetLength(1); j++) band[d, j] += scale * addend[d, j];
    }

    // Fisher-scoring IRLS for the canonical (log-link) Poisson GLM: working response z = η + (y−μ)/μ,
    // weight w = μ, solving the banded weighted normal equations (BᵀWB + Λ) β⁺ = BᵀW z + Λβ₀ each iteration, with
    // N(β₀, Λ⁻¹) the g-prior, Λ = (ȳ/m)·BᵀB and β₀ = log ȳ·(BᵀB)⁻¹Bᵀ1, or Λ = 0 under BIC. The g-prior's iteration
    // starts from β₀, BIC's from 0. The rejection is ZeroResponse when the g-prior has no counts to scale it,
    // Inadmissible when BᵀB or a weighted Gram fails DesignAdmission, and NotConverged when the iteration cap is
    // reached before the step settles.
    private (double[] Beta, DesignRejection Rejection, GPriorTerms? Prior) FitGlm(double[,] design, double[] y)
    {
        int nu = design.GetLength(1);
        var bd = new BandedDesign(design);
        int m = bd.Rows;
        GPriorTerms? prior = null;
        if (Marginal == PoissonMarginal.GPrior)
        {
            prior = NullPrior(bd, y, out DesignRejection refused);
            if (prior is null) return (new double[nu], refused, null);
        }
        double[] beta = prior is null ? new double[nu] : (double[])prior.Mean.Clone();
        var w = new double[m];
        var z = new double[m];
        var band = new double[bd.Bandwidth + 1, nu];
        var b = new double[nu];
        var chol = new BandCholesky(nu, bd.Bandwidth, BandFactorization.Ldlt, strict: true);

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
            if (prior is not null)
            {
                AddScaled(band, prior.GramBand, prior.Scale);
                for (int p = 0; p < nu; p++) b[p] += prior.PrecisionMean[p];
            }

            chol.DecomposeBanded(band);
            if (!DesignAdmission.PivotsAdmit(chol.Status))
                return (beta, DesignRejection.Inadmissible, prior);
            double[] sol = chol.Solve(b);

            double shift = 0.0;
            for (int p = 0; p < nu; p++)
            {
                double d = sol[p] - beta[p];
                shift += d * d;
                beta[p] = sol[p];
            }
            if (shift < Tolerance)   // the final system: the full gate
                return (beta, DesignAdmission.Admits(chol) ? DesignRejection.None : DesignRejection.Inadmissible, prior);
        }
        return (beta, DesignRejection.NotConverged, prior);
    }

    // The g-prior from the null fit μ = ȳ: Scale = ȳ/g with g = m, and β₀ the projection of log ȳ·1 onto the basis
    // (β₀ = log ȳ·1 for a B-spline basis, whose functions sum to one). Null, with the reason, when ȳ = 0 or BᵀB is
    // inadmissible.
    private static GPriorTerms? NullPrior(BandedDesign bd, double[] y, out DesignRejection refused)
    {
        int m = bd.Rows, nu = bd.Dimension;
        double mean = 0.0;
        foreach (double v in y) mean += v;
        mean /= m;
        refused = DesignRejection.ZeroResponse;
        if (!(mean > 0.0)) return null;

        var ones = new double[m];
        Array.Fill(ones, 1.0);
        var gramBand = new double[bd.Bandwidth + 1, nu];
        var columnSums = new double[nu];
        bd.Accumulate(null, gramBand, ones, columnSums);
        var gram = new BandCholesky(nu, bd.Bandwidth, BandFactorization.Ldlt, strict: true);
        gram.DecomposeBanded(gramBand);
        refused = DesignRejection.Inadmissible;
        if (!DesignAdmission.Admits(gram)) return null;

        double level = Math.Log(mean), scale = mean / m;
        double[] priorMean = gram.Solve(columnSums);
        var precisionMean = new double[nu];
        for (int j = 0; j < nu; j++)
        {
            priorMean[j] *= level;
            precisionMean[j] = scale * level * columnSums[j];   // Λβ₀ = Scale·BᵀBβ₀ = Scale·log ȳ·Bᵀ1
        }
        var nullFit = new double[m];
        bd.Apply(priorMean, nullFit);
        refused = DesignRejection.None;
        return new GPriorTerms(scale, gramBand, gram.LogDet, priorMean, nullFit, precisionMean);
    }

    private static double LinearPredictor(double[,] design, double[] beta, int row)
    {
        double eta = 0.0;
        for (int j = 0; j < beta.Length; j++) eta += design[row, j] * beta[j];
        return eta;
    }

    private static double Clamp(double eta) => eta > EtaClamp ? EtaClamp : eta < -EtaClamp ? -EtaClamp : eta;
}
