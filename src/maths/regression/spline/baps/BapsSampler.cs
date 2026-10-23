using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;
using Maths.Samplers.Ensemble;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Baps;

/// <summary>How the smoothing parameter λ is sampled in <see cref="BapsSampler"/>.</summary>
public enum BapsLambdaUpdate
{
    /// <summary>
    /// Conjugate Gibbs over the variance components: <c>β | σ²,τ²</c> from its banded Gaussian full conditional,
    /// <c>σ² | β</c> and <c>τ² | β</c> from inverse-gamma full conditionals (λ = σ²/τ² derived). Full Bayesian.
    /// </summary>
    Gibbs,

    /// <summary>
    /// Metropolis on <c>log λ</c> against the profiled REML marginal evidence (β integrated out) under a uniform prior
    /// on the effective degrees of freedom (<see cref="PenalizedSpline.EdfSlope"/>); β is Rao-Blackwellized to its
    /// posterior mean at each draw. Marginalizes the spline. The evidence levels off as λ → ∞, so a flat prior on
    /// <c>log λ</c> would be improper; the edf prior decays like 1/λ and makes the target proper.
    /// </summary>
    MarginalEvidence,

    /// <summary>
    /// Locally-adaptive Gibbs (the "A" in BAPS): as <see cref="Gibbs"/> but each difference carries its own mean-1
    /// multiplier <c>η_i ~ Gamma(φ, φ)</c>, so the smoothing <c>λ_i = λ_global·η_i</c> varies across the domain
    /// — heavy where the fit is flat, light where it must track a feature. Every conditional stays conjugate
    /// (the prior normalizer factorizes as <c>(Πη_i)|DDᵀ|</c>, giving each <c>η_i</c> a Gamma full conditional).
    /// </summary>
    AdaptiveGibbs
}

/// <summary>Posterior summary of a BAPS run: the Bayes fit coefficients and the smoothing / noise posteriors.</summary>
/// <remarks>
/// <paramref name="LambdaMedian"/>/<paramref name="LambdaLo"/>/<paramref name="LambdaHi"/> are the 50%, 2.5% and 97.5%
/// quantiles of the pooled draws of the global smoothing (in the adaptive mode, λ_global). A median, not a mean,
/// because λ can be heavy-tailed enough that its posterior mean does not exist. <paramref name="LocalSmoothing"/> is
/// the posterior-mean local smoothing field λ_i (length ν−r) — non-null only for <see cref="BapsLambdaUpdate.AdaptiveGibbs"/>; its dips
/// mark where the fit was allowed to track features. <paramref name="RHatLogLambda"/> is R̂ on log λ across chains
/// with its <see cref="DiagnosticState"/> (unavailable for a single chain). <paramref name="PosteriorDraws"/> keeps
/// every retained draw (<see cref="BapsDraw"/>), per chain in draw order; every summary above is read from it.
/// </remarks>
public sealed record BapsResult(
    double[] Coefficients,
    double LambdaMedian,
    double LambdaLo,
    double LambdaHi,
    double SigmaMean,
    RHatDiagnostic RHatLogLambda,
    int Draws,
    IReadOnlyList<IReadOnlyList<BapsDraw>> PosteriorDraws,
    double[]? LocalSmoothing = null);

/// <summary>
/// One retained BAPS draw: the coefficients, the global smoothing λ, the noise variance σ², the penalty variance τ²,
/// and the adaptive mode's local multipliers η.
/// <list type="bullet">
///   <item>In the Gibbs modes λ = σ²/τ², all from the end of the sweep.</item>
///   <item>The marginal-evidence mode samples λ alone. Its coefficients are the conditional mean and σ² the profiled
///   variance at λ, not draws, and τ² is NaN.</item>
///   <item><see cref="Eta"/> is set in the adaptive mode only, where the local smoothing is λᵢ = λ·ηᵢ.</item>
/// </list>
/// </summary>
public sealed record BapsDraw(double[] Coefficients, double Lambda, double Sigma2, double Tau2, double[]? Eta);

/// <summary>
/// The BAPS (Bayesian adaptive P-spline) sampler — the measure-side dual of free-knot BARS. A rich fixed B-spline
/// basis with a difference penalty whose strength λ is inferred, sampled either by conjugate Gibbs over the
/// variance components or by Metropolis against the REML marginal evidence (<see cref="BapsLambdaUpdate"/>). Every
/// inner solve is banded (<c>A = ZᵀZ + λP</c> via <see cref="BandCholesky"/>'s root-free LDLᵀ, whose factor also
/// supplies the Gaussian β-draw); chains run in parallel from one seed tree and retain every draw, from which the Bayes
/// fit, the λ posterior (median + 95% interval), σ, and R̂ on log λ are read.
/// </summary>
public sealed class BapsSampler
{
    private readonly double[,] _z;
    private readonly double[] _y;
    private readonly DifferencePenalty _penalty;
    private readonly BandedDesign _design;
    private readonly PenalizedSpline _model;
    private readonly int _n;
    private readonly int _nu;
    private readonly int _bw;
    private readonly int _r;
    private readonly double[,] _gram0;   // ZᵀZ band — constant across draws
    private readonly double[] _zty;      // Zᵀy — constant
    private readonly double _priorShape; // inverse-gamma hyperprior on each variance component
    private readonly double _priorScale;
    private readonly double _mhLogStep;  // Metropolis log-λ proposal sd (marginal-evidence mode)
    private readonly double _phi;        // adaptivity: η_i ~ Gamma(φ, φ); φ→∞ recovers global smoothing

    public BapsSampler(double[,] design, double[] y, DifferencePenalty penalty,
                       double priorShape = 1e-3, double priorScale = 1e-3, double mhLogStep = 0.6, double phi = 1.0)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(penalty);
        _z = design;
        _y = y;
        _penalty = penalty;
        _design = new BandedDesign(design);
        _model = new PenalizedSpline(design, y, penalty);
        _n = _design.Rows;
        _nu = _design.Dimension;
        _r = penalty.Order;
        if (y.Length != _n)
            throw new ArgumentException("Response length must match the design row count.", nameof(y));
        _bw = Math.Max(_design.Bandwidth, penalty.Order);
        _gram0 = new double[_bw + 1, _nu];
        _zty = new double[_nu];
        _design.Accumulate(null, _gram0, y, _zty);   // ZᵀZ and Zᵀy, once
        _priorShape = priorShape;
        _priorScale = priorScale;
        _mhLogStep = mhLogStep;
        _phi = phi;
    }

    public BapsResult Run(BapsLambdaUpdate mode, int chains = 4, int masterSeed = 1, int burn = 500, int samples = 1500)
    {
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));
        if (samples < 2 || burn < 0) throw new ArgumentOutOfRangeException(nameof(samples));

        int[] seeds = SeedTree.Derive(masterSeed, chains);
        var draws = new List<BapsDraw>[chains];
        Parallel.For(0, chains, c =>
        {
            var rng = new Xoshiro256PlusPlus(seeds[c]);
            draws[c] = mode switch
            {
                BapsLambdaUpdate.Gibbs => RunGibbs(rng, burn, samples),
                BapsLambdaUpdate.MarginalEvidence => RunMarginal(rng, burn, samples),
                BapsLambdaUpdate.AdaptiveGibbs => RunAdaptiveGibbs(rng, burn, samples),
                _ => throw new ArgumentOutOfRangeException(nameof(mode)),
            };
        });

        // Every summary is read from the retained draws: each chain's sums in draw order, then over chains; the log-λ
        // moments per chain for R̂; and the pooled λ sample for its quantiles.
        double total = chains * (double)samples;
        var coef = new double[_nu];
        var logLam = new RunningMoments[chains];
        double sigSum = 0.0;
        var allLam = new double[(int)total];
        int nLoc = _nu - _r;
        double[]? localField = mode == BapsLambdaUpdate.AdaptiveGibbs && nLoc > 0 ? new double[nLoc] : null;
        var posteriorDraws = new IReadOnlyList<BapsDraw>[chains];
        int idx = 0;
        for (int c = 0; c < chains; c++)
        {
            posteriorDraws[c] = draws[c].ToArray();
            var coefSum = new double[_nu];
            var moments = new RunningMoments();
            double sigmaSum = 0.0;
            double[]? localSum = localField is null ? null : new double[nLoc];
            foreach (BapsDraw d in draws[c])
            {
                for (int p = 0; p < _nu; p++) coefSum[p] += d.Coefficients[p];
                allLam[idx++] = d.Lambda;
                moments.Add(Math.Log(d.Lambda));
                sigmaSum += Math.Sqrt(d.Sigma2);
                if (localSum is not null)
                    for (int k = 0; k < nLoc; k++) localSum[k] += d.Lambda * d.Eta![k];
            }
            for (int p = 0; p < _nu; p++) coef[p] += coefSum[p];
            logLam[c] = moments;
            sigSum += sigmaSum;
            if (localSum is not null)
                for (int k = 0; k < nLoc; k++) localField![k] += localSum[k];
        }
        for (int p = 0; p < _nu; p++) coef[p] /= total;
        if (localField is not null)
            for (int k = 0; k < nLoc; k++) localField[k] /= total;

        Array.Sort(allLam);

        return new BapsResult(
            coef, Quantile(allLam, 0.5), Quantile(allLam, 0.025), Quantile(allLam, 0.975),
            sigSum / total, ChainDiagnostics.RHat(logLam), (int)total, posteriorDraws, localField);
    }

    // ── Gibbs: β | σ²,τ²  then  σ² | β  then  τ² | β ───────────────────────────
    private List<BapsDraw> RunGibbs(Xoshiro256PlusPlus rng, int burn, int samples)
    {
        var draws = new List<BapsDraw>(samples);
        var band = new double[_bw + 1, _nu];
        var z = new double[_nu];
        var innov = new double[_nu];
        var beta = new double[_nu];
        var chol = new BandCholesky(_nu, _bw, BandFactorization.Ldlt);

        double sigma2 = 1.0, tau2 = 1.0;
        int steps = burn + samples;
        for (int it = 0; it < steps; it++)
        {
            double lambda = sigma2 / tau2;
            Array.Copy(_gram0, band, _gram0.Length);        // A = ZᵀZ + λP
            _penalty.AccumulateInto(band, _nu, lambda);
            chol.DecomposeBanded(band);
            if (chol.Status is not { PositiveDefinite: true, FlooredPivots: 0 })
                throw new InvalidOperationException($"ZᵀZ + λP is not numerically positive definite at λ = {lambda} ({chol.Status}).");
            double[] betaHat = chol.Solve(_zty);            // posterior mean A⁻¹Zᵀy

            for (int p = 0; p < _nu; p++) z[p] = Normal.SampleStandard(rng);
            chol.SampleInnovation(z, innov);                // v ~ N(0, A⁻¹)
            double sigma = Math.Sqrt(sigma2);
            for (int p = 0; p < _nu; p++) beta[p] = betaHat[p] + sigma * innov[p];

            sigma2 = InverseGamma.ShapeScale(_priorShape + 0.5 * _n, _priorScale + 0.5 * ResidualSS(beta)).Sample(rng);
            tau2 = InverseGamma.ShapeScale(_priorShape + 0.5 * (_nu - _r), _priorScale + 0.5 * _penalty.Roughness(beta)).Sample(rng);

            if (it >= burn) draws.Add(new BapsDraw((double[])beta.Clone(), sigma2 / tau2, sigma2, tau2, null));
        }
        return draws;
    }

    // ── Marginal evidence: Metropolis on log λ against REML evidence × uniform-edf prior; β Rao-Blackwellized ──
    private List<BapsDraw> RunMarginal(Xoshiro256PlusPlus rng, int burn, int samples)
    {
        var draws = new List<BapsDraw>(samples);
        double logLam = 0.0;
        double ev = MarginalLogTarget(1.0);
        int steps = burn + samples;
        for (int it = 0; it < steps; it++)
        {
            double prop = logLam + _mhLogStep * Normal.SampleStandard(rng);
            double evProp = MarginalLogTarget(Math.Exp(prop));
            if (Math.Log(1.0 - rng.NextDouble()) < evProp - ev)
            {
                logLam = prop;
                ev = evProp;
            }
            if (it >= burn)
            {
                double lambda = Math.Exp(logLam);
                draws.Add(new BapsDraw((double[])_model.Coefficients(lambda).Clone(), lambda, _model.ProfiledVariance(lambda),
                                       double.NaN, null));
            }
        }
        return draws;
    }

    // log p(u | y) up to a constant on u = ln λ: the profiled evidence plus the uniform-edf prior density |d edf/du|.
    private double MarginalLogTarget(double lambda) => _model.RemlLogEvidence(lambda) + Math.Log(_model.EdfSlope(lambda));

    // ── Adaptive Gibbs: as RunGibbs, plus per-difference multipliers η_i ~ Gamma(φ, φ) ──
    private List<BapsDraw> RunAdaptiveGibbs(Xoshiro256PlusPlus rng, int burn, int samples)
    {
        int nLoc = _nu - _r;
        var draws = new List<BapsDraw>(samples);
        var band = new double[_bw + 1, _nu];
        var z = new double[_nu];
        var innov = new double[_nu];
        var beta = new double[_nu];
        var eta = new double[nLoc];
        var sqDiff = new double[nLoc];
        var weights = new double[nLoc];        // λ_i = (σ²/τ²)·η_i
        var chol = new BandCholesky(_nu, _bw, BandFactorization.Ldlt);

        double sigma2 = 1.0, tau2 = 1.0;
        for (int i = 0; i < nLoc; i++) eta[i] = 1.0;   // start unadapted (= global)
        int steps = burn + samples;
        for (int it = 0; it < steps; it++)
        {
            double lamGlobal = sigma2 / tau2;
            for (int i = 0; i < nLoc; i++) weights[i] = lamGlobal * eta[i];
            Array.Copy(_gram0, band, _gram0.Length);            // A = ZᵀZ + Dᵀdiag(λ_i)D
            _penalty.AccumulateInto(band, _nu, weights);
            chol.DecomposeBanded(band);
            if (chol.Status is not { PositiveDefinite: true, FlooredPivots: 0 })
                throw new InvalidOperationException($"ZᵀZ + Dᵀdiag(λᵢ)D is not numerically positive definite at global λ = {lamGlobal} ({chol.Status}).");
            double[] betaHat = chol.Solve(_zty);

            for (int p = 0; p < _nu; p++) z[p] = Normal.SampleStandard(rng);
            chol.SampleInnovation(z, innov);
            double sigma = Math.Sqrt(sigma2);
            for (int p = 0; p < _nu; p++) beta[p] = betaHat[p] + sigma * innov[p];

            _penalty.SquaredDifferencesInto(beta, sqDiff);       // (Δ^rβ)_i²

            sigma2 = InverseGamma.ShapeScale(_priorShape + 0.5 * _n, _priorScale + 0.5 * ResidualSS(beta)).Sample(rng);

            double weightedRough = 0.0;                          // Σ η_i (Δβ)_i²
            for (int i = 0; i < nLoc; i++) weightedRough += eta[i] * sqDiff[i];
            tau2 = InverseGamma.ShapeScale(_priorShape + 0.5 * nLoc, _priorScale + 0.5 * weightedRough).Sample(rng);

            // η_i | β, τ² ~ Gamma(φ + ½, φ + (Δβ)_i²/(2τ²)) — the +½ is the prior-normalizer Jacobian.
            for (int i = 0; i < nLoc; i++)
                eta[i] = Gamma.ShapeRate(_phi + 0.5, _phi + 0.5 * sqDiff[i] / tau2).Sample(rng);

            if (it >= burn) draws.Add(new BapsDraw((double[])beta.Clone(), sigma2 / tau2, sigma2, tau2, (double[])eta.Clone()));
        }
        return draws;
    }

    private double ResidualSS(double[] beta)
    {
        double ss = 0.0;
        for (int i = 0; i < _n; i++)
        {
            double fit = 0.0;
            for (int j = 0; j < _nu; j++) fit += _z[i, j] * beta[j];
            double d = _y[i] - fit;
            ss += d * d;
        }
        return ss;
    }

    private static double Quantile(double[] sorted, double p)
    {
        double pos = p * (sorted.Length - 1);
        int lo = (int)Math.Floor(pos);
        int hi = Math.Min(lo + 1, sorted.Length - 1);
        return sorted[lo] + (pos - lo) * (sorted[hi] - sorted[lo]);
    }
}
