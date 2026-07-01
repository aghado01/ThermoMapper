using System;
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
    /// Metropolis on <c>log λ</c> against the profiled REML marginal evidence (β integrated out), with a flat
    /// prior on <c>log λ</c>; β is Rao-Blackwellized to its posterior mean at each draw. Marginalizes the spline.
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
/// <paramref name="LambdaMean"/>/<paramref name="LambdaLo"/>/<paramref name="LambdaHi"/> summarize the global
/// smoothing (in the adaptive mode, λ_global). <paramref name="LocalSmoothing"/> is the posterior-mean local
/// smoothing field λ_i (length ν−r) — non-null only for <see cref="BapsLambdaUpdate.AdaptiveGibbs"/>; its dips
/// mark where the fit was allowed to track features.
/// </remarks>
public sealed record BapsResult(
    double[] Coefficients,
    double LambdaMean,
    double LambdaLo,
    double LambdaHi,
    double SigmaMean,
    double RHatLogLambda,
    int Draws,
    double[]? LocalSmoothing = null);

/// <summary>
/// The BAPS (Bayesian adaptive P-spline) sampler — the measure-side dual of free-knot BARS. A rich fixed B-spline
/// basis with a difference penalty whose strength λ is inferred, sampled either by conjugate Gibbs over the
/// variance components or by Metropolis against the REML marginal evidence (<see cref="BapsLambdaUpdate"/>). Every
/// inner solve is banded (<c>A = ZᵀZ + λP</c> via <see cref="BandCholesky"/>'s root-free LDLᵀ, whose factor also
/// supplies the Gaussian β-draw); chains run in parallel from one seed tree and reduce per draw to the Bayes fit,
/// the λ posterior (mean + 95% interval), σ, and R̂ on log λ.
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
        var accums = new ChainAccum[chains];
        Parallel.For(0, chains, c =>
        {
            var rng = new Xoshiro256PlusPlus(seeds[c]);
            accums[c] = mode switch
            {
                BapsLambdaUpdate.Gibbs => RunGibbs(rng, burn, samples),
                BapsLambdaUpdate.MarginalEvidence => RunMarginal(rng, burn, samples),
                BapsLambdaUpdate.AdaptiveGibbs => RunAdaptiveGibbs(rng, burn, samples),
                _ => throw new ArgumentOutOfRangeException(nameof(mode)),
            };
        });

        double total = chains * (double)samples;
        var coef = new double[_nu];
        var logLamSum = new double[chains];
        var logLamSumSq = new double[chains];
        double sigSum = 0.0;
        var allLam = new double[(int)total];
        int idx = 0;
        for (int c = 0; c < chains; c++)
        {
            ChainAccum a = accums[c];
            for (int p = 0; p < _nu; p++) coef[p] += a.CoefSum[p];
            logLamSum[c] = a.LogLamSum;
            logLamSumSq[c] = a.LogLamSumSq;
            sigSum += a.SigmaSum;
            Array.Copy(a.Lambdas, 0, allLam, idx, samples);
            idx += samples;
        }
        for (int p = 0; p < _nu; p++) coef[p] /= total;

        double[]? localField = null;
        if (accums[0].LocalLambdaSum is not null)
        {
            int nLoc = accums[0].LocalLambdaSum!.Length;
            localField = new double[nLoc];
            for (int c = 0; c < chains; c++)
                for (int k = 0; k < nLoc; k++) localField[k] += accums[c].LocalLambdaSum![k];
            for (int k = 0; k < nLoc; k++) localField[k] /= total;
        }

        double lamMean = 0.0;
        foreach (double l in allLam) lamMean += l;
        lamMean /= total;
        Array.Sort(allLam);

        return new BapsResult(
            coef, lamMean, Quantile(allLam, 0.025), Quantile(allLam, 0.975),
            sigSum / total, ChainDiagnostics.RHat(logLamSum, logLamSumSq, samples), (int)total, localField);
    }

    // ── Gibbs: β | σ²,τ²  then  σ² | β  then  τ² | β ───────────────────────────
    private ChainAccum RunGibbs(Xoshiro256PlusPlus rng, int burn, int samples)
    {
        var acc = new ChainAccum(_nu, samples, 0);
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
            double[] betaHat = chol.Solve(_zty);            // posterior mean A⁻¹Zᵀy

            for (int p = 0; p < _nu; p++) z[p] = StandardNormal(rng);
            chol.SampleInnovation(z, innov);                // v ~ N(0, A⁻¹)
            double sigma = Math.Sqrt(sigma2);
            for (int p = 0; p < _nu; p++) beta[p] = betaHat[p] + sigma * innov[p];

            sigma2 = InvGamma(rng, _priorShape + 0.5 * _n, _priorScale + 0.5 * ResidualSS(beta));
            tau2 = InvGamma(rng, _priorShape + 0.5 * (_nu - _r), _priorScale + 0.5 * _penalty.Roughness(beta));

            if (it >= burn) acc.Record(it - burn, beta, sigma2 / tau2, Math.Sqrt(sigma2));
        }
        return acc;
    }

    // ── Marginal evidence: Metropolis on log λ against the REML evidence; β Rao-Blackwellized ──
    private ChainAccum RunMarginal(Xoshiro256PlusPlus rng, int burn, int samples)
    {
        var acc = new ChainAccum(_nu, samples, 0);
        double logLam = 0.0;
        double ev = _model.RemlLogEvidence(1.0);
        int steps = burn + samples;
        for (int it = 0; it < steps; it++)
        {
            double prop = logLam + _mhLogStep * StandardNormal(rng);
            double evProp = _model.RemlLogEvidence(Math.Exp(prop));
            if (Math.Log(1.0 - rng.NextDouble()) < evProp - ev)   // flat prior on log λ ⇒ pure evidence ratio
            {
                logLam = prop;
                ev = evProp;
            }
            if (it >= burn)
            {
                double lambda = Math.Exp(logLam);
                acc.Record(it - burn, _model.Coefficients(lambda), lambda, Math.Sqrt(_model.ProfiledVariance(lambda)));
            }
        }
        return acc;
    }

    // ── Adaptive Gibbs: as RunGibbs, plus per-difference multipliers η_i ~ Gamma(φ, φ) ──
    private ChainAccum RunAdaptiveGibbs(Xoshiro256PlusPlus rng, int burn, int samples)
    {
        int nLoc = _nu - _r;
        var acc = new ChainAccum(_nu, samples, nLoc);
        var band = new double[_bw + 1, _nu];
        var z = new double[_nu];
        var innov = new double[_nu];
        var beta = new double[_nu];
        var eta = new double[nLoc];
        var sqDiff = new double[nLoc];
        var weights = new double[nLoc];        // λ_i = (σ²/τ²)·η_i
        var localLambda = new double[nLoc];
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
            double[] betaHat = chol.Solve(_zty);

            for (int p = 0; p < _nu; p++) z[p] = StandardNormal(rng);
            chol.SampleInnovation(z, innov);
            double sigma = Math.Sqrt(sigma2);
            for (int p = 0; p < _nu; p++) beta[p] = betaHat[p] + sigma * innov[p];

            _penalty.SquaredDifferencesInto(beta, sqDiff);       // (Δ^rβ)_i²

            sigma2 = InvGamma(rng, _priorShape + 0.5 * _n, _priorScale + 0.5 * ResidualSS(beta));

            double weightedRough = 0.0;                          // Σ η_i (Δβ)_i²
            for (int i = 0; i < nLoc; i++) weightedRough += eta[i] * sqDiff[i];
            tau2 = InvGamma(rng, _priorShape + 0.5 * nLoc, _priorScale + 0.5 * weightedRough);

            // η_i | β, τ² ~ Gamma(φ + ½, φ + (Δβ)_i²/(2τ²)) — the +½ is the prior-normalizer Jacobian.
            for (int i = 0; i < nLoc; i++)
                eta[i] = Gamma.Sample(rng, _phi + 0.5, 1.0 / (_phi + 0.5 * sqDiff[i] / tau2));

            if (it >= burn)
            {
                double lg = sigma2 / tau2;
                for (int i = 0; i < nLoc; i++) localLambda[i] = lg * eta[i];
                acc.Record(it - burn, beta, lg, Math.Sqrt(sigma2), localLambda);
            }
        }
        return acc;
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

    private static double InvGamma(Xoshiro256PlusPlus rng, double shape, double scale)
        => 1.0 / Gamma.Sample(rng, shape, 1.0 / scale);

    private static double StandardNormal(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double Quantile(double[] sorted, double p)
    {
        double pos = p * (sorted.Length - 1);
        int lo = (int)Math.Floor(pos);
        int hi = Math.Min(lo + 1, sorted.Length - 1);
        return sorted[lo] + (pos - lo) * (sorted[hi] - sorted[lo]);
    }

    // Per-chain reduction: running Σβ, the λ draws, and the log-λ moments for R̂.
    private sealed class ChainAccum
    {
        public readonly double[] CoefSum;
        public readonly double[] Lambdas;
        public readonly double[]? LocalLambdaSum;   // running Σ λ_i field (adaptive mode only)
        public double LogLamSum;
        public double LogLamSumSq;
        public double SigmaSum;

        public ChainAccum(int nu, int samples, int nLocal)
        {
            CoefSum = new double[nu];
            Lambdas = new double[samples];
            LocalLambdaSum = nLocal > 0 ? new double[nLocal] : null;
        }

        public void Record(int i, double[] beta, double lambda, double sigma, double[]? localLambda = null)
        {
            for (int p = 0; p < CoefSum.Length; p++) CoefSum[p] += beta[p];
            Lambdas[i] = lambda;
            double ll = Math.Log(lambda);
            LogLamSum += ll;
            LogLamSumSq += ll * ll;
            SigmaSum += sigma;
            if (localLambda is not null && LocalLambdaSum is not null)
                for (int k = 0; k < LocalLambdaSum.Length; k++) LocalLambdaSum[k] += localLambda[k];
        }
    }
}
