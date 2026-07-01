using System;
using System.Collections.Generic;
using Maths.Rng;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Bgp;

/// <summary>The pooled posterior of a Bayesian-GP run: the Bayes-averaged conditional mean at the requested points
/// plus the bandwidth posterior (the data-chosen length scale, with a credible interval) and the sampler's mixing
/// and convergence.</summary>
/// <param name="PosteriorMean">f̂(x*) at each test point — averaged over the bandwidth posterior (model averaging).</param>
/// <param name="BandwidthMean">t̂ = posterior-mean bandwidth.</param>
/// <param name="BandwidthLo">2.5% bandwidth quantile.</param>
/// <param name="BandwidthHi">97.5% bandwidth quantile.</param>
/// <param name="AcceptanceRate">Metropolis acceptance over the sampled steps (pooled across chains).</param>
/// <param name="Draws">Post-burn draws averaged (pooled = chains × per-chain draws).</param>
/// <param name="BandwidthRHat">Gelman–Rubin R̂ on the bandwidth across chains (NaN for a single chain) — ≈1 at convergence.</param>
/// <param name="BandwidthEss">Effective sample size of the bandwidth across the pooled chains.</param>
/// <param name="ChainSeeds">The resolved per-chain seeds (reproducibility provenance).</param>
public sealed record BgpResult(
    double[] PosteriorMean,
    double BandwidthMean, double BandwidthLo, double BandwidthHi,
    double AcceptanceRate, int Draws,
    double BandwidthRHat, double BandwidthEss, int[] ChainSeeds);

/// <summary>
/// The top of the Bayesian-GP instrument (Tang, Wu, Cheng &amp; Dunson 2025): with the regression function
/// marginalized out, this samples the one remaining parameter — the kernel bandwidth t — from its posterior
/// ∝ p(Y|t)·p(t) by random-walk Metropolis–Hastings on log t (the marginal evidence from <see cref="GpRegression"/>
/// times the empirical-Bayes prior from <see cref="EmpiricalBayesBandwidthPrior"/>), and reduces the t-draws to
/// the Bayes-averaged conditional mean. Sampling in log t carries the change-of-variables term so the stationary
/// law is the true t-posterior; the readout averages each draw's GP conditional mean, so the bandwidth uncertainty
/// is integrated rather than plugged in. Structurally the marginalize-coefficients / sample-the-hyperparameter
/// twin of the BAPS λ-sampler — here the hyperparameter is a length scale, not a penalty.
/// As a client of the kernel-agnostic ensemble shell (<see cref="ChainEnsemble"/>) it now runs multiple
/// warm-started replicas and gains the cross-chain R̂ + ESS it previously lacked; the target is unimodal /
/// log-concave, so the chains warm-start at the prior midpoint rather than overdisperse.
/// </summary>
public sealed class BgpSampler : IEnsembleModel<double>
{
    private readonly GpRegression _gp;
    private readonly EmpiricalBayesBandwidthPrior _prior;

    // Per-run context (set at the top of Run, read by StartChain during the synchronous fan-out).
    private double[,] _xTest = null!;
    private double _t0;
    private double _proposalSd;
    private int _testCount;

    public BgpSampler(GpRegression gp, EmpiricalBayesBandwidthPrior prior)
    {
        ArgumentNullException.ThrowIfNull(gp);
        ArgumentNullException.ThrowIfNull(prior);
        _gp = gp; _prior = prior;
    }

    /// <summary>The single R̂/ESS functional: the bandwidth t.</summary>
    public int FunctionalDim => 1;

    /// <summary>ESS on the bandwidth.</summary>
    public int EssDim => 1;

    /// <summary>Warm-start one bandwidth chain at the prior midpoint (<paramref name="readoutSeed"/> unused — the readout is deterministic).</summary>
    public IChainRun<double> StartChain(int chainSeed, int readoutSeed)
    {
        var chain = new BgpBandwidthChain(_gp, _prior, _t0, _proposalSd, new Xoshiro256PlusPlus(chainSeed));
        return new BgpChainRun(chain, _gp, _xTest, _testCount);
    }

    /// <param name="xTest">Points at which to report the posterior mean (the training inputs for the in-sample fit).</param>
    /// <param name="draws">Post-burn draws to average, per chain.</param>
    /// <param name="burn">Burn-in steps per chain.</param>
    /// <param name="proposalSd">Random-walk standard deviation on log t.</param>
    /// <param name="seed">Master RNG seed (the ensemble fans out per-chain streams from it).</param>
    /// <param name="tInit">Optional starting bandwidth (default: the geometric midpoint of the prior's support); shared by every chain.</param>
    /// <param name="chains">Number of warm-started replicas (≥ 2 makes R̂ defined; default 2).</param>
    public BgpResult Run(double[,] xTest, int draws = 600, int burn = 200, double proposalSd = 0.4,
                         int seed = 0, double? tInit = null, int chains = 2)
    {
        ArgumentNullException.ThrowIfNull(xTest);
        if (draws < 1) throw new ArgumentOutOfRangeException(nameof(draws));
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));

        double t0 = tInit ?? Math.Sqrt(Math.Max(_prior.LowerBound, 1e-12));   // geometric mid of (γ₁T_n², 1]
        t0 = Math.Min(0.999, Math.Max(t0, _prior.LowerBound * 1.000001));

        _xTest = xTest;
        _t0 = t0;
        _proposalSd = proposalSd;
        _testCount = xTest.GetLength(0);

        (EnsembleRun summary, IReadOnlyList<IChainRun<double>> runs) =
            new ChainEnsemble().Run<double>(this, chains, seed, burn, draws);

        int m = _testCount;
        var meanAcc = new double[m];
        var allT = new List<double>(chains * summary.SamplesUsed);
        foreach (IChainRun<double> r in runs)
        {
            var br = (BgpChainRun)r;
            for (int i = 0; i < m; i++) meanAcc[i] += br.MeanAcc[i];
            allT.AddRange(br.TDraws);
        }
        double denom = chains * (double)summary.SamplesUsed;
        for (int i = 0; i < m; i++) meanAcc[i] /= denom;
        allT.Sort();

        return new BgpResult(meanAcc, summary.FunctionalMean[0], Quantile(allT, 0.025), Quantile(allT, 0.975),
                             summary.AcceptanceRate, chains * summary.SamplesUsed,
                             summary.FunctionalRHat[0], summary.FunctionalEss[0], summary.ChainSeeds);
    }

    private static double Quantile(List<double> sorted, double q)
    {
        int n = sorted.Count;
        if (n == 0) return double.NaN;
        if (n == 1) return sorted[0];
        double pos = q * (n - 1);
        int lo = (int)Math.Floor(pos);
        double frac = pos - lo;
        return lo + 1 < n ? sorted[lo] * (1.0 - frac) + sorted[lo + 1] * frac : sorted[lo];
    }

    /// <summary>
    /// One random-walk Metropolis chain on log t (the bandwidth kernel): each <see cref="Step"/> proposes
    /// u' = u + σ·N(0,1), accepts under the marginal-evidence × prior ratio with the +u log-t Jacobian preserved,
    /// and returns the current bandwidth t = exp(u). <see cref="CurrentFit"/> is the accepted GP fit the readout reuses.
    /// </summary>
    private sealed class BgpBandwidthChain : IChain<double>
    {
        private readonly GpRegression _gp;
        private readonly EmpiricalBayesBandwidthPrior _prior;
        private readonly Xoshiro256PlusPlus _rng;
        private readonly double _proposalSd;
        private double _u;          // log t
        private double _cur;        // log target at the current state
        private GpFit _fit;

        public long Accepted { get; private set; }
        public long Attempts { get; private set; }

        /// <summary>The GP fit at the current accepted bandwidth — the readout's conditional-mean source.</summary>
        public GpFit CurrentFit => _fit;

        public BgpBandwidthChain(GpRegression gp, EmpiricalBayesBandwidthPrior prior, double t0, double proposalSd, Xoshiro256PlusPlus rng)
        {
            _gp = gp; _prior = prior; _proposalSd = proposalSd; _rng = rng;
            _u = Math.Log(t0);
            _fit = gp.Fit(t0);
            _cur = _fit.LogMarginal + prior.LogDensity(t0) + _u;   // +u: log t Jacobian
        }

        public double Step()
        {
            double uProp = _u + _proposalSd * Gaussian(_rng);
            double tProp = Math.Exp(uProp);
            double priorProp = _prior.LogDensity(tProp);
            Attempts++;

            GpFit fitProp = null;
            double prop = double.NegativeInfinity;
            if (!double.IsNegativeInfinity(priorProp))
            {
                fitProp = _gp.Fit(tProp);
                prop = fitProp.LogMarginal + priorProp + uProp;
            }

            if (prop - _cur >= 0.0 || _rng.NextDouble() < Math.Exp(prop - _cur))
            {
                _u = uProp; _cur = prop; _fit = fitProp; Accepted++;
            }
            return Math.Exp(_u);
        }

        private static double Gaussian(Xoshiro256PlusPlus rng)
        {
            double u1 = 1.0 - rng.NextDouble();
            double u2 = 1.0 - rng.NextDouble();
            return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }
    }

    /// <summary>
    /// One bandwidth chain + its readout sinks: the per-draw bandwidth t (the R̂/ESS functional) plus the model's
    /// own accumulators — the running sum of the GP conditional mean (pooled → the Bayes-averaged f̂) and the
    /// t-draw cloud (pooled → the bandwidth credible interval).
    /// </summary>
    private sealed class BgpChainRun : IChainRun<double>
    {
        private readonly BgpBandwidthChain _chain;
        private readonly GpRegression _gp;
        private readonly double[,] _xTest;

        public readonly double[] MeanAcc;       // Σ over draws of the GP conditional mean at the test points
        public readonly List<double> TDraws;    // the post-burn bandwidth cloud

        public BgpChainRun(BgpBandwidthChain chain, GpRegression gp, double[,] xTest, int testCount)
        {
            _chain = chain; _gp = gp; _xTest = xTest;
            MeanAcc = new double[testCount];
            TDraws = new List<double>();
        }

        public IChain<double> Chain => _chain;

        public void Burn(int steps)
        {
            for (int i = 0; i < steps; i++) _chain.Step();
        }

        public void Accumulate(in double draw, Span<double> functionals)
        {
            functionals[0] = draw;
            TDraws.Add(draw);
            double[] pred = _gp.PredictMean(_chain.CurrentFit, _xTest);
            for (int i = 0; i < MeanAcc.Length; i++) MeanAcc[i] += pred[i];
        }
    }
}
