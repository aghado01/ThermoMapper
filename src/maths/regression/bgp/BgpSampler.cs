using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Bgp;

/// <summary>A bandwidth posterior in one frame: the posterior mean and the type-7 2.5% and 97.5% quantiles.</summary>
public sealed record BandwidthSummary(double Mean, double Lo, double Hi);

/// <summary>
/// One retained bandwidth draw: the bandwidth t in the mapped frame, the GP fit the chain held there (its log marginal,
/// α and mean), and that fit's conditional variance of f at the run's test points
/// (<see cref="GpRegression.PredictVariance"/>). Given t the GP posterior is fixed, so this is the draw's whole state.
/// A rejected proposal keeps the previous fit, so neighbouring draws share one <see cref="GpFit"/> and one variance
/// array.
/// </summary>
public sealed record BgpDraw(double Bandwidth, GpFit Fit, double[] ConditionalVariance);

/// <summary>The pooled posterior of a Bayesian-GP run: the Bayes-averaged conditional mean at the requested points
/// plus the bandwidth posterior (the data-chosen length scale, with a credible interval) and the sampler's mixing
/// and convergence.</summary>
/// <param name="PosteriorMean">f̂(x*) at each test point — averaged over the bandwidth posterior (model averaging).</param>
/// <param name="PosteriorSd">The posterior sd of f(x*) at each test point, by the law of total variance: the spread of the
/// draws' conditional means plus their mean conditional variance (<see cref="BgpDraw.ConditionalVariance"/>). It is the
/// uncertainty of the regression function; a new observation adds σ².</param>
/// <param name="Bandwidth">The bandwidth posterior in the GP's mapped frame, where the sampler and the prior work.</param>
/// <param name="BandwidthOriginal">The same posterior in the original input units: every value is s²·t (exact for the mean
/// and, the map being monotone, for the quantiles).</param>
/// <param name="AcceptanceRate">Metropolis acceptance over the sampled steps (pooled across chains).</param>
/// <param name="Draws">Post-burn draws averaged (pooled = chains × per-chain draws).</param>
/// <param name="BandwidthConvergence">The bandwidth's rank-normalized diagnostics across chains: split R̂ (≈1 at
/// convergence) and bulk and tail effective sample sizes.</param>
/// <param name="ChainSeeds">The resolved per-chain seeds (reproducibility provenance).</param>
/// <param name="InputScale">The GP's fitted input scale s (<see cref="GpTransform.InputScale"/>): a bandwidth t in the mapped frame is s²·t in the original units.</param>
/// <param name="BandwidthDraws">Every retained draw (<see cref="BgpDraw"/>), per chain in draw order. The posterior mean
/// and the bandwidth summary are read from these.</param>
/// <param name="Checkpoint">The run's restart state when it was requested, and null otherwise.</param>
public sealed record BgpResult(
    double[] PosteriorMean,
    double[] PosteriorSd,
    BandwidthSummary Bandwidth, BandwidthSummary BandwidthOriginal,
    double AcceptanceRate, int Draws,
    ConvergenceDiagnostic BandwidthConvergence, int[] ChainSeeds, double InputScale,
    IReadOnlyList<IReadOnlyList<BgpDraw>> BandwidthDraws,
    EnsembleCheckpoint? Checkpoint);

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
    public int MonitoredDim => 1;

    /// <summary>Warm-start one bandwidth chain at the prior midpoint (<paramref name="readoutSeed"/> unused — the readout is deterministic).</summary>
    public IChainRun<double> StartChain(int chainSeed, int readoutSeed)
    {
        var chain = new BgpBandwidthChain(_gp, _prior, _t0, _proposalSd, new Xoshiro256PlusPlus(chainSeed));
        return new BgpChainRun(chain, _gp, _xTest, new List<BgpDraw>());
    }

    /// <summary>The run's fingerprint for checkpoints: the GP's data, noise variance, kernel and input map, the bandwidth
    /// prior, the test points, the starting bandwidth and the proposal scale.</summary>
    public string Fingerprint
    {
        get
        {
            var f = new RunFingerprint().Add("BGP");
            _gp.AddTo(f);
            _prior.AddTo(f);
            foreach (double v in _xTest) f.Add(v);
            return f.Add(_testCount).Add(_t0).Add(_proposalSd).ToString();
        }
    }

    /// <summary>Rebuild one chain from its restart state: its stream, its log bandwidth, the log target there, its
    /// counters and its retained draws. The current fit is recomputed at the bandwidth it was made at, which gives the
    /// same fit.</summary>
    public IChainRun<double> RestoreChain(IChainState state)
    {
        if (state is not BgpChainState s) throw new ArgumentException("Not a BGP chain's restart state.", nameof(state));
        var chain = new BgpBandwidthChain(_gp, _prior, _proposalSd, s);
        return new BgpChainRun(chain, _gp, _xTest, new List<BgpDraw>(s.Draws));
    }

    /// <param name="xTest">Points (original units) at which to report the posterior mean (the training inputs for the in-sample fit).</param>
    /// <param name="draws">Post-burn draws to average, per chain.</param>
    /// <param name="burn">Burn-in steps per chain.</param>
    /// <param name="proposalSd">Random-walk standard deviation on log t.</param>
    /// <param name="seed">Master RNG seed (the ensemble fans out per-chain streams from it).</param>
    /// <param name="tInit">Optional starting bandwidth inside the prior's support (default: the support's geometric midpoint); shared by every chain.</param>
    /// <param name="chains">Number of warm-started replicas (≥ 2 makes R̂ defined; default 2).</param>
    /// <param name="rHatTarget">When &gt; 0, stop after the first round in which the bandwidth meets this R̂ target and its
    /// bulk and tail ESS reach <paramref name="minEssPerChain"/> per chain (<see cref="ChainEnsemble"/>'s rule). It
    /// needs at least two chains.</param>
    /// <param name="batchSize">Draws per round between stop checks; ≤ 0 is a single round.</param>
    /// <param name="minEssPerChain">The early stop's ESS gate, per chain.</param>
    /// <param name="checkpoint">When true, the result carries the run's restart state.</param>
    /// <param name="resumeFrom">Continue that checkpoint's chains for <paramref name="draws"/> more draws, exactly as the
    /// uninterrupted run would have; the result covers every draw since burn-in. It needs the same GP, prior, test
    /// points, settings, master seed and chain count.</param>
    public BgpResult Run(double[,] xTest, int draws = 600, int burn = 200, double proposalSd = 0.4,
                         int seed = 0, double? tInit = null, int chains = 2,
                         double rHatTarget = 0.0, int batchSize = 0, double minEssPerChain = ChainEnsemble.MinEssPerChain,
                         bool checkpoint = false, EnsembleCheckpoint? resumeFrom = null)
    {
        ArgumentNullException.ThrowIfNull(xTest);
        if (draws < 1) throw new ArgumentOutOfRangeException(nameof(draws));
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));

        double t0 = tInit ?? Math.Sqrt(_prior.LowerBound);   // geometric mid of (γ₁T_n², 1], inside any valid support
        if (!(t0 > _prior.LowerBound && t0 <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(tInit), t0, "The starting bandwidth must lie in the prior support (γ₁T_n², 1].");

        _xTest = _gp.Transform.MapInputs(xTest);   // mapped once; every draw predicts in the fitted frame
        _t0 = t0;
        _proposalSd = proposalSd;
        _testCount = xTest.GetLength(0);

        (EnsembleRun summary, IReadOnlyList<IChainRun<double>> runs) =
            new ChainEnsemble().Run<double>(this, chains, seed, burn, draws, rHatTarget, batchSize, minEssPerChain,
                                            checkpoint, resumeFrom);

        // Every summary is read from the retained draws. The posterior mean sums each draw's conditional mean at the test
        // points per chain in draw order, then over chains; a fit shared by neighbouring draws is predicted once. The
        // posterior sd adds, by the law of total variance, the spread of those conditional means about the mean (a
        // second pass) to the draws' mean conditional variance. The bandwidth mean averages the chains' running means,
        // and its interval reads the pooled sample.
        int m = _testCount;
        var meanAcc = new double[m];
        var varianceAcc = new double[m];
        var allT = new List<double>(chains * summary.SamplesUsed);
        var preds = new List<double[]>(chains * summary.SamplesUsed);
        var bandwidthDraws = new IReadOnlyList<BgpDraw>[chains];
        double chainMeans = 0.0;
        for (int ci = 0; ci < chains; ci++)
        {
            List<BgpDraw> chainDraws = ((BgpChainRun)runs[ci]).Draws;
            bandwidthDraws[ci] = chainDraws.ToArray();
            var chainSum = new double[m];
            var moments = new RunningMoments();
            GpFit? predicted = null;
            double[] pred = [];
            foreach (BgpDraw d in chainDraws)
            {
                if (!ReferenceEquals(d.Fit, predicted))
                {
                    pred = _gp.PredictMeanInFrame(d.Fit, _xTest);
                    predicted = d.Fit;
                }
                for (int i = 0; i < m; i++) chainSum[i] += pred[i];
                for (int i = 0; i < m; i++) varianceAcc[i] += d.ConditionalVariance[i];
                preds.Add(pred);
                moments.Add(d.Bandwidth);
                allT.Add(d.Bandwidth);
            }
            for (int i = 0; i < m; i++) meanAcc[i] += chainSum[i];
            chainMeans += moments.Mean;
        }
        double denom = chains * (double)summary.SamplesUsed;
        for (int i = 0; i < m; i++) meanAcc[i] /= denom;
        allT.Sort();

        var spread = new double[m];
        foreach (double[] pred in preds)
            for (int i = 0; i < m; i++)
            {
                double dev = pred[i] - meanAcc[i];
                spread[i] += dev * dev;
            }
        var posteriorSd = new double[m];
        for (int i = 0; i < m; i++) posteriorSd[i] = Math.Sqrt((spread[i] + varianceAcc[i]) / denom);

        var fitted = new BandwidthSummary(chainMeans / chains, Quantile(allT, 0.025), Quantile(allT, 0.975));
        double s2 = _gp.Transform.InputScale * _gp.Transform.InputScale;   // t_original = s²·t
        var original = new BandwidthSummary(s2 * fitted.Mean, s2 * fitted.Lo, s2 * fitted.Hi);

        return new BgpResult(meanAcc, posteriorSd, fitted, original,
                             summary.AcceptanceRate, chains * summary.SamplesUsed,
                             summary.Monitored[0], summary.ChainSeeds,
                             _gp.Transform.InputScale, bandwidthDraws, summary.Checkpoint);
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

    /// <summary>One bandwidth chain's restart state: its stream's words, its log bandwidth and the log target there, the
    /// bandwidth its current fit was made at, its counters and its retained draws.</summary>
    private sealed record BgpChainState(
        (ulong, ulong, ulong, ulong) Rng, double LogBandwidth, double LogTarget, double FitBandwidth, long Accepted,
        long Attempts, List<BgpDraw> Draws) : IChainState;

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
        private GpFactoredFit _factored;   // the current fit with its factor, held only while current

        public long Accepted { get; private set; }
        public long Attempts { get; private set; }

        /// <summary>The GP fit at the current accepted bandwidth — the readout's conditional-mean source.</summary>
        public GpFit CurrentFit => _factored.Fit;

        /// <summary>The current fit with its factor — the readout's conditional-variance source.</summary>
        public GpFactoredFit CurrentFactored => _factored;

        public BgpBandwidthChain(GpRegression gp, EmpiricalBayesBandwidthPrior prior, double t0, double proposalSd, Xoshiro256PlusPlus rng)
        {
            _gp = gp; _prior = prior; _proposalSd = proposalSd; _rng = rng;
            _u = Math.Log(t0);
            _factored = gp.FitFactored(t0);
            _cur = _factored.Fit.LogMarginal + prior.LogDensity(t0) + _u;   // +u: log t Jacobian
        }

        /// <summary>A chain continuing from its restart state. The fit is recomputed at the bandwidth it was made at —
        /// the start's t0 can differ from exp(log t0) in the last bit — so it is the same fit.</summary>
        public BgpBandwidthChain(GpRegression gp, EmpiricalBayesBandwidthPrior prior, double proposalSd, BgpChainState s)
        {
            _gp = gp; _prior = prior; _proposalSd = proposalSd;
            _rng = Xoshiro256PlusPlus.Restore(s.Rng.Item1, s.Rng.Item2, s.Rng.Item3, s.Rng.Item4);
            _u = s.LogBandwidth;
            _cur = s.LogTarget;
            _factored = gp.FitFactored(s.FitBandwidth);
            Accepted = s.Accepted;
            Attempts = s.Attempts;
        }

        /// <summary>The chain's part of its restart state.</summary>
        public BgpChainState Save(List<BgpDraw> draws) =>
            new(_rng.SaveState(), _u, _cur, _factored.Fit.Bandwidth, Accepted, Attempts, new List<BgpDraw>(draws));

        public double Step()
        {
            double uProp = _u + _proposalSd * Normal.SampleStandard(_rng);
            double tProp = Math.Exp(uProp);
            double priorProp = _prior.LogDensity(tProp);
            Attempts++;

            GpFactoredFit? fitProp = null;
            double prop = double.NegativeInfinity;
            if (!double.IsNegativeInfinity(priorProp))
            {
                fitProp = _gp.FitFactored(tProp);
                prop = fitProp.Fit.LogMarginal + priorProp + uProp;
            }

            if (prop - _cur >= 0.0 || _rng.NextDouble() < Math.Exp(prop - _cur))
            {
                _u = uProp; _cur = prop; _factored = fitProp!; Accepted++;   // prop is −∞ without a fit, so never accepted
            }
            return Math.Exp(_u);
        }

    }

    /// <summary>
    /// One bandwidth chain and its retained draws: the per-draw bandwidth t (the R̂/ESS functional) and, per draw, the
    /// bandwidth with the GP fit the chain held there and its conditional variance at the test points, computed once
    /// per distinct fit while the fit's factor exists. The run's summaries are read from these.
    /// </summary>
    private sealed class BgpChainRun : IResumableChainRun<double>
    {
        private readonly BgpBandwidthChain _chain;
        private readonly GpRegression _gp;
        private readonly double[,] _xTest;
        private GpFit? _varianceFit;          // the fit _variance belongs to
        private double[] _variance = [];

        public readonly List<BgpDraw> Draws;   // every post-burn draw, in order

        public BgpChainRun(BgpBandwidthChain chain, GpRegression gp, double[,] xTest, List<BgpDraw> draws)
        {
            _chain = chain; _gp = gp; _xTest = xTest; Draws = draws;
        }

        public IChainState SaveState() => _chain.Save(Draws);

        public IChain<double> Chain => _chain;

        public void Burn(int steps)
        {
            for (int i = 0; i < steps; i++) _chain.Step();
        }

        public void Accumulate(in double draw, Span<double> functionals)
        {
            functionals[0] = draw;
            GpFit fit = _chain.CurrentFit;
            if (!ReferenceEquals(fit, _varianceFit))
            {
                _variance = _gp.ConditionalVarianceInFrame(_chain.CurrentFactored, _xTest);
                _varianceFit = fit;
            }
            Draws.Add(new BgpDraw(draw, fit, _variance));
        }
    }
}
