using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Bgp;

/// <summary>A bandwidth posterior in one frame: the posterior mean and the type-7 2.5% and 97.5% quantiles.</summary>
public sealed record BandwidthSummary(double Mean, double Lo, double Hi);

/// <summary>
/// One retained draw: the bandwidth t in the mapped frame, the GP fit the chain held there (its log marginal, α, mean,
/// noise ratio and residual quadratic), the signal variance s², and the conditional variance of f at the run's test
/// points given all of them (<see cref="GpRegression.PredictVariance"/> times s²).
/// <list type="bullet">
///   <item>For a known σ², s² = 1 and the noise ratio is σ². Given t the GP posterior is fixed, so (t, fit) is the
///   draw's whole state. A rejected proposal keeps the previous fit, so neighbouring draws share one fit and one
///   variance array.</item>
///   <item>For an unknown scale, the state is (t, η, s²). s² is drawn exactly given (t, η) for each draw, so each draw
///   has its own variance array.</item>
/// </list>
/// </summary>
public sealed record BgpDraw(double Bandwidth, GpFit Fit, double[] ConditionalVariance, double SignalVariance)
{
    /// <summary>The noise-to-signal ratio η = σ²/s² (σ² itself for a known σ²).</summary>
    public double NoiseRatio => Fit.NoiseRatio;

    /// <summary>The noise variance σ² = s²η.</summary>
    public double NoiseVariance => SignalVariance * Fit.NoiseRatio;
}

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
/// <param name="NoiseRatioConvergence">For an unknown scale, the noise ratio η's rank-normalized diagnostics (it is
/// monitored beside t and gates an early stop); null for a known σ².</param>
public sealed record BgpResult(
    double[] PosteriorMean,
    double[] PosteriorSd,
    BandwidthSummary Bandwidth, BandwidthSummary BandwidthOriginal,
    double AcceptanceRate, int Draws,
    ConvergenceDiagnostic BandwidthConvergence, int[] ChainSeeds, double InputScale,
    IReadOnlyList<IReadOnlyList<BgpDraw>> BandwidthDraws,
    EnsembleCheckpoint? Checkpoint,
    ConvergenceDiagnostic? NoiseRatioConvergence)
{
    /// <summary>The noise variance σ²'s posterior median and 95% interval (a point for a known σ²).</summary>
    public PosteriorInterval NoiseVarianceInterval => Interval(d => d.NoiseVariance);

    /// <summary>The signal variance s²'s posterior median and 95% interval (1 for a known σ²).</summary>
    public PosteriorInterval SignalVarianceInterval => Interval(d => d.SignalVariance);

    /// <summary>The noise ratio η's posterior median and 95% interval (σ² for a known σ²).</summary>
    public PosteriorInterval NoiseRatioInterval => Interval(d => d.NoiseRatio);

    private PosteriorInterval Interval(Func<BgpDraw, double> value)
    {
        var sample = new double[Draws];
        int s = 0;
        foreach (IReadOnlyList<BgpDraw> chain in BandwidthDraws)
            foreach (BgpDraw d in chain) sample[s++] = value(d);
        return PosteriorInterval.Of(sample);
    }
}

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
    private readonly double _logEtaLower;   // log-uniform prior on the noise ratio η (unknown scale only)
    private readonly double _logEtaUpper;
    private readonly bool _unknownScale;

    // Per-run context (set at the top of Run, read by StartChain during the synchronous fan-out).
    private double[,] _xTest = null!;
    private double _t0;
    private double _proposalSd;
    private double _eta0;
    private double _etaProposalSd;
    private int _testCount;

    /// <param name="gp">The GP: a known σ², or an unknown scale (<see cref="GpRegression.WithUnknownScale"/>).</param>
    /// <param name="prior">The empirical-Bayes bandwidth prior.</param>
    /// <param name="noiseRatioLower">For an unknown scale, the lower end of the log-uniform prior on η = σ²/s². It also
    /// floors the nugget, which keeps Σ_t + ηI well conditioned.</param>
    /// <param name="noiseRatioUpper">For an unknown scale, the upper end of η's prior.</param>
    public BgpSampler(GpRegression gp, EmpiricalBayesBandwidthPrior prior, double noiseRatioLower = 1e-6, double noiseRatioUpper = 1e2)
    {
        ArgumentNullException.ThrowIfNull(gp);
        ArgumentNullException.ThrowIfNull(prior);
        if (!(noiseRatioLower > 0.0 && noiseRatioUpper > noiseRatioLower && double.IsFinite(noiseRatioUpper)))
            throw new ArgumentOutOfRangeException(nameof(noiseRatioLower), "The noise-ratio prior needs 0 < lower < upper < ∞.");
        _gp = gp; _prior = prior;
        _logEtaLower = Math.Log(noiseRatioLower);
        _logEtaUpper = Math.Log(noiseRatioUpper);
        _unknownScale = gp.NoiseModel == GpNoiseModel.UnknownScale;
    }

    /// <summary>The R̂/ESS functionals: the bandwidth t, and for an unknown scale the noise ratio η.</summary>
    public int FunctionalDim => _unknownScale ? 2 : 1;

    /// <summary>ESS on each functional.</summary>
    public int MonitoredDim => FunctionalDim;

    /// <summary>Warm-start one chain at the starting bandwidth (and noise ratio). <paramref name="readoutSeed"/> drives the
    /// signal variance's exact draws for an unknown scale; a known σ² never uses it.</summary>
    public IChainRun<double> StartChain(int chainSeed, int readoutSeed)
    {
        var chain = new BgpBandwidthChain(this, new Xoshiro256PlusPlus(chainSeed));
        return new BgpChainRun(chain, this, new Xoshiro256PlusPlus(readoutSeed), new List<BgpDraw>());
    }

    /// <summary>The run's fingerprint for checkpoints: the GP's data, noise model, noise variance, kernel and input map,
    /// the bandwidth prior, the noise-ratio prior, the test points, the starting values and the proposal scales.</summary>
    public string Fingerprint
    {
        get
        {
            var f = new RunFingerprint().Add("BGP");
            _gp.AddTo(f);
            _prior.AddTo(f);
            foreach (double v in _xTest) f.Add(v);
            return f.Add(_testCount).Add(_t0).Add(_proposalSd).Add(_logEtaLower).Add(_logEtaUpper).Add(_eta0)
                    .Add(_etaProposalSd).ToString();
        }
    }

    /// <summary>Rebuild one chain from its restart state: its streams, its log bandwidth (and log noise ratio), the log
    /// target there, its counters and its retained draws. The current fit is recomputed at the bandwidth and nugget it
    /// was made at, which gives the same fit.</summary>
    public IChainRun<double> RestoreChain(IChainState state)
    {
        if (state is not BgpChainState s) throw new ArgumentException("Not a BGP chain's restart state.", nameof(state));
        var chain = new BgpBandwidthChain(this, s);
        var readout = Xoshiro256PlusPlus.Restore(s.Readout.Item1, s.Readout.Item2, s.Readout.Item3, s.Readout.Item4);
        return new BgpChainRun(chain, this, readout, new List<BgpDraw>(s.Draws));
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
    /// <param name="noiseRatioProposalSd">For an unknown scale, the random-walk standard deviation on log η, proposed
    /// jointly with log t.</param>
    /// <param name="noiseRatioInit">For an unknown scale, the starting η inside its prior (default: the geometric middle of
    /// the prior's range); shared by every chain.</param>
    public BgpResult Run(double[,] xTest, int draws = 600, int burn = 200, double proposalSd = 0.4,
                         int seed = 0, double? tInit = null, int chains = 2,
                         double rHatTarget = 0.0, int batchSize = 0, double minEssPerChain = ChainEnsemble.MinEssPerChain,
                         bool checkpoint = false, EnsembleCheckpoint? resumeFrom = null,
                         double noiseRatioProposalSd = 0.4, double? noiseRatioInit = null)
    {
        ArgumentNullException.ThrowIfNull(xTest);
        if (draws < 1) throw new ArgumentOutOfRangeException(nameof(draws));
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));

        double t0 = tInit ?? Math.Sqrt(_prior.LowerBound);   // geometric mid of (γ₁T_n², 1], inside any valid support
        if (!(t0 > _prior.LowerBound && t0 <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(tInit), t0, "The starting bandwidth must lie in the prior support (γ₁T_n², 1].");
        double eta0 = double.NaN;
        if (_unknownScale)
        {
            eta0 = noiseRatioInit ?? Math.Exp(0.5 * (_logEtaLower + _logEtaUpper));
            if (!(Math.Log(eta0) >= _logEtaLower && Math.Log(eta0) <= _logEtaUpper))
                throw new ArgumentOutOfRangeException(nameof(noiseRatioInit), eta0, "The starting noise ratio must lie in its prior's range.");
        }

        _xTest = _gp.Transform.MapInputs(xTest);   // mapped once; every draw predicts in the fitted frame
        _t0 = t0;
        _proposalSd = proposalSd;
        _eta0 = eta0;
        _etaProposalSd = _unknownScale ? noiseRatioProposalSd : double.NaN;
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
                             _gp.Transform.InputScale, bandwidthDraws, summary.Checkpoint,
                             _unknownScale ? summary.Monitored[1] : null);
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

    /// <summary>One chain's restart state: both streams' words, its log bandwidth (and log noise ratio) with the log target
    /// there, the bandwidth and nugget its current fit was made at, its counters and its retained draws.</summary>
    private sealed record BgpChainState(
        (ulong, ulong, ulong, ulong) Rng, (ulong, ulong, ulong, ulong) Readout, double LogBandwidth, double LogNoiseRatio,
        double LogTarget, double FitBandwidth, double FitNoiseRatio, long Accepted, long Attempts, List<BgpDraw> Draws)
        : IChainState;

    /// <summary>
    /// One random-walk Metropolis chain. For a known σ² it walks log t: each <see cref="Step"/> proposes
    /// u' = u + σ·N(0,1) and accepts under the marginal-evidence × prior ratio with the +u log-t Jacobian preserved. For
    /// an unknown scale it walks (log t, log η) jointly, one refactor per step: η's log-uniform prior is flat in log η,
    /// so it adds no Jacobian, and a proposal outside its range is refused. It returns the current bandwidth;
    /// <see cref="CurrentFit"/> is the accepted GP fit the readout reuses.
    /// </summary>
    private sealed class BgpBandwidthChain : IChain<double>
    {
        private readonly BgpSampler _o;
        private readonly Xoshiro256PlusPlus _rng;
        private double _u;          // log t
        private double _v;          // log η (unknown scale only)
        private double _cur;        // log target at the current state
        private GpFactoredFit _factored;   // the current fit with its factor, held only while current

        public long Accepted { get; private set; }
        public long Attempts { get; private set; }

        /// <summary>The GP fit at the current accepted state — the readout's conditional-mean source.</summary>
        public GpFit CurrentFit => _factored.Fit;

        /// <summary>The current fit with its factor — the readout's conditional-variance source.</summary>
        public GpFactoredFit CurrentFactored => _factored;

        public BgpBandwidthChain(BgpSampler owner, Xoshiro256PlusPlus rng)
        {
            _o = owner; _rng = rng;
            double t0 = owner._t0;
            _u = Math.Log(t0);
            if (owner._unknownScale)
            {
                _v = Math.Log(owner._eta0);
                _factored = owner._gp.FitFactored(t0, owner._eta0);
            }
            else
            {
                _factored = owner._gp.FitFactored(t0);
            }
            _cur = _factored.Fit.LogMarginal + owner._prior.LogDensity(t0) + _u;   // +u: log t Jacobian
        }

        /// <summary>A chain continuing from its restart state. The fit is recomputed at the bandwidth and nugget it was made
        /// at — the start's t0 can differ from exp(log t0) in the last bit — so it is the same fit.</summary>
        public BgpBandwidthChain(BgpSampler owner, BgpChainState s)
        {
            _o = owner;
            _rng = Xoshiro256PlusPlus.Restore(s.Rng.Item1, s.Rng.Item2, s.Rng.Item3, s.Rng.Item4);
            _u = s.LogBandwidth;
            _v = s.LogNoiseRatio;
            _cur = s.LogTarget;
            _factored = owner._unknownScale ? owner._gp.FitFactored(s.FitBandwidth, s.FitNoiseRatio)
                                            : owner._gp.FitFactored(s.FitBandwidth);
            Accepted = s.Accepted;
            Attempts = s.Attempts;
        }

        /// <summary>The chain's part of its restart state, with the readout stream's words and the retained draws.</summary>
        public BgpChainState Save((ulong, ulong, ulong, ulong) readout, List<BgpDraw> draws) =>
            new(_rng.SaveState(), readout, _u, _v, _cur, _factored.Fit.Bandwidth, _factored.Fit.NoiseRatio, Accepted, Attempts,
                new List<BgpDraw>(draws));

        public double Step()
        {
            double uProp = _u + _o._proposalSd * Normal.SampleStandard(_rng);
            double vProp = _o._unknownScale ? _v + _o._etaProposalSd * Normal.SampleStandard(_rng) : _v;
            double tProp = Math.Exp(uProp);
            double priorProp = _o._prior.LogDensity(tProp);
            if (_o._unknownScale && !(vProp >= _o._logEtaLower && vProp <= _o._logEtaUpper)) priorProp = double.NegativeInfinity;
            Attempts++;

            GpFactoredFit? fitProp = null;
            double prop = double.NegativeInfinity;
            if (!double.IsNegativeInfinity(priorProp))
            {
                fitProp = _o._unknownScale ? _o._gp.FitFactored(tProp, Math.Exp(vProp)) : _o._gp.FitFactored(tProp);
                prop = fitProp.Fit.LogMarginal + priorProp + uProp;
            }

            if (prop - _cur >= 0.0 || _rng.NextDouble() < Math.Exp(prop - _cur))
            {
                _u = uProp; _v = vProp; _cur = prop; _factored = fitProp!; Accepted++;   // prop is −∞ without a fit, so never accepted
            }
            return Math.Exp(_u);
        }
    }

    /// <summary>
    /// One chain and its retained draws. Per draw it keeps the bandwidth, the GP fit the chain held there, the signal
    /// variance and the conditional variance at the test points. The unit-scale variance is computed once per distinct
    /// fit while the fit's factor exists. For an unknown scale, s² | t, η, y ~ IG((n − 1)/2, Q/2) is drawn per draw on the
    /// readout stream, which scales that variance, and η is the second functional. The run's summaries are read from
    /// these draws.
    /// </summary>
    private sealed class BgpChainRun : IResumableChainRun<double>
    {
        private readonly BgpBandwidthChain _chain;
        private readonly BgpSampler _o;
        private readonly Xoshiro256PlusPlus _readout;   // the signal variance's draws (unknown scale only)
        private GpFit? _varianceFit;          // the fit _variance belongs to
        private double[] _variance = [];

        public readonly List<BgpDraw> Draws;   // every post-burn draw, in order

        public BgpChainRun(BgpBandwidthChain chain, BgpSampler owner, Xoshiro256PlusPlus readout, List<BgpDraw> draws)
        {
            _chain = chain; _o = owner; _readout = readout; Draws = draws;
        }

        public IChainState SaveState() => _chain.Save(_readout.SaveState(), Draws);

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
                _variance = _o._gp.ConditionalVarianceInFrame(_chain.CurrentFactored, _o._xTest);
                _varianceFit = fit;
            }
            if (!_o._unknownScale)
            {
                Draws.Add(new BgpDraw(draw, fit, _variance, 1.0));
                return;
            }
            functionals[1] = fit.NoiseRatio;
            double s2 = InverseGamma.ShapeScale(0.5 * (_o._gp.Count - 1), 0.5 * fit.Quadratic).Sample(_readout);
            var variance = new double[_variance.Length];
            for (int i = 0; i < variance.Length; i++) variance[i] = s2 * _variance[i];
            Draws.Add(new BgpDraw(draw, fit, variance, s2));
        }
    }
}
