using System;
using System.Collections.Generic;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// The curve's peak (argmax) as a posterior, not a point estimate — the downstream T_c readout. Each draw
/// contributes its exact closed-form peak (reduce-per-draw, never argmax-of-the-mean), pooled into a location
/// with a 95% credible interval and its own rank-normalized diagnostics (<see cref="LocationConvergence"/>).
/// </summary>
public sealed record PeakPosterior(
    double LocationMean, double LocationLo, double LocationHi,
    double HeightMean, ConvergenceDiagnostic LocationConvergence);

/// <summary>
/// One transition resolved from the pooled fields — the per-mode generalization of <see cref="PeakPosterior"/>
/// (the dominant mode, k = 1). <see cref="Location"/> and its credible interval come from the peak-intensity λ(T)
/// mass over the mode's support (epistemic — where the transition is); <see cref="SpanLeft"/>/<see cref="SpanRight"/>
/// come from the coverage π(T) (structural — how wide the feature is); <see cref="Mass"/> = ∫λ over the support
/// (≈ posterior probability of a transition there). <see cref="LocationRHat"/> is the cross-chain R̂ on the mode
/// location, computed matching-free from the per-chain λ histograms (each chain's λ over the support IS its own
/// peak-location distribution there) — ≈1 when the chains agree, &gt;1 a disagreement / multimodality flag,
/// unavailable if fewer than two chains resolved the mode. Modes are ordered by mass, dominant first.
/// </summary>
public sealed record PeakMode(
    double Location, double LocationLo, double LocationHi,
    double SpanLeft, double SpanRight, double Mass, RHatDiagnostic LocationRHat);

/// <summary>
/// What the tempering ladders did while sampling (burn-in excluded), averaged over chains: the round-trip and
/// annealed-restart rates per scan, each adjacent pair's mean rejection rate (cold first), their sum the barrier
/// estimate Λ̂, the round-trip rate Theorem 1 of Syed et al. (2022) would give under ideal local exploration (NaN for
/// the sequential sweep), and each chain's final rungs. Compare round trips only between even-odd schedules on one
/// path; the cold chain's ESS is the neutral measure.
/// </summary>
public sealed record TemperingReadout(
    TemperingPath Path, SwapSchedule Schedule, double RoundTripsPerScan, double AnnealedRestartsPerScan,
    double[] PairRejection, double Barrier, double IdealRoundTripRate, double[][] Ladders);

/// <summary>
/// The pooled posterior fit of a BARS run plus its convergence diagnostics, all on the reported functionals
/// (knot count, fitted grid values) — the quantities you read, monitored for convergence.
/// <see cref="Fit"/> averages each draw's conditional mean curve (exact for this affine readout, and lower in variance
/// than averaging coefficient draws). <see cref="FitSd"/> is the pointwise posterior sd by the law of total variance:
/// the spread of the conditional means across draws plus the mean conditional variance the model reports. Readouts
/// that are nonlinear in the curve (<see cref="Peak"/>, <see cref="PeakCountMean"/>, λ, π, <see cref="PeakModes"/>)
/// use each draw's coefficient draw instead; <see cref="Uncertainty"/> says how exact those draws, and the
/// conditional variance, are.
/// <see cref="RHatFit"/> is the per-grid-point consensus map: where it is elevated the chains disagree on the
/// fit (multimodality / under-determination), where it is ≈1 the fit is settled; <see cref="RHatFitMax"/> is its
/// least favourable reading (<see cref="RHatDiagnostic.Worst"/>). <see cref="Knots"/> holds the knot count's
/// rank-normalized diagnostics (split R̂, bulk and tail ESS), and <see cref="AcceptanceRate"/> counts the sampling
/// phase only. <see cref="Tempering"/> reports the ladders when the run was tempered, and is null otherwise.
/// <see cref="PeakIntensity"/> is the grid-aligned peak-intensity λ(T): the matching-free posterior density
/// of significant transitions (∑ = <see cref="PeakCountMean"/>), its modes the persistent transitions.
/// <see cref="SpanCoverage"/> is the grid-aligned span-coverage π(T): the posterior probability each T lies
/// within some significant peak's FWHM-esque span, counted once per draw where nested spans overlap — its plateaus
/// the structural feature widths, its complement the stable regimes.
/// <see cref="PeakModes"/> resolves each transition from those fields — location + credible interval from λ, span
/// from π — with <see cref="Peak"/> the dominant mode (the k = 1 case).
/// </summary>
public sealed record BarsResult(
    double[] Grid,
    double[] Fit,
    double[] FitSd,
    RHatDiagnostic[] RHatFit,
    double MeanKnots,
    double PeakCountMean,
    double[] PeakIntensity,
    double[] SpanCoverage,
    ConvergenceDiagnostic Knots,
    RHatDiagnostic RHatFitMax,
    PeakPosterior Peak,
    IReadOnlyList<PeakMode> PeakModes,
    int SamplesUsed,
    double AcceptanceRate,
    int[] ChainSeeds,
    IReadOnlyList<MoveStat> MoveStats,
    CoefficientUncertainty Uncertainty,
    TemperingReadout? Tempering);

/// <summary>
/// The BARS client of the kernel-agnostic ensemble shell (<see cref="ChainEnsemble"/>): it knows how to start one
/// reversible-jump replica (overdispersed start, optional adaptive-τ kernel or parallel-tempering ladder, optional
/// scale-mixture weight resampler) and how to reduce each draw to the readout functionals (knot count + fitted
/// grid values). The shell owns the fan-out, the burn → rounds → R̂ early-stop, and the pooled R̂/ESS; this owns
/// the BARS semantics and assembles <see cref="BarsResult"/> from the pooled functionals + its own pooled sinks
/// (peak-intensity λ(T), span-coverage π(T), the per-draw argmax cloud). Parallelism is across chains, never data.
/// </summary>
public sealed class BarsEnsemble : IEnsembleModel<KnotConfig>
{
    private readonly IBasis _basis;
    private readonly IObservationModel _model;
    private readonly IComplexityPrior _prior;
    private readonly IKnotKernel _kernel;
    private readonly IWeightResampler? _resampler;
    private readonly double? _dmgkConstant;     // null ⇒ equal birth/death/relocate weights
    private readonly double? _adaptKernelTau;    // set ⇒ each chain tunes an AdaptiveLocalBetaKernel in burn-in
    private readonly int? _temperLevels;        // set ⇒ each chain is a parallel-tempering ladder (cold = the sample)
    private readonly double _temperBetaMin;
    private readonly TemperingPath _temperPath;
    private readonly SwapSchedule _temperSchedule;
    private readonly bool _temperAdaptLadder;

    // Per-run context, set at the top of Run and read (only) by StartChain/the chain runs during that call. Run
    // drives the shell synchronously, so these are valid for the whole run; the fan-out reads them, never writes.
    private double[] _x = null!;
    private double[] _y = null!;
    private double[] _grid = null!;
    private int _weightEvery;
    private int _startDispersion;
    private double _peakProminence;
    private double _spanDropFraction;

    /// <param name="kernel">The knot proposal kernel, shared by every chain and tempering replica, so it must be
    /// immutable. To tune τ, pass <paramref name="adaptKernelTau"/>, which builds a private adaptive kernel per chain.</param>
    /// <param name="dmgkConstant">When set (DMGK uses 0.4), birth/death/relocate use the prior-aware DMGK
    /// schedule (<see cref="DmgkSchedule"/>) instead of constant weights. The prior cancels from the acceptance; dimension mixing is not faster.</param>
    /// <param name="adaptKernelTau">When set, each chain proposes with an <see cref="AdaptiveLocalBetaKernel"/>
    /// started at this τ and tunes it to the relocate target during burn-in (then freezes) — robust to a
    /// mis-scaled proposal. Overrides <paramref name="kernel"/> for proposals.</param>
    /// <param name="temperLevels">When set (≥ 2), each "chain" becomes a <see cref="ParallelTempering{T}"/> ladder of
    /// this many replicas, and the cold (β = 1) replica supplies the draws — so the readout is unchanged but the
    /// chain can cross modal barriers (a multi-peak T_c posterior). Mutually exclusive with the resampler and with
    /// adaptive-τ in this version (a plain fixed-kernel tempered ladder).</param>
    /// <param name="temperBetaMin">The hot end of the geometric β-ladder (default 0.05), in (0, 1).</param>
    /// <param name="temperPath">The ladder's tempering path. The default, the likelihood-tempered path (Q5), keeps the
    /// knot prior whole, so hot replicas approach the prior itself; it needs a model that
    /// <see cref="IObservationModel.CanTemperLikelihood"/>. The powered path tempers prior and likelihood together and
    /// serves every model.</param>
    /// <param name="temperSchedule">Which adjacent pairs each scan attempts to swap; deterministic even-odd by
    /// default (Syed et al. 2022).</param>
    /// <param name="temperAdaptLadder">When true (the default), each ladder re-spaces its interior rungs during burn-in
    /// to equalize adjacent rejection rates, in rounds of doubling length, and then holds them fixed for sampling.
    /// The geometric ladder from <paramref name="temperBetaMin"/> is the start; the endpoints never move.</param>
    public BarsEnsemble(IBasis basis, IObservationModel model, IComplexityPrior prior, IKnotKernel kernel,
                        IWeightResampler? resampler = null, double? dmgkConstant = null, double? adaptKernelTau = null,
                        int? temperLevels = null, double temperBetaMin = 0.05,
                        TemperingPath temperPath = TemperingPath.Likelihood,
                        SwapSchedule temperSchedule = SwapSchedule.DeterministicEvenOdd, bool temperAdaptLadder = true)
    {
        ArgumentNullException.ThrowIfNull(basis);
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(prior);
        ArgumentNullException.ThrowIfNull(kernel);
        if (kernel is AdaptiveLocalBetaKernel)
            throw new ArgumentException(
                "An adaptive kernel holds per-chain state and cannot be shared; request it with adaptKernelTau.", nameof(kernel));
        if (resampler is not null && model is not WeightedNormalModel)
            throw new ArgumentException(
                "Robust weights need the Normal model: only WeightedNormalModel conditions on weights and draws σ².", nameof(resampler));
        if (temperLevels is int lv)
        {
            if (lv < 2) throw new ArgumentOutOfRangeException(nameof(temperLevels), "Tempering needs ≥ 2 levels.");
            if (resampler is not null || adaptKernelTau is not null)
                throw new ArgumentException("Tempering is mutually exclusive with the resampler and adaptive-τ.", nameof(temperLevels));
            if (!(temperBetaMin > 0.0 && temperBetaMin < 1.0))
                throw new ArgumentOutOfRangeException(nameof(temperBetaMin), "The hot end of the ladder must lie in (0, 1).");
            if (temperPath is not (TemperingPath.Powered or TemperingPath.Likelihood))
                throw new ArgumentOutOfRangeException(nameof(temperPath));
            if (temperPath == TemperingPath.Likelihood && !model.CanTemperLikelihood)
                throw new ArgumentException(
                    $"{model.GetType().Name} has no tempered marginal; temper it on the powered path.", nameof(temperPath));
            if (temperSchedule is not (SwapSchedule.DeterministicEvenOdd or SwapSchedule.StochasticEvenOdd or SwapSchedule.Sequential))
                throw new ArgumentOutOfRangeException(nameof(temperSchedule));
        }
        _basis = basis;
        _model = model;
        _prior = prior;
        _kernel = kernel;
        _resampler = resampler;
        _dmgkConstant = dmgkConstant;
        _adaptKernelTau = adaptKernelTau;
        _temperLevels = temperLevels;
        _temperBetaMin = temperBetaMin;
        _temperPath = temperPath;
        _temperSchedule = temperSchedule;
        _temperAdaptLadder = temperAdaptLadder;
    }

    /// <summary>Functionals R̂ runs on: knot count (index 0) + each fitted grid value (the per-grid consensus map).</summary>
    public int FunctionalDim => 1 + _grid.Length;

    /// <summary>Only the knot count (index 0) gets ESS — the per-grid sequences would be a needless O(n²)·g blow-up.</summary>
    public int MonitoredDim => 1;

    /// <summary>
    /// Start one reversible-jump replica from the shell's seeds. <paramref name="chainSeed"/> fans out into three
    /// child seeds: the overdispersed start and the untempered chain, the tempering sub-ladder's own fan-out, and the
    /// weight block. <paramref name="readoutSeed"/> drives the per-sample coefficient draws, decoupled from the chain.
    /// </summary>
    public IChainRun<KnotConfig> StartChain(int chainSeed, int readoutSeed)
    {
        // One child per consumer, so no two generators share a stream. Seeding them all from chainSeed directly
        // would not: Xoshiro and the ladder's SeedTree read the same SplitMix64 sequence, so a 4+-level ladder
        // would duplicate the weight stream.
        int[] seeds = SeedTree.Derive(chainSeed, 3);
        var rng = new Xoshiro256PlusPlus(seeds[0]);
        var run = new BarsChainRun(this, rng, OverdispersedStart(rng, _startDispersion), _grid.Length)
        {
            ReadoutRng = new Xoshiro256PlusPlus(readoutSeed),
            WeightRng = new Xoshiro256PlusPlus(seeds[2]),
            Kernel = _adaptKernelTau is double t0 ? new AdaptiveLocalBetaKernel(t0) : _kernel,
        };
        if (_temperLevels is int levels)
        {
            run.Ladder = NewLadder(run, levels, seeds[1]);
            run.SetChain(run.Ladder.ColdChain);   // the cold (β = 1) replica is the sampler; its draws are the posterior
        }
        else
        {
            run.SetChain(NewChain(run));
        }
        return run;
    }

    /// <remarks>
    /// <paramref name="rHatTarget"/> &gt; 0 stops sampling early once R̂ on knot count (rank-normalized) and every
    /// grid point (classic) meets it (<see cref="RHatDiagnostic.MeetsTarget"/>) and the knot count's bulk and tail
    /// ESS reach <paramref name="minEssPerChain"/> per chain (<see cref="ChainEnsemble.MinEssPerChain"/> by default);
    /// checked after each round of
    /// <paramref name="batchSize"/> samples (a non-positive batch is a single round of <paramref name="samples"/>).
    /// It needs at least two chains.
    /// </remarks>
    public BarsResult Run(double[] x, double[] y, double[] grid,
                          int chains = 4, int masterSeed = 0, int burn = 2000, int samples = 2000,
                          int startDispersion = 4, int weightEvery = 25, double rHatTarget = 0.0, int batchSize = 0,
                          double peakProminence = 0.1, double spanDropFraction = 0.5,
                          double minEssPerChain = ChainEnsemble.MinEssPerChain)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(grid);
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));
        if (samples < 1) throw new ArgumentOutOfRangeException(nameof(samples));

        _x = x;
        _y = y;
        _grid = grid;
        _weightEvery = weightEvery;
        _startDispersion = startDispersion;
        _peakProminence = peakProminence;
        _spanDropFraction = spanDropFraction;

        (EnsembleRun summary, IReadOnlyList<IChainRun<KnotConfig>> chainRuns) =
            new ChainEnsemble().Run<KnotConfig>(this, chains, masterSeed, burn, samples, rHatTarget, batchSize, minEssPerChain);

        int g = grid.Length;
        int collected = summary.SamplesUsed;
        double denom = chains * (double)collected;

        var bars = new BarsChainRun[chains];
        for (int ci = 0; ci < chains; ci++) bars[ci] = (BarsChainRun)chainRuns[ci];

        // Pooled functionals: index 0 = knot count, indices 1..g = fitted grid values.
        var fit = new double[g];
        for (int gi = 0; gi < g; gi++) fit[gi] = summary.FunctionalMean[1 + gi];
        double meanK = summary.FunctionalMean[0];

        // Law of total variance: Var[f(T)] = Var_draws(conditional mean) + E_draws(conditional variance).
        var fitSd = new double[g];
        for (int gi = 0; gi < g; gi++)
        {
            double sv = 0.0;
            for (int ci = 0; ci < chains; ci++) sv += bars[ci].ConditionalVarianceSum[gi];
            fitSd[gi] = Math.Sqrt(summary.FunctionalVariance[1 + gi] + sv / denom);
        }

        // Consensus map: each grid point's R̂ as the shell computed it, with its state; the summary is the worst point.
        var rHatFit = summary.FunctionalRHat[1..];
        RHatDiagnostic rHatFitMax = RHatDiagnostic.Worst(rHatFit);
        ConvergenceDiagnostic knots = summary.Monitored[0];

        // BARS's own pooled sinks (the orchestrator does not know what these mean).
        var peakIntensity = new double[g];
        var spanCoverage = new double[g];
        for (int gi = 0; gi < g; gi++)
        {
            double si = 0.0, sc = 0.0;
            for (int ci = 0; ci < chains; ci++) { si += bars[ci].PeakIntensitySum[gi]; sc += bars[ci].SpanCoverageSum[gi]; }
            peakIntensity[gi] = si / denom;
            spanCoverage[gi] = sc / denom;
        }

        double peakCountMean = 0.0;
        for (int ci = 0; ci < chains; ci++) peakCountMean += bars[ci].PeakCountSum;
        peakCountMean /= denom;

        var moveAtt = new Dictionary<string, long>();
        var moveAcc = new Dictionary<string, long>();
        for (int ci = 0; ci < chains; ci++)
            foreach (MoveStat m in bars[ci].RjChain.MoveStats())
            {
                moveAtt[m.Key] = moveAtt.GetValueOrDefault(m.Key) + m.Attempts;
                moveAcc[m.Key] = moveAcc.GetValueOrDefault(m.Key) + m.Accepted;
            }
        var moveStats = new List<MoveStat>(moveAtt.Count);
        foreach (string key in moveAtt.Keys)
            moveStats.Add(new MoveStat(key, moveAtt[key], moveAcc.GetValueOrDefault(key)));

        PeakPosterior peak = SummarisePeak(bars, denom);
        var chainLambda = new double[chains][];
        for (int ci = 0; ci < chains; ci++) chainLambda[ci] = bars[ci].PeakIntensitySum;
        var peakModes = ExtractModes(grid, peakIntensity, spanCoverage, chainLambda);
        return new BarsResult(grid, fit, fitSd, rHatFit, meanK, peakCountMean, peakIntensity, spanCoverage, knots, rHatFitMax, peak, peakModes, collected, summary.AcceptanceRate, summary.ChainSeeds, moveStats, _model.Uncertainty, Tempering(bars));
    }

    // The ladders' sampling-phase statistics, averaged over chains; null for an untempered run.
    private TemperingReadout? Tempering(BarsChainRun[] states)
    {
        if (_temperLevels is not int levels) return null;
        int pairs = levels - 1, c = states.Length;
        double trips = 0.0, restarts = 0.0;
        var rejection = new double[pairs];
        var ladders = new double[c][];
        foreach (BarsChainRun st in states)
        {
            ParallelTempering<KnotConfig> ladder = st.Ladder!;
            double scans = ladder.Scans - st.BurnScans;
            trips += (ladder.RoundTrips - st.BurnRoundTrips) / scans;
            restarts += (ladder.AnnealedRestarts - st.BurnRestarts) / scans;
            IReadOnlyList<SwapPairStat> stats = ladder.PairStats;
            for (int i = 0; i < pairs; i++) rejection[i] += (stats[i].MeanRejection * ladder.Scans - st.BurnRejection![i]) / scans;
        }
        for (int i = 0; i < pairs; i++) rejection[i] /= c;
        for (int ci = 0; ci < c; ci++) ladders[ci] = states[ci].Ladder!.Betas;
        double ideal = _temperSchedule == SwapSchedule.Sequential ? double.NaN
                     : TemperingLadder.IdealRoundTripRate(rejection, _temperSchedule);
        return new TemperingReadout(_temperPath, _temperSchedule, trips / c, restarts / c, rejection,
                                    TemperingLadder.Barrier(rejection), ideal, ladders);
    }

    private static PeakPosterior SummarisePeak(BarsChainRun[] states, double denom)
    {
        int c = states.Length;
        var traces = new ChainTrace[c];
        var all = new List<double>();
        double mean = 0.0, height = 0.0;
        for (int i = 0; i < c; i++)
        {
            traces[i] = states[i].PeakTrace;
            all.AddRange(states[i].PeakDraws);
            mean += states[i].PeakLocSum;
            height += states[i].PeakHeightSum;
        }
        mean /= denom;
        height /= denom;
        all.Sort();
        double lo = Quantile(all, 0.025);
        double hi = Quantile(all, 0.975);
        return new PeakPosterior(mean, lo, hi, height, ChainDiagnostics.Diagnose(traces, 0));
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

    // Advance the chain one step (the IChain.Step body): step the lone chain or the tempered ladder (returning the
    // cold state), then run the optional scale-mixture weight block. Every weightEvery steps it draws (σ², β) from
    // their exact conditional given the knots and w, redraws w given the residuals and σ², and hands the RJ knot chain
    // the target | w, keeping its state and cumulative move counters. The block has its own seeded stream.
    private KnotConfig StepOnce(BarsChainRun st)
    {
        st.Config = st.Ladder is not null ? st.Ladder.Step() : st.RjChain.Step();

        if (_resampler is not null)
        {
            st.SinceWeight++;
            if (st.SinceWeight >= _weightEvery)
            {
                double[,] z = _basis.Design(st.Config, _x);
                (double[] beta, double sigma2) = ((WeightedNormalModel)_model).SampleConditional(z, _y, st.W, st.WeightRng);
                if (sigma2 > 0.0)
                {
                    var r = new double[_y.Length];
                    for (int i = 0; i < r.Length; i++)
                    {
                        double f = 0.0;
                        for (int j = 0; j < beta.Length; j++) f += z[i, j] * beta[j];
                        r[i] = _y[i] - f;
                    }
                    st.W = _resampler.Resample(r, sigma2, st.WeightRng);
                    st.RjChain.ReplaceTarget(new SplineTarget(_basis, _model, _prior, _x, _y, st.W));
                }
                st.SinceWeight = 0;
            }
        }
        return st.Config;
    }

    // Reduce one post-burn draw. Readouts affine in the curve use the conditional mean: the grid values into the
    // functionals (index 0 = knot count), and the conditional variance into its sink for FitSd. Readouts nonlinear in
    // the curve use a coefficient draw: the closed-form argmax cloud, and per significant peak its location into the
    // grid-aligned intensity accumulator (pooled → λ(T), ∑_grid λ = PeakCountMean) and its FWHM-esque span into the
    // coverage accumulator (pooled → π(T)).
    private void AccumulateDraw(BarsChainRun st, KnotConfig draw, Span<double> functionals)
    {
        int g = _grid.Length;
        double[,] zTrain = _basis.Design(draw, _x);
        double[] mean = _model.PosteriorMeanCoefficients(zTrain, _y, st.W);
        double[,] zGrid = _basis.Design(draw, _grid);
        for (int gi = 0; gi < g; gi++)
        {
            double f = 0.0;
            for (int j = 0; j < mean.Length; j++) f += zGrid[gi, j] * mean[j];
            functionals[1 + gi] = f;
        }
        functionals[0] = draw.Count;

        if (_model.ConditionalVariance(zTrain, _y, st.W, zGrid) is double[] v)
            for (int gi = 0; gi < g; gi++) st.ConditionalVarianceSum[gi] += v[gi];

        if (_basis is SplineBasis splineBasis)
        {
            double[] coef = _model.SampleCoefficients(zTrain, _y, st.W, st.ReadoutRng);
            (double peakX, double peakF) = SplineExtrema.Argmax(draw, coef, splineBasis);
            st.PeakLocSum += peakX;
            st.PeakHeightSum += peakF;
            st.PeakDraws.Add(peakX);
            st.PeakTrace.Add(peakX);

            var spans = SplineExtrema.SignificantPeakSpans(draw, coef, splineBasis, _peakProminence, _spanDropFraction);
            st.PeakCountSum += spans.Count;
            for (int p = 0; p < spans.Count; p++)
                st.PeakIntensitySum[NearestGridIndex(_grid, spans[p].Location)] += 1.0;
            AddSpanCoverage(spans, _grid, st.SpanCoverageSum);
        }
    }

    // π(T) is P(T lies in some significant span), so a draw adds 1 to a grid point at most once, however many spans
    // cover it: a nested peak's span lies inside its neighbour's (R14). λ(T) stays a count, one per span.
    internal static void AddSpanCoverage(IReadOnlyList<PeakSpan> spans, double[] grid, double[] sum)
    {
        for (int gi = 0; gi < grid.Length; gi++)
            for (int p = 0; p < spans.Count; p++)
                if (grid[gi] >= spans[p].Left && grid[gi] <= spans[p].Right) { sum[gi] += 1.0; break; }
    }

    private ReversibleJumpChain<KnotConfig> NewChain(BarsChainRun st)
        => new(Moves(st.Kernel), new SplineTarget(_basis, _model, _prior, _x, _y, st.W), st.Config, st.Rng);

    private ParallelTempering<KnotConfig> NewLadder(BarsChainRun st, int levels, int seed)
    {
        var target = new SplineTarget(_basis, _model, _prior, _x, _y, st.W);
        double[] betas = ParallelTempering<KnotConfig>.GeometricLadder(levels, _temperBetaMin);
        return new ParallelTempering<KnotConfig>(_ => Moves(st.Kernel), target, st.Config, betas, _temperPath, seed,
                                                 _temperSchedule);
    }

    private IRjMove<KnotConfig>[] Moves(IKnotKernel kernel)
    {
        var birth = new KnotBirthMove(kernel);
        var death = new KnotDeathMove(kernel);
        var relocate = new KnotRelocateMove(kernel);
        return _dmgkConstant is double c
            ? DmgkSchedule.Wrap(birth, death, relocate, _prior, c)
            : new IRjMove<KnotConfig>[] { birth, death, relocate };
    }

    // Burn-in with optional adaptive proposal-scale tuning: when the chain's kernel is adaptive, nudge τ toward
    // the relocate target every adaptEvery steps (vanishing Robbins–Monro step), then leave it frozen for the
    // sampling phase (the diminishing-adaptation condition for ergodicity).
    private void AdaptiveBurn(BarsChainRun st, int burn)
    {
        if (st.Ladder is ParallelTempering<KnotConfig> ladder && _temperAdaptLadder)
        {
            LadderBurn(st, ladder, burn);
            return;
        }
        if (st.Kernel is not AdaptiveLocalBetaKernel ak)
        {
            for (int k = 0; k < burn; k++) StepOnce(st);
            return;
        }

        const int adaptEvery = 200;
        long prevAtt = 0, prevAcc = 0;
        int done = 0, round = 0;
        while (done < burn)
        {
            int step = Math.Min(adaptEvery, burn - done);
            for (int k = 0; k < step; k++) StepOnce(st);
            done += step;
            (long att, long acc) = RelocateStat(st.RjChain);   // cumulative across weight resamples
            long dA = att - prevAtt, dC = acc - prevAcc;         // dA = 0 ⇒ no relocate attempted this window
            prevAtt = att; prevAcc = acc;
            if (dA > 0) ak.Adapt((double)dC / dA, stepSize: 3.0 / Math.Sqrt(round + 1));
            round++;
        }
    }

    // Burn-in that tunes the ladder (Syed et al. 2022, Algorithm 4): rounds of 1, 2, 4, … scans, each followed by a
    // re-spacing of the interior rungs from that round's rejection rates, while the burn budget lasts; the rest of the
    // burn runs on the last ladder, which sampling then keeps.
    private void LadderBurn(BarsChainRun st, ParallelTempering<KnotConfig> ladder, int burn)
    {
        int done = 0;
        for (int round = 1; done + round <= burn; round *= 2)
        {
            for (int k = 0; k < round; k++) StepOnce(st);
            done += round;
            ladder.AdaptLadder();
        }
        for (; done < burn; done++) StepOnce(st);
    }

    private static (long Attempts, long Accepted) RelocateStat(ReversibleJumpChain<KnotConfig> chain)
    {
        foreach (MoveStat m in chain.MoveStats())
            if (m.Key == "relocate") return (m.Attempts, m.Accepted);
        return (0, 0);
    }

    // Nearest grid index to x for the peak-intensity deposit (grid ascending — it is the evaluation grid).
    private static int NearestGridIndex(double[] grid, double x)
    {
        int idx = Array.BinarySearch(grid, x);
        if (idx >= 0) return idx;
        int ins = ~idx;                              // first index with grid[index] > x
        if (ins <= 0) return 0;
        if (ins >= grid.Length) return grid.Length - 1;
        return (x - grid[ins - 1]) <= (grid[ins] - x) ? ins - 1 : ins;
    }

    // Resolve transitions from the pooled fields: each local maximum of λ(T) is a mode; its location + credible
    // interval come from λ's mass over the mode's support, its span from the π ≥ ½ region around it. The dominant
    // (highest-mass) mode is the k = 1 case carried by PeakPosterior; this generalizes it to every transition.
    private static List<PeakMode> ExtractModes(double[] grid, double[] lambda, double[] pi, double[][] chainLambda)
    {
        int g = grid.Length;
        var modes = new List<PeakMode>();
        for (int i = 1; i < g - 1; i++)
        {
            if (lambda[i] <= 0.0 || !(lambda[i] > lambda[i - 1] && lambda[i] >= lambda[i + 1])) continue;

            int lo = i;
            while (lo > 0 && lambda[lo - 1] > 0.0 && lambda[lo - 1] <= lambda[lo]) lo--;
            int hi = i;
            while (hi < g - 1 && lambda[hi + 1] > 0.0 && lambda[hi + 1] <= lambda[hi]) hi++;

            double mass = 0.0, wsum = 0.0;
            for (int j = lo; j <= hi; j++) { mass += lambda[j]; wsum += lambda[j] * grid[j]; }
            if (mass <= 0.0) continue;

            double location = wsum / mass;
            double loCi = MassQuantile(grid, lambda, lo, hi, mass, 0.025);
            double hiCi = MassQuantile(grid, lambda, lo, hi, mass, 0.975);

            int sl = i, sr = i;
            while (sl > 0 && pi[sl - 1] >= 0.5) sl--;
            while (sr < g - 1 && pi[sr + 1] >= 0.5) sr++;
            bool covered = pi[i] >= 0.5;

            RHatDiagnostic rhat = ModeRHat(chainLambda, grid, lo, hi);
            modes.Add(new PeakMode(location, loCi, hiCi,
                covered ? grid[sl] : grid[i], covered ? grid[sr] : grid[i], mass, rhat));
        }
        modes.Sort((u, v) => v.Mass.CompareTo(u.Mass));   // dominant first
        return modes;
    }

    // Matching-free between-chain R̂ on a mode's location: each chain's λ over the support [lo,hi] is that chain's
    // own peak-location histogram there, giving its mean and within-chain variance directly. Gelman–Rubin on those
    // (n = mean per-chain peak count) → ≈1 when chains agree. The variance is two-pass and exactly zero for a chain
    // whose mass sits on one grid point, so agreeing point masses read as constant and disagreeing ones as frozen.
    private static RHatDiagnostic ModeRHat(double[][] chainLambda, double[] grid, int lo, int hi)
    {
        int m = chainLambda.Length;
        var mean = new double[m];
        var varc = new double[m];
        double nbar = 0.0;
        int seen = 0;
        for (int c = 0; c < m; c++)
        {
            double cnt = 0.0, s = 0.0;
            int support = 0, only = -1;
            for (int j = lo; j <= hi; j++)
            {
                double w = chainLambda[c][j];
                if (w <= 0.0) continue;
                cnt += w; s += w * grid[j]; support++; only = j;
            }
            if (cnt <= 0.0) continue;
            if (support == 1)
            {
                mean[seen] = grid[only];
                varc[seen] = 0.0;
            }
            else
            {
                double mu = s / cnt, ss = 0.0;
                for (int j = lo; j <= hi; j++)
                {
                    double w = chainLambda[c][j];
                    if (w > 0.0) { double d = grid[j] - mu; ss += w * d * d; }
                }
                mean[seen] = mu;
                varc[seen] = ss / cnt;
            }
            nbar += cnt;
            seen++;
        }
        if (seen < 2) return RHatDiagnostic.Unavailable;   // can't assess convergence from one chain
        nbar = Math.Max(1.0, nbar / seen);
        return ChainDiagnostics.RHat(mean.AsSpan(0, seen), varc.AsSpan(0, seen), nbar);
    }

    // Step-CDF quantile of λ-mass over the support [lo, hi].
    private static double MassQuantile(double[] grid, double[] w, int lo, int hi, double total, double q)
    {
        double target = q * total, cum = 0.0;
        for (int j = lo; j <= hi; j++) { cum += w[j]; if (cum >= target) return grid[j]; }
        return grid[hi];
    }

    private static KnotConfig OverdispersedStart(Xoshiro256PlusPlus rng, int dispersion)
    {
        int k = dispersion <= 0 ? 0 : rng.NextInt(dispersion + 1);
        var knots = new double[k];
        for (int i = 0; i < k; i++) knots[i] = rng.NextDouble();
        Array.Sort(knots);
        return new KnotConfig(knots);
    }

    /// <summary>
    /// One BARS reversible-jump replica + its thread-local accumulators. As an <see cref="IChain{KnotConfig}"/> its
    /// <see cref="Step"/> advances the chain (or tempered ladder) and runs the resampler block, and its acceptance
    /// is the RJ chain's cumulative count; as an <see cref="IChainRun{KnotConfig}"/>
    /// it burns adaptively, reduces a draw to functionals, and carries the BARS sinks (peak cloud, λ/π fields).
    /// </summary>
    private sealed class BarsChainRun : IChain<KnotConfig>, IChainRun<KnotConfig>
    {
        private readonly BarsEnsemble _owner;
        public readonly Xoshiro256PlusPlus Rng;
        public readonly double[] PeakIntensitySum;
        public readonly double[] SpanCoverageSum;
        public readonly double[] ConditionalVarianceSum;
        public readonly List<double> PeakDraws = new();
        public readonly ChainTrace PeakTrace = new(1);   // the same draws, as the diagnostics' ordered trace
        public KnotConfig Config;
        public double[]? W;
        public IKnotKernel Kernel = null!;
        public Xoshiro256PlusPlus ReadoutRng = null!;   // per-draw coefficient draws, decoupled from the chain RNG
        public Xoshiro256PlusPlus WeightRng = null!;    // the robust weight block's (σ², β, w) draws, derived from the chain seed
        public ParallelTempering<KnotConfig>? Ladder;   // set ⇒ the RJ chain is this ladder's cold replica
        public double PeakLocSum;
        public double PeakHeightSum;
        public double PeakCountSum;
        public int SinceWeight;
        public long BurnScans, BurnRoundTrips, BurnRestarts;   // the ladder's counts at the end of burn-in
        public double[]? BurnRejection;                        // its per-pair rejection sums then

        private ReversibleJumpChain<KnotConfig> _rjChain = null!;

        public BarsChainRun(BarsEnsemble owner, Xoshiro256PlusPlus rng, KnotConfig start, int gridSize)
        {
            _owner = owner;
            Rng = rng;
            Config = start;
            PeakIntensitySum = new double[gridSize];
            SpanCoverageSum = new double[gridSize];
            ConditionalVarianceSum = new double[gridSize];
        }

        /// <summary>The RJ chain — for the tempered case, the ladder's cold replica. Set once, at start.</summary>
        public ReversibleJumpChain<KnotConfig> RjChain => _rjChain;
        public void SetChain(ReversibleJumpChain<KnotConfig> chain) => _rjChain = chain;

        // IChain
        public KnotConfig Step() => _owner.StepOnce(this);
        public long Accepted => _rjChain.Accepted;
        public long Attempts => _rjChain.Attempts;

        // IChainRun
        public IChain<KnotConfig> Chain => this;
        public void Burn(int steps)
        {
            _owner.AdaptiveBurn(this, steps);
            if (Ladder is not ParallelTempering<KnotConfig> ladder) return;
            BurnScans = ladder.Scans;
            BurnRoundTrips = ladder.RoundTrips;
            BurnRestarts = ladder.AnnealedRestarts;
            IReadOnlyList<SwapPairStat> stats = ladder.PairStats;
            BurnRejection = new double[stats.Count];
            for (int i = 0; i < stats.Count; i++) BurnRejection[i] = ladder.Scans > 0 ? stats[i].MeanRejection * ladder.Scans : 0.0;
        }
        public void Accumulate(in KnotConfig draw, Span<double> functionals) => _owner.AccumulateDraw(this, draw, functionals);
    }
}
