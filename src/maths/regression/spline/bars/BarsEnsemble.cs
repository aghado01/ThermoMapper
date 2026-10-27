using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
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
/// One retained draw's exact peak readout, read off its coefficient draw: the closed-form global argmax (location and
/// height) and the significant peaks' spans in location order (<see cref="PeakSpan"/>: location, height, the two
/// crossings and both clipping flags). Nothing is snapped to the grid. A basis without closed-form peaks gives a NaN
/// argmax and no spans.
/// </summary>
public sealed record PeakDraw(double ArgmaxLocation, double ArgmaxHeight, IReadOnlyList<PeakSpan> Spans);

/// <summary>
/// One retained BARS draw's state: the knot configuration the chain held, the coefficients and noise variance σ² drawn
/// from their conditional posterior given those knots (and the robust weights), the robust weights in force at the
/// draw, and the draw's exact peak readout, read off its coefficients.
/// <list type="bullet">
///   <item><see cref="Coefficients"/> and <see cref="Sigma2"/> are one joint draw on the readout stream
///   (<see cref="IObservationModel.SampleConditional"/>): exact for the Normal model, and otherwise as exact as
///   <see cref="BarsResult.Uncertainty"/> says. σ² is NaN for a model without a noise scale.</item>
///   <item><see cref="Weights"/> is null without the robust weight block. The block replaces the array at each update,
///   so the draws between two updates share one array.</item>
///   <item>The conditional mean curve behind <see cref="BarsResult.Fit"/> is recomputable from <see cref="Knots"/>
///   and <see cref="Weights"/> (<see cref="IObservationModel.PosteriorMeanCoefficients"/>), so it is not stored.</item>
/// </list>
/// </summary>
public sealed record BarsDraw(KnotConfig Knots, double[] Coefficients, double Sigma2, double[]? Weights, PeakDraw Peak);

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
/// <see cref="Evaluations"/> counts every evaluated knot configuration of every chain and ladder rung, burn-in
/// included, by admission outcome (<see cref="DesignRejection.None"/> for admitted), so the truncation the support
/// rules and admission impose stays visible. <see cref="Checkpoint"/> is the run's restart state when it was requested
/// (<see cref="BarsEnsemble.Run"/>'s <c>checkpoint</c>), and null otherwise.
/// <see cref="PeakIntensity"/> is the grid-aligned peak-intensity λ(T): the matching-free posterior density
/// of significant transitions (∑ = <see cref="PeakCountMean"/>), its modes the persistent transitions.
/// <see cref="SpanCoverage"/> is the grid-aligned span-coverage π(T): the posterior probability each T lies
/// within some significant peak's FWHM-esque span, counted once per draw where nested spans overlap — its plateaus
/// the structural feature widths, its complement the stable regimes.
/// <see cref="PeakModes"/> resolves each transition from those fields — location + credible interval from λ, span
/// from π — with <see cref="Peak"/> the dominant mode (the k = 1 case).
/// <see cref="Draws"/> keeps every retained draw's state (<see cref="BarsDraw"/>), per chain in draw order: its knots,
/// its coefficient and σ² draw, its robust weights, and its exact peak readout (the closed-form locations, heights,
/// crossings and clipping flags that λ and π count on the grid).
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
    IReadOnlyList<IReadOnlyList<BarsDraw>> Draws,
    int SamplesUsed,
    double AcceptanceRate,
    int[] ChainSeeds,
    IReadOnlyList<MoveStat> MoveStats,
    CoefficientUncertainty Uncertainty,
    TemperingReadout? Tempering,
    Tally<DesignRejection> Evaluations,
    EnsembleCheckpoint? Checkpoint);

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
    private readonly bool _temperReferenceRung;
    private readonly double _minColumnMass;
    private readonly double _birthExploration;
    private readonly BirthBase _birthBase;
    private IKnotBase _knotBase = UniformKnotBase.Instance;

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
    /// <param name="temperBetaMin">The hot end of the geometric β-ladder (default 0.05), in (0, 1); with
    /// <paramref name="temperReferenceRung"/>, the lowest rung above the reference.</param>
    /// <param name="temperPath">The ladder's tempering path. The default, the likelihood-tempered path (Q5), keeps the
    /// knot prior whole, so hot replicas approach the prior itself; it needs a model that
    /// <see cref="IObservationModel.CanTemperLikelihood"/>. The powered path tempers prior and likelihood together and
    /// serves every model.</param>
    /// <param name="temperSchedule">Which adjacent pairs each scan attempts to swap; deterministic even-odd by
    /// default (Syed et al. 2022).</param>
    /// <param name="temperAdaptLadder">When true (the default), each ladder re-spaces its interior rungs during burn-in
    /// to equalize adjacent rejection rates, in rounds of doubling length, and then holds them fixed for sampling.
    /// The geometric ladder from <paramref name="temperBetaMin"/> is the start; the endpoints never move.</param>
    /// <param name="minColumnMass">The support gate: every basis function needs column mass Σᵢ Bⱼ(xᵢ) at least this
    /// (one data point's worth by default; 0 turns it off). A declared restriction of the knot prior's support
    /// (<see cref="SplineTarget.MinColumnMass"/>).</param>
    /// <param name="birthExploration">ε, the share of births drawn from <paramref name="birthBase"/> rather than near
    /// an existing knot (<see cref="KnotBirthDensity"/>). The default 1 draws every birth from the base, which step 3.9's
    /// benchmark found best or tied on every fit; values below 1 restore DMGK's local births in part, and 0 entirely.</param>
    /// <param name="birthBase">The base measure of the exploratory births and of the first knot.</param>
    /// <param name="temperReferenceRung">When true, each ladder ends at β = 0 with an exact reference rung
    /// (<see cref="SplineTarget.Reference"/>, OQ2): its hot replica draws independently from the path's b = 0 member every
    /// scan. It needs at least three levels, and on the likelihood path a model that
    /// <see cref="IObservationModel.CanTemperToPrior"/>. Off by default: at equal work it did not pay at 4 or 8 levels on
    /// unimodal fits or on separated modes, where a ladder that spans to the prior is too short to carry its
    /// regenerations to the cold chain (OQ2, preparation brief).</param>
    public BarsEnsemble(IBasis basis, IObservationModel model, IComplexityPrior prior, IKnotKernel kernel,
                        IWeightResampler? resampler = null, double? dmgkConstant = null, double? adaptKernelTau = null,
                        int? temperLevels = null, double temperBetaMin = 0.05,
                        TemperingPath temperPath = TemperingPath.Likelihood,
                        SwapSchedule temperSchedule = SwapSchedule.DeterministicEvenOdd, bool temperAdaptLadder = true,
                        double minColumnMass = SplineTarget.DefaultMinColumnMass,
                        double birthExploration = 1.0, BirthBase birthBase = BirthBase.Uniform, bool temperReferenceRung = false)
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
            if (temperReferenceRung && lv < 3)
                throw new ArgumentOutOfRangeException(nameof(temperLevels), "A reference rung needs at least three levels.");
            if (temperReferenceRung && temperPath == TemperingPath.Likelihood && !model.CanTemperToPrior)
                throw new ArgumentException(
                    $"{model.GetType().Name} does not serve b = 0 on the likelihood path; use its linear form or the powered path.",
                    nameof(temperReferenceRung));
        }
        else if (temperReferenceRung)
            throw new ArgumentException("A reference rung needs a tempered ladder.", nameof(temperReferenceRung));
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
        _temperReferenceRung = temperReferenceRung;
        if (!(minColumnMass >= 0.0 && double.IsFinite(minColumnMass)))
            throw new ArgumentOutOfRangeException(nameof(minColumnMass), "The minimum column mass must be finite and non-negative.");
        _minColumnMass = minColumnMass;
        if (!(birthExploration >= 0.0 && birthExploration <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(birthExploration), "The exploratory share must lie in [0, 1].");
        if (birthBase is not (BirthBase.Uniform or BirthBase.Design)) throw new ArgumentOutOfRangeException(nameof(birthBase));
        _birthExploration = birthExploration;
        _birthBase = birthBase;
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
            Kernel = _adaptKernelTau is double t0 ? new AdaptiveLocalBetaKernel(t0) : _kernel,
        };
        run.Config = OnSupport(run.Config, rng);
        if (_temperLevels is int levels)
        {
            run.Ladder = NewLadder(run, levels, seeds[1]);
            run.SetChain(run.Ladder.ColdChain);   // the cold (β = 1) replica is the sampler; its draws are the posterior
        }
        else
        {
            run.SetChain(NewChain(run));
        }
        run.Scan = new BlockScan<BarsChainRun>(StepBlocks(), run, new Xoshiro256PlusPlus(seeds[2]), x => x);   // the weight stream
        return run;
    }

    /// <summary>
    /// The run's fingerprint for checkpoints: the data (design points, responses, grid), the run settings, and each
    /// component's type and the ensemble's own settings. The components' internal parameters (a prior's mean, a
    /// kernel's spread) are the caller's to keep the same.
    /// </summary>
    public string Fingerprint => new RunFingerprint().Add("BARS").Add(_x).Add(_y).Add(_grid)
        .Add(_weightEvery).Add(_startDispersion).Add(_peakProminence).Add(_spanDropFraction)
        .Add(_basis.GetType().FullName!).Add(_model.GetType().FullName!).Add(_prior.GetType().FullName!)
        .Add(_kernel.GetType().FullName!).Add(_resampler?.GetType().FullName ?? "")
        .Add(_dmgkConstant ?? double.NaN).Add(_adaptKernelTau ?? double.NaN).Add(_temperLevels ?? 0).Add(_temperBetaMin)
        .Add((int)_temperPath).Add((int)_temperSchedule).Add(_temperAdaptLadder ? 1 : 0).Add(_minColumnMass).Add(_birthExploration).Add((int)_birthBase)
        .Add(_temperReferenceRung ? 1 : 0).Add((_model as WeightedNormalModel)?.Form.ToString() ?? "").ToString();

    /// <summary>
    /// Rebuild one chain from its saved restart state (<see cref="BarsChainState"/>): the chain or ladder, every RNG
    /// stream, the kernel's tuning, the robust weights, and the chain's readouts with its retained draws. Rebuilding evaluates
    /// the current states once; those evaluations are not counted.
    /// </summary>
    public IChainRun<KnotConfig> RestoreChain(IChainState state)
    {
        if (state is not BarsChainState s) throw new ArgumentException("Not a BARS chain's restart state.", nameof(state));
        static Xoshiro256PlusPlus Rng((ulong, ulong, ulong, ulong) w) => Xoshiro256PlusPlus.Restore(w.Item1, w.Item2, w.Item3, w.Item4);

        var run = new BarsChainRun(this, Rng(s.StartRng), s.Config, s.Readout.Copy())
        {
            W = s.Weights is null ? null : (double[])s.Weights.Clone(),
            ReadoutRng = Rng(s.ReadoutRng),
            Kernel = s.KernelLogTau is double logTau ? AdaptiveLocalBetaKernel.FromLogTau(logTau) : _kernel,
            BurnScans = s.BurnScans,
            BurnRoundTrips = s.BurnRoundTrips,
            BurnRestarts = s.BurnRestarts,
            BurnRejection = s.BurnRejection is null ? null : (double[])s.BurnRejection.Clone(),
        };
        var target = new SplineTarget(_basis, _model, _prior, _x, _y, run.W, run.Readout.Evaluations, _minColumnMass);
        if (s.Ladder is LadderSnapshot<KnotConfig> ladder)
        {
            run.Ladder = ParallelTempering<KnotConfig>.Restore(_ => Moves(run.Kernel), target, ladder,
                _temperReferenceRung ? target.Reference(_temperPath) : null);
            run.SetChain(run.Ladder.ColdChain);
        }
        else
        {
            run.SetChain(ReversibleJumpChain<KnotConfig>.Restore(Moves(run.Kernel), target, s.Chain!));
        }
        run.Scan = BlockScan<BarsChainRun>.Restore(StepBlocks(), new BlockScanSnapshot<BarsChainRun>(
            run, s.WeightRng.Item1, s.WeightRng.Item2, s.WeightRng.Item3, s.WeightRng.Item4, s.ScanSteps, s.BlockAttempts,
            s.BlockAccepts), x => x);
        run.Readout.Evaluations.Assign(s.Readout.Evaluations);
        return run;
    }

    // One chain's restart state. The untempered chain shares the start stream, so its snapshot carries those words. The
    // step scan's snapshot holds the live run as its state, so only its stream (the weight stream), step count and block
    // counters are kept.
    private IChainState Save(BarsChainRun st)
    {
        BlockScanSnapshot<BarsChainRun> scan = st.Scan.Snapshot();
        return new BarsChainState(
            st.Config, st.W is null ? null : (double[])st.W.Clone(), scan.Steps, scan.BlockAttempts, scan.BlockAccepts,
            st.Rng.SaveState(), st.ReadoutRng.SaveState(), (scan.Rng0, scan.Rng1, scan.Rng2, scan.Rng3),
            st.Kernel is AdaptiveLocalBetaKernel adaptive ? adaptive.LogTau : null,
            st.Ladder is null ? st.RjChain.Snapshot() : null, st.Ladder?.Snapshot(),
            st.Readout.Copy(), st.BurnScans, st.BurnRoundTrips, st.BurnRestarts,
            st.BurnRejection is null ? null : (double[])st.BurnRejection.Clone());
    }

    /// <summary>
    /// One BARS chain's restart state, beside its readouts: the knot configuration, the robust weights, the step scan's
    /// count and block counters (the weight block's cadence), every RNG stream's words, the adaptive kernel's log τ, the chain's or ladder's snapshot (rungs,
    /// replica assignments, swap stream, round-trip and swap statistics), and the ladder's counts at the end of
    /// burn-in. Taken after burn-in, so burn-in adaptation is never in flight.
    /// </summary>
    private sealed record BarsChainState(
        KnotConfig Config, double[]? Weights, long ScanSteps, long[] BlockAttempts, long[] BlockAccepts,
        (ulong, ulong, ulong, ulong) StartRng, (ulong, ulong, ulong, ulong) ReadoutRng, (ulong, ulong, ulong, ulong) WeightRng,
        double? KernelLogTau, ChainSnapshot<KnotConfig>? Chain, LadderSnapshot<KnotConfig>? Ladder,
        BarsReadout Readout, long BurnScans, long BurnRoundTrips, long BurnRestarts,
        double[]? BurnRejection) : IChainState;

    /// <remarks>
    /// <paramref name="rHatTarget"/> &gt; 0 stops sampling early once R̂ on knot count (rank-normalized) and every
    /// grid point (classic) meets it (<see cref="RHatDiagnostic.MeetsTarget"/>) and the knot count's bulk and tail
    /// ESS reach <paramref name="minEssPerChain"/> per chain (<see cref="ChainEnsemble.MinEssPerChain"/> by default);
    /// checked after each round of
    /// <paramref name="batchSize"/> samples (a non-positive batch is a single round of <paramref name="samples"/>).
    /// It needs at least two chains.
    /// <para>With <paramref name="checkpoint"/>, the result carries the run's restart state; with
    /// <paramref name="resumeFrom"/>, the run continues that checkpoint's chains for <paramref name="samples"/> more
    /// draws, exactly as the uninterrupted run would have, and its result covers every draw since burn-in. Resuming
    /// needs the same data and settings (the run's fingerprint) and the same master seed and chain count.</para>
    /// </remarks>
    public BarsResult Run(double[] x, double[] y, double[] grid,
                          int chains = 4, int masterSeed = 0, int burn = 2000, int samples = 2000,
                          int startDispersion = 4, int weightEvery = 25, double rHatTarget = 0.0, int batchSize = 0,
                          double peakProminence = 0.1, double spanDropFraction = 0.5,
                          double minEssPerChain = ChainEnsemble.MinEssPerChain,
                          bool checkpoint = false, EnsembleCheckpoint? resumeFrom = null)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(grid);
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));
        if (samples < 1) throw new ArgumentOutOfRangeException(nameof(samples));

        _x = x;
        _knotBase = _birthBase == BirthBase.Design ? new DesignHistogramBase(x) : UniformKnotBase.Instance;
        _y = y;
        _grid = grid;
        _weightEvery = weightEvery;
        _startDispersion = startDispersion;
        _peakProminence = peakProminence;
        _spanDropFraction = spanDropFraction;
        RequireAdmissibleZeroKnots();

        (EnsembleRun summary, IReadOnlyList<IChainRun<KnotConfig>> chainRuns) =
            new ChainEnsemble().Run<KnotConfig>(this, chains, masterSeed, burn, samples, rHatTarget, batchSize, minEssPerChain,
                                                checkpoint, resumeFrom);

        int g = grid.Length;
        int collected = summary.SamplesUsed;

        var bars = new BarsChainRun[chains];
        for (int ci = 0; ci < chains; ci++) bars[ci] = (BarsChainRun)chainRuns[ci];

        // Merge the chains' readouts once; the reductions below read the merged one.
        var pooled = new BarsReadout(g);
        foreach (BarsChainRun st in bars) pooled.MergeFrom(st.Readout);

        // Pooled functionals: index 0 = knot count, indices 1..g = fitted grid values.
        var fit = new double[g];
        for (int gi = 0; gi < g; gi++) fit[gi] = summary.FunctionalMean[1 + gi];
        double meanK = summary.FunctionalMean[0];

        // Law of total variance: Var[f(T)] = Var_draws(conditional mean) + E_draws(conditional variance).
        var fitSd = new double[g];
        WeightedSums cv = pooled.ConditionalVariance;
        for (int gi = 0; gi < g; gi++)
            fitSd[gi] = Math.Sqrt(summary.FunctionalVariance[1 + gi] + (cv.TotalWeight > 0.0 ? cv.Mean(gi) : 0.0));

        // Consensus map: each grid point's R̂ as the shell computed it, with its state; the summary is the worst point.
        var rHatFit = summary.FunctionalRHat[1..];
        RHatDiagnostic rHatFitMax = RHatDiagnostic.Worst(rHatFit);
        ConvergenceDiagnostic knots = summary.Monitored[0];

        // λ(T) and π(T) are views over the retained draws: each chain's grid counts, pooled over chains and divided by
        // the draw count. The per-chain λ counts are each chain's peak-location histogram, which ModeRHat reads.
        var chainLambda = new double[chains][];
        var peakIntensity = new double[g];
        var spanCoverage = new double[g];
        double draws = 0.0, spans = 0.0;
        for (int ci = 0; ci < chains; ci++)
        {
            (double[] intensity, double[] coverage) = PeakCounts(bars[ci].Readout.Draws, grid);
            chainLambda[ci] = intensity;
            for (int gi = 0; gi < g; gi++)
            {
                peakIntensity[gi] += intensity[gi];
                spanCoverage[gi] += coverage[gi];
            }
            draws += bars[ci].Readout.Draws.Count;
            foreach (BarsDraw d in bars[ci].Readout.Draws) spans += d.Peak.Spans.Count;
        }
        for (int gi = 0; gi < g; gi++)
        {
            peakIntensity[gi] /= draws;
            spanCoverage[gi] /= draws;
        }
        double peakCountMean = spans / draws;

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

        PeakPosterior peak = SummarisePeak(bars, _basis is SplineBasis);
        var peakModes = ExtractModes(grid, peakIntensity, spanCoverage, chainLambda);
        var retained = new IReadOnlyList<BarsDraw>[chains];
        for (int ci = 0; ci < chains; ci++) retained[ci] = bars[ci].Readout.Draws.ToArray();
        return new BarsResult(grid, fit, fitSd, rHatFit, meanK, peakCountMean, peakIntensity, spanCoverage, knots, rHatFitMax, peak, peakModes, retained, collected, summary.AcceptanceRate, summary.ChainSeeds, moveStats, _model.Uncertainty, Tempering(bars), pooled.Evaluations, summary.Checkpoint);
    }

    /// <summary>
    /// log p(y_unit | draw) for every retained draw of a result this ensemble produced and every unit (OQ17), on the data it
    /// fitted. A draw's knots, coefficients and σ² are global, so every unit conditions on them. A robust weight belongs to
    /// one observation, so it is integrated out under its prior: the resampler's marginal (Student t, or the
    /// contaminated-normal mixture) replaces the Normal term N(yᵢ; zᵢβ, σ²). The Poisson model's term is
    /// Poisson(yᵢ; exp(zᵢβ)); the Whittle model's is the exponential density of a periodogram ordinate with mean exp(zᵢβ),
    /// under Whittle's approximation of independent ordinates. The Gibbs-loss model declines: its likelihood is a
    /// pseudo-likelihood, exp(−η·loss), not a density of a new observation.
    /// <para>
    /// <see cref="PointwiseContent.CollapsedPredictive"/> keeps each draw's knots and integrates its coefficients and σ²
    /// out given the observations outside the unit: ln p(y_U | y_rest, knots), the g-prior's leave-unit-out Student t
    /// (<see cref="IConditionalFit.LeaveOutLogDensities"/>). The Normal model has it without robust weights. With them it
    /// would condition on the unit's own fitted weights, which biases leave-one-out, and the other models have no closed
    /// form; all of these decline.
    /// </para>
    /// </summary>
    /// <param name="units">The partition of the n observations; null takes each observation as its own unit.</param>
    /// <param name="content">The draw content (<see cref="PointwiseContent"/>); the retained parameters by default.</param>
    public PointwiseLogDensity PointwiseLogDensity(BarsResult result, double[] x, double[] y, UnitPartition? units = null,
                                                   PointwiseContent content = PointwiseContent.RetainedParameters)
    {
        ArgumentNullException.ThrowIfNull(result);
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(y);
        if (x.Length != y.Length) throw new ArgumentException("One response per design point is required.", nameof(y));
        UnitPartition partition = units ?? UnitPartition.EachObservation(y.Length);
        if (partition.Observations != y.Length)
            throw new ArgumentException("The partition must cover the fitted observations.", nameof(units));
        if (content is not (PointwiseContent.RetainedParameters or PointwiseContent.CollapsedPredictive))
            throw new ArgumentOutOfRangeException(nameof(content));
        if (content == PointwiseContent.CollapsedPredictive) return CollapsedPointwise(result, x, y, partition);
        Func<double, double, double, double> term = _model switch   // (yᵢ, ηᵢ = zᵢβ, σ²) → log density
        {
            WeightedNormalModel when _resampler is IWeightResampler resampler
                => (yi, eta, s2) => resampler.MarginalLogDensity(yi - eta, s2),
            WeightedNormalModel => (yi, eta, s2) => -0.5 * Math.Log(2.0 * Math.PI * s2) - 0.5 * (yi - eta) * (yi - eta) / s2,
            PoissonModel => (yi, eta, _) => yi * eta - Math.Exp(eta) - SpecialFunctions.LogGamma(yi + 1.0),
            WhittleModel => (yi, eta, _) => -eta - yi * Math.Exp(-eta),
            _ => throw new NotSupportedException(
                $"{_model.GetType().Name} has no pointwise density: its likelihood is a pseudo-likelihood, not a density of a new observation."),
        };

        int n = y.Length, draws = 0;
        foreach (IReadOnlyList<BarsDraw> chain in result.Draws) draws += chain.Count;
        var values = new double[draws, partition.Units];
        var terms = new double[n];
        int row = 0;
        foreach (IReadOnlyList<BarsDraw> chain in result.Draws)
            foreach (BarsDraw d in chain)
            {
                double[,] z = _basis.Design(d.Knots, x);
                for (int i = 0; i < n; i++)
                {
                    double eta = 0.0;
                    for (int j = 0; j < d.Coefficients.Length; j++) eta += z[i, j] * d.Coefficients[j];
                    terms[i] = term(y[i], eta, d.Sigma2);
                }
                for (int u = 0; u < partition.Units; u++)
                {
                    double sum = 0.0;
                    foreach (int i in partition.Members[u]) sum += terms[i];
                    values[row, u] = sum;
                }
                row++;
            }

        var latents = new List<LatentReport>
        {
            new("knot configuration", LatentScope.Global, LatentTreatment.Conditioned),
            new("coefficients β", LatentScope.Global, LatentTreatment.Conditioned),
        };
        if (_model is WeightedNormalModel)
        {
            latents.Add(new LatentReport("noise variance σ²", LatentScope.Global, LatentTreatment.Conditioned));
            if (_resampler is not null) latents.Add(new LatentReport("robust weights wᵢ", LatentScope.Observation, LatentTreatment.Integrated));
        }
        return new PointwiseLogDensity(values, partition, latents, ChainLengths.Of(result.Draws));
    }

    // The collapsed content: each draw's knots refitted, and its fit's leave-unit-out predictive per unit.
    private PointwiseLogDensity CollapsedPointwise(BarsResult result, double[] x, double[] y, UnitPartition partition)
    {
        if (_model is not WeightedNormalModel)
            throw new NotSupportedException($"{_model.GetType().Name} has no closed-form collapsed predictive; use the retained parameters.");
        if (_resampler is not null)
            throw new NotSupportedException(
                "With robust weights the collapsed predictive would condition on each unit's own fitted weights; use the retained parameters, which integrate them.");
        int draws = 0;
        foreach (IReadOnlyList<BarsDraw> chain in result.Draws) draws += chain.Count;
        var values = new double[draws, partition.Units];
        int row = 0;
        foreach (IReadOnlyList<BarsDraw> chain in result.Draws)
            foreach (BarsDraw d in chain)
            {
                double[,] z = _basis.Design(d.Knots, x);
                double[] terms = _model.FitConditional(z, y, null).LeaveOutLogDensities(partition)
                                 ?? throw new NotSupportedException($"{_model.GetType().Name} returned no collapsed predictive.");
                for (int u = 0; u < partition.Units; u++) values[row, u] = terms[u];
                row++;
            }
        var latents = new List<LatentReport>
        {
            new("knot configuration", LatentScope.Global, LatentTreatment.Conditioned),
            new("coefficients β", LatentScope.Global, LatentTreatment.IntegratedGivenOthers),
            new("noise variance σ²", LatentScope.Global, LatentTreatment.IntegratedGivenOthers),
        };
        return new PointwiseLogDensity(values, partition, latents, ChainLengths.Of(result.Draws));
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

    // The peak location's mean, height and credible interval, and its diagnostics, read from the retained draws' argmax:
    // each chain's draws in order form its trace, and the pooled ones the sample. A basis without closed-form peaks
    // (closedForm false) has no argmax, so every field is NaN or unavailable rather than a location of 0.
    private static PeakPosterior SummarisePeak(BarsChainRun[] states, bool closedForm)
    {
        var traces = new ChainTrace[states.Length];
        var sample = new DrawSample();
        double location = 0.0, height = 0.0, count = 0.0;
        for (int i = 0; i < states.Length; i++)
        {
            traces[i] = new ChainTrace(1);
            if (!closedForm) continue;
            double chainLocation = 0.0, chainHeight = 0.0;
            foreach (BarsDraw d in states[i].Readout.Draws)
            {
                chainLocation += d.Peak.ArgmaxLocation;
                chainHeight += d.Peak.ArgmaxHeight;
                sample.Add(d.Peak.ArgmaxLocation);
                traces[i].Add(d.Peak.ArgmaxLocation);
            }
            location += chainLocation;
            height += chainHeight;
            count += states[i].Readout.Draws.Count;
        }
        return new PeakPosterior(count > 0.0 ? location / count : double.NaN, sample.Quantile(0.025), sample.Quantile(0.975),
                                 count > 0.0 ? height / count : double.NaN, ChainDiagnostics.Diagnose(traces, 0));
    }

    // Advance the chain one step (the IChain.Step body): the step scan runs the knot step, the lone chain or the tempered
    // ladder (returning the cold state), and every weightEvery steps the scale-mixture weight block.
    private KnotConfig StepOnce(BarsChainRun st)
    {
        st.Scan.Step();
        return st.Config;
    }

    // The step scan's blocks: the knot step on every step, and the weight block on its cadence when the run resamples
    // weights. The chain's own stream drives the knot step; the scan's stream is the weight stream.
    private IBlock<BarsChainRun>[] StepBlocks() => _resampler is null
        ? new IBlock<BarsChainRun>[] { new KnotStepBlock() }
        : new IBlock<BarsChainRun>[] { new KnotStepBlock(), new WeightBlock(this) };

    // The reversible-jump step, or the tempered ladder's scan returning its cold state, each on its own stream; accepted
    // when the cold chain's move was.
    private sealed class KnotStepBlock : IBlock<BarsChainRun>
    {
        public string Name => "knots";

        public bool Update(BarsChainRun st, Xoshiro256PlusPlus rng)
        {
            long accepted = st.RjChain.Accepted;
            st.Config = st.Ladder is not null ? st.Ladder.Step() : st.RjChain.Step();
            return st.RjChain.Accepted > accepted;
        }
    }

    // The scale-mixture weight block, every weightEvery steps on the weight stream: it draws (σ², β) from their exact
    // conditional given the knots and w, redraws w given the residuals and σ², and hands the RJ knot chain the target | w,
    // keeping its state and cumulative move counters.
    private sealed class WeightBlock(BarsEnsemble owner) : IBlock<BarsChainRun>
    {
        public string Name => "robust weights";

        public int Every => Math.Max(1, owner._weightEvery);

        public bool Update(BarsChainRun st, Xoshiro256PlusPlus rng)
        {
            double[,] z = owner._basis.Design(st.Config, owner._x);
            (double[] beta, double sigma2) = ((WeightedNormalModel)owner._model).SampleConditional(z, owner._y, st.W, rng);
            if (sigma2 > 0.0)
            {
                var r = new double[owner._y.Length];
                for (int i = 0; i < r.Length; i++)
                {
                    double f = 0.0;
                    for (int j = 0; j < beta.Length; j++) f += z[i, j] * beta[j];
                    r[i] = owner._y[i] - f;
                }
                st.W = owner._resampler!.Resample(r, sigma2, rng);
                st.RjChain.ReplaceTarget(new SplineTarget(owner._basis, owner._model, owner._prior, owner._x, owner._y, st.W,
                                                          st.Readout.Evaluations, owner._minColumnMass));
            }
            return true;
        }
    }

    // Reduce one post-burn draw and retain its state. Readouts affine in the curve use the conditional mean: the grid
    // values into the functionals (index 0 = knot count), and the conditional variance into its sink for FitSd. The
    // draw's coefficients and σ² come from one conditional draw on the readout stream; readouts nonlinear in the curve
    // use those coefficients: the draw's exact argmax and significant spans, retained as its PeakDraw (the peak summary,
    // λ(T) and π(T) are read from those after the run).
    private void AccumulateDraw(BarsChainRun st, KnotConfig draw, Span<double> functionals)
    {
        // Chains start on the support and never leave it; a draw off it has no conditional posterior to read, and its
        // NaN or unconverged coefficients would pass into the curve and the peaks unnoticed.
        if (!(st.RjChain.CurrentLogDensity > double.NegativeInfinity))
            throw new InvalidOperationException(
                $"A retained draw ({draw.Count} knots) lies off the posterior's support, so it has no conditional posterior.");
        int g = _grid.Length;
        double[,] zTrain = _basis.Design(draw, _x);
        IConditionalFit conditional = _model.FitConditional(zTrain, _y, st.W);   // one fit for all three readouts
        double[] mean = conditional.MeanCoefficients();
        double[,] zGrid = _basis.Design(draw, _grid);
        for (int gi = 0; gi < g; gi++)
        {
            double f = 0.0;
            for (int j = 0; j < mean.Length; j++) f += zGrid[gi, j] * mean[j];
            functionals[1 + gi] = f;
        }
        functionals[0] = draw.Count;

        BarsReadout r = st.Readout;
        if (conditional.ConditionalVariance(zGrid) is double[] v)
            r.ConditionalVariance.Add(v);

        // Every draw is retained, so the draw count stays every peak readout's denominator.
        ConditionalDraw c = conditional.Sample(st.ReadoutRng);
        r.Draws.Add(new BarsDraw(draw, c.Coefficients, c.Sigma2, st.W, PeakReadout(draw, c.Coefficients)));
    }

    // A draw's exact peak readout from its coefficients; a basis without closed-form peaks gives a NaN argmax and no
    // spans.
    private PeakDraw PeakReadout(KnotConfig draw, double[] coef)
    {
        if (_basis is not SplineBasis splineBasis) return new PeakDraw(double.NaN, double.NaN, Array.Empty<PeakSpan>());
        (double peakX, double peakF) = SplineExtrema.Argmax(draw, coef, splineBasis);
        var spans = SplineExtrema.SignificantPeakSpans(draw, coef, splineBasis, _peakProminence, _spanDropFraction);
        return new PeakDraw(peakX, peakF, spans.ToArray());
    }

    // A chain's λ and π grid counts, read from its retained draws: λ adds one per significant span at the grid point
    // nearest its location, and π adds one per draw at each grid point some span covers. Each draw still counts in the
    // denominator when it has no spans.
    private static (double[] Intensity, double[] Coverage) PeakCounts(IReadOnlyList<BarsDraw> draws, double[] grid)
    {
        var intensity = new double[grid.Length];
        var coverage = new double[grid.Length];
        foreach (BarsDraw d in draws)
        {
            foreach (PeakSpan s in d.Peak.Spans) intensity[NearestGridIndex(grid, s.Location)]++;
            foreach (int gi in CoveredGridPoints(d.Peak.Spans, grid)) coverage[gi]++;
        }
        return (intensity, coverage);
    }

    // π(T) is P(T lies in some significant span), so a draw adds 1 to a grid point at most once, however many spans
    // cover it: a nested peak's span lies inside its neighbour's (R14). λ(T) stays a count, one per span.
    internal static int[] CoveredGridPoints(IReadOnlyList<PeakSpan> spans, double[] grid)
    {
        var covered = new List<int>();
        for (int gi = 0; gi < grid.Length; gi++)
            for (int p = 0; p < spans.Count; p++)
                if (grid[gi] >= spans[p].Left && grid[gi] <= spans[p].Right) { covered.Add(gi); break; }
        return covered.ToArray();
    }

    private ReversibleJumpChain<KnotConfig> NewChain(BarsChainRun st)
        => new(Moves(st.Kernel), new SplineTarget(_basis, _model, _prior, _x, _y, st.W, st.Readout.Evaluations, _minColumnMass), st.Config, st.Rng);

    private ParallelTempering<KnotConfig> NewLadder(BarsChainRun st, int levels, int seed)
    {
        var target = new SplineTarget(_basis, _model, _prior, _x, _y, st.W, st.Readout.Evaluations, _minColumnMass);
        double[] betas = _temperReferenceRung
            ? ParallelTempering<KnotConfig>.ReferenceLadder(levels, _temperBetaMin)
            : ParallelTempering<KnotConfig>.GeometricLadder(levels, _temperBetaMin);
        return new ParallelTempering<KnotConfig>(_ => Moves(st.Kernel), target, st.Config, betas, _temperPath, seed,
                                                 _temperSchedule, _temperReferenceRung ? target.Reference(_temperPath) : null);
    }

    private IRjMove<KnotConfig>[] Moves(IKnotKernel kernel)
    {
        var density = new KnotBirthDensity(kernel, _birthExploration, _knotBase);
        var birth = new KnotBirthMove(density);
        var death = new KnotDeathMove(density);
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

    // Nearest grid index to x for the peak-intensity count (grid ascending — it is the evaluation grid).
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

    // Every count prior admits zero knots, the configuration with the fewest basis functions. If the data refuse even that,
    // no chain can start on the posterior's support, and the run fails here rather than return draws from off it.
    private void RequireAdmissibleZeroKnots()
    {
        var outcomes = new Tally<DesignRejection>();
        var target = new SplineTarget(_basis, _model, _prior, _x, _y, null, outcomes, _minColumnMass);
        if (target.LogPosterior(new KnotConfig([])) > double.NegativeInfinity) return;
        DesignRejection reason = DesignRejection.None;
        foreach ((DesignRejection category, long count) in outcomes.Entries)
            if (count > 0 && category != DesignRejection.None) reason = category;
        throw new InvalidOperationException(
            $"BARS cannot start on these {_x.Length} observations: the zero-knot design, the smallest, is refused ({reason}). " +
            "More observations, a lower-order basis or a smaller minColumnMass would admit it.");
    }

    // The start moved onto the posterior's support. A chain started off it leaves only by a proposal onto it, and where
    // every neighbour is off it too, as with more knots than a few observations admit, it never moves and its draws are
    // not posterior draws. So knots are dropped one at a time, at random, until the configuration is admissible: at worst
    // zero knots, which Run has checked. A start already on the support is kept, and draws nothing.
    private KnotConfig OnSupport(KnotConfig start, Xoshiro256PlusPlus rng)
    {
        var target = new SplineTarget(_basis, _model, _prior, _x, _y, null, null, _minColumnMass);
        KnotConfig config = start;
        while (!(target.LogPosterior(config) > double.NegativeInfinity))
        {
            double[] knots = config.InteriorKnots;
            if (knots.Length == 0) throw new InvalidOperationException("The zero-knot configuration is off the support.");
            int drop = rng.NextInt(knots.Length);
            var next = new double[knots.Length - 1];
            Array.Copy(knots, next, drop);
            Array.Copy(knots, drop + 1, next, drop, knots.Length - drop - 1);
            config = new KnotConfig(next);
        }
        return config;
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
    /// <see cref="Step"/> runs the step scan (the chain or tempered ladder, and the weight block), and its acceptance
    /// is the RJ chain's cumulative count; as an <see cref="IChainRun{KnotConfig}"/>
    /// it burns adaptively, reduces a draw to functionals, and carries the BARS readout (its retained draws).
    /// </summary>
    private sealed class BarsChainRun : IChain<KnotConfig>, IResumableChainRun<KnotConfig>
    {
        private readonly BarsEnsemble _owner;
        public readonly Xoshiro256PlusPlus Rng;
        public readonly BarsReadout Readout;              // this chain's mergeable readout statistics
        public KnotConfig Config;
        public double[]? W;
        public IKnotKernel Kernel = null!;
        public Xoshiro256PlusPlus ReadoutRng = null!;   // per-draw coefficient draws, decoupled from the chain RNG
        public BlockScan<BarsChainRun> Scan = null!;     // the step: the knot step, and the weight block on its cadence
        public ParallelTempering<KnotConfig>? Ladder;   // set ⇒ the RJ chain is this ladder's cold replica
        public long BurnScans, BurnRoundTrips, BurnRestarts;   // the ladder's counts at the end of burn-in
        public double[]? BurnRejection;                        // its per-pair rejection sums then

        private ReversibleJumpChain<KnotConfig> _rjChain = null!;

        public BarsChainRun(BarsEnsemble owner, Xoshiro256PlusPlus rng, KnotConfig start, int gridSize)
            : this(owner, rng, start, new BarsReadout(gridSize)) { }

        public BarsChainRun(BarsEnsemble owner, Xoshiro256PlusPlus rng, KnotConfig start, BarsReadout readout)
        {
            _owner = owner;
            Rng = rng;
            Config = start;
            Readout = readout;
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
        public IChainState SaveState() => _owner.Save(this);
    }
}
