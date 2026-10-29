using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass.Calibration;

/// <summary>
/// The conditional sampler's settings (D09): the run's lengths, R's decorrelation cadence and the optional tempering
/// ladder on D05's likelihood path. Values only.
/// </summary>
public sealed record CalibrationSampling
{
    /// <summary>The seed every stream derives from.</summary>
    public int MasterSeed { get; init; }

    /// <summary>Scans discarded, during which the θ walk and the variance steps adapt; then they are held.</summary>
    public int BurnIn { get; init; } = 2000;

    /// <summary>Draws kept.</summary>
    public int Samples { get; init; } = 2000;

    /// <summary>Scans per kept draw.</summary>
    public int Thin { get; init; } = 1;

    /// <summary>Scans between R's decorrelation steps (<c>decor.step.every</c>, 100): a uniform independence proposal on
    /// each continuous coordinate in turn. 0 runs none.</summary>
    public int DecorrelationEvery { get; init; } = 100;

    /// <summary>
    /// The inverse temperatures of a tempered run, cold first: 1, then strictly descending to at least 0, on the
    /// likelihood path through <see cref="ParallelTempering{TState}"/>; null runs one untempered chain. A rung at 0
    /// samples the prior, which every prior must make proper.
    /// </summary>
    public double[]? Tempering { get; init; }

    /// <summary>Checks the fields, naming each in the reason.</summary>
    public void Validate()
    {
        if (BurnIn < 0) throw Bad("BurnIn", "must not be negative");
        if (Samples < 1) throw Bad("Samples", "must be at least 1");
        if (Thin < 1) throw Bad("Thin", "must be at least 1");
        if (DecorrelationEvery < 0) throw Bad("DecorrelationEvery", "must not be negative");
        if (Tempering is not { } b) return;
        if (b.Length < 2) throw Bad("Tempering", "needs at least two rungs");
        if (b[0] != 1.0) throw Bad("Tempering", "starts at the cold rung, exactly 1");
        for (int i = 1; i < b.Length; i++)
            if (!(b[i] < b[i - 1] && b[i] >= 0.0)) throw Bad("Tempering", "must descend strictly to at least 0");
    }

    private static ArgumentException Bad(string field, string reason) => new($"CalibrationConfig.Sampling.{field}: {reason}");
}

/// <summary>
/// θ's random walk in [0, 1] coordinates, the reversible-jump chain's one move: u′ = u + L·z with L the factor of the
/// proposal covariance, symmetric; a proposal outside the box leaves the prior's support, so the chain rejects it (R's
/// −9999). Its covariance is set by <see cref="ThetaAdaptation"/> during burn-in, then held. Not available without a
/// continuous parameter.
/// </summary>
public sealed class ThetaWalk : IRjMove<CalibrationState>
{
    private double[,] _factor;

    public ThetaWalk(int dimension)
    {
        if (dimension < 0) throw new ArgumentOutOfRangeException(nameof(dimension));
        Dimension = dimension;
        _factor = new double[dimension, dimension];
        for (int i = 0; i < dimension; i++) _factor[i, i] = 0.1;
    }

    /// <summary>The continuous parameters p.</summary>
    public int Dimension { get; }

    /// <inheritdoc/>
    public string Key => "theta";

    /// <inheritdoc/>
    public string ReverseKey => "theta";

    /// <inheritdoc/>
    public double Weight(CalibrationState state) => Dimension > 0 ? 1.0 : 0.0;

    /// <summary>The proposal covariance's lower factor L.</summary>
    public double[,] Factor => (double[,])_factor.Clone();

    /// <summary>Sets the proposal covariance to <paramref name="covariance"/>; one that does not factor leaves it as it
    /// was. Returns whether it was set.</summary>
    public bool SetCovariance(double[,] covariance)
    {
        ArgumentNullException.ThrowIfNull(covariance);
        int p = Dimension;
        var l = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = covariance[i, j];
                for (int k = 0; k < j; k++) s -= l[i, k] * l[j, k];
                if (i == j)
                {
                    if (!(s > 0.0 && double.IsFinite(s))) return false;
                    l[i, i] = Math.Sqrt(s);
                }
                else l[i, j] = s / l[j, j];
            }
        _factor = l;
        return true;
    }

    /// <inheritdoc/>
    public Proposal<CalibrationState>? Propose(CalibrationState current, Xoshiro256PlusPlus rng)
    {
        if (Dimension == 0) return null;
        Span<double> z = stackalloc double[Dimension];
        for (int i = 0; i < Dimension; i++) z[i] = Normal.SampleStandard(rng);
        var next = current.Unit.ToArray();
        for (int i = 0; i < Dimension; i++)
        {
            double step = 0.0;
            for (int j = 0; j <= i; j++) step += _factor[i, j] * z[j];
            next[i] += step;
        }
        return new Proposal<CalibrationState>(current.WithUnit(next), 0.0, 0.0);
    }
}

/// <summary>
/// The θ walk's adaptation during burn-in, R's scheme (Haario et al.'s adaptive Metropolis): the running covariance of θ's
/// draws, scaled by 2.4²/p and by e^τ, with τ moved every 100 scans by min(0.1, 5/√scans) towards a 23% acceptance; the
/// covariance enters from scan 300, before which the walk keeps its start, 0.1² on each coordinate, scaled by e^τ. R adapts
/// for the whole run; here <see cref="Adapting"/> stops at the end of burn-in, so the kept draws come from one fixed
/// kernel.
/// </summary>
public sealed class ThetaAdaptation : IBlock<ReversibleJumpChain<CalibrationState>>
{
    private readonly ThetaWalk _walk;
    private readonly double[] _mean;
    private readonly double[,] _scatter;
    private long _scans, _lastAttempts, _lastAccepts;
    private double _tau;
    private const double Epsilon = 1e-13;

    public ThetaAdaptation(ThetaWalk walk)
    {
        _walk = walk ?? throw new ArgumentNullException(nameof(walk));
        _mean = new double[walk.Dimension];
        _scatter = new double[walk.Dimension, walk.Dimension];
    }

    /// <inheritdoc/>
    public string Name => "theta.adapt";

    /// <summary>Whether the walk adapts (burn-in).</summary>
    public bool Adapting { get; set; } = true;

    /// <summary>The log scale τ.</summary>
    public double LogScale => _tau;

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<CalibrationState> state, Xoshiro256PlusPlus rng)
    {
        int p = _walk.Dimension;
        if (!Adapting || p == 0) return false;
        _scans++;
        ReadOnlySpan<double> u = state.Current.Unit;
        for (int i = 0; i < p; i++)
        {
            double delta = u[i] - _mean[i];
            _mean[i] += delta / _scans;
            for (int j = 0; j < p; j++) _scatter[i, j] += delta * (u[j] - _mean[j]);
        }
        if (_scans % 100 == 0)
        {
            MoveStat walk = state.MoveStats()[0];
            long attempts = walk.Attempts - _lastAttempts, accepts = walk.Accepted - _lastAccepts;
            (_lastAttempts, _lastAccepts) = (walk.Attempts, walk.Accepted);
            double step = Math.Min(0.1, 5.0 / Math.Sqrt(_scans));
            if (attempts > 0) _tau += accepts / (double)attempts > 0.23 ? step : -step;
            var covariance = new double[p, p];
            double scale = Math.Exp(_tau);
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++)
                    covariance[i, j] = _scans >= 300
                        ? scale * 2.4 * 2.4 / p * (_scatter[i, j] / (_scans - 1) + (i == j ? Epsilon : 0.0))
                        : i == j ? scale * 0.01 : 0.0;
            _walk.SetCovariance(covariance);
        }
        return false;
    }
}

/// <summary>
/// R's decorrelation step: on its cadence, each continuous coordinate in turn takes a uniform independence proposal on
/// [0, 1], accepted by the tempered density's ratio (the uniform prior and proposal cancel), "to decorrelate the thetas
/// that don't change anything".
/// </summary>
public sealed class ThetaDecorrelation : CalibrationBlock, IBlock<ReversibleJumpChain<CalibrationState>>
{
    public ThetaDecorrelation(CalibrationTarget target, int every) : base(target)
    {
        if (every < 1) throw new ArgumentOutOfRangeException(nameof(every));
        Every = every;
    }

    /// <inheritdoc/>
    public override string Name => "theta.decorrelate";

    /// <inheritdoc/>
    public int Every { get; }

    /// <inheritdoc/>
    public override bool Update(ReversibleJumpChain<CalibrationState> state, Xoshiro256PlusPlus rng)
    {
        _ = Rung(state);
        bool any = false;
        for (int k = 0; k < Target.Model.Layout.Continuous.Count; k++)
        {
            CalibrationState current = state.Current;
            double[] unit = current.Unit.ToArray();
            unit[k] = rng.NextDouble();
            CalibrationState candidate = current.WithUnit(unit);
            CalibrationEvaluation evaluation = Target.EvaluateState(candidate);
            double logAccept = state.Temperature.LogDensity(evaluation) - state.CurrentLogDensity;
            if (Math.Log(1.0 - rng.NextDouble()) < logAccept)
            {
                Hand(state, candidate, evaluation);
                any = true;
            }
        }
        return any;
    }
}

/// <summary>
/// The conditional sampler for p(θ, σ², v | y, e_s) (D09, stage 7), one emulator draw: a scan of θ's walk, its
/// adaptation, R's decorrelation step, the categorical block, the variance groups on ln σ² and the discrepancy's Gaussian
/// draw. Untempered it is a <see cref="BlockScan{TState}"/> over a reversible-jump chain; tempered, every rung of a
/// <see cref="ParallelTempering{TState}"/> on the likelihood path runs the same blocks between its moves and the swaps,
/// whose ratio is (b_i − b_j)(ℓ(x_j) − ℓ(x_i)) (C-O1's corrected swap). As an <see cref="IChain{TDraw}"/> it presents
/// the cold state.
/// </summary>
public sealed class CalibrationChain : IChain<CalibrationState>
{
    private readonly List<ThetaAdaptation> _adaptations = [];
    private readonly List<VarianceBlock> _variances = [];
    private readonly ReversibleJumpChain<CalibrationState>? _lone;
    private readonly BlockScan<ReversibleJumpChain<CalibrationState>>? _scan;

    /// <param name="target">The conditional target.</param>
    /// <param name="sampling">The settings; only the decorrelation cadence and the ladder are read here.</param>
    /// <param name="start">The start, inside the prior's support.</param>
    /// <param name="seed">Every stream's seed.</param>
    public CalibrationChain(CalibrationTarget target, CalibrationSampling sampling, CalibrationState start, int seed)
    {
        ArgumentNullException.ThrowIfNull(target);
        ArgumentNullException.ThrowIfNull(sampling);
        ArgumentNullException.ThrowIfNull(start);
        sampling.Validate();
        Target = target;
        CalibrationModel model = target.Model;
        if (double.IsNegativeInfinity(model.LogPrior(start))) throw new ArgumentException("The start lies outside the prior's support.", nameof(start));
        int p = model.Layout.Continuous.Count;
        if (sampling.Tempering is not { } ladder)
        {
            var walk = new ThetaWalk(p);
            int[] seeds = SeedTree.Derive(seed, 2);
            _lone = new ReversibleJumpChain<CalibrationState>([walk], target, start, new Xoshiro256PlusPlus(seeds[0]));
            List<IBlock<ReversibleJumpChain<CalibrationState>>> blocks = [new MoveStep<CalibrationState>(), .. Blocks(walk, sampling)];
            _scan = new BlockScan<ReversibleJumpChain<CalibrationState>>(blocks, _lone, new Xoshiro256PlusPlus(seeds[1]), chain => chain);
            return;
        }
        if (ladder[^1] == 0.0 && !model.PriorsProper)
            throw new ArgumentException("A rung at b = 0 samples the prior, and a variance group's prior is improper.", nameof(sampling));
        var walks = new ThetaWalk[ladder.Length];
        Ladder = new ParallelTempering<CalibrationState>(
            i => [walks[i] = new ThetaWalk(p)], target, start, (double[])ladder.Clone(), TemperingPath.Likelihood, seed,
            blocksForReplica: i => Blocks(walks[i], sampling));
    }

    /// <summary>The target.</summary>
    public CalibrationTarget Target { get; }

    /// <summary>The ladder of a tempered run; null untempered.</summary>
    public ParallelTempering<CalibrationState>? Ladder { get; }

    /// <summary>The cold chain: the lone chain, or the ladder's b = 1 replica.</summary>
    public ReversibleJumpChain<CalibrationState> Cold => _lone ?? Ladder!.ColdChain;

    /// <summary>The cold state.</summary>
    public CalibrationState Current => Cold.Current;

    /// <inheritdoc/>
    public long Accepted => Cold.Accepted;

    /// <inheritdoc/>
    public long Attempts => Cold.Attempts;

    /// <summary>The cold chain's blocks' attempts and acceptances, the walk's first; for a ladder, its cold replica's.</summary>
    public IReadOnlyList<MoveStat> Statistics()
    {
        var stats = new List<MoveStat>(Cold.MoveStats());
        stats.AddRange(_scan is not null ? _scan.BlockStats().Skip(1) : Ladder!.BlockStats(0));
        return stats;
    }

    /// <inheritdoc/>
    public CalibrationState Step()
    {
        if (_scan is not null) _scan.Step();
        else Ladder!.Step();
        return Current;
    }

    /// <summary>Whether the θ walks and the variance steps adapt; burn-in sets it, and the kept draws need it off.</summary>
    public bool Adapting
    {
        get => _adaptations.Count > 0 ? _adaptations[0].Adapting : _variances.Count > 0 && _variances[0].Adapting;
        set
        {
            foreach (ThetaAdaptation a in _adaptations) a.Adapting = value;
            foreach (VarianceBlock v in _variances) v.Adapting = value;
        }
    }

    /// <summary><paramref name="scans"/> scans with adaptation on, then held.</summary>
    public void Burn(int scans)
    {
        if (scans < 0) throw new ArgumentOutOfRangeException(nameof(scans));
        Adapting = true;
        for (int i = 0; i < scans; i++) Step();
        Adapting = false;
    }

    // One replica's blocks after its θ walk: adaptation, decorrelation, levels, variances, discrepancy.
    private List<IBlock<ReversibleJumpChain<CalibrationState>>> Blocks(ThetaWalk walk, CalibrationSampling sampling)
    {
        CalibrationModel model = Target.Model;
        var adaptation = new ThetaAdaptation(walk) { Adapting = false };
        var variances = new VarianceBlock(Target);
        lock (_adaptations)
        {
            _adaptations.Add(adaptation);
            _variances.Add(variances);
        }
        List<IBlock<ReversibleJumpChain<CalibrationState>>> blocks = [adaptation];
        if (sampling.DecorrelationEvery > 0 && walk.Dimension > 0) blocks.Add(new ThetaDecorrelation(Target, sampling.DecorrelationEvery));
        if (model.Layout.Categorical.Count > 0) blocks.Add(new CategoricalBlock(Target, model.Config.EnumerationCap));
        blocks.Add(variances);
        if (model.Coefficients > 0) blocks.Add(new DiscrepancyBlock(Target));
        return blocks;
    }
}

/// <summary>
/// One kept draw of the conditional sampler (D09's outputs): the emulator draw it conditions on, θ in physical
/// coordinates through the recorded transforms, the levels (indices into the declared levels), the group variances, the
/// discrepancy coefficients v, and the untempered log-likelihood.
/// </summary>
public sealed record CalibrationDraw(int EmulatorDraw, double[] Theta, int[] Levels, double[] Variances, double[] Coefficients, double LogLikelihood)
{
    /// <summary>The draw of <paramref name="state"/>, evaluated as <paramref name="evaluation"/> under <paramref name="target"/>.</summary>
    public static CalibrationDraw Of(CalibrationTarget target, CalibrationState state, ITargetEvaluation evaluation)
        => new(target.Draw, target.Model.Physical(state), state.Levels.ToArray(), state.Variances.ToArray(), state.Coefficients.ToArray(),
               ((CalibrationEvaluation)evaluation).LogLikelihoodValue);
}

/// <summary>
/// A conditional run on one emulator draw: its kept draws, with the burn-in and thinning that produced them, the cold
/// chain's move and block statistics, and a ladder's swap acceptance (NaN untempered).
/// </summary>
public sealed record ConditionalRun(int EmulatorDraw, IReadOnlyList<CalibrationDraw> Draws, int BurnIn, int Thin,
                                    IReadOnlyList<MoveStat> Statistics, double SwapRate);

/// <summary>The conditional sampler's entry point (D09, stage 7).</summary>
public static class CalibrationSampler
{
    /// <summary>
    /// A run of <see cref="CalibrationChain"/> on emulator draw <paramref name="draw"/> under the model's
    /// <see cref="CalibrationSampling"/> (or <paramref name="sampling"/>), from <paramref name="start"/> or a draw from the
    /// prior: burn-in with adaptation, then the kept draws, one every <see cref="CalibrationSampling.Thin"/> scans.
    /// </summary>
    public static ConditionalRun Run(CalibrationModel model, int draw, CalibrationSampling sampling, CalibrationState? start = null)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(sampling);
        sampling.Validate();
        CalibrationTarget target = model.Target(draw);
        int[] seeds = SeedTree.Derive(sampling.MasterSeed, 2);
        CalibrationState initial = start ?? model.Start(new Xoshiro256PlusPlus(seeds[0]));
        var chain = new CalibrationChain(target, sampling, initial, seeds[1]);
        chain.Burn(sampling.BurnIn);
        var draws = new List<CalibrationDraw>(sampling.Samples);
        for (int s = 0; s < sampling.Samples; s++)
        {
            for (int t = 0; t < sampling.Thin; t++) chain.Step();
            draws.Add(CalibrationDraw.Of(target, chain.Current, chain.Cold.CurrentEvaluation));
        }
        double swaps = chain.Ladder is { } ladder && ladder.SwapAttempts > 0 ? ladder.SwapAccepts / (double)ladder.SwapAttempts : double.NaN;
        return new ConditionalRun(draw, draws, sampling.BurnIn, sampling.Thin, chain.Statistics(), swaps);
    }
}
