using System;
using System.Collections.Generic;
using System.Threading;
using Maths.LinAlg;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Smc;

namespace Maths.Regression.Spline.Bass.Calibration;

/// <summary>
/// One stage's block acceptances in a cut SMC run's mutation: the θ walk's, the variance steps' and, above the
/// enumeration cap, the single-site level steps' attempts and acceptances, over every particle and sweep of the stage.
/// The levels' enumerated draw and v's Gaussian draw always accept.
/// </summary>
public sealed record CutSmcBlockRates(long ThetaAttempts, long ThetaAccepts, long VarianceAttempts, long VarianceAccepts,
                                      long SiteAttempts, long SiteAccepts)
{
    /// <summary>The θ walk's acceptance; NaN without a continuous parameter.</summary>
    public double Theta => ThetaAttempts > 0 ? ThetaAccepts / (double)ThetaAttempts : double.NaN;

    /// <summary>The variance steps' acceptance, over the groups.</summary>
    public double Variances => VarianceAttempts > 0 ? VarianceAccepts / (double)VarianceAttempts : double.NaN;

    /// <summary>The single-site level steps' acceptance; NaN when the levels are enumerated or absent.</summary>
    public double Sites => SiteAttempts > 0 ? SiteAccepts / (double)SiteAttempts : double.NaN;
}

/// <summary>
/// The conditional sampler's block scan (stage 7) aimed at a stage of the cut SMC (D09, stage 8), as the
/// <see cref="ISmcKernel{TState}"/> that mutates its particles. A stage's density is log π + Σ_k w_k·ℓ_k: one emulator
/// draw at rung b on D05's likelihood path (the first leg, prior to μ_1), or two at 1 − t and t
/// (a <see cref="CalibrationBridge"/>). One sweep runs θ's walk, the levels, the variance groups on ln σ² and the
/// discrepancy's draw, each leaving that density invariant:
/// <list type="bullet">
/// <item>the walk and the variance steps are Metropolis on the stage's density, the latter with its Jacobian;</item>
/// <item>the levels' conditional is enumerated under it while the combinations number at most the cap, else each level
/// takes a single-site Metropolis step;</item>
/// <item>v is drawn from the bridged Gaussian, whose precision and mean add the terms', each quadratic in v.</item>
/// </list>
/// The walk's covariance, (2.38²/p) times θ's weighted covariance, and each variance step, 2.38 times ln σ²_g's weighted
/// spread, are read from the population once per stage (<see cref="Prepare"/>), adaptation from the past (D10), and held
/// for the stage's moves. R's decorrelation step is left out: its cadence of scans has no counterpart in a few sweeps per
/// stage, and the population already spans θ. Every change is evaluated afresh; the states' caches keep two draws.
/// </summary>
public sealed class CutSmcKernel : ISmcKernel<CalibrationState>
{
    private readonly CalibrationModel _model;
    private readonly int[][]? _combinations;   // the level combinations within the cap; null above it or without levels
    private readonly List<CutSmcBlockRates> _rates = [];
    private ThetaWalk _walk;
    private double[] _varianceSteps;
    private Counters _counters = new();
    private const double Spread = 2.38;

    /// <summary>A kernel for <paramref name="model"/>'s stages.</summary>
    public CutSmcKernel(CalibrationModel model)
    {
        ArgumentNullException.ThrowIfNull(model);
        _model = model;
        ParameterLayout layout = model.Layout;
        if (layout.Categorical.Count > 0 && layout.Combinations <= model.Config.EnumerationCap)
        {
            _combinations = new int[layout.Combinations][];
            for (long c = 0; c < layout.Combinations; c++)
            {
                _combinations[c] = new int[layout.Categorical.Count];
                layout.Levels(c, _combinations[c]);
            }
        }
        _walk = new ThetaWalk(layout.Continuous.Count);
        _varianceSteps = new double[model.Measurement.Groups];
        Array.Fill(_varianceSteps, 0.5);
    }

    /// <summary>Whether the levels' conditional is enumerated (the combinations within the cap).</summary>
    public bool Enumerated => _combinations is not null;

    /// <summary>The walk's proposal factor (lower, unit coordinates) set by the last <see cref="Prepare"/>.</summary>
    public double[,] WalkFactor => _walk.Factor;

    /// <summary>Each variance group's step on ln σ² set by the last <see cref="Prepare"/>.</summary>
    public IReadOnlyList<double> VarianceSteps => (double[])_varianceSteps.Clone();

    /// <summary>One entry per prepared stage, in order: its blocks' attempts and acceptances.</summary>
    public IReadOnlyList<CutSmcBlockRates> StageRates
    {
        get
        {
            lock (_rates)
            {
                var all = new List<CutSmcBlockRates>(_rates);
                if (_rates.Count > 0) all[^1] = _counters.Rates();
                return all;
            }
        }
    }

    /// <inheritdoc/>
    public void Prepare(SmcStage<CalibrationState> stage, IReadOnlyList<CalibrationState> particles, ReadOnlySpan<double> weights)
    {
        _ = Terms(stage);
        ArgumentNullException.ThrowIfNull(particles);
        if (weights.Length != particles.Count) throw new ArgumentException("One weight per particle.", nameof(weights));
        int p = _model.Layout.Continuous.Count, n = particles.Count;
        if (p > 0)
        {
            var mean = new double[p];
            for (int i = 0; i < n; i++)
                for (int a = 0; a < p; a++) mean[a] += weights[i] * particles[i].Unit[a];
            var covariance = new double[p, p];
            for (int i = 0; i < n; i++)
                for (int a = 0; a < p; a++)
                    for (int c = 0; c <= a; c++) covariance[a, c] += weights[i] * (particles[i].Unit[a] - mean[a]) * (particles[i].Unit[c] - mean[c]);
            double scale = Spread * Spread / p;
            for (int a = 0; a < p; a++)
                for (int c = 0; c <= a; c++)
                    covariance[a, c] = covariance[c, a] = scale * (covariance[a, c] + (a == c ? Floor : 0.0));
            var walk = new ThetaWalk(p);
            if (!walk.SetCovariance(covariance))
            {
                var diagonal = new double[p, p];
                for (int a = 0; a < p; a++) diagonal[a, a] = covariance[a, a];
                walk.SetCovariance(diagonal);
            }
            _walk = walk;
        }
        var steps = new double[_varianceSteps.Length];
        for (int g = 0; g < steps.Length; g++)
        {
            double mean = 0.0, variance = 0.0;
            for (int i = 0; i < n; i++) mean += weights[i] * Math.Log(particles[i].Variances[g]);
            for (int i = 0; i < n; i++)
            {
                double e = Math.Log(particles[i].Variances[g]) - mean;
                variance += weights[i] * e * e;
            }
            steps[g] = Spread * Math.Sqrt(variance + Floor);
        }
        _varianceSteps = steps;
        lock (_rates)
        {
            if (_rates.Count > 0) _rates[^1] = _counters.Rates();
            _counters = new Counters();
            _rates.Add(_counters.Rates());
        }
    }

    // Added to each variance read from the population, so a collapsed one still moves.
    private const double Floor = 1e-12;

    /// <summary>One sweep of the scan at <paramref name="stage"/>: θ's walk, the levels, the variances, v. Accepted reports
    /// the walk's acceptance, or the variance steps' without a continuous parameter.</summary>
    public SmcMove<CalibrationState> Move(CalibrationState state, ITargetEvaluation evaluation, SmcStage<CalibrationState> stage, Xoshiro256PlusPlus rng)
    {
        LikelihoodTerm[] terms = Terms(stage);
        var sweep = new Sweep(state, evaluation, stage);
        bool walked = _model.Layout.Continuous.Count > 0 && Walk(ref sweep, rng);
        if (_model.Layout.Categorical.Count > 0) Levels(ref sweep, terms, rng);
        bool varied = Variances(ref sweep, rng);
        if (_model.Coefficients > 0) Discrepancy(ref sweep, terms, rng);
        return new SmcMove<CalibrationState>(sweep.State, sweep.Evaluation, _model.Layout.Continuous.Count > 0 ? walked : varied);
    }

    // The state, its evaluation and density under the stage, as the sweep's blocks hand them on.
    internal struct Sweep(CalibrationState state, ITargetEvaluation evaluation, SmcStage<CalibrationState> stage)
    {
        public CalibrationState State = state;
        public ITargetEvaluation Evaluation = evaluation;
        public double Density = stage.LogDensity(evaluation);
        public readonly SmcStage<CalibrationState> Stage = stage;

        public void Take(CalibrationState state, ITargetEvaluation evaluation, double density)
        {
            State = state;
            Evaluation = evaluation;
            Density = density;
        }

        public void Take(CalibrationState state)
        {
            ITargetEvaluation evaluation = Stage.Target.Evaluate(state);
            Take(state, evaluation, Stage.LogDensity(evaluation));
        }
    }

    // θ's walk: a symmetric proposal from the stage's factor, Metropolis on the stage's density.
    internal bool Walk(ref Sweep sweep, Xoshiro256PlusPlus rng)
    {
        CalibrationState candidate = _walk.Propose(sweep.State, rng)!.Value.Candidate;
        ITargetEvaluation evaluation = sweep.Stage.Target.Evaluate(candidate);
        double density = sweep.Stage.LogDensity(evaluation);
        bool accepted = Math.Log(1.0 - rng.NextDouble()) < density - sweep.Density;
        if (accepted) sweep.Take(candidate, evaluation, density);
        Interlocked.Increment(ref _counters.ThetaAttempts);
        if (accepted) Interlocked.Increment(ref _counters.ThetaAccepts);
        return accepted;
    }

    // The levels: their conditional under the stage's terms enumerated and drawn from, each term's prediction kept with
    // the new state; above the cap, a Metropolis step per categorical parameter to a uniformly chosen other level.
    internal void Levels(ref Sweep sweep, LikelihoodTerm[] terms, Xoshiro256PlusPlus rng)
    {
        if (_combinations is null)
        {
            ParameterLayout layout = _model.Layout;
            for (int j = 0; j < layout.Categorical.Count; j++)
            {
                int[] levels = sweep.State.Levels.ToArray();
                int other = rng.NextInt(layout.Categorical[j].Levels.Length - 1);
                levels[j] = other >= levels[j] ? other + 1 : other;
                CalibrationState candidate = sweep.State.WithLevels(levels);
                ITargetEvaluation evaluation = sweep.Stage.Target.Evaluate(candidate);
                double density = sweep.Stage.LogDensity(evaluation);
                bool accepted = Math.Log(1.0 - rng.NextDouble()) < density - sweep.Density;
                if (accepted) sweep.Take(candidate, evaluation, density);
                Interlocked.Increment(ref _counters.SiteAttempts);
                if (accepted) Interlocked.Increment(ref _counters.SiteAccepts);
            }
            return;
        }
        (double[] probabilities, double[][,] predictions) = CategoricalBlock.Enumerate(_model, sweep.State, terms, _combinations);
        double u = rng.NextDouble(), cumulative = 0.0;
        int chosen = probabilities.Length - 1;
        for (int c = 0; c < probabilities.Length; c++)
            if (u < (cumulative += probabilities[c])) { chosen = c; break; }
        PredictionCache? cache = null;
        for (int t = 0; t < terms.Length; t++)
        {
            var row = new double[_model.Outputs];
            for (int j = 0; j < row.Length; j++) row[j] = predictions[t][chosen, j];
            cache = PredictionCache.Add(cache, terms[t].Target, row);
        }
        sweep.Take(cache is null ? sweep.State.WithLevels(_combinations[chosen]) : sweep.State.WithLevels(_combinations[chosen], cache));
    }

    // Each group's σ²_g by a random-walk step on u = ln σ²_g, whose target carries the Jacobian e^u.
    internal bool Variances(ref Sweep sweep, Xoshiro256PlusPlus rng)
    {
        bool any = false;
        for (int g = 0; g < _varianceSteps.Length; g++)
        {
            CalibrationState current = sweep.State;
            SmcStage<CalibrationState> stage = sweep.Stage;
            double[] variances = current.Variances.ToArray();
            double u = Math.Log(variances[g]);
            double logTarget = sweep.Density + u;
            CalibrationState? candidate = null;
            ITargetEvaluation? evaluation = null;
            double density = double.NegativeInfinity;
            double LogDensity(double proposal)
            {
                double[] next = (double[])variances.Clone();
                next[g] = Math.Exp(proposal);
                candidate = current.WithVariances(next);
                evaluation = stage.Target.Evaluate(candidate);
                density = stage.LogDensity(evaluation);
                return double.IsNegativeInfinity(density) ? double.NegativeInfinity : density + proposal;
            }
            bool accepted = RandomWalkMetropolis.Step(ref u, ref logTarget, _varianceSteps[g], LogDensity, rng);
            if (accepted) sweep.Take(candidate!, evaluation!, density);
            Interlocked.Increment(ref _counters.VarianceAttempts);
            if (accepted) Interlocked.Increment(ref _counters.VarianceAccepts);
            any |= accepted;
        }
        return any;
    }

    // v from its Gaussian conditional under the stage's terms.
    internal void Discrepancy(ref Sweep sweep, LikelihoodTerm[] terms, Xoshiro256PlusPlus rng)
    {
        (DenseCholesky precision, double[] rhs) = DiscrepancyBlock.Conditional(_model, sweep.State, terms);
        var v = new double[_model.Coefficients];
        new GaussianConditional(v.Length).Draw(precision, rhs, rng, v);
        sweep.Take(sweep.State.WithCoefficients(v));
    }

    /// <summary>
    /// The weighted likelihoods of <paramref name="stage"/>'s density log π + Σ_k w_k·ℓ_k: [(e, b)] for a draw's
    /// conditional at rung b, [(e_a, 1 − t), (e_b, t)] for a bridge at t, a zero weight dropped and a bridge from a draw
    /// to itself taken whole. Refuses another path, target or model.
    /// </summary>
    internal LikelihoodTerm[] Terms(SmcStage<CalibrationState> stage)
    {
        if (stage.Temperature.Path != TemperingPath.Likelihood)
            throw new InvalidOperationException("The cut SMC tempers on D05's likelihood path; its blocks draw from that path's conditionals.");
        double b = stage.Temperature.InverseTemperature;
        switch (stage.Target)
        {
            case CalibrationTarget target when ReferenceEquals(target.Model, _model):
                return b > 0.0 ? [new LikelihoodTerm(target, b)] : [];
            case CalibrationBridge bridge when ReferenceEquals(bridge.Model, _model):
                if (ReferenceEquals(bridge.From, bridge.To)) return [new LikelihoodTerm(bridge.To, 1.0)];
                if (b == 0.0) return [new LikelihoodTerm(bridge.From, 1.0)];
                if (b == 1.0) return [new LikelihoodTerm(bridge.To, 1.0)];
                return [new LikelihoodTerm(bridge.From, 1.0 - b), new LikelihoodTerm(bridge.To, b)];
            default:
                throw new InvalidOperationException("The kernel moves its model's conditional targets and the bridges between them.");
        }
    }

    private sealed class Counters
    {
        public long ThetaAttempts, ThetaAccepts, VarianceAttempts, VarianceAccepts, SiteAttempts, SiteAccepts;

        public CutSmcBlockRates Rates() => new(Interlocked.Read(ref ThetaAttempts), Interlocked.Read(ref ThetaAccepts),
            Interlocked.Read(ref VarianceAttempts), Interlocked.Read(ref VarianceAccepts),
            Interlocked.Read(ref SiteAttempts), Interlocked.Read(ref SiteAccepts));
    }
}

/// <summary>The cut SMC's particles: each drawn exactly from the prior (<see cref="CalibrationModel.DrawPrior"/>), the
/// first leg's normalized first stage; an offspring a distinct state with the same values.</summary>
public sealed class CutSmcParticles : ISmcParticles<CalibrationState>
{
    private readonly CalibrationModel _model;

    /// <summary>Particles for <paramref name="model"/>, whose priors must all be proper.</summary>
    public CutSmcParticles(CalibrationModel model)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (!model.PriorsProper) throw new ArgumentException("The first leg starts from the prior, and a variance group's prior is improper.", nameof(model));
        _model = model;
    }

    /// <inheritdoc/>
    public CalibrationState Draw(Xoshiro256PlusPlus rng) => _model.DrawPrior(rng);

    /// <inheritdoc/>
    public CalibrationState Fork(CalibrationState state)
    {
        ArgumentNullException.ThrowIfNull(state);
        return state.Fork();
    }
}
