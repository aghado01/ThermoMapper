using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Maths.Rng;

namespace Maths.Samplers.Smc;

/// <summary>Settings for an SMC run.</summary>
public sealed record SmcOptions
{
    /// <summary>The population size N.</summary>
    public required int Particles { get; init; }

    /// <summary>How the population is resampled.</summary>
    public ResamplingScheme Resampling { get; init; } = ResamplingScheme.Systematic;

    /// <summary>Resample when ESS/N falls below this fraction: 1 resamples at every stage, 0 never.</summary>
    public double ResampleBelow { get; init; } = 0.5;

    /// <summary>Kernel moves per particle at each stage; 0 makes the run sequential importance sampling.</summary>
    public int MutationSteps { get; init; } = 1;

    /// <summary>The most stages a run may take after its first: the adaptive schedule's cap, and the length of the
    /// stage seeds derived up front.</summary>
    public int MaxStages { get; init; } = 1000;

    /// <summary>Particles moved concurrently. Every particle has its own stream, so the result does not depend on
    /// this.</summary>
    public int MaxDegreeOfParallelism { get; init; } = 1;
}

/// <summary>What one stage did. Stage 0 is the initial draw.</summary>
/// <param name="Index">The stage's position in the run.</param>
/// <param name="Temperature">The stage's temperature.</param>
/// <param name="EffectiveSampleSize">ESS of the reweighted population, before any resampling.</param>
/// <param name="ConditionalEss">The reweighting step's conditional ESS.</param>
/// <param name="Resampled">Whether the population was resampled.</param>
/// <param name="DistinctAncestors">Particles with at least one offspring (N when not resampled): the step's loss of
/// diversity.</param>
/// <param name="Families">Distinct initial particles among the current particles' lineages: the run's cumulative
/// loss of diversity. Equal weights do not make the particles independent; this says how many lineages remain.</param>
/// <param name="AcceptanceRate">The kernel's acceptance over the stage's moves; NaN at stage 0 or with no moves.</param>
/// <param name="LogNormalizerIncrement">log(Z_stage/Z_previous) as the reweighting estimates it; 0 at stage 0.</param>
public readonly record struct SmcStageReport(
    int Index, Temperature Temperature, double EffectiveSampleSize, double ConditionalEss, bool Resampled,
    int DistinctAncestors, int Families, double AcceptanceRate, double LogNormalizerIncrement);

/// <summary>A run's population and its account.</summary>
/// <param name="Particles">The current particles.</param>
/// <param name="Weights">Their normalized weights.</param>
/// <param name="LogNormalizerRatio">log(Z_current/Z_first), the sum of the stages' increments, for the stage densities
/// as the targets supply them. It is a model evidence only when those densities are complete and the first stage is a
/// normalized prior.</param>
/// <param name="Stages">One report per stage, the initial draw first.</param>
public sealed record SmcResult<TState>(
    IReadOnlyList<TState> Particles, double[] Weights, double LogNormalizerRatio, IReadOnlyList<SmcStageReport> Stages);

/// <summary>Everything a run needs to continue: restarting from it reproduces the uninterrupted run bit for bit.</summary>
public sealed record SmcCheckpoint<TState>(
    int Seed, int StagesDone, SmcStage<TState> Stage, TState[] Particles, ITargetEvaluation[] Evaluations,
    double[] Weights, int[] Families, double LogNormalizerRatio, SmcStageReport[] Reports);

/// <summary>
/// A sequential Monte Carlo population carried through a sequence of stages (Del Moral, Doucet and Jasra 2006; the
/// invitation of Dai, Heng, Jacob and Whiteley 2022): reweight by the ratio of consecutive stage densities, resample
/// when the weights concentrate, then mutate every particle with a kernel that leaves the new stage invariant. It
/// reweights through <see cref="ITargetEvaluation"/>s, so moving to another temperature of the same target queries
/// each particle's evaluation instead of repeating its work.
/// <para>Every random draw comes from a stream derived from the run's seed by position (<see cref="SeedTree"/>): one
/// per initial particle, then, at each stage, one for resampling and one per particle. The result therefore does not
/// depend on how many particles move at once, and a checkpoint restarts the uninterrupted run exactly.</para>
/// </summary>
public sealed class SmcRun<TState>
{
    private readonly ISmcParticles<TState> _particles;
    private readonly ISmcKernel<TState> _kernel;
    private readonly SmcOptions _options;
    private readonly int _seed;
    private readonly int[] _stageSeeds;
    private readonly List<SmcStageReport> _reports;
    private TState[] _state;
    private ITargetEvaluation[] _evaluations;
    private double[] _logDensity;
    private double[] _weights;
    private int[] _families;
    private double _logNormalizer;

    private SmcRun(ISmcParticles<TState> particles, ISmcKernel<TState> kernel, SmcOptions options, int seed,
                   SmcStage<TState> stage, TState[] state, ITargetEvaluation[] evaluations, double[] weights,
                   int[] families, int stagesDone, double logNormalizer, List<SmcStageReport> reports)
    {
        _particles = particles;
        _kernel = kernel;
        _options = options;
        _seed = seed;
        _stageSeeds = SeedTree.Derive(SeedTree.Derive(seed, 2)[1], options.MaxStages + 1);
        Stage = stage;
        _state = state;
        _evaluations = evaluations;
        _weights = weights;
        _families = families;
        StagesDone = stagesDone;
        _logNormalizer = logNormalizer;
        _reports = reports;
        _logDensity = new double[state.Length];
        for (int i = 0; i < state.Length; i++) _logDensity[i] = stage.LogDensity(evaluations[i]);
    }

    /// <summary>The current stage.</summary>
    public SmcStage<TState> Stage { get; private set; }

    /// <summary>Stages taken after the initial draw.</summary>
    public int StagesDone { get; private set; }

    /// <summary>log(Z_current/Z_first).</summary>
    public double LogNormalizerRatio => _logNormalizer;

    /// <summary>One report per stage, the initial draw first.</summary>
    public IReadOnlyList<SmcStageReport> Reports => _reports;

    /// <summary>The population size.</summary>
    public int Count => _state.Length;

    /// <summary>Draws the initial population from <paramref name="first"/>, with equal weights.</summary>
    public static SmcRun<TState> Start(SmcStage<TState> first, ISmcParticles<TState> particles, ISmcKernel<TState> kernel, SmcOptions options, int seed)
    {
        ArgumentNullException.ThrowIfNull(first.Target);
        ArgumentNullException.ThrowIfNull(particles);
        ArgumentNullException.ThrowIfNull(kernel);
        Validate(options);
        int n = options.Particles;
        int[] seeds = SeedTree.Derive(SeedTree.Derive(seed, 2)[0], n);
        var state = new TState[n];
        var evaluations = new ITargetEvaluation[n];
        Parallel.For(0, n, new ParallelOptions { MaxDegreeOfParallelism = options.MaxDegreeOfParallelism }, i =>
        {
            state[i] = particles.Draw(new Xoshiro256PlusPlus(seeds[i]));
            evaluations[i] = first.Target.Evaluate(state[i]);
        });
        for (int i = 0; i < n; i++)
            if (double.IsNegativeInfinity(first.LogDensity(evaluations[i])))
                throw new InvalidOperationException($"Initial particle {i} lies outside the first stage's support.");
        var weights = new double[n];
        Array.Fill(weights, 1.0 / n);
        var families = new int[n];
        for (int i = 0; i < n; i++) families[i] = i;
        var reports = new List<SmcStageReport> { new(0, first.Temperature, n, n, false, n, n, double.NaN, 0.0) };
        return new SmcRun<TState>(particles, kernel, options, seed, first, state, evaluations, weights, families, 0, 0.0, reports);
    }

    /// <summary>Continues a checkpointed run; the checkpoint itself stays reusable.</summary>
    public static SmcRun<TState> Resume(SmcCheckpoint<TState> checkpoint, ISmcParticles<TState> particles, ISmcKernel<TState> kernel, SmcOptions options)
    {
        ArgumentNullException.ThrowIfNull(checkpoint);
        ArgumentNullException.ThrowIfNull(particles);
        ArgumentNullException.ThrowIfNull(kernel);
        Validate(options);
        if (checkpoint.Particles.Length != options.Particles)
            throw new ArgumentException("The checkpoint's population size differs from the options'.", nameof(options));
        var state = new TState[checkpoint.Particles.Length];
        for (int i = 0; i < state.Length; i++) state[i] = particles.Fork(checkpoint.Particles[i]);
        return new SmcRun<TState>(particles, kernel, options, checkpoint.Seed, checkpoint.Stage, state,
            (ITargetEvaluation[])checkpoint.Evaluations.Clone(), (double[])checkpoint.Weights.Clone(),
            (int[])checkpoint.Families.Clone(), checkpoint.StagesDone, checkpoint.LogNormalizerRatio,
            new List<SmcStageReport>(checkpoint.Reports));
    }

    /// <summary>
    /// The conditional ESS a move to <paramref name="next"/> would have, without moving: what an adaptive schedule
    /// tunes. When <paramref name="next"/> shares the current target, the particles' evaluations answer it with no new
    /// model work.
    /// </summary>
    public double ConditionalEss(SmcStage<TState> next)
    {
        ITargetEvaluation[] evaluations = EvaluationsFor(next);
        var increments = new double[_state.Length];
        for (int i = 0; i < increments.Length; i++) increments[i] = Increment(next.LogDensity(evaluations[i]), _logDensity[i]);
        return SmcWeights.ConditionalEss(_weights, increments);
    }

    /// <summary>Reweights to <paramref name="next"/>, resamples if the weights have concentrated, then mutates.</summary>
    public void Advance(SmcStage<TState> next)
    {
        ArgumentNullException.ThrowIfNull(next.Target);
        if (StagesDone == _options.MaxStages) throw new InvalidOperationException($"The run has taken its {_options.MaxStages} stages.");
        int n = _state.Length;

        // Reweight.
        ITargetEvaluation[] evaluations = EvaluationsFor(next);
        var density = new double[n];
        var increments = new double[n];
        for (int i = 0; i < n; i++)
        {
            density[i] = next.LogDensity(evaluations[i]);
            increments[i] = Increment(density[i], _logDensity[i]);
        }
        double logIncrement = SmcWeights.LogNormalizerIncrement(_weights, increments);
        double conditionalEss = SmcWeights.ConditionalEss(_weights, increments);
        var logWeights = new double[n];
        for (int i = 0; i < n; i++) logWeights[i] = Math.Log(_weights[i]) + increments[i];
        SmcWeights.Normalize(logWeights, _weights);
        double ess = SmcWeights.EffectiveSampleSize(_weights);
        _evaluations = evaluations;
        _logDensity = density;
        Stage = next;
        StagesDone++;
        _logNormalizer += logIncrement;
        int[] seeds = SeedTree.Derive(_stageSeeds[StagesDone], n + 1);

        // Resample: the first offspring of a particle takes it, later ones fork it.
        bool resampled = _options.ResampleBelow >= 1.0 || ess < _options.ResampleBelow * n;
        int distinct = n;
        if (resampled)
        {
            var ancestors = new int[n];
            Resampling.Ancestors(_options.Resampling, _weights, new Xoshiro256PlusPlus(seeds[0]), ancestors);
            var state = new TState[n];
            var evals = new ITargetEvaluation[n];
            var logDensity = new double[n];
            var families = new int[n];
            var taken = new bool[n];
            distinct = 0;
            for (int k = 0; k < n; k++)
            {
                int a = ancestors[k];
                if (taken[a]) state[k] = _particles.Fork(_state[a]);
                else
                {
                    state[k] = _state[a];
                    taken[a] = true;
                    distinct++;
                }
                evals[k] = _evaluations[a];
                logDensity[k] = _logDensity[a];
                families[k] = _families[a];
            }
            _state = state;
            _evaluations = evals;
            _logDensity = logDensity;
            _families = families;
            Array.Fill(_weights, 1.0 / n);
        }

        // Mutate, each particle on its own stream; particles of zero weight are left as they are.
        _kernel.Prepare(next, _state, _weights);
        var accepted = new int[n];
        int steps = _options.MutationSteps;
        Parallel.For(0, n, new ParallelOptions { MaxDegreeOfParallelism = _options.MaxDegreeOfParallelism }, i =>
        {
            if (_weights[i] == 0.0 || steps == 0) return;
            var rng = new Xoshiro256PlusPlus(seeds[i + 1]);
            TState x = _state[i];
            ITargetEvaluation e = _evaluations[i];
            for (int s = 0; s < steps; s++)
            {
                SmcMove<TState> move = _kernel.Move(x, e, next, rng);
                x = move.State;
                e = move.Evaluation;
                if (move.Accepted) accepted[i]++;
            }
            _state[i] = x;
            _evaluations[i] = e;
            _logDensity[i] = next.LogDensity(e);
        });

        long acceptedTotal = 0, attempted = 0;
        for (int i = 0; i < n; i++)
        {
            acceptedTotal += accepted[i];
            if (_weights[i] > 0.0) attempted += steps;
        }
        var lineages = new HashSet<int>(_families);
        _reports.Add(new SmcStageReport(StagesDone, next.Temperature, ess, conditionalEss, resampled, distinct, lineages.Count,
            attempted > 0 ? (double)acceptedTotal / attempted : double.NaN, logIncrement));
    }

    /// <summary>An independent copy of the run's state, forking every particle.</summary>
    public SmcCheckpoint<TState> Checkpoint()
    {
        var particles = new TState[_state.Length];
        for (int i = 0; i < particles.Length; i++) particles[i] = _particles.Fork(_state[i]);
        return new SmcCheckpoint<TState>(_seed, StagesDone, Stage, particles, (ITargetEvaluation[])_evaluations.Clone(),
            (double[])_weights.Clone(), (int[])_families.Clone(), _logNormalizer, _reports.ToArray());
    }

    /// <summary>The population as it stands.</summary>
    public SmcResult<TState> Result() => new((TState[])_state.Clone(), (double[])_weights.Clone(), _logNormalizer, _reports.ToArray());

    // The particles' evaluations under the next stage's target: reused when the target is the current one.
    private ITargetEvaluation[] EvaluationsFor(SmcStage<TState> next)
    {
        if (ReferenceEquals(next.Target, Stage.Target)) return (ITargetEvaluation[])_evaluations.Clone();
        var evaluations = new ITargetEvaluation[_state.Length];
        Parallel.For(0, _state.Length, new ParallelOptions { MaxDegreeOfParallelism = _options.MaxDegreeOfParallelism },
            i => evaluations[i] = next.Target.Evaluate(_state[i]));
        return evaluations;
    }

    // log π_next(x) − log π_current(x); a state outside the next stage's support gets −∞.
    private static double Increment(double next, double current) =>
        double.IsNegativeInfinity(next) ? double.NegativeInfinity : next - current;

    private static void Validate(SmcOptions options)
    {
        ArgumentNullException.ThrowIfNull(options);
        if (options.Particles < 1) throw new ArgumentOutOfRangeException(nameof(options), "The population needs at least one particle.");
        if (!(options.ResampleBelow >= 0.0)) throw new ArgumentOutOfRangeException(nameof(options), "ResampleBelow must be at least 0.");
        if (options.MutationSteps < 0) throw new ArgumentOutOfRangeException(nameof(options), "MutationSteps must be at least 0.");
        if (options.MaxStages < 0) throw new ArgumentOutOfRangeException(nameof(options), "MaxStages must be at least 0.");
        if (options.MaxDegreeOfParallelism < 1) throw new ArgumentOutOfRangeException(nameof(options), "MaxDegreeOfParallelism must be at least 1.");
    }
}

/// <summary>Whole runs: through a fixed stage sequence, or tempered adaptively from b = 0 to 1 on one target.</summary>
public static class SmcSampler
{
    /// <summary>Starts at <paramref name="stages"/>[0] and advances through the rest.</summary>
    public static SmcResult<TState> Run<TState>(IReadOnlyList<SmcStage<TState>> stages, ISmcParticles<TState> particles,
                                                ISmcKernel<TState> kernel, SmcOptions options, int seed)
    {
        ArgumentNullException.ThrowIfNull(stages);
        if (stages.Count == 0) throw new ArgumentException("There are no stages.", nameof(stages));
        SmcRun<TState> run = SmcRun<TState>.Start(stages[0], particles, kernel, options, seed);
        for (int s = 1; s < stages.Count; s++) run.Advance(stages[s]);
        return run.Result();
    }

    /// <summary>
    /// Tempers <paramref name="target"/> along <paramref name="path"/> from b = 0, whose distribution
    /// <paramref name="particles"/> draws from (on the likelihood path, the prior), to b = 1, choosing each next b so that
    /// the reweighting's conditional ESS is <paramref name="targetConditionalEss"/>·N (Zhou, Johansen and Aston 2016).
    /// </summary>
    public static SmcResult<TState> RunTempered<TState>(ISamplingTarget<TState> target, TemperingPath path, double targetConditionalEss,
                                                        ISmcParticles<TState> particles, ISmcKernel<TState> kernel, SmcOptions options, int seed)
    {
        ArgumentNullException.ThrowIfNull(target);
        SmcRun<TState> run = SmcRun<TState>.Start(new SmcStage<TState>(target, new Temperature(0.0, path)), particles, kernel, options, seed);
        return TemperToPosterior(run, targetConditionalEss);
    }

    /// <summary>Continues a run on its current target to b = 1 with the adaptive schedule; it can follow a resume.</summary>
    public static SmcResult<TState> TemperToPosterior<TState>(SmcRun<TState> run, double targetConditionalEss)
    {
        ArgumentNullException.ThrowIfNull(run);
        if (!(targetConditionalEss > 0.0 && targetConditionalEss < 1.0))
            throw new ArgumentOutOfRangeException(nameof(targetConditionalEss), "The target fraction must lie in (0, 1).");
        while (run.Stage.Temperature.InverseTemperature < 1.0)
            run.Advance(new SmcStage<TState>(run.Stage.Target, new Temperature(NextInverseTemperature(run, targetConditionalEss), run.Stage.Temperature.Path)));
        return run.Result();
    }

    /// <summary>
    /// The largest b′ in (b, 1] whose reweighting keeps the conditional ESS at or above the target fraction, by
    /// bisection; the conditional ESS falls as b′ moves away from b. When even the smallest step falls below it, as
    /// when part of the population leaves the support at any b′ &gt; b, the smallest step bisection reaches.
    /// </summary>
    public static double NextInverseTemperature<TState>(SmcRun<TState> run, double targetConditionalEss)
    {
        ArgumentNullException.ThrowIfNull(run);
        double b = run.Stage.Temperature.InverseTemperature;
        TemperingPath path = run.Stage.Temperature.Path;
        ISamplingTarget<TState> target = run.Stage.Target;
        double goal = targetConditionalEss * run.Count;
        double Cess(double candidate) => run.ConditionalEss(new SmcStage<TState>(target, new Temperature(candidate, path)));
        if (Cess(1.0) >= goal) return 1.0;
        double lo = b, hi = 1.0;
        for (int iteration = 0; iteration < 60 && hi - lo > 1e-12; iteration++)
        {
            double mid = 0.5 * (lo + hi);
            if (Cess(mid) >= goal) lo = mid;
            else hi = mid;
        }
        return lo > b ? lo : hi;
    }
}
