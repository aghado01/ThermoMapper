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

    /// <summary>How the population is resampled. Multinomial by default: the scheme whose offspring the normalizer's
    /// genealogy error describes (<see cref="SmcStageReport.NormalizerRelativeVariance"/>). The lower-variance schemes
    /// stay selectable for runs that need only posterior expectations, and leave that error NaN.</summary>
    public ResamplingScheme Resampling { get; init; } = ResamplingScheme.Multinomial;

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

    /// <summary>Record each stage's population (<see cref="SmcStageTrace"/>), for an evidence over the stages. Off by
    /// default.</summary>
    public bool TraceStages { get; init; }
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
/// <param name="NormalizerRelativeVariance">The genealogy estimate of Var(Ẑ)/Z² for log(Z_stage/Z_first), read after
/// the reweighting and before any resampling (<see cref="SmcWeights.GenealogyRelativeVariance"/>); 0 at stage 0. NaN
/// from the stage after the population was first resampled by a scheme other than multinomial, whose offspring the
/// estimate's factors do not describe.</param>
public readonly record struct SmcStageReport(
    int Index, Temperature Temperature, double EffectiveSampleSize, double ConditionalEss, bool Resampled,
    int DistinctAncestors, int Families, double AcceptanceRate, double LogNormalizerIncrement,
    double NormalizerRelativeVariance);

/// <summary>
/// One stage's population, recorded when the run traces its stages (<see cref="SmcOptions.TraceStages"/>), as the stage
/// left it, after its mutation: the per-stage draws a ladder evidence reads (<see cref="Evidence.LadderEvidence"/>).
/// </summary>
/// <param name="Temperature">The stage's temperature.</param>
/// <param name="NewTarget">Whether the stage's target differs from the previous stage's.</param>
/// <param name="Weights">The particles' normalized weights.</param>
/// <param name="Slopes">Each particle's ∂_b log q_b at the stage (<see cref="Temperature.LogDensitySlope"/>); NaN
/// throughout when the target supplies no slope.</param>
/// <param name="LogRatios">Each particle's log q_next − log q_stage, the increments the next reweighting reads; null
/// until the run advances past the stage.</param>
public sealed record SmcStageTrace(Temperature Temperature, bool NewTarget, double[] Weights, double[] Slopes, double[]? LogRatios);

/// <summary>A run's population and its account.</summary>
/// <param name="Particles">The current particles.</param>
/// <param name="Weights">Their normalized weights.</param>
/// <param name="LogNormalizerRatio">log(Z_current/Z_first), the sum of the stages' increments, for the stage densities
/// as the targets supply them. It is a model evidence only when those densities are complete and the first stage is a
/// normalized prior.</param>
/// <param name="Stages">One report per stage, the initial draw first.</param>
/// <param name="StageTraces">One population per stage, the initial draw first, when the run traced them; null
/// otherwise.</param>
public sealed record SmcResult<TState>(
    IReadOnlyList<TState> Particles, double[] Weights, double LogNormalizerRatio, IReadOnlyList<SmcStageReport> Stages,
    IReadOnlyList<SmcStageTrace>? StageTraces = null)
{
    /// <summary>The standard error of <see cref="LogNormalizerRatio"/>, from the last stage's genealogy estimate
    /// (<see cref="SmcStageReport.NormalizerRelativeVariance"/>, <see cref="SmcWeights.LogNormalizerError"/>): NaN unless
    /// every resampling was multinomial.</summary>
    public double LogNormalizerError => SmcWeights.LogNormalizerError(Stages[^1].NormalizerRelativeVariance);
}

/// <summary>Everything a run needs to continue: restarting from it reproduces the uninterrupted run bit for bit.
/// <paramref name="GenealogyHolds"/> says no resampling so far was by a scheme other than multinomial, so the genealogy
/// estimate continues.</summary>
public sealed record SmcCheckpoint<TState>(
    int Seed, int StagesDone, SmcStage<TState> Stage, TState[] Particles, ITargetEvaluation[] Evaluations,
    double[] Weights, int[] Families, double LogNormalizerRatio, SmcStageReport[] Reports, bool GenealogyHolds,
    SmcStageTrace[]? StageTraces = null);

/// <summary>
/// A sequential Monte Carlo population carried through a sequence of stages (Del Moral, Doucet and Jasra 2006; the
/// invitation of Dai, Heng, Jacob and Whiteley 2022): reweight by the ratio of consecutive stage densities, resample
/// when the weights concentrate, then mutate every particle with a kernel that leaves the new stage invariant. It
/// reweights through <see cref="ITargetEvaluation"/>s, so moving to another temperature of the same target queries
/// each particle's evaluation instead of repeating its work. Each stage reports the normalizer's relative variance
/// from the particles' genealogy (Lee and Whiteley 2018), which sees the correlation that shared ancestors leave,
/// while every resampling has been multinomial.
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
    private readonly List<SmcStageTrace>? _traces;   // when the run traces its stages
    private TState[] _state;
    private ITargetEvaluation[] _evaluations;
    private double[] _logDensity;
    private double[] _weights;
    private int[] _families;
    private int _resamplings;   // resampling events so far, the genealogy estimate's R
    private bool _genealogyHolds;   // no resampling so far by a scheme other than multinomial
    private double _logNormalizer;

    private SmcRun(ISmcParticles<TState> particles, ISmcKernel<TState> kernel, SmcOptions options, int seed,
                   SmcStage<TState> stage, TState[] state, ITargetEvaluation[] evaluations, double[] weights,
                   int[] families, int stagesDone, double logNormalizer, List<SmcStageReport> reports, bool genealogyHolds,
                   List<SmcStageTrace>? traces)
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
        foreach (SmcStageReport report in reports)
            if (report.Resampled) _resamplings++;
        _genealogyHolds = genealogyHolds;
        _traces = traces;
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

    /// <summary>One population per stage, the initial draw first, when the run traces its stages; null otherwise.</summary>
    public IReadOnlyList<SmcStageTrace>? StageTraces => _traces;

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
        var reports = new List<SmcStageReport> { new(0, first.Temperature, n, n, false, n, n, double.NaN, 0.0, 0.0) };
        var run = new SmcRun<TState>(particles, kernel, options, seed, first, state, evaluations, weights, families, 0, 0.0, reports,
                                     genealogyHolds: true, options.TraceStages ? [] : null);
        run.RecordStage(newTarget: false);
        return run;
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
        if (options.TraceStages != checkpoint.StageTraces is not null)
            throw new ArgumentException("The checkpoint's run " + (options.TraceStages ? "did not trace" : "traced") + " its stages.", nameof(options));
        var state = new TState[checkpoint.Particles.Length];
        for (int i = 0; i < state.Length; i++) state[i] = particles.Fork(checkpoint.Particles[i]);
        return new SmcRun<TState>(particles, kernel, options, checkpoint.Seed, checkpoint.Stage, state,
            (ITargetEvaluation[])checkpoint.Evaluations.Clone(), (double[])checkpoint.Weights.Clone(),
            (int[])checkpoint.Families.Clone(), checkpoint.StagesDone, checkpoint.LogNormalizerRatio,
            new List<SmcStageReport>(checkpoint.Reports), checkpoint.GenealogyHolds,
            checkpoint.StageTraces is null ? null : new List<SmcStageTrace>(checkpoint.StageTraces));
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
        bool newTarget = !ReferenceEquals(next.Target, Stage.Target);
        if (_traces is not null) _traces[^1] = _traces[^1] with { LogRatios = (double[])increments.Clone() };
        var logWeights = new double[n];
        for (int i = 0; i < n; i++) logWeights[i] = Math.Log(_weights[i]) + increments[i];
        SmcWeights.Normalize(logWeights, _weights);
        double ess = SmcWeights.EffectiveSampleSize(_weights);
        double relativeVariance = _genealogyHolds ? SmcWeights.GenealogyRelativeVariance(_weights, _families, _resamplings) : double.NaN;
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
            _resamplings++;
            if (_options.Resampling != ResamplingScheme.Multinomial) _genealogyHolds = false;
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
            attempted > 0 ? (double)acceptedTotal / attempted : double.NaN, logIncrement, relativeVariance));
        RecordStage(newTarget);
    }

    /// <summary>An independent copy of the run's state, forking every particle.</summary>
    public SmcCheckpoint<TState> Checkpoint()
    {
        var particles = new TState[_state.Length];
        for (int i = 0; i < particles.Length; i++) particles[i] = _particles.Fork(_state[i]);
        return new SmcCheckpoint<TState>(_seed, StagesDone, Stage, particles, (ITargetEvaluation[])_evaluations.Clone(),
            (double[])_weights.Clone(), (int[])_families.Clone(), _logNormalizer, _reports.ToArray(), _genealogyHolds,
            _traces?.ToArray());
    }

    /// <summary>The population as it stands.</summary>
    public SmcResult<TState> Result() => new((TState[])_state.Clone(), (double[])_weights.Clone(), _logNormalizer, _reports.ToArray(),
                                             _traces?.ToArray());

    // The stage's population as it now stands, when the run traces its stages. Its arrays are never written again.
    private void RecordStage(bool newTarget)
    {
        if (_traces is null) return;
        int n = _state.Length;
        Temperature t = Stage.Temperature;
        var slopes = new double[n];
        if (SuppliesSlope(t, _evaluations[0])) for (int i = 0; i < n; i++) slopes[i] = t.LogDensitySlope(_evaluations[i]);
        else Array.Fill(slopes, double.NaN);
        _traces.Add(new SmcStageTrace(t, newTarget, (double[])_weights.Clone(), slopes, null));
    }

    private static bool SuppliesSlope(Temperature t, ITargetEvaluation evaluation)
    {
        try
        {
            _ = t.LogDensitySlope(evaluation);
            return true;
        }
        catch (NotSupportedException)
        {
            return false;
        }
    }

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
