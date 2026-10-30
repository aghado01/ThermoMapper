using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Smc;

namespace Maths.Regression.Spline.Bass.Calibration;

/// <summary>
/// The cut SMC's settings (D09, stage 8), values only. Its defaults are provisional: read in the checks, not tuned.
/// </summary>
public sealed record CutSmcSpec
{
    /// <summary>The population size N.</summary>
    public int Particles { get; init; } = 1000;

    /// <summary>Each next stage, on the first leg and on every bridge, keeps the reweighting's conditional ESS at this
    /// fraction of N (Zhou, Johansen and Aston 2016).</summary>
    public double TargetConditionalEss { get; init; } = 0.99;

    /// <summary>Sweeps of the block scan (<see cref="CutSmcKernel"/>) per particle at every stage.</summary>
    public int MutationSteps { get; init; } = 5;

    /// <summary>Resample when ESS/N falls below this fraction.</summary>
    public double ResampleBelow { get; init; } = 0.5;

    /// <summary>The resampling scheme: multinomial, the one whose offspring the normalizer's genealogy error describes.</summary>
    public ResamplingScheme Resampling { get; init; } = ResamplingScheme.Multinomial;

    /// <summary>The emulator draws S the run visits, as the nested sampler's (<see cref="ReferenceSampling.OuterDraws"/>):
    /// null takes every retained draw once; weighted draws are resampled by weight.</summary>
    public int? OuterDraws { get; init; }

    /// <summary>The order the draws are visited in: the tour, or a random order, its comparison.</summary>
    public CutSmcOrder Order { get; init; } = CutSmcOrder.Tour;

    /// <summary>The stages each leg may take, the first leg's and each bridge's: the run's cap is this times S + 1.</summary>
    public int MaxStagesPerLeg { get; init; } = 1000;

    /// <summary>Checks the fields, naming each in the reason.</summary>
    public void Validate()
    {
        if (Particles < 2) throw Bad("Particles", "must be at least 2");
        if (!(TargetConditionalEss > 0.0 && TargetConditionalEss < 1.0)) throw Bad("TargetConditionalEss", "must lie in (0, 1)");
        if (MutationSteps < 1) throw Bad("MutationSteps", "must be at least 1: the bridges need their kernel");
        if (!(ResampleBelow >= 0.0 && ResampleBelow <= 1.0)) throw Bad("ResampleBelow", "must lie in [0, 1]");
        if (!Enum.IsDefined(Resampling)) throw Bad("Resampling", "is not a known scheme");
        if (OuterDraws is < 1) throw Bad("OuterDraws", "must be at least 1");
        if (!Enum.IsDefined(Order)) throw Bad("Order", "is not a known order");
        if (MaxStagesPerLeg < 1) throw Bad("MaxStagesPerLeg", "must be at least 1");
    }

    private static ArgumentException Bad(string field, string reason) => new($"CalibrationConfig.Smc.{field}: {reason}");
}

/// <summary>
/// One leg of a cut SMC run: the first, from the prior to μ_1 on D05's likelihood path, or a bridge from one emulator
/// draw's conditional to the next. Its stages' reports (ESS, conditional ESS, families, the walk's acceptance) and its
/// blocks' acceptances, kept apart from any chain's R̂: equal weights after resampling do not make particles
/// independent, and <see cref="Families"/> says how many initial lineages remain.
/// </summary>
/// <param name="From">The draw the leg leaves; null for the first leg, which leaves the prior.</param>
/// <param name="To">The draw whose conditional the leg reaches.</param>
/// <param name="Distance">The ordering's distance between the two draws (<see cref="CutSmcOrdering.Distances"/>); NaN
/// for the first leg and when the order was not a tour.</param>
/// <param name="FirstStep">A bridge's first t, bisected on the cached increments; NaN for the first leg.</param>
/// <param name="Stages">The leg's stage reports, in order.</param>
/// <param name="Blocks">The leg's block acceptances, summed over its stages.</param>
public sealed record CutSmcLeg(int? From, int To, double Distance, double FirstStep, IReadOnlyList<SmcStageReport> Stages, CutSmcBlockRates Blocks)
{
    /// <summary>The stages the leg took.</summary>
    public int StageCount => Stages.Count;

    /// <summary>The least ESS over the leg's stages, before any resampling.</summary>
    public double MinEss => Stages.Min(s => s.EffectiveSampleSize);

    /// <summary>The least conditional ESS over the leg's stages.</summary>
    public double MinConditionalEss => Stages.Min(s => s.ConditionalEss);

    /// <summary>The initial lineages left at the leg's end.</summary>
    public int Families => Stages[^1].Families;

    /// <summary>The θ walk's acceptance over the leg's moves (the variance steps' without a continuous parameter).</summary>
    public double Acceptance
    {
        get
        {
            double[] rates = [.. Stages.Select(s => s.AcceptanceRate).Where(double.IsFinite)];
            return rates.Length > 0 ? rates.Average() : double.NaN;
        }
    }
}

/// <summary>
/// The population at one emulator draw's conditional μ_s, read when the leg reaching it ends: the draws in physical
/// coordinates (D09's outputs) with their normalized weights. <see cref="LogEvidence"/> is log p(y | e_s) as the run
/// estimates it, the first leg's evidence and each bridge's ratio; the cut never reads it.
/// </summary>
public sealed record CutSmcPopulation(int EmulatorDraw, int Visit, IReadOnlyList<CalibrationDraw> Draws, double[] Weights,
                                      double LogEvidence, double LogEvidenceError, CutSmcLeg Leg);

/// <summary>
/// A cut SMC run's result (D09's outputs): the outer draws, the order they were visited in, one population per visit
/// with its leg's diagnostics, and the pooled cut posterior.
/// </summary>
/// <param name="Parameters">The calibration parameters, in declaration order.</param>
/// <param name="OuterDraws">The emulator draws S the run took, in the outer sample's order.</param>
/// <param name="Order">The visits, as positions in <see cref="OuterDraws"/>; the first is the first leg's.</param>
/// <param name="Distances">The ordering's distances between the positions, for a tour; null for a random order.</param>
/// <param name="Populations">One population per visit, in the order visited.</param>
/// <param name="Stages">Every stage's report, the initial draw first.</param>
public sealed record CutSmcCalibration(IReadOnlyList<CalibrationParameter> Parameters, IReadOnlyList<int> OuterDraws, IReadOnlyList<int> Order,
                                       double[,]? Distances, IReadOnlyList<CutSmcPopulation> Populations, IReadOnlyList<SmcStageReport> Stages)
{
    /// <summary>The tour's length under <see cref="Distances"/>; NaN without them.</summary>
    public double TourLength => Distances is null ? double.NaN : CutSmcOrdering.Length(Distances, [.. Order]);

    /// <summary>The stages taken after the initial draw, over every leg.</summary>
    public int TotalStages => Stages.Count - 1;

    /// <summary>The initial lineages left at the run's end.</summary>
    public int Families => Stages[^1].Families;

    /// <summary>
    /// The cut posterior's draws, p_cut(θ | y) = (1/S)·Σ_s p(θ | y, e_s): each population's draws, weighted by its
    /// normalized weights over S. The weighted emulator draws entered through the outer sample, resampled by weight, so
    /// every visit counts 1/S; no evidence reweights a draw, which keeps the cut.
    /// </summary>
    public (IReadOnlyList<CalibrationDraw> Draws, double[] Weights) Pooled()
    {
        var draws = new List<CalibrationDraw>();
        var weights = new List<double>();
        int s = Populations.Count;
        foreach (CutSmcPopulation population in Populations)
            for (int i = 0; i < population.Draws.Count; i++)
            {
                draws.Add(population.Draws[i]);
                weights.Add(population.Weights[i] / s);
            }
        return (draws, [.. weights]);
    }
}

/// <summary>
/// Cut-calibration SMC (D09, stage 8; Mathews, Gopalan, Gattiker, Smith and Francom 2024): one population carried
/// through every emulator draw's conditional μ_s = p(θ, σ², v | y, e_s) in turn, the draws visited in a tour that keeps
/// successive conditionals close (<see cref="CutSmcOrdering"/>). The first leg tempers the prior to μ_1 on D05's
/// likelihood path, every prior proper; each next draw is reached by a geometric bridge through the draws' induced
/// densities (<see cref="CalibrationBridge"/>), t placed by conditional ESS. A bridge's first step moves to a new target,
/// so e_s is evaluated on the particles once, t₁ bisected on the cached t·Δ, and the run then advanced and tempered to
/// t = 1 on evaluations it reuses. One <see cref="SmcRun{TState}"/> spans the whole tour, so its families track the
/// initial lineages throughout; the mutation is the block scan aimed at each stage (<see cref="CutSmcKernel"/>). Each
/// μ_s's population is kept, and the cut posterior pools them at 1/S each: no stage reweights the emulator, which is
/// the cut (ordinary joint SMC would).
/// </summary>
public static class CutSmcSampler
{
    /// <summary>
    /// Runs the cut SMC on <paramref name="model"/> under its config's <see cref="CalibrationConfig.Smc"/> and the seed
    /// of its <see cref="CalibrationConfig.Sampling"/>. The outer draws are the nested sampler's for the same seed.
    /// Particles move <paramref name="maxDegreeOfParallelism"/> at a time (−1, no bound); each has its own stream, so the
    /// result does not depend on it.
    /// </summary>
    public static CutSmcCalibration Run(CalibrationModel model, int maxDegreeOfParallelism = -1)
    {
        ArgumentNullException.ThrowIfNull(model);
        CutSmcSpec spec = model.Config.Smc ?? throw new ArgumentException(
            "CalibrationConfig.Smc is null: the cut posterior is then the nested reference sampler's (NestedCalibrationSampler.Run).", nameof(model));
        spec.Validate();
        CalibrationSampling sampling = model.Config.Sampling;
        sampling.Validate();
        if (!model.PriorsProper) throw new ArgumentException("The first leg starts from the prior, and a variance group's prior is improper.", nameof(model));
        if (maxDegreeOfParallelism is 0 or < -1) throw new ArgumentOutOfRangeException(nameof(maxDegreeOfParallelism), "−1, no bound, or at least 1.");
        int parallelism = maxDegreeOfParallelism == -1 ? Environment.ProcessorCount : maxDegreeOfParallelism;

        // The outer draws from the nested sampler's stream for the same seed; the run's own streams beside its inner runs'.
        int[] seeds = SeedTree.Derive(sampling.MasterSeed, 3);
        int[] outer = NestedCalibrationSampler.OuterDraws(model.Emulator, spec.OuterDraws, new Xoshiro256PlusPlus(seeds[0]));
        int[] streams = SeedTree.Derive(seeds[2], 2);
        int s = outer.Length;
        var targets = new Dictionary<int, CalibrationTarget>();   // one per draw, so the states' caches carry across legs
        CalibrationTarget Target(int draw) => targets.TryGetValue(draw, out CalibrationTarget? t) ? t : targets[draw] = model.Target(draw);

        var options = new SmcOptions
        {
            Particles = spec.Particles,
            Resampling = spec.Resampling,
            ResampleBelow = spec.ResampleBelow,
            MutationSteps = spec.MutationSteps,
            MaxStages = checked(spec.MaxStagesPerLeg * (s + 1)),
            MaxDegreeOfParallelism = parallelism,
        };
        var kernel = new CutSmcKernel(model);
        CalibrationTarget first = Target(outer[0]);
        SmcRun<CalibrationState> run = SmcRun<CalibrationState>.Start(
            new SmcStage<CalibrationState>(first, new Temperature(0.0, TemperingPath.Likelihood)), new CutSmcParticles(model), kernel, options, streams[0]);
        var populations = new List<CutSmcPopulation>(s);
        Temper(run, spec, first, options.MaxStages);
        populations.Add(Population(run, kernel, first, 0, new CutSmcLeg(null, outer[0], double.NaN, double.NaN, [], Empty), 0, parallelism));

        int[] order;
        double[,]? distances = null;
        if (s == 1) order = [0];
        else if (spec.Order == CutSmcOrder.Random) order = CutSmcOrdering.Random(s, 0, new Xoshiro256PlusPlus(streams[1]));
        else
        {
            SmcResult<CalibrationState> population = run.Result();
            int[] distinct = [.. outer.Distinct().Order()];
            double[,] byDraw = CutSmcOrdering.Distances(
                CutSmcOrdering.LogLikelihoods([.. distinct.Select(Target)], population.Particles, parallelism), population.Weights);
            distances = new double[s, s];
            for (int a = 0; a < s; a++)
                for (int b = 0; b < s; b++) distances[a, b] = byDraw[Array.IndexOf(distinct, outer[a]), Array.IndexOf(distinct, outer[b])];
            order = CutSmcOrdering.Tour(distances, 0);
        }

        for (int k = 1; k < s; k++)
        {
            CalibrationTarget from = Target(outer[order[k - 1]]), to = Target(outer[order[k]]);
            var bridge = new CalibrationBridge(from, to);
            int stagesBefore = run.StagesDone;
            double t1 = FirstStep(run, bridge, spec.TargetConditionalEss, parallelism);
            Guard(run, () => run.Advance(new SmcStage<CalibrationState>(bridge, new Temperature(t1, TemperingPath.Likelihood))), to, options.MaxStages);
            Temper(run, spec, to, options.MaxStages);
            var leg = new CutSmcLeg(from.Draw, to.Draw, distances?[order[k - 1], order[k]] ?? double.NaN, t1, [], Empty);
            populations.Add(Population(run, kernel, to, k, leg, stagesBefore, parallelism));
        }
        return new CutSmcCalibration(model.Config.Parameters, outer, order, distances, populations, run.Reports);
    }

    /// <summary>
    /// A bridge's first t: e_s evaluated on the particles once, then the largest t₁ in (0, 1] whose reweighting by t₁·Δ
    /// keeps the conditional ESS at or above the target fraction, bisected on those cached increments as
    /// <see cref="SmcSampler.NextInverseTemperature"/> bisects; the smallest step it reaches when even that falls below.
    /// The evaluations stay in the particles' caches, so the advance to t₁ reads them again for nothing.
    /// </summary>
    internal static double FirstStep(SmcRun<CalibrationState> run, CalibrationBridge bridge, double targetConditionalEss, int parallelism)
    {
        SmcResult<CalibrationState> population = run.Result();
        int n = population.Particles.Count;
        var difference = new double[n];
        Parallel.For(0, n, new ParallelOptions { MaxDegreeOfParallelism = parallelism },
            i => difference[i] = bridge.EvaluateState(population.Particles[i]).Difference);
        return FirstStep(population.Weights, difference, targetConditionalEss);
    }

    internal static double FirstStep(double[] weights, double[] difference, double targetConditionalEss)
    {
        double goal = targetConditionalEss * weights.Length;
        var increments = new double[difference.Length];
        double Cess(double t)
        {
            for (int i = 0; i < increments.Length; i++) increments[i] = double.IsNegativeInfinity(difference[i]) ? double.NegativeInfinity : t * difference[i];
            return SmcWeights.ConditionalEss(weights, increments);
        }
        if (Cess(1.0) >= goal) return 1.0;
        double lo = 0.0, hi = 1.0;
        for (int iteration = 0; iteration < 60 && hi - lo > 1e-12; iteration++)
        {
            double mid = 0.5 * (lo + hi);
            if (Cess(mid) >= goal) lo = mid;
            else hi = mid;
        }
        return lo > 0.0 ? lo : hi;
    }

    private static readonly CutSmcBlockRates Empty = new(0, 0, 0, 0, 0, 0);

    // Tempers the run's current target to 1, naming the leg if the run's stages run out.
    private static void Temper(SmcRun<CalibrationState> run, CutSmcSpec spec, CalibrationTarget to, int maxStages)
        => Guard(run, () => SmcSampler.TemperToPosterior(run, spec.TargetConditionalEss), to, maxStages);

    private static void Guard(SmcRun<CalibrationState> run, Action step, CalibrationTarget to, int maxStages)
    {
        try
        {
            step();
        }
        catch (InvalidOperationException e) when (run.StagesDone == maxStages)
        {
            throw new InvalidOperationException(
                $"The cut SMC run took all its stages on the leg to emulator draw {to.Draw}: raise CalibrationConfig.Smc.MaxStagesPerLeg, or read the leg's diagnostics.", e);
        }
    }

    // The population at `to`'s conditional as the leg leaves it, with the leg's stages from `stagesBefore` on.
    private static CutSmcPopulation Population(SmcRun<CalibrationState> run, CutSmcKernel kernel, CalibrationTarget to, int visit, CutSmcLeg leg,
                                               int stagesBefore, int parallelism)
    {
        SmcResult<CalibrationState> result = run.Result();
        int n = result.Particles.Count;
        var draws = new CalibrationDraw[n];
        Parallel.For(0, n, new ParallelOptions { MaxDegreeOfParallelism = parallelism }, i =>
        {
            CalibrationState state = result.Particles[i];
            draws[i] = CalibrationDraw.Of(to, state, to.EvaluateState(state));
        });
        SmcStageReport[] stages = [.. run.Reports.Skip(stagesBefore + 1)];
        IReadOnlyList<CutSmcBlockRates> rates = kernel.StageRates;
        CutSmcBlockRates blocks = Empty;
        for (int k = stagesBefore; k < run.StagesDone; k++) blocks = Sum(blocks, rates[k]);
        return new CutSmcPopulation(to.Draw, visit, draws, result.Weights, run.LogNormalizerRatio, result.LogNormalizerError,
                                    leg with { Stages = stages, Blocks = blocks });
    }

    private static CutSmcBlockRates Sum(CutSmcBlockRates a, CutSmcBlockRates b)
        => new(a.ThetaAttempts + b.ThetaAttempts, a.ThetaAccepts + b.ThetaAccepts, a.VarianceAttempts + b.VarianceAttempts,
               a.VarianceAccepts + b.VarianceAccepts, a.SiteAttempts + b.SiteAttempts, a.SiteAccepts + b.SiteAccepts);
}
