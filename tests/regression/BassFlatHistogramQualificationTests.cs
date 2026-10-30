using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using Maths.Distributions;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Smc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// N17 Phase 2's Friedman qualification (the flat-histogram Phase 2 brief, step 4), behind <c>BASS_QUALIFICATION</c> = 1:
/// the flat-histogram sampler on S10's Friedman fixture (n = 500, <see cref="BassMixing.Friedman"/>), each arm from the
/// WHAM warm start of its own self-sized, linear-form tempering pilot. It reads; it sets no default, and asserts only its
/// own soundness (the weights normalized, the basins' fractions summing to one, the curve's ends).
/// <list type="bullet">
/// <item><b>Arm (i):</b> the window from b_lo = 0, the pilot's reference rung, to b = 1, so the full evidence; from the
/// default gain start (the pilot's draws over the chains) and, on the same warm start, from the gain counted from 0.</item>
/// <item><b>Arm (ii):</b> the window from b_lo = 0.10, just below the lower split at b ≈ 0.12–0.15, its pilot's hot end there
/// with no reference rung: the posterior and ln Z(1) − ln Z(0.10). Fitting SMC from the prior stops at b = 0.10 on the
/// linear form with no production change (the SMC core, as the reference harness drives it), so its ln Ẑ(0.10) reads
/// ln Z(1) through the ratio.</item>
/// </list>
/// <para><b>What would count, written before the run (2026-10-30)</b>, so the reading is not fitted to its numbers:</para>
/// <list type="number">
/// <item><b>The chains tunnel:</b> every production chain completes round trips between the window's end bins. The runs are
/// sized for at least 20 per chain; a chain with none has not crossed the window, and its arm's agreement and evidence then
/// read nothing about Friedman's posterior.</item>
/// <item><b>They agree:</b> the chains' own ln Ẑ_c(1) (arm i; arm ii's ratios) lie within about a nat of one another, read
/// as an SD over chains of at most 1 nat; and for each labelled basin holding 5% or more of the pooled weight, the chains'
/// reweighted fractions pass a χ² test of equality against their own batch-means errors (20 batches per chain) at
/// p ≥ 0.01.</item>
/// <item><b>They reach the bound:</b> the pooled ln Ẑ(1) sits at or above −864.6, the 99.9% lower bound that R2's
/// from-prior normalizers put on ln p(y) by Markov's inequality, against the −887 of the pilots' own structures. Arm (ii)'s
/// ln Ẑ(1) through SMC's ln Ẑ(0.10) is read against the same bound, the errors added in quadrature.</item>
/// </list>
/// Agreement below the bound would repeat the reference path's reading, a sampler held on part of the posterior; agreement
/// at or above it, with every chain tunnelling, would be the first sampler that computes Friedman's structural posterior.
/// </summary>
public sealed class BassFlatHistogramQualificationTests
{
    // Master seeds no pilot or earlier harness uses: the band measurement and the reference path took 53 and its tree, 2026
    // and 4101–4202; this chip's checks 7101–7804, its sizing 8101.
    private const int ArmOneSeed = 8017, ArmTwoSeed = 8027, SmcSeed = 8039;
    private const double LowerSplit = 0.10;
    private const double Bound = -864.6, PilotsOwn = -887.15;
    private const int Chains = 4;

    // The lengths, the same for every arm. Sized on master seed 8101, which no arm uses: there a chain swept about 30,000
    // times a second, learning's 1M sweeps a chain made one or two round trips, and production's 200,000 none; so about 30
    // minutes an arm buys 45M sweeps a chain, 20M learning and 25M production, kept every 200th (memory) for 125,000 draws.
    private const long LearningIterations = 100_000;
    private const int ProductionDraws = 125_000, Thin = 200;

    private static readonly string[] Labels = ["A", "B", "G", "S", "T", "a1", "a2", "a5"];
    private static readonly Dictionary<string, string> Basins = new()
    {
        ["A"] = "x1 + 2 × x1·x2 + x2 + 2 × x3 + x4 + x5",
        ["B"] = "x1 + 3 × x1·x2 + x2 + 2 × x3 + x4 + x5",
        ["G"] = "x1 + 5 × x1·x2 + 2 × x3 + x4 + x5",
        ["S"] = "2 × x1 + 3 × x1·x2 + 2 × x2 + 2 × x3 + x4 + x5",
        ["T"] = "x1 + 3 × x1·x2 + 2 × x2 + 2 × x3 + x4 + x5",
        ["a1"] = "4 × x1·x2 + 2 × x3 + x4 + x5",
        ["a2"] = "4 × x1·x2 + x2 + 2 × x3 + x4 + x5",
        ["a5"] = "x1 + 2 × x1·x2 + 2 × x2 + 2 × x3 + x4 + x5",
    };

    private readonly ITestOutputHelper _out;
    public BassFlatHistogramQualificationTests(ITestOutputHelper output) => _out = output;

    private static readonly Lazy<BassMixing.Fixture> Fixture = new(() => BassMixing.Get(BassMixing.Friedman), LazyThreadSafetyMode.ExecutionAndPublication);

    // The pilot: four linear-form ladders, self-sizing, burned 5,000 scans and traced 2,000, as R1 burned S10's ladders;
    // with the reference rung (arm i), or with the hot end at b_lo and none (arm ii).
    private static BassFlatHistogramPilot Pilot(bool reference) => new()
    {
        Chains = 4, BurnIn = 5_000, TracedScans = 2_000,
        Tempering = new TemperingSpec { Form = TemperedForm.Linear, ReferenceRung = reference, HotEnd = reference ? null : LowerSplit },
    };

    private static BassConfig ArmConfig(bool reference, long? gainStart, int seed)
    {
        BassMixing.Fixture f = Fixture.Value;
        var spec = new BassFlatHistogramSpec
        {
            Run = new FlatHistogramSpec
            {
                Chains = Chains, LearningIterations = LearningIterations, ProductionDraws = ProductionDraws, Start = FlatHistogramStart.Wham,
                GainStart = gainStart,
            },
            Pilot = Pilot(reference),
        };
        return f.Config with { Sampling = f.Config.Sampling with { MasterSeed = seed, Thin = Thin, FlatHistogram = spec } };
    }

    [Theory]
    [InlineData("i", false)]
    [InlineData("i", true)]
    [InlineData("ii", false)]
    public void Friedman_Arm(string arm, bool gainFromZero)
    {
        if (!BassQualification.Enabled) return;   // heavy: about half an hour an arm
        bool one = arm == "i";
        Assert.False(gainFromZero && !one, "only arm (i) runs the gain from 0");
        BassConfig config = ArmConfig(one, gainFromZero ? 0L : null, one ? ArmOneSeed : ArmTwoSeed);
        string label = $"Arm ({arm}), {(gainFromZero ? "the gain from 0" : "the default gain start")}, master seed {config.Sampling.MasterSeed}";
        RunArm(label, config, one ? null : SmcAtLowerSplit.Value);
    }

    // ---- A second master seed (the main session, 2026-10-30) ----

    // Arm (i) at the default gain start on a second master seed, added after the first run, which read one seed per arm and
    // 4–17 round trips a chain against the 20 sized for. Its sweeps are doubled through Thin (400), so the learning's
    // iterations, its gain sequence and the draws kept are the first run's, and only the mixing per step grows.
    // What would count, written before its run (2026-10-30):
    // 1. the first run's three criteria, read on this seed alone: every chain tunnels; the chains' ln Ẑ_c(1) SD at most 1
    //    nat, with a1 (or any labelled basin at 5% or more) passing the χ² at p ≥ 0.01; the pooled ln Ẑ(1) at or above −864.6;
    // 2. agreement with seed 8017's arm (i): ln Ẑ(1) within 1 nat plus both runs' errors of −841.36 ± 0.17, and a1's pooled
    //    share within 5 points of 65.9%.
    private const int SecondSeed = 8217;

    [Fact]
    public void Friedman_ArmOne_SecondSeed_TwiceTheSweeps()
    {
        if (!BassQualification.Enabled) return;   // heavy: about half an hour
        BassConfig first = ArmConfig(true, null, SecondSeed);
        BassConfig config = first with { Sampling = first.Sampling with { Thin = 2 * Thin } };
        RunArm($"Arm (i), the default gain start, master seed {SecondSeed}, twice the sweeps (Thin {2 * Thin})", config, null);
    }

    // The run above was refused at its start (2026-10-30, 48 s): SAMC's fixed point could not flatten this seed's pilot ĝ
    // over the window [−1552.56, −773.559], K = 128, b_lo = 0, no node empty, a node's occupancy 125 from flat. So the WHAM
    // warm start is not robust across seeds on Friedman, the refusal the Phase 2 report foresaw. The claim under test is the
    // posterior and its evidence, which the start only speeds, so the seed is rerun from a flat start over the window, K and
    // b_lo its own pilot resolved, with the gain from 0 and the same lengths. What would count is unchanged from above,
    // written before this rerun.
    [Fact]
    public void Friedman_ArmOne_SecondSeed_FlatStart()
    {
        if (!BassQualification.Enabled) return;   // heavy: about half an hour
        BassConfig first = ArmConfig(true, null, SecondSeed);
        BassFlatHistogramSpec spec = first.Sampling.FlatHistogram! with
        {
            Run = first.Sampling.FlatHistogram!.Run with
            {
                Start = FlatHistogramStart.Flat, Window = new FlatHistogramWindow(-1552.56, -773.559), LowSlope = 0.0, Nodes = 128,
            },
            Pilot = null,
        };
        BassConfig config = first with { Sampling = first.Sampling with { Thin = 2 * Thin, FlatHistogram = spec } };
        RunArm($"Arm (i), a flat start over seed {SecondSeed}'s pilot window, master seed {SecondSeed}, twice the sweeps (Thin {2 * Thin})", config, null);
    }

    private void RunArm(string label, BassConfig config, SmcAt? smc)
    {
        BassMixing.Fixture fixture = Fixture.Value;
        BassFlatHistogramPilot? pilot = config.Sampling.FlatHistogram!.Pilot;
        FlatHistogramSpec run = config.Sampling.FlatHistogram.Run;
        long pilotScans = pilot is null ? 0 : (long)pilot.Chains * (pilot.BurnIn + pilot.TracedScans);
        long learningSweeps = run.Chains * run.LearningIterations * config.Sampling.Thin;
        var clock = Stopwatch.StartNew();
        double pilotEnd = double.NaN, learningEnd = double.NaN;
        var progress = new Progress(p =>
        {
            double now = clock.Elapsed.TotalSeconds;
            if (double.IsNaN(pilotEnd) && p.Sweeps >= pilotScans) pilotEnd = now;
            if (double.IsNaN(learningEnd) && p.Sweeps >= pilotScans + learningSweeps) learningEnd = now;
        });
        BassFit fit;
        try
        {
            fit = Bass.Fit(config, fixture.Data, new BassRunOptions { Progress = progress, ProgressInterval = 250_000 });
        }
        catch (InvalidOperationException e) when (e.Message.Contains("WHAM start", StringComparison.Ordinal))
        {
            _out.WriteLine($"## {label}: SAMC's fixed point refused the pilot's ĝ after {clock.Elapsed.TotalSeconds:F0} s");
            _out.WriteLine(e.Message);
            return;
        }
        double total = clock.Elapsed.TotalSeconds;
        Read(label, fit, fixture, (pilotEnd, learningEnd - pilotEnd, total - learningEnd, total), smc);
    }

    [Fact]
    public void Friedman_SmcFromThePrior_ToTheLowerSplit()
    {
        if (!BassQualification.Enabled) return;
        SmcAt smc = SmcAtLowerSplit.Value;
        _out.WriteLine($"## Fitting SMC from the prior to b = {LowerSplit} on the linear form: {SmcParticles:N0} particles, {smc.Replicates.Length} replicates " +
                       $"from master seed {SmcSeed}, C4's other defaults; {smc.Seconds:F0} s");
        _out.WriteLine("replicate | ln Ẑ(0.10) ± genealogy error | stages | families at the end");
        _out.WriteLine("--- | --- | --- | ---");
        for (int r = 0; r < smc.Replicates.Length; r++)
            _out.WriteLine($"{r} | {smc.Replicates[r]:F2} ± {smc.Errors[r]:F2} | {smc.Stages[r]} | {smc.Families[r]}");
        _out.WriteLine($"pooled {smc.LogEvidence:F2} ± {smc.Error:F2} (within {smc.Within:F2}, between {smc.Between:F2}); sd {Sd(smc.Replicates):F2}");
        Assert.All(smc.Replicates, v => Assert.True(double.IsFinite(v)));
    }

    // ---- SMC from the prior to b_lo, through the SMC core as the reference harness drives it (no production change) ----

    private const int SmcParticles = 4_000, SmcReplicates = 4;

    private sealed record SmcAt(double LogEvidence, double Error, double Within, double Between, double[] Replicates, double[] Errors, int[] Stages,
                                int[] Families, double Seconds, double[][]? Marks = null);

    private static readonly Lazy<SmcAt> SmcAtLowerSplit = new(() => SmcTo(LowerSplit), LazyThreadSafetyMode.ExecutionAndPublication);

    // ln Ẑ(b) = ln E_prior[e^{b·ℓ₁}] on the prior on F by fitting SMC's loop (C4's defaults) from exact prior draws to b on
    // the linear form, each next b the conditional-ESS rule's, capped at b: each replicate's normalizer, pooled by the larger
    // of the within-run and between-run errors. With `marks`, a stage is forced at each, and each replicate's ln Ẑ there kept.
    private static SmcAt SmcTo(double b, double[]? marks = null)
    {
        BassMixing.Fixture f = Fixture.Value;
        BassModel linear = BassModel.Create(f.Config, f.Data, TemperedForm.Linear);
        var defaults = new SmcSpec();
        int[] seeds = SeedTree.Derive(SmcSeed, SmcReplicates);
        var values = new double[SmcReplicates];
        var errors = new double[SmcReplicates];
        var stages = new int[SmcReplicates];
        var families = new int[SmcReplicates];
        double[][]? atMarks = marks is null ? null : [.. marks.Select(_ => new double[SmcReplicates])];
        var clock = Stopwatch.StartNew();
        int inner = Math.Max(1, Environment.ProcessorCount / SmcReplicates);
        Parallel.For(0, SmcReplicates, r =>
        {
            int[] streams = SeedTree.Derive(seeds[r], 2);   // the run's, as BassSmcFit derives it
            var target = new BassParticleTarget(linear.Target);
            var options = new SmcOptions
            {
                Particles = SmcParticles, Resampling = defaults.Resampling, ResampleBelow = defaults.ResampleBelow, MutationSteps = defaults.MutationSteps,
                MaxStages = defaults.MaxStages, MaxDegreeOfParallelism = inner,
            };
            SmcRun<BassParticle> run = SmcRun<BassParticle>.Start(new SmcStage<BassParticle>(target, new Temperature(0.0, TemperingPath.Likelihood)),
                                                                  new BassParticles(linear), new BassSmcKernel(linear), options, streams[0]);
            while (run.Stage.Temperature.InverseTemperature < b)
            {
                double now = run.Stage.Temperature.InverseTemperature;
                double next = Math.Min(SmcSampler.NextInverseTemperature(run, defaults.TargetConditionalEss), b);
                if (marks is not null && marks.FirstOrDefault(m => m > now, double.NaN) is double mark && !double.IsNaN(mark)) next = Math.Min(next, mark);
                run.Advance(new SmcStage<BassParticle>(target, new Temperature(next, TemperingPath.Likelihood)));
                int at = marks is null ? -1 : Array.IndexOf(marks, next);
                if (at >= 0) atMarks![at][r] = run.Result().LogNormalizerRatio;
            }
            SmcResult<BassParticle> result = run.Result();
            values[r] = result.LogNormalizerRatio;
            errors[r] = result.LogNormalizerError;
            stages[r] = result.Stages.Count - 1;
            families[r] = result.Stages[^1].Families;
        });
        double within = Math.Sqrt(errors.Sum(e => e * e)) / SmcReplicates, between = Sd(values) / Math.Sqrt(SmcReplicates);
        return new SmcAt(values.Average(), Math.Max(within, between), within, between, values, errors, stages, families, clock.Elapsed.TotalSeconds,
                         atMarks);
    }

    // ---- Added after the run (2026-10-30): where the arms' ln Ẑ(0.10) parts ----

    private static readonly double[] Marks = [0.01, 0.02, 0.03, 0.05, 0.07, 0.10];

    /// <summary>
    /// Added after the run, to locate a difference it found: arm (i) read ln Ẑ(0.10) ≈ −125.9 and SMC from the prior
    /// −131.35 ± 0.21, while the arms agree on ln Z(1) − ln Z(0.10). Below b = 0.10 three estimates of ln Z(b) on the linear
    /// form are read side by side: arm (i)'s curve (the same run as <see cref="Friedman_Arm"/>'s, bit for bit), SMC from
    /// the prior with a stage forced at each mark, and plain importance sampling from exact prior draws, reliable where b
    /// times ℓ₁'s spread under the prior is small.
    /// </summary>
    [Fact]
    public void Friedman_ArmOneAgainstSmcAndThePrior_BelowTheLowerSplit()
    {
        if (!BassQualification.Enabled) return;   // heavy: arm (i) again
        BassMixing.Fixture fixture = Fixture.Value;
        BassModel linear = BassModel.Create(fixture.Config, fixture.Data, TemperedForm.Linear);

        // Plain importance sampling from the prior on F: ln mean e^{b·ℓ₁} over exact draws, its error from 8 streams' spread.
        const int streams = 8, perStream = 25_000;
        var levels = new double[streams][];
        int[] seeds = SeedTree.Derive(SmcSeed + 1, streams);
        Parallel.For(0, streams, s =>
        {
            var sampler = new BassPriorSampler(linear);
            var rng = new Xoshiro256PlusPlus(seeds[s]);
            levels[s] = new double[perStream];
            for (int i = 0; i < perStream; i++) levels[s][i] = sampler.DrawState(rng).Evaluation.LogMarginal;
        });
        double[] all = [.. levels.SelectMany(l => l)];
        double intercept = linear.Target.Evaluate(BassMixing.Loaded(linear, [], 0.0, out _)!).LogMarginal;
        _out.WriteLine($"## Below the lower split: prior draws' ℓ₁ (n = {all.Length:N0}): min {all.Min():F1}, 1% {BassMixing.Quantile(all, 0.01):F1}, " +
                       $"50% {BassMixing.Quantile(all, 0.5):F1}, 99% {BassMixing.Quantile(all, 0.99):F1}, max {all.Max():F1}; the intercept alone {intercept:F1}");

        SmcAt smc = SmcTo(LowerSplit, Marks);
        var clock = Stopwatch.StartNew();
        BassFit fit = Bass.Fit(ArmConfig(true, null, ArmOneSeed), fixture.Data);
        double armSeconds = clock.Elapsed.TotalSeconds;
        BassFlatHistogramReadout fh = fit.Diagnostics.FlatHistogram!;
        LogEvidenceCurve curve = fh.Evidence!.Curve;
        Assert.Equal(-841.36, fh.Evidence.LogRatio, 2);   // the qualification's run
        _out.WriteLine($"arm (i) rerun in {armSeconds:F0} s; window [{fh.Window.Low:F1}, {fh.Window.High:F1}]; the denominator's Kish ESS " +
                       $"{curve.NormalizerEffectiveSampleSize:F0} of {fit.DrawWeights!.Length:N0}; production below the window " +
                       $"{fh.Production!.BelowWindow:N0}, in its lowest bin's chains {string.Join(", ", fh.ProductionChains.Select(c => c.Flatness[0].ToString("F2")))} of flat");
        _out.WriteLine("b | prior IS ± error | SMC (forced marks) per replicate, pooled ± error | arm (i) pooled ± error | arm (i) per chain | arm (i) Kish ESS");
        _out.WriteLine("--- | --- | --- | --- | --- | ---");
        for (int m = 0; m < Marks.Length; m++)
        {
            double b = Marks[m];
            double[] perStreamZ = [.. levels.Select(l => LogMeanExp(l, b))];
            double prior = LogMeanExp(all, b), priorError = Sd(perStreamZ) / Math.Sqrt(streams);
            int k = (int)Math.Round(b * 100);
            Assert.Equal(b, curve.InverseTemperatures[k], 12);
            double[] chains = [.. curve.ChainLogEvidence.Select(c => c[k])];
            double[] smcAt = smc.Marks![m];
            _out.WriteLine($"{b:F2} | {prior:F2} ± {priorError:F2} | {string.Join(", ", smcAt.Select(v => v.ToString("F2")))}; {smcAt.Average():F2} ± " +
                           $"{Sd(smcAt) / Math.Sqrt(smcAt.Length):F2} | {curve.LogEvidence[k]:F2} ± {curve.MonteCarloError[k]:F2} | " +
                           $"{string.Join(", ", chains.Select(v => v.ToString("F2")))} | {curve.EffectiveSampleSize[k]:F0}");
        }
        _out.WriteLine($"SMC with the marks forced: ln Ẑ(0.10) {smc.LogEvidence:F2} ± {smc.Error:F2}; stages {string.Join(", ", smc.Stages)}; families " +
                       $"{string.Join(", ", smc.Families)}; {smc.Seconds:F0} s");
    }

    private static double LogMeanExp(IReadOnlyList<double> levels, double b)
    {
        double max = levels.Max(l => b * l);
        return max + Math.Log(levels.Average(l => Math.Exp(b * l - max)));
    }


    // ---- The readings ----

    private void Read(string label, BassFit fit, BassMixing.Fixture fixture, (double Pilot, double Learning, double Production, double Total) seconds, SmcAt? smc)
    {
        BassFlatHistogramReadout fh = fit.Diagnostics.FlatHistogram!;
        BassFlatHistogramPilotReadout? pilot = fh.Pilot;
        BassFlatHistogramEvidence curve = fh.Evidence!;
        double[] w = fit.DrawWeights!;
        int chains = fit.Draws.Count, perChain = fit.Draws[0].Count;
        FlatHistogramSpec run = fit.Config.Sampling.FlatHistogram!.Run;
        long learning = run.LearningIterations;
        int production = run.ProductionDraws, thin = fit.Config.Sampling.Thin;
        Assert.True(fh.Complete);
        Assert.Equal(1.0, w.Sum(), 9);

        // The run.
        _out.WriteLine($"## {label}");
        _out.WriteLine($"- window [{fh.Window.Low:F1}, {fh.Window.High:F1}] ({fh.Window.High - fh.Window.Low:F0} nats), K = {fh.Nodes}, b_lo = {fh.LowSlope}, " +
                       $"t₀ = {fh.GainScale}, the gain counted from {fh.GainStart:N0}");
        _out.WriteLine($"- budget: {chains} chains × ({learning:N0} learning + {production:N0} production iterations) × {thin} sweeps " +
                       $"= {(learning + production) * thin:N0} sweeps a chain; " +
                       (pilot is null ? "no pilot (a flat start over a declared window)" :
                           $"pilot {pilot.Ladders.Levels.Length} ladders × ({Pilot(true).BurnIn:N0} + {Pilot(true).TracedScans:N0}) scans"));
        _out.WriteLine($"- time: pilot {seconds.Pilot:F0} s, learning {seconds.Learning:F0} s, production {seconds.Production:F0} s, all {seconds.Total:F0} s " +
                       $"({(learning + production) * thin / (seconds.Learning + seconds.Production):N0} sweeps a second a chain)");

        // The pilot, and what WHAM's start had to fill.
        if (pilot is not null) ReadPilot(fh, pilot);

        // The learning: its trace, and whether the weights settled.
        ReadRest(label, fit, fixture, fh, curve, w, chains, perChain, learning, smc);
    }

    private void ReadPilot(BassFlatHistogramReadout fh, BassFlatHistogramPilotReadout pilot)
    {
        _out.WriteLine($"- pilot: ladders of {string.Join(", ", pilot.Ladders.Levels)} rungs; {pilot.Ladders.RoundTripsPerScan:F4} round trips per traced scan; " +
                       $"barrier {pilot.Ladders.Barrier:F2}, worst pair rejection {pilot.Ladders.MaxPairRejection:F3}; {pilot.InverseTemperatures.Length} pooled " +
                       $"rungs, {pilot.Draws:N0} draws, inefficiencies {pilot.Inefficiencies.Min():F1}–{pilot.Inefficiencies.Max():F1}");
        _out.WriteLine($"- SAMC's fixed point accepted the pilot's ĝ; empty nodes [{string.Join(", ", pilot.EmptyNodes)}]; nodes under 10 pilot draws " +
                       $"[{string.Join(", ", Enumerable.Range(0, pilot.NodeDraws.Length).Where(k => pilot.NodeDraws[k] < 10).Select(k => $"{k} ({fh.Weights.Node(k):F0})"))}]");
        _out.WriteLine($"- pilot draws per node: {string.Join(" ", pilot.NodeDraws.Select(v => v.ToString("F0")))}");
    }

    private void ReadRest(string label, BassFit fit, BassMixing.Fixture fixture, BassFlatHistogramReadout fh, BassFlatHistogramEvidence curve,
                          double[] w, int chains, int perChain, long learning, SmcAt? smc)
    {
        IReadOnlyList<LearningSnapshot> trace = fh.WeightTrace;
        LearningSnapshot At(double share) => trace.Last(s => s.Iteration <= share * learning);
        _out.WriteLine($"- learning: round trips per chain {string.Join(", ", fh.Learning.RoundTrips)}; flatness {fh.Learning.MaxFlatnessDeviation:F2}; " +
                       $"visits in/below/above {fh.Learning.InWindow:N0}/{fh.Learning.BelowWindow:N0}/{fh.Learning.AboveWindow:N0}");
        _out.WriteLine($"- the weights' centred RMS change: start to 10% {Change(trace[0], At(0.1)):F2}, 10–50% {Change(At(0.1), At(0.5)):F2}, " +
                       $"50–90% {Change(At(0.5), At(0.9)):F2}, 90–100% {Change(At(0.9), trace[^1]):F3}; start to end {Change(trace[0], trace[^1]):F2}");
        _out.WriteLine($"- ln w at the nodes, start: {string.Join(" ", trace[0].Weights.LogWeights.Select(v => v.ToString("F0")))}");
        _out.WriteLine($"- ln w at the nodes, end: {string.Join(" ", fh.Weights.LogWeights.Select(v => v.ToString("F0")))}");

        // Production, per chain.
        double[] chainEvidence = [.. curve.Curve.ChainLogEvidence.Select(c => c[^1] - c[0])];
        double[] logPosterior = [.. Enumerable.Range(0, chains).Select(c => Reweighting.Expectation(fh.Draws!, 1.0, fit.Config.Coefficients.Scale.IsHeld ? 1 : 2).ChainMeans[c])];
        _out.WriteLine($"- production: round trips {fh.Production!.TotalRoundTrips:N0}, chains without one [{string.Join(", ", fh.Production.ChainsWithoutRoundTrip)}]; " +
                       $"flatness {fh.Production.MaxFlatnessDeviation:F2}; b = 1 Kish ESS {fh.EffectiveSampleSize:F0} of {chains * perChain:N0}");
        _out.WriteLine("chain | round trips | flatness | in / below / above | b = 1 Kish ESS | " + (curve.IsEvidence ? "ln Ẑ_c(1)" : $"ln Ẑ_c(1) − ln Ẑ_c({curve.From})") +
                       " | reweighted log posterior");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | ---");
        for (int c = 0; c < chains; c++)
        {
            FlatHistogramDiagnostics d = fh.ProductionChains[c];
            _out.WriteLine($"{c} | {d.RoundTrips[0]} | {d.MaxFlatnessDeviation:F2} | {d.InWindow:N0} / {d.BelowWindow:N0} / {d.AboveWindow:N0} | " +
                           $"{fh.ChainEffectiveSampleSizes[c]:F0} | {chainEvidence[c]:F2} | {logPosterior[c]:F1}");
        }
        for (int c = 0; c < chains; c++)
            _out.WriteLine($"  occupancy over flat, chain {c}: {string.Join(" ", fh.ProductionChains[c].Flatness.Select(v => v.ToString("F1")))}");

        // The evidence or the ratio, and through SMC's ln Ẑ(b_lo).
        double logZ = curve.IsEvidence ? curve.LogRatio : smc is null ? double.NaN : smc.LogEvidence + curve.LogRatio;
        double logZError = curve.IsEvidence ? curve.MonteCarloError : smc is null ? double.NaN : Math.Sqrt(curve.MonteCarloError * curve.MonteCarloError + smc.Error * smc.Error);
        _out.WriteLine(curve.IsEvidence
            ? $"- ln Ẑ(1) pooled {curve.LogRatio:F2} ± {curve.MonteCarloError:F2}; chains' sd {Sd(chainEvidence):F2}, range {chainEvidence.Max() - chainEvidence.Min():F2}; " +
              $"against the bound {Bound} and the pilots' own structures' {PilotsOwn}"
            : $"- ln Ẑ(1) − ln Ẑ({curve.From}) pooled {curve.LogRatio:F2} ± {curve.MonteCarloError:F2} (a ratio); chains' sd {Sd(chainEvidence):F2}, range " +
              $"{chainEvidence.Max() - chainEvidence.Min():F2}; through SMC's ln Ẑ({curve.From}) {smc?.LogEvidence:F2} ± {smc?.Error:F2}: ln Ẑ(1) {logZ:F2} ± {logZError:F2}, " +
              $"against the bound {Bound}");
        _out.WriteLine($"- ln Ẑ(b) on the grid: {string.Join(", ", new[] { 0, 10, 20, 25, 30, 40, 60, 80, 100 }.Select(k => $"{curve.Curve.InverseTemperatures[k]:F3}: {curve.Curve.LogEvidence[k]:F1}"))}");

        // The basins, reweighted: the labels' fractions and log masses, each chain's, and the χ² of their agreement.
        string[][] basins = [.. fit.Draws.Select(store => store.Select(d => fixture.Basin(d.Structure.TermIds.Select(id => fit.Terms[id]))).ToArray())];
        var pooledShare = new Dictionary<string, double>();
        int i = 0;
        foreach (string[] chain in basins)
            foreach (string basin in chain) pooledShare[basin] = pooledShare.GetValueOrDefault(basin) + w[i++];
        Assert.Equal(1.0, pooledShare.Values.Sum(), 9);
        _out.WriteLine("basin | pooled | log mass | " + string.Join(" | ", Enumerable.Range(0, chains).Select(c => $"chain {c}")) + " | χ² p (5%+) | draws visiting");
        _out.WriteLine("--- | --- | --- | " + string.Join(" | ", Enumerable.Repeat("---", chains)) + " | --- | ---");
        IEnumerable<(string Label, string Basin)> rows = Labels.Select(l => (l, Basins[l]))
            .Concat(pooledShare.Where(x => !Basins.ContainsValue(x.Key)).OrderByDescending(x => x.Value).Take(6).Select(x => ("—", x.Key)));
        foreach ((string name, string basin) in rows)
        {
            double share = pooledShare.GetValueOrDefault(basin);
            (double[] fractions, double[] errors) = ChainFractions(basins, fit.Draws.Select(s => s.Count).ToArray(), w, basin);
            string p = share >= 0.05 ? AgreementP(fractions, errors).ToString("F3") : "—";
            int visits = basins.Sum(chain => chain.Count(b => b == basin));
            _out.WriteLine($"{name} {(name == "—" ? basin : "")} | {share:P1} | {(share > 0 ? $"{logZ + Math.Log(share):F1}" : "—")} | " +
                           string.Join(" | ", fractions.Select((f, c) => $"{f:P1} ± {errors[c]:P1}")) + $" | {p} | {visits:N0}");
        }
        var byPairs = new SortedDictionary<int, double>();
        foreach ((string basin, double share) in pooledShare) byPairs[Pairs(basin)] = byPairs.GetValueOrDefault(Pairs(basin)) + share;
        _out.WriteLine($"- share by the count of x1·x2 terms: {string.Join(", ", byPairs.Select(x => $"{x.Key}: {x.Value:P1}"))}");
        for (int c = 0; c < chains; c++)
        {
            var own = new SortedDictionary<int, double>();
            double mass = 0.0;
            for (int s = 0; s < perChain; s++)
            {
                double weight = w[c * perChain + s];
                own[Pairs(basins[c][s])] = own.GetValueOrDefault(Pairs(basins[c][s])) + weight;
                mass += weight;
            }
            _out.WriteLine($"  chain {c}: {string.Join(", ", own.Select(x => $"{x.Key}: {x.Value / mass:P1}"))}");
        }
        _out.WriteLine($"- a1 or a2 visited: {(basins.Any(chain => chain.Any(b => b == Basins["a1"] || b == Basins["a2"])) ? "yes" : "no")}, " +
                       $"chains [{string.Join(", ", Enumerable.Range(0, chains).Where(c => basins[c].Any(b => b == Basins["a1"] || b == Basins["a2"])))}]");

        // σ, and D10's verdict on π_w's functionals.
        double[] sigma = [.. fit.Draws.SelectMany(s => s).Select(d => Math.Sqrt(d.Sigma2))];
        PosteriorSummary interval = PosteriorSummary.Of(sigma, w, [0.025, 0.975]);
        _out.WriteLine($"- σ's 95% interval [{interval.Quantiles[0]:F4}, {interval.Quantiles[1]:F4}] (the truth is 1)");
        ConvergenceVerdict verdict = fit.Diagnostics.Verdict;
        _out.WriteLine($"- D10's verdict on π_w's functionals: {verdict.Outcome}; " + string.Join("; ", verdict.Groups.Select(g =>
            $"{g.Name} {g.Outcome}{(g.Failing.Count > 0 ? " (" + string.Join(", ", g.Failing.Select(r => $"{r.Name} R̂ {r.RHat.Value:F2}")) + ")" : "")}")));
        _out.WriteLine($"- reading against what would count: every chain tunnels {(fh.Production.ChainsWithoutRoundTrip.Length == 0 ? "yes" : "no")}" +
                       $" (fewest {fh.ProductionChains.Min(d => d.RoundTrips[0])}); chains' ln Ẑ_c sd {Sd(chainEvidence):F2} nats; " +
                       $"ln Ẑ(1) {logZ:F2} {(logZ >= Bound ? "at or above" : "below")} {Bound}");
    }

    // The centred RMS change of ln w between two snapshots.
    private static double Change(LearningSnapshot a, LearningSnapshot b)
    {
        int k = a.Weights.Nodes;
        double shift = Enumerable.Range(0, k).Average(j => b.Weights.LogWeights[j] - a.Weights.LogWeights[j]);
        return Math.Sqrt(Enumerable.Range(0, k).Average(j => Math.Pow(b.Weights.LogWeights[j] - a.Weights.LogWeights[j] - shift, 2)));
    }

    // Each chain's reweighted fraction of `basin` and its batch-means error, 20 batches of its draws in order.
    private static (double[] Fractions, double[] Errors) ChainFractions(string[][] basins, int[] counts, double[] w, string basin)
    {
        const int batches = 20;
        int chains = basins.Length, offset = 0;
        var fractions = new double[chains];
        var errors = new double[chains];
        for (int c = 0; c < chains; c++)
        {
            int n = counts[c], size = n / batches;
            double hit = 0.0, mass = 0.0;
            var batch = new double[batches];
            for (int b = 0; b < batches; b++)
            {
                double bh = 0.0, bm = 0.0;
                for (int s = b * size; s < (b + 1) * size; s++)
                {
                    double weight = w[offset + s];
                    bm += weight;
                    if (basins[c][s] == basin) bh += weight;
                }
                batch[b] = bm > 0.0 ? bh / bm : 0.0;
            }
            for (int s = 0; s < n; s++)
            {
                mass += w[offset + s];
                if (basins[c][s] == basin) hit += w[offset + s];
            }
            fractions[c] = hit / mass;
            errors[c] = Sd(batch) / Math.Sqrt(batches);
            offset += n;
        }
        return (fractions, errors);
    }

    // The χ² test of the chains' fractions' equality against their own errors, C − 1 degrees of freedom: its p-value; 0
    // where a chain's error is 0 and its fraction departs.
    private static double AgreementP(double[] fractions, double[] errors)
    {
        double precision = 0.0, sum = 0.0;
        for (int c = 0; c < fractions.Length; c++)
        {
            if (errors[c] == 0.0) continue;
            precision += 1.0 / (errors[c] * errors[c]);
            sum += fractions[c] / (errors[c] * errors[c]);
        }
        if (precision == 0.0) return fractions.All(f => f == fractions[0]) ? 1.0 : 0.0;
        double mean = sum / precision, chi2 = 0.0;
        for (int c = 0; c < fractions.Length; c++)
        {
            double d = fractions[c] - mean;
            if (errors[c] == 0.0)
            {
                if (Math.Abs(d) > 1e-12) return 0.0;
                continue;
            }
            chi2 += d * d / (errors[c] * errors[c]);
        }
        return SpecialFunctions.RegularizedGammaQ((fractions.Length - 1) / 2.0, chi2 / 2.0);
    }

    // The basin's count of x1·x2 terms.
    private static int Pairs(string basin)
        => basin.Split(" + ").Sum(t => t == "x1·x2" ? 1 : t.EndsWith(" × x1·x2", StringComparison.Ordinal) ? int.Parse(t[..t.IndexOf(' ')]) : 0);

    private static double Sd(IReadOnlyList<double> values)
    {
        double mean = values.Average();
        return Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / (values.Count - 1));
    }

    private sealed class Progress(Action<BassProgress> report) : IProgress<BassProgress>
    {
        public void Report(BassProgress value) => report(value);
    }
}
