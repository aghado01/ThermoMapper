using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Evidence;
using Maths.Samplers.Smc;
using Repo.TestHarness;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// KHAOS's SMC trap (<c>smc-adaptive-mutation-report.md</c>, "Open"), measured for the user's choice of remedy, behind
/// <see cref="BassQualification"/>, on that measurement's Ishigami fixture (<see cref="BassSmcMutationQualificationTests.Khaos"/>:
/// the raw ridge at τ held 10⁻⁵, degree 8, order 3, n 120, σ² ~ IG(3, 0.5)) and its seeds, so the engine arms' first replicates
/// are that measurement's, ln Ẑ for ln Ẑ.
/// <list type="bullet">
/// <item><b>The engine's schedule</b> (conditional ESS 0.99) at 5 and 20 sweeps per stage, each stage read after its mutation: b, the
/// population's weighted mean M, its distinct structures, its weighted mean collapsed log posterior at b = 1 (log marginal and
/// structure prior at the particle's τ) and the stage's moved fraction. Then two check stages at b = 1, which leave ln Ẑ alone
/// and show whether the population still drifts there.</item>
/// <item><b>A schedule that cannot jump:</b> the same at 5 sweeps, with each b′ at most r·b, r = 1.25 and 1.1.</item>
/// <item><b>KHAOS's g-prior</b> (ζ = 1, g₀² ~ IG(0.001, 1,000) drawn) under the engine's schedule.</item>
/// <item><b>The tempering ladder</b> with its reference rung, self-sizing, four chains: the cold chains' M and σ², the stepping-stone
/// and path-sampling ln Ẑ, the pair rejections and the evidence's saturation flag.</item>
/// </list>
/// A replicate is found when its weighted σ² is within twice the long chains' (0.2245 on the raw ridge; the same cut serves the g-prior,
/// whose chains read about 0.20). Settings:
/// <c>KHAOS_TRAP_PARTS</c> (engine,capped,gprior,ladder), <c>KHAOS_TRAP_REPLICATES</c> (32 for the engine at 5 sweeps, half that
/// elsewhere), <c>KHAOS_TRAP_THREADS</c> (4), and the ladder's <c>KHAOS_TRAP_LADDER_BURN</c> and <c>KHAOS_TRAP_LADDER_SAMPLES</c> (5,000 each). Tables: <c>khaos-trap-replicates.csv</c> and <c>khaos-trap-stages.csv</c> in the run's
/// artifact directory.
/// </summary>
public sealed class KhaosSmcTrapQualificationTests(ITestOutputHelper output)
{
    private const int Particles = 500, CheckStages = 2, MasterSeed = 20261108;   // the adaptive-mutation measurement's
    private const double TargetCess = 0.99, FoundBelow = 2 * 0.2245;

    private sealed record Arm(string Name, bool GPrior, int Rate, double? Cap, int Replicates);

    private sealed record Stage(double B, double MeanM, int Distinct, double MeanLogPi, double SdLogPi, double Ess, double Moved);

    private sealed record Replicate(Arm Arm, int Index, double LogZ, double LogZError, double M, double Sigma2, double Seconds, Stage[] Stages, Stage[] Checks)
    {
        public bool Found => Sigma2 < FoundBelow;

        // The flags a run could raise on itself, each read on the replicate.
        public double LastFrom => Stages.Length > 1 ? Stages[^2].B : 0.0;
        public bool Jumped => LastFrom < 0.5;
        public double LeastLate => Stages.TakeLast(3).Min(s => s.Moved);
        public double Drift => Checks[^1].MeanLogPi - Stages[^1].MeanLogPi;
        public double DriftError => Stages[^1].SdLogPi / Math.Sqrt(Stages[^1].Ess);
        public bool Drifted => Drift > Math.Max(1.0, 3.0 * DriftError);
        public int StagesAbove(double b) => Stages.Count(s => s.B >= b);
    }

    private static int Setting(string name, int fallback) => int.TryParse(Environment.GetEnvironmentVariable(name), out int v) ? v : fallback;

    private static BassSmcMutationQualificationTests.Fixture GPrior(BassSmcMutationQualificationTests.Fixture f)
        => f with { Config = f.Config with { Coefficients = BassCoefficientPrior.ModifiedG(TauPrior.Gamma(0.001, 1000.0)) } };

    private static string Key(IEnumerable<BassTerm> terms)
        => string.Join("+", terms.Select(t => string.Concat(t.Legendre.Select(f => $"{f.Variable}^{f.Degree}"))).Order(StringComparer.Ordinal));

    // The population as the stage left it.
    private static Stage Read(BassModel model, SmcRun<BassParticle> run)
    {
        SmcResult<BassParticle> r = run.Result();
        double m = 0.0, lp = 0.0, lp2 = 0.0, w2 = 0.0;
        var keys = new HashSet<string>();
        for (int i = 0; i < r.Particles.Count; i++)
        {
            double w = r.Weights[i];
            BassState s = r.Particles[i].State;
            double l = BassMixing.LogPosterior(model, s);
            m += w * s.Terms.Count;
            lp += w * l;
            lp2 += w * l * l;
            w2 += w * w;
            keys.Add(Key(s.Terms));
        }
        return new Stage(run.Stage.Temperature.InverseTemperature, m, keys.Count, lp, Math.Sqrt(Math.Max(0.0, lp2 - lp * lp)), 1.0 / w2,
                         run.Reports[^1].MovedFraction);
    }

    // One replicate, as BassSmcFit runs it (BassSmcMutationQualificationTests.Replicate), its schedule capped when the arm says.
    private static Replicate Run(Arm arm, BassModel model, int index, int seed)
    {
        int[] streams = SeedTree.Derive(seed, 2);   // the run's, the composition's
        var target = new BassParticleTarget(model.Target);
        var options = new SmcOptions
        {
            Particles = Particles, Resampling = ResamplingScheme.Multinomial, ResampleBelow = 0.5, MutationSteps = arm.Rate,
            MaxStages = 5000, MaxDegreeOfParallelism = 1,
        };
        SmcStage<BassParticle> At(double b) => new(target, new Temperature(b, TemperingPath.Likelihood));
        long t0 = BassCost.Tick();
        SmcRun<BassParticle> run = SmcRun<BassParticle>.Start(At(0.0), new BassParticles(model), new BassSmcKernel(model), options, streams[0]);
        var stages = new List<Stage>();
        while (run.Stage.Temperature.InverseTemperature < 1.0)
        {
            double b = run.Stage.Temperature.InverseTemperature, next = SmcSampler.NextInverseTemperature(run, TargetCess);
            if (arm.Cap is double r && b > 0.0) next = Math.Min(next, Math.Min(1.0, r * b));
            run.Advance(At(next));
            stages.Add(Read(model, run));
        }
        double seconds = (BassCost.Tick() - t0) / BassCost.TicksPerSecond;
        SmcResult<BassParticle> result = run.Result();

        int[] composition = SeedTree.Derive(streams[1], Particles);
        var terms = new BassTermTable();
        double m = 0.0, sigma2 = 0.0;
        for (int i = 0; i < Particles; i++)
        {
            BassParticle p = result.Particles[i];
            BassDraw d = p.Moves.Compose(p.State, new Xoshiro256PlusPlus(composition[i]), terms);
            m += result.Weights[i] * d.Structure.TermIds.Length;
            sigma2 += result.Weights[i] * d.Sigma2;
        }
        var checks = new List<Stage>();
        for (int c = 0; c < CheckStages; c++)
        {
            run.Advance(At(1.0));
            checks.Add(Read(model, run));
        }
        return new Replicate(arm, index, result.LogNormalizerRatio, result.LogNormalizerError, m, sigma2, seconds, [.. stages], [.. checks]);
    }

    [Fact]
    public void Trap_ItsMechanism_ItsFlags_AndTheLadder()
    {
        if (!BassQualification.Enabled) return;   // heavy
        int replicates = Setting("KHAOS_TRAP_REPLICATES", 32), threads = Setting("KHAOS_TRAP_THREADS", 4), half = Math.Max(1, replicates / 2);
        string[] parts = (Environment.GetEnvironmentVariable("KHAOS_TRAP_PARTS") ?? "engine,capped,gprior,ladder").Split(',', StringSplitOptions.RemoveEmptyEntries);
        ArtifactRun artifacts = HarnessArtifacts.Create("test-runs", "KhaosSmcTrap", "Maths.Regression.Tests", new Dictionary<string, object?>
        {
            ["Particles"] = Particles, ["Replicates"] = replicates, ["Threads"] = threads, ["Parts"] = string.Join(",", parts),
            ["Clock"] = BassCost.Cycles ? "cycles" : "wall",
        });
        output.WriteLine($"artifacts: {artifacts.RunDirectory}");
        BassSmcMutationQualificationTests.Fixture raw = BassSmcMutationQualificationTests.Khaos(), g = GPrior(raw);
        BassModel rawModel = BassModel.Create(raw.Config, raw.Data), gModel = BassModel.Create(g.Config, g.Data);
        Assert.True(rawModel.Target.HasCompleteValues && gModel.Target.HasCompleteValues);

        var arms = new List<Arm>();
        if (parts.Contains("engine")) arms.AddRange([new("engine, 5 sweeps", false, 5, null, replicates), new("engine, 20 sweeps", false, 20, null, half)]);
        if (parts.Contains("capped")) arms.AddRange([new("capped r 1.25, 5 sweeps", false, 5, 1.25, half), new("capped r 1.1, 5 sweeps", false, 5, 1.1, half)]);
        if (parts.Contains("gprior")) arms.Add(new("g-prior, engine, 5 sweeps", true, 5, null, half));

        using var rowsCsv = new StreamWriter(Path.Combine(artifacts.RunDirectory, "khaos-trap-replicates.csv")) { NewLine = "\n" };
        using var stagesCsv = new StreamWriter(Path.Combine(artifacts.RunDirectory, "khaos-trap-stages.csv")) { NewLine = "\n" };
        rowsCsv.WriteLine("arm,replicate,found,log_z,log_z_error,m,sigma2,seconds,stages,stages_above_0.01,last_from,least_late,drift,drift_error");
        stagesCsv.WriteLine("arm,replicate,stage,check,b,mean_m,distinct,mean_log_pi,sd_log_pi,ess,moved");
        int[] seeds = SeedTree.Derive(MasterSeed, replicates);
        var jobs = arms.SelectMany(a => Enumerable.Range(0, a.Replicates).Select(r => (a, r))).ToArray();
        var done = new Replicate[jobs.Length];
        Parallel.For(0, jobs.Length, new ParallelOptions { MaxDegreeOfParallelism = threads },
                     k => done[k] = Run(jobs[k].a, jobs[k].a.GPrior ? gModel : rawModel, jobs[k].r, seeds[jobs[k].r]));

        foreach (Replicate x in done)
        {
            rowsCsv.WriteLine(string.Join(",", Csv(x.Arm.Name), x.Index, x.Found ? 1 : 0, F(x.LogZ), F(x.LogZError), F(x.M), F(x.Sigma2), F(x.Seconds),
                                          x.Stages.Length, x.StagesAbove(0.01), F(x.LastFrom), F(x.LeastLate), F(x.Drift), F(x.DriftError)));
            foreach ((Stage s, int i, bool check) in x.Stages.Select((s, i) => (s, i, false)).Concat(x.Checks.Select((s, i) => (s, x.Stages.Length + i, true))))
                stagesCsv.WriteLine(string.Join(",", Csv(x.Arm.Name), x.Index, i, check ? 1 : 0, F(s.B), F(s.MeanM), s.Distinct, F(s.MeanLogPi), F(s.SdLogPi), F(s.Ess), F(s.Moved)));
        }

        output.WriteLine("arm | found | ln Ẑ found: mean, sd, range | ln Ẑ trapped: mean, sd | stages found, trapped | s found, trapped | " +
                         "jumped (found, trapped) | late least < 0.1 (found, trapped) | drifted (found, trapped) | drift found median, trapped median");
        foreach (Arm arm in arms)
        {
            Replicate[] all = [.. done.Where(x => x.Arm == arm)], f = [.. all.Where(x => x.Found)], t = [.. all.Where(x => !x.Found)];
            string Z(Replicate[] xs) => xs.Length == 0 ? "-" : $"{xs.Average(x => x.LogZ):F2}, {Sd(xs.Select(x => x.LogZ)):F2}" + (xs == f ? $", [{xs.Min(x => x.LogZ):F2}, {xs.Max(x => x.LogZ):F2}]" : "");
            string Mean(Replicate[] xs, Func<Replicate, double> v) => xs.Length == 0 ? "-" : xs.Average(v).ToString("F1", CultureInfo.InvariantCulture);
            string Count(Func<Replicate, bool> flag) => $"{f.Count(flag)}/{f.Length}, {t.Count(flag)}/{t.Length}";
            string Median(Replicate[] xs) => xs.Length == 0 ? "-" : BassMixing.Quantile([.. xs.Select(x => x.Drift)], 0.5).ToString("F2", CultureInfo.InvariantCulture);
            output.WriteLine($"{arm.Name} | {f.Length}/{all.Length} | {Z(f)} | {Z(t)} | {Mean(f, x => x.Stages.Length)}, {Mean(t, x => x.Stages.Length)} | " +
                             $"{Mean(f, x => x.Seconds)}, {Mean(t, x => x.Seconds)} | {Count(x => x.Jumped)} | {Count(x => x.LeastLate < 0.1)} | {Count(x => x.Drifted)} | " +
                             $"{Median(f)}, {Median(t)}");
        }
        foreach (Replicate x in done.Where(x => x.Index < 2 || x.Found))
            output.WriteLine($"  {x.Arm.Name} r{x.Index}: ln Ẑ {x.LogZ:F2} ± {x.LogZError:F2}, M {x.M:F2}, σ² {x.Sigma2:F3}, {x.Stages.Length} stages " +
                             $"({x.StagesAbove(0.01)} at b ≥ 0.01), last from b {x.LastFrom:G3}, drift {x.Drift:F2} ± {x.DriftError:F2}; mean M at the last stages " +
                             string.Join(" ", x.Stages.TakeLast(4).Concat(x.Checks).Select(s => $"{s.MeanM:F2}@{s.B:G3}({s.Distinct})")));

        if (parts.Contains("ladder")) Ladder(raw, threads);
    }

    // ---- the flags as built (the main session's decision, 2026-11-08) ------------------------------------------------------------

    /// <summary>
    /// SMC's flags as BASS's fit raises them (<see cref="BassSmcReplicate.JumpedOverACollapse"/>, <see cref="BassSmcReadout.Overdispersed"/>),
    /// on production fits of N 500 at 5 sweeps a stage: the raw-ridge trap fixture with 32 replicates and KHAOS's g-prior with 16, each
    /// replicate read trapped or found by its own weighted σ²; and, for false alarms, GBASS slice 5's t(5) and skewed GH evidence fixtures
    /// at the adaptive-mutation measurement's size, 32 replicates. The disagreement is read over each whole fit and over its groups of
    /// four, a default fit's replicates. Gates: no replicate that reached the posterior jumps; every fit and group on the trap fixtures is
    /// flagged one way or the other; the healthy fits raise nothing.
    /// </summary>
    [Fact]
    public void Flags_OnTheTrap_AndNoneOnHealthyFits()
    {
        if (!BassQualification.Enabled) return;   // heavy
        var options = new BassRunOptions { MaxDegreeOfParallelism = Setting("KHAOS_TRAP_THREADS", 4) };
        (string Name, BassSmcMutationQualificationTests.Fixture F, int Replicates, bool Healthy)[] fits =
        [
            ("KHAOS raw ridge", BassSmcMutationQualificationTests.Khaos(), 32, false),
            ("KHAOS g-prior", GPrior(BassSmcMutationQualificationTests.Khaos()), 16, false),
            ("GBASS t(5)", BassSmcMutationQualificationTests.Gbass("t"), 32, true),
            ("GBASS skewed GH", BassSmcMutationQualificationTests.Gbass("gh"), 32, true),
        ];
        output.WriteLine("fit | reached | jumped (reached, trapped) | between-to-within, whole | groups of four flagged | flags");
        output.WriteLine("--- | --- | --- | --- | --- | ---");
        foreach ((string name, BassSmcMutationQualificationTests.Fixture f, int replicates, bool healthy) in fits)
        {
            BassConfig config = f.Config with
            {
                Sampling = f.Config.Sampling with
                {
                    MasterSeed = 20261109, ReportEvidence = true, Smc = new SmcSpec { Particles = Particles, Replicates = replicates, MutationSteps = 5 },
                },
            };
            BassSmcReadout readout = Bass.Fit(config, f.Data, options).Diagnostics.Smc!;
            Assert.True(readout.ReachedPosterior);
            double[] sigma2 = readout.Quantity("sigma2").ReplicateMeans;
            bool[] reached = [.. sigma2.Select(v => healthy || v < FoundBelow)];
            int jumpedReached = readout.Replicates.Where((r, c) => reached[c] && r.JumpedOverACollapse == true).Count();
            int jumpedTrapped = readout.Replicates.Where((r, c) => !reached[c] && r.JumpedOverACollapse == true).Count();
            int groups = replicates / 4, flaggedGroups = 0;
            for (int g = 0; g < groups; g++)
            {
                var group = new BassSmcReadout([.. readout.Replicates.Skip(4 * g).Take(4)], []);
                if (group.JumpedReplicates > 0 || group.Overdispersed == true) flaggedGroups++;
            }
            output.WriteLine($"{name} | {reached.Count(r => r)}/{replicates} | {jumpedReached}, {jumpedTrapped}/{reached.Count(r => !r)} | " +
                             $"{readout.NormalizerBetweenToWithin:F2} | {flaggedGroups}/{groups} | {readout.Flags.Count}");
            foreach (BassSmcReplicate r in readout.Replicates.Take(4))
                output.WriteLine($"  last step from b {r.LastStepFrom:G3} over {r.StructuresBeforeLastStep} of {r.Particles}, ln Ẑ {r.LogNormalizer:F2} ± {r.LogNormalizerError:F2}");
            Assert.Equal(0, jumpedReached);
            if (healthy)
            {
                Assert.Empty(readout.Flags);
                Assert.Equal(0, flaggedGroups);
            }
            else
            {
                Assert.True(readout.JumpedReplicates > 0 || readout.Overdispersed == true, $"{name}: unflagged");
                Assert.Equal(groups, flaggedGroups);
            }
        }
    }

    // The tempering ladder on the raw-ridge fixture: the cold chains against the long chains, and its evidence.
    private void Ladder(BassSmcMutationQualificationTests.Fixture f, int threads)
    {
        int burn = Setting("KHAOS_TRAP_LADDER_BURN", 5000), samples = Setting("KHAOS_TRAP_LADDER_SAMPLES", 5000);
        BassConfig config = f.Config with
        {
            Sampling = f.Config.Sampling with
            {
                Chains = 4, MasterSeed = 7, BurnIn = burn, Samples = samples, ReportEvidence = true,
                Tempering = new TemperingSpec { ReferenceRung = true },
            },
        };
        long t0 = BassCost.Tick();
        BassFit fit = Bass.Fit(config, f.Data, new BassRunOptions { MaxDegreeOfParallelism = threads });
        double seconds = (BassCost.Tick() - t0) / BassCost.TicksPerSecond;
        for (int c = 0; c < fit.Draws.Count; c++)
        {
            BassDraw[] draws = [.. fit.Draws[c]];
            output.WriteLine($"ladder chain {c}: cold M {draws.Average(d => (double)d.Structure.TermIds.Length):F2}, σ² {draws.Average(d => d.Sigma2):F4}");
        }
        BassTemperingReadout ladder = fit.Diagnostics.Tempering!;
        for (int c = 0; c < ladder.Levels.Length; c++)
        {
            int worst = Array.IndexOf(ladder.PairRejection[c], ladder.PairRejection[c].Max());
            output.WriteLine($"ladder chain {c}: {ladder.Levels[c]} rungs, b {string.Join(" ", ladder.Ladders[c].Select(b => b.ToString("G3", CultureInfo.InvariantCulture)))}; " +
                             $"worst pair {worst} ({ladder.Ladders[c][worst]:G3}, {ladder.Ladders[c][worst + 1]:G3}) rejects {ladder.PairRejection[c][worst]:F3}");
        }
        BassEvidence e = fit.Evidence!;
        output.WriteLine($"ladder: {seconds:F0} s; stepping stone {e.SteppingStone?.LogEvidence:F2} ± {e.SteppingStone?.MonteCarloError:F2} " +
                         $"(between/within {e.SteppingStone?.BetweenToWithin:F2}); path sampling {e.PathSampling?.LogEvidence:F2} ± {e.PathSampling?.MonteCarloError:F2}; " +
                         $"max pair rejection {e.MaxPairRejection:F3}, saturated {e.Saturated}, at max levels {e.AtMaxLevels}; barrier {ladder.Barrier:F2}, " +
                         $"round trips per scan {ladder.RoundTripsPerScan:G3}");
        foreach (LadderEstimate chain in e.SteppingStone?.Chains ?? [])
            output.WriteLine($"  stepping stone, a chain: {chain.LogRatio:F2} ± {chain.MonteCarloError:F2}");
    }

    private static double Sd(IEnumerable<double> values)
    {
        double[] v = [.. values];
        if (v.Length < 2) return double.NaN;
        double mean = v.Average();
        return Math.Sqrt(v.Sum(x => (x - mean) * (x - mean)) / (v.Length - 1));
    }

    private static string F(double v) => v.ToString("R", CultureInfo.InvariantCulture);

    private static string Csv(string s) => s.Contains(',') ? $"\"{s}\"" : s;
}
