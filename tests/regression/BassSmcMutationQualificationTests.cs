using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Smc;
using Repo.TestHarness;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// SMC's adaptive count of moves per stage (<see cref="AdaptiveMutation"/>) against the fixed count, at equal cost, behind
/// <see cref="BassQualification"/>: the measurement that sets <see cref="SmcSpec"/>'s default. On GBASS slice 5's evidence fixtures
/// (t(5) and the skewed generalized hyperbolic, n 12) and on KHAOS's Ishigami fit under the raw ridge (degree 8, order 3, n 120,
/// σ² made proper so the normalizer is an evidence), each arm runs R replicates of BASS's SMC, each replicate the engine's own run
/// (<see cref="BassParticles"/>, <see cref="BassSmcKernel"/>, the schedule's conditional ESS 0.99), timed alone on its thread.
/// <list type="bullet">
/// <item><b>The arms:</b> the fixed count at 2, 5, 10 and 20 sweeps per stage, the cost frontier; the adaptive count on a budget of
/// 5 per stage at c = 0.1 and 0.01, and of 20 at c = 0.01. Each adaptive arm is read against the fixed arm at its budget.</item>
/// <item><b>The references:</b> ln p(y) by plain Monte Carlo over the prior's exact draws (slice 5's); on KHAOS, which has none, the
/// pooled ln of the mean Ẑ over every replicate of every arm, which leans on the replicates that found the posterior since Ẑ is
/// unbiased in Z (a floor where none did). The posterior means of M, σ² and the skew by eight long MCMC chains.</item>
/// <item><b>The readings:</b> each arm's ln Ẑ spread and error against the reference, the posterior means' errors, the replicates
/// that reached the posterior (weighted σ² within twice the chains'), and its cost in sweeps and in seconds on the thread's own
/// clock (<c>BASS_COST_CLOCK</c> = cycles; indicative on a shared machine). Equal cost is read through the work-normalized error,
/// the mean squared error times the cost, which a change of N leaves alone; its ratio to the fixed arm at the same budget takes a
/// 95% interval from a bootstrap over the replicates, paired since every arm's replicate r starts from the same seed.</item>
/// <item><b>The check that the harness is the engine:</b> a replicate here gives the ln Ẑ <see cref="Bass.Fit"/>'s replicate gives, bit
/// for bit, on the fixed count.</item>
/// </list>
/// Settings: <c>SMC_MUTATION_PARTICLES</c> (500), <c>SMC_MUTATION_REPLICATES</c> (32), <c>SMC_MUTATION_THREADS</c> (4),
/// <c>SMC_MUTATION_FIXTURES</c> (t,gh,khaos), and the references' sizes, <c>SMC_MUTATION_PLAIN_DRAWS</c> (2,000,000) and
/// <c>SMC_MUTATION_CHAIN_SWEEPS</c> (each chain's burn-in and samples, 20,000; KHAOS 30,000). Tables: <c>smc-mutation-replicates.csv</c>, <c>smc-mutation-stages.csv</c> and
/// <c>smc-mutation-results.csv</c> in the run's artifact directory.
/// </summary>
public sealed class BassSmcMutationQualificationTests(ITestOutputHelper output)
{
    private const double TargetCess = 0.99;   // SmcSpec's default
    private const int MasterSeed = 20261108;

    private sealed record Arm(string Name, int Rate, AdaptiveMutation? Adaptive);

    private static readonly Arm[] Arms =
    [
        new("fixed 2", 2, null), new("fixed 5", 5, null), new("fixed 10", 10, null), new("fixed 20", 20, null),
        new("adaptive 5, c 0.1", 5, new AdaptiveMutation { StayProbability = 0.1 }),
        new("adaptive 5, c 0.01", 5, new AdaptiveMutation { StayProbability = 0.01 }),
        new("adaptive 20, c 0.01", 20, new AdaptiveMutation { StayProbability = 0.01 }),
    ];

    internal sealed record Fixture(string Name, BassConfig Config, BassData Data);

    private sealed record Reference(double LogZ, double LogZError, double M, double MError, double Sigma2, double Sigma2Error, double Skew, double SkewError);

    private sealed record Row(string Fixture, string Arm, int Replicate, int Seed, int Stages, long Moves, double Seconds, double LogZ, double LogZError,
                              double M, double Sigma2, double Skew, int LeastMovedStage, double LeastMoved, int Families, SmcStageReadout Readout);

    private static int Setting(string name, int fallback) => int.TryParse(Environment.GetEnvironmentVariable(name), out int v) ? v : fallback;

    // ---- the fixtures ---------------------------------------------------------------------------------------------------------------

    // GBASS slice 5's evidence fixture (BassNoiseReachQualificationTests): twelve cases of sin(7x) + 0.4·N(0, 1), M_max 2, σ² ~ IG(8, 6).
    internal static Fixture Gbass(string form)
    {
        double[] x = [.. Enumerable.Range(0, 12).Select(i => i / 11.0)];
        var noiseRng = new Xoshiro256PlusPlus(20261109);
        BassData data = new([new ContinuousPredictor("x", x)], new ScalarResponse([.. x.Select(v => Math.Sin(7 * v) + 0.4 * Normal.SampleStandard(noiseRng))]));
        BassNoise noise = form == "t"
            ? BassNoise.StudentT(5, reference: BassNoiseReference.None)
            : BassNoise.GeneralizedHyperbolic(new GigParameters(0.7, 1.2, 0.8), new BassNoisePrior(0.3, 0.5), reference: BassNoiseReference.None);
        BassConfig config = new()
        {
            Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 2, Support = new SupportRule(SupportKind.Count, 2) },
            Coefficients = BassCoefficientPrior.G(TauPrior.UnitInformationGamma(5)),
            Intercept = GPriorCentre.Zero,
            Sigma2 = Sigma2Prior.InverseGamma(8.0, 6.0),
            Noise = noise,
            Sampling = new BassSampling { Chains = 1, MasterSeed = 1, BurnIn = 0, Samples = 1 },
        };
        return new Fixture(form == "t" ? "t(5)" : "skewed GH", config, data);
    }

    // KHAOS slice 0's Ishigami readout fixture on its raw inputs with [0, 1] declared (KhaosRichStartTests), under the raw ridge at
    // khaos's τ² = 10⁵, h₁ = 4, h₂ = 40/n and the zero centre; σ² ~ IG(3, 0.5) in place of Jeffreys, so the run starts from a proper
    // prior in the integrated form and its normalizer is the evidence.
    internal static Fixture Khaos()
    {
        const int inputs = 4, cases = 120;
        var rng = new Xoshiro256PlusPlus(20261107);
        var x = new double[inputs][];
        for (int j = 0; j < inputs; j++) x[j] = new double[cases];
        var y = new double[cases];
        for (int i = 0; i < cases; i++)
        {
            for (int j = 0; j < inputs; j++) x[j][i] = rng.NextDouble();
            double z1 = Math.PI * (2 * x[0][i] - 1), z2 = Math.PI * (2 * x[1][i] - 1), z3 = Math.PI * (2 * x[2][i] - 1);
            y[i] = Math.Sin(z1) + 7 * Math.Pow(Math.Sin(z2), 2) + 0.1 * Math.Pow(z3, 4) * Math.Sin(z1) + 0.5 * Normal.SampleStandard(rng);
        }
        var domain = new InputDomain(0.0, 1.0);
        BassData data = new([.. Enumerable.Range(0, inputs).Select(j => new ContinuousPredictor($"x{j}", x[j], domain))], new ScalarResponse(y));
        BassConfig config = new()
        {
            Structure = new BassStructureSpec { MaxInteraction = 3, Legendre = new LegendreBasis { MaxDegree = 8 } },
            Count = new BassCountPrior(4, 40.0 / cases),
            Coefficients = BassCoefficientPrior.RawRidge(TauPrior.Held(1e-5)),
            Intercept = GPriorCentre.Zero,
            Sigma2 = Sigma2Prior.InverseGamma(3.0, 0.5),
            Sampling = new BassSampling { Chains = 1, MasterSeed = 1 },
        };
        return new Fixture("KHAOS Ishigami", config, data);
    }

    private static BassConfig SmcConfig(Fixture f, SmcSpec spec) => f.Config with { Sampling = f.Config.Sampling with { MasterSeed = MasterSeed, Smc = spec, ReportEvidence = true } };

    // ---- one replicate, as BassSmcFit runs it ----------------------------------------------------------------------------------------

    private static Row Replicate(Fixture f, BassModel model, Arm arm, int particles, int replicate, int seed)
    {
        int[] streams = SeedTree.Derive(seed, 2);   // the run's, the composition's
        var target = new BassParticleTarget(model.Target);
        var options = new SmcOptions
        {
            Particles = particles, Resampling = ResamplingScheme.Multinomial, ResampleBelow = 0.5, MutationSteps = arm.Rate,
            AdaptiveMutation = arm.Adaptive, MaxStages = 1000, MaxDegreeOfParallelism = 1,
        };
        long t0 = BassCost.Tick();
        SmcRun<BassParticle> run = SmcRun<BassParticle>.Start(new SmcStage<BassParticle>(target, new Temperature(0.0, TemperingPath.Likelihood)),
                                                             new BassParticles(model), new BassSmcKernel(model), options, streams[0]);
        SmcResult<BassParticle> result = SmcSampler.TemperToPosterior(run, TargetCess);
        double seconds = (BassCost.Tick() - t0) / BassCost.TicksPerSecond;
        if (!result.ReachedPosterior) throw new InvalidOperationException($"{f.Name}, {arm.Name}: replicate {replicate} stopped short of b = 1.");

        int[] composition = SeedTree.Derive(streams[1], particles);
        var terms = new BassTermTable();
        double m = 0.0, sigma2 = 0.0, skew = 0.0;
        for (int i = 0; i < particles; i++)
        {
            BassParticle p = result.Particles[i];
            BassDraw d = p.Moves.Compose(p.State, new Xoshiro256PlusPlus(composition[i]), terms);
            double w = result.Weights[i];
            m += w * d.Structure.TermIds.Length;
            sigma2 += w * d.Sigma2;
            skew += w * (d.Latents?.Skew ?? 0.0);
        }
        SmcStageReadout readout = SmcStageReadout.Of(result.Stages);
        (int leastStage, double least) = readout.LeastMoved;
        return new Row(f.Name, arm.Name, replicate, seed, result.Stages.Count - 1, readout.Moves, seconds, result.LogNormalizerRatio,
                       result.LogNormalizerError, m, sigma2, skew, leastStage, least, result.Stages[^1].Families, readout);
    }

    // ---- the references ------------------------------------------------------------------------------------------------------------

    // ln p(y) = ln E[N(y | f + σβv, σ²c·v)] over exact draws of the whole prior on F, the latents with them (slice 5's), in blocks
    // on their own streams.
    private static (double LogZ, double Error) PlainMonteCarlo(Fixture f, int draws, int threads)
    {
        BassModel model = BassModel.Create(f.Config, f.Data);
        double c = f.Config.Noise!.VarianceScale;
        double[] y = f.Data.Response.Values.ToArray();
        const int blocks = 16;
        int[] seeds = SeedTree.Derive(20261110, blocks);
        var logs = new double[draws];
        int per = draws / blocks;
        Parallel.For(0, blocks, new ParallelOptions { MaxDegreeOfParallelism = threads }, b =>
        {
            var sampler = new BassPriorSampler(model);
            var rng = new Xoshiro256PlusPlus(seeds[b]);
            for (int k = b * per; k < (b + 1) * per; k++)
            {
                BassState state = sampler.DrawPrior(rng).State;
                var moves = new BassStructureMoves(model, state, new BassMovesSnapshot(sampler.LastDesign(), new Tally<BassRefusal>(), 0,
                                                                                       new BassSupportTally(model.Prior), new BassMoveTally()));
                BassSweepFit fit = moves.DrawSweepFit(state, rng, model.Target, power: 0.0);
                ObservationLatents latents = state.Fit.Working!.Latents!;
                double s2 = latents.Sigma2 ?? fit.Sigma2, sigma = Math.Sqrt(s2), log = 0.0;
                for (int i = 0; i < y.Length; i++)
                {
                    double vi = latents.Values[i], v = s2 * c * vi, r = y[i] - fit.Fitted[i] - sigma * latents.Skew * vi;
                    log += -0.5 * Math.Log(2.0 * Math.PI * v) - 0.5 * r * r / v;
                }
                logs[k] = log;
            }
        });
        int n = per * blocks;
        double top = logs.Take(n).Max(), mean = logs.Take(n).Average(l => Math.Exp(l - top)), second = logs.Take(n).Average(l => Math.Exp(2.0 * (l - top)));
        return (top + Math.Log(mean), Math.Sqrt((second - mean * mean) / n) / mean);
    }

    // The posterior means of M, σ² and the skew over eight long chains, each with its error from the chains' spread.
    private static (double M, double MError, double Sigma2, double Sigma2Error, double Skew, double SkewError) LongChains(Fixture f, int burnIn, int samples, int threads)
    {
        const int chains = 8;
        BassConfig config = f.Config with { Sampling = f.Config.Sampling with { Chains = chains, MasterSeed = 11, BurnIn = burnIn, Samples = samples } };
        BassFit fit = Bass.Fit(config, f.Data, new BassRunOptions { MaxDegreeOfParallelism = threads });
        double[][] means = [.. fit.Draws.Select(chain =>
        {
            BassDraw[] draws = [.. chain];
            return new[] { draws.Average(d => (double)d.Structure.TermIds.Length), draws.Average(d => d.Sigma2), draws.Average(d => d.Latents?.Skew ?? 0.0) };
        })];
        (double, double) Pooled(int q)
        {
            double[] v = [.. means.Select(m => m[q])];
            double mean = v.Average();
            return (mean, Math.Sqrt(v.Sum(x => (x - mean) * (x - mean)) / (v.Length - 1) / v.Length));
        }
        (double m, double me) = Pooled(0);
        (double s, double se) = Pooled(1);
        (double k, double ke) = Pooled(2);
        return (m, me, s, se, k, ke);
    }

    // ---- the measurement ------------------------------------------------------------------------------------------------------------

    [Fact]
    public void AdaptiveAgainstFixed_AtEqualCost()
    {
        if (!BassQualification.Enabled) return;   // heavy
        int particles = Setting("SMC_MUTATION_PARTICLES", 500), replicates = Setting("SMC_MUTATION_REPLICATES", 32), threads = Setting("SMC_MUTATION_THREADS", 4);
        string[] chosen = (Environment.GetEnvironmentVariable("SMC_MUTATION_FIXTURES") ?? "t,gh,khaos").Split(',', StringSplitOptions.RemoveEmptyEntries);
        ArtifactRun artifacts = HarnessArtifacts.Create("test-runs", "BassSmcMutation", "Maths.Regression.Tests", new Dictionary<string, object?>
        {
            ["Particles"] = particles, ["Replicates"] = replicates, ["Threads"] = threads, ["Fixtures"] = string.Join(",", chosen),
            ["Clock"] = BassCost.Cycles ? "cycles" : "wall", ["Processors"] = Environment.ProcessorCount,
        });
        using var rowsCsv = new StreamWriter(Path.Combine(artifacts.RunDirectory, "smc-mutation-replicates.csv")) { NewLine = "\n" };
        using var stagesCsv = new StreamWriter(Path.Combine(artifacts.RunDirectory, "smc-mutation-stages.csv")) { NewLine = "\n" };
        using var resultsCsv = new StreamWriter(Path.Combine(artifacts.RunDirectory, "smc-mutation-results.csv")) { NewLine = "\n" };
        rowsCsv.WriteLine("fixture,arm,replicate,seed,stages,sweeps,seconds,log_z,log_z_error,m,sigma2,skew,least_moved_stage,least_moved,families");
        stagesCsv.WriteLine("fixture,arm,replicate,stage,b,steps,acceptance,moved,ess_fraction,resampled,families");
        resultsCsv.WriteLine("fixture,arm,baseline,replicates,particles,stages,sweeps,seconds,reached,log_z_reference,log_z_mean,log_z_sd,log_z_bias," +
                             "log_z_rmse,m_rmse,sigma2_rmse,skew_rmse,work_ratio_log_z,work_ratio_low,work_ratio_high,work_ratio_m,work_ratio_sigma2," +
                             "work_ratio_skew,least_moved_median");
        output.WriteLine($"artifacts: {artifacts.RunDirectory}; N {particles}, R {replicates}, {threads} threads, clock {(BassCost.Cycles ? "cycles" : "wall")}");

        foreach (string name in chosen)
        {
            Fixture f = name == "khaos" ? Khaos() : Gbass(name);
            bool khaos = name == "khaos";
            BassModel model = BassModel.Create(SmcConfig(f, new SmcSpec { Particles = particles, Replicates = replicates }), f.Data);
            Assert.True(model.Target.HasCompleteValues, $"{f.Name}: the target's values are not complete, so its normalizer is no evidence");

            // The harness is the engine: a replicate here gives Bass.Fit's replicate's ln Ẑ bit for bit.
            const int checkParticles = 60, checkSteps = 2;
            BassFit engine = Bass.Fit(SmcConfig(f, new SmcSpec { Particles = checkParticles, Replicates = 2, MutationSteps = checkSteps }), f.Data,
                                      new BassRunOptions { MaxDegreeOfParallelism = 2 });
            BassModel checkModel = BassModel.Create(SmcConfig(f, new SmcSpec { Particles = checkParticles, Replicates = 2, MutationSteps = checkSteps }), f.Data);
            int[] checkSeeds = SeedTree.Derive(MasterSeed, 2);
            for (int c = 0; c < 2; c++)
                Assert.Equal(engine.Diagnostics.Smc!.Replicates[c].LogNormalizer,
                             Replicate(f, checkModel, new Arm("check", checkSteps, null), checkParticles, c, checkSeeds[c]).LogZ);

            // The references.
            int plain = Setting("SMC_MUTATION_PLAIN_DRAWS", 2_000_000), sweeps = Setting("SMC_MUTATION_CHAIN_SWEEPS", khaos ? 30_000 : 20_000);
            (double logZ, double logZError) = khaos ? (double.NaN, double.NaN) : PlainMonteCarlo(f, plain, threads);
            var chains = LongChains(f, sweeps, sweeps, threads);
            var reference = new Reference(logZ, logZError, chains.M, chains.MError, chains.Sigma2, chains.Sigma2Error, chains.Skew, chains.SkewError);
            output.WriteLine($"{f.Name}: reference ln p(y) {reference.LogZ:F4} ± {reference.LogZError:F4}; M {reference.M:F4} ± {reference.MError:F4}, " +
                             $"σ² {reference.Sigma2:F4} ± {reference.Sigma2Error:F4}, skew {reference.Skew:F4} ± {reference.SkewError:F4}");

            // A warm-up replicate per arm, discarded, so no timed replicate pays a first compilation.
            foreach (Arm arm in Arms) Replicate(f, model, arm, Math.Min(particles, 50), -1, 1);

            // Every arm's replicates on the same seeds, interleaved so the machine's load falls on all of them alike.
            int[] seeds = SeedTree.Derive(MasterSeed, replicates);
            var rows = new Row[Arms.Length, replicates];
            Parallel.For(0, Arms.Length * replicates, new ParallelOptions { MaxDegreeOfParallelism = threads }, k =>
            {
                int r = k / Arms.Length, a = k % Arms.Length;
                rows[a, r] = Replicate(f, model, Arms[a], particles, r, seeds[r]);
            });
            for (int a = 0; a < Arms.Length; a++)
                for (int r = 0; r < replicates; r++)
                {
                    Row x = rows[a, r];
                    rowsCsv.WriteLine(string.Join(",", Csv(x.Fixture), Csv(x.Arm), x.Replicate, x.Seed, x.Stages, x.Moves, F(x.Seconds), F(x.LogZ), F(x.LogZError),
                                                  F(x.M), F(x.Sigma2), F(x.Skew), x.LeastMovedStage, F(x.LeastMoved), x.Families));
                    if (r < 2)
                        for (int s = 0; s < x.Readout.Steps.Length; s++)
                            stagesCsv.WriteLine(string.Join(",", Csv(x.Fixture), Csv(x.Arm), r, s, F(x.Readout.InverseTemperature[s]), x.Readout.Steps[s],
                                                            F(x.Readout.AcceptanceRate[s]), F(x.Readout.MovedFraction[s]), F(x.Readout.EssFraction[s]),
                                                            x.Readout.Resampled[s] ? 1 : 0, x.Readout.Families[s]));
                }
            Summarize(f, reference, rows, particles, khaos, resultsCsv);
            rowsCsv.Flush();
            stagesCsv.Flush();
            resultsCsv.Flush();
        }
    }

    // Each arm's readings and its work-normalized errors (MSE times mean seconds) against the fixed arm at its budget, with a paired
    // bootstrap's 95% interval on the ln Ẑ ratio. Without a plain Monte Carlo ln p(y), the reference is the pooled ln mean Ẑ.
    private void Summarize(Fixture f, Reference reference, Row[,] rows, int particles, bool noEvidence, StreamWriter csv)
    {
        int arms = rows.GetLength(0), r = rows.GetLength(1);
        Row[] Of(int a) => [.. Enumerable.Range(0, r).Select(k => rows[a, k])];
        double logZ = reference.LogZ;
        if (noEvidence)
        {
            double[] all = [.. rows.Cast<Row>().Select(x => x.LogZ)];
            double top = all.Max();
            logZ = top + Math.Log(all.Average(v => Math.Exp(v - top)));
        }
        var boot = new Xoshiro256PlusPlus(20261111);
        const int draws = 4000;
        int[][] picks = [.. Enumerable.Range(0, draws).Select(_ => Enumerable.Range(0, r).Select(_ => (int)(boot.NextDouble() * r)).ToArray())];
        output.WriteLine($"{f.Name}: ln p(y) reference {logZ:F4}{(noEvidence ? " (pooled ln mean Ẑ)" : "")}");
        output.WriteLine($"{f.Name}: arm | stages | sweeps | s | reached | ln Ẑ mean ± sd (bias) | rmse ln Ẑ, M, σ², skew | " +
                         "work ratio to the fixed arm at its budget: ln Ẑ [95%], M, σ², skew");
        for (int a = 0; a < arms; a++)
        {
            Row[] x = Of(a);
            int b = Array.FindIndex(Arms, v => v.Adaptive is null && v.Rate == Arms[a].Rate);
            Row[] y = Of(b);
            double seconds = x.Average(v => v.Seconds), baseSeconds = y.Average(v => v.Seconds);
            double Ratio(Func<Row, double> value, double centre, int[]? pick = null)
            {
                double denominator = Mse(y, value, centre, pick) * baseSeconds;
                return denominator > 0.0 ? Mse(x, value, centre, pick) * seconds / denominator : double.NaN;
            }
            double ratio = Ratio(v => v.LogZ, logZ), ratioM = Ratio(v => v.M, reference.M), ratioS = Ratio(v => v.Sigma2, reference.Sigma2);
            double ratioK = Ratio(v => v.Skew, reference.Skew);
            double[] boots = [.. picks.Select(pick => Ratio(v => v.LogZ, logZ, pick)).Where(v => !double.IsNaN(v)).Order()];
            double low = boots.Length > 0 ? boots[(int)(0.025 * (boots.Length - 1))] : double.NaN;
            double high = boots.Length > 0 ? boots[(int)(0.975 * (boots.Length - 1))] : double.NaN;
            double mean = x.Average(v => v.LogZ), sd = Math.Sqrt(x.Sum(v => (v.LogZ - mean) * (v.LogZ - mean)) / (r - 1));
            double rmse = Math.Sqrt(Mse(x, v => v.LogZ, logZ)), rmseM = Math.Sqrt(Mse(x, v => v.M, reference.M));
            double rmseS = Math.Sqrt(Mse(x, v => v.Sigma2, reference.Sigma2)), rmseK = Math.Sqrt(Mse(x, v => v.Skew, reference.Skew));
            int reached = x.Count(v => v.Sigma2 < 2.0 * reference.Sigma2);
            double[] least = [.. x.Select(v => v.LeastMoved).Order()];
            double stages = x.Average(v => v.Stages), sweeps = x.Average(v => (double)v.Moves);
            output.WriteLine($"  {Arms[a].Name,-20} {stages,6:F1} {sweeps,7:F1} {seconds,8:F3} {reached,3}/{r} | {mean,10:F4} ± {sd:F4} ({mean - logZ:F4}) | " +
                             $"{rmse:F4} {rmseM:F3} {rmseS:F4} {rmseK:F4} | {ratio:F2} [{low:F2}, {high:F2}] {ratioM:F2} {ratioS:F2} {ratioK:F2}; " +
                             $"least moved (median) {least[least.Length / 2]:F3}");
            csv.WriteLine(string.Join(",", Csv(f.Name), Csv(Arms[a].Name), Csv(Arms[b].Name), r, particles, F(stages), F(sweeps), F(seconds), reached, F(logZ),
                                      F(mean), F(sd), F(mean - logZ), F(rmse), F(rmseM), F(rmseS), F(rmseK), F(ratio), F(low), F(high), F(ratioM), F(ratioS),
                                      F(ratioK), F(least[least.Length / 2])));
        }
    }

    // The mean squared error of `value` about `centre` over the replicates, or over a bootstrap's pick of them.
    private static double Mse(Row[] x, Func<Row, double> value, double centre, int[]? pick = null)
    {
        double sum = 0.0;
        int n = pick?.Length ?? x.Length;
        for (int k = 0; k < n; k++)
        {
            double e = value(x[pick?[k] ?? k]) - centre;
            sum += e * e;
        }
        return sum / n;
    }

    private static string F(double v) => v.ToString("R", CultureInfo.InvariantCulture);

    private static string Csv(string s) => s.Contains(',') ? $"\"{s}\"" : s;
}
