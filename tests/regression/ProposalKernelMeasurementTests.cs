using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Repo.TestHarness;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The proposal-kernels brief's pre-registered measurement: mixing per second of each new birth against the current kernels, on hinge
/// BASS (Friedman) and on KHAOS (Ishigami), at equal sweeps.
/// <list type="bullet">
/// <item><b>Fixtures.</b> S10's Friedman (n 500, p 10, five inert; the g-form); the same with 40 more inert inputs (p 50), the history
/// birth's own regime; KHAOS's Ishigami (n 120, p 4, d_max 8, q_max 3) under the raw ridge and under the modified g-prior.</item>
/// <item><b>Arms.</b> The fit's own birth, and the birth mixture with the history birth at ω_H = ¼, δ = 0.01, a refresh every 100
/// sweeps: settings fixed before any run.</item>
/// <item><b>The design.</b> Four master seeds of four chains, 5,000 burn-in sweeps and 5,000 kept, each fit on one thread, the arms
/// interleaved by seed so they share the machine's conditions, their order alternating by seed. A short fit of each arm per
/// fixture runs first and is discarded, so no timed fit carries the process's compilation of its paths.</item>
/// <item><b>Readings,</b> per fit: the process's CPU seconds and the wall clock; for M, τ where drawn, the log posterior, σ² and the
/// monitored fitted values, R̂ and the bulk and tail ESS; births attempted and accepted in the sampling phase.</item>
/// </list>
/// Rows go to <c>proposal-kernels-measure.csv</c> in the run's directory as each fit ends; with <c>BASS_KERNEL_CSV</c> naming an
/// earlier table, its fits are skipped and new rows appended to it. Settings: <c>BASS_KERNEL_MEASURE=1</c> to run,
/// <c>BASS_KERNEL_FIXTURES</c> (a comma list; all by default), <c>BASS_KERNEL_ARMS</c> (base,history; residual for slice 2's birth at
/// ω_R = ¼), <c>BASS_KERNEL_SEEDS</c> (4), <c>BASS_KERNEL_BURN</c> and
/// <c>BASS_KERNEL_SAMPLES</c> (5,000 each).
/// </summary>
public sealed class ProposalKernelMeasurementTests(ITestOutputHelper output)
{
    private static bool Enabled => Environment.GetEnvironmentVariable("BASS_KERNEL_MEASURE") is "1";

    private static int Setting(string name, int fallback) => int.TryParse(Environment.GetEnvironmentVariable(name), out int v) ? v : fallback;

    internal static readonly HistoryBirth History = new(Share: 0.25, Floor: 0.01, Refresh: 100);

    internal static readonly ResidualBirth Residual = new(Share: 0.25);

    internal sealed record Fixture(string Name, BassConfig Config, BassData Data);

    // S10's Friedman with `extra` inert inputs appended, each uniform from its own stream.
    private static BassData Friedman(int extra)
    {
        BassData d = BassQualificationRecoveryTests.Friedman(500, 51).Data;
        if (extra == 0) return d;
        var rng = new Xoshiro256PlusPlus(52);
        int n = d.Response.Cases;
        var added = new List<BassPredictor>();
        for (int v = 0; v < extra; v++)
            added.Add(new ContinuousPredictor($"x{d.Predictors.Count + v + 1}", [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())]));
        return new BassData([.. d.Predictors, .. added], d.Response);
    }

    internal static Fixture[] Fixtures()
    {
        BassSmcMutationQualificationTests.Fixture khaos = BassSmcMutationQualificationTests.Khaos();
        return
        [
            new("friedman p10", new BassConfig(), Friedman(0)),
            new("friedman p50", new BassConfig(), Friedman(40)),
            new("ishigami raw ridge", khaos.Config, khaos.Data),
            new("ishigami modified g", khaos.Config with { Coefficients = BassCoefficientPrior.ModifiedG(TauPrior.Gamma(0.001, 1000.0)) }, khaos.Data),
        ];
    }

    private static BassConfig Arm(Fixture f, string arm, int seed, int burn, int samples)
    {
        BirthMixture? births = arm switch
        {
            "base" => null,
            "history" => new BirthMixture { History = History },
            "residual" => new BirthMixture { Residual = Residual },
            _ => throw new ArgumentOutOfRangeException(nameof(arm)),
        };
        BassSampling s = f.Config.Sampling with { Chains = 4, MasterSeed = seed, BurnIn = burn, Samples = samples, Births = births };
        return f.Config with { Sampling = s };
    }

    private static readonly string Header = "fixture,arm,seed,chains,burn,samples,cpu_s,wall_s,functional,rhat,bulk_ess,tail_ess,births_attempted,births_accepted";

    [Fact]
    public void Kernels_MixingPerSecond()
    {
        if (!Enabled) return;   // heavy, and its own decision
        string[] wanted = (Environment.GetEnvironmentVariable("BASS_KERNEL_FIXTURES") ?? "").Split(',', StringSplitOptions.RemoveEmptyEntries);
        string[] arms = (Environment.GetEnvironmentVariable("BASS_KERNEL_ARMS") ?? "base,history").Split(',', StringSplitOptions.RemoveEmptyEntries);
        string[] reversed = [.. arms.Reverse()];
        int seeds = Setting("BASS_KERNEL_SEEDS", 4), burn = Setting("BASS_KERNEL_BURN", 5_000), samples = Setting("BASS_KERNEL_SAMPLES", 5_000);
        ArtifactRun artifacts = HarnessArtifacts.Create("test-runs", "ProposalKernels", "Maths.Regression.Tests", new Dictionary<string, object?>
        {
            ["Seeds"] = seeds, ["BurnIn"] = burn, ["Samples"] = samples, ["Share"] = History.Share, ["Floor"] = History.Floor, ["Refresh"] = History.Refresh,
        });
        string csv = Environment.GetEnvironmentVariable("BASS_KERNEL_CSV") is { Length: > 0 } given ? given : Path.Combine(artifacts.RunDirectory, "proposal-kernels-measure.csv");
        var done = new HashSet<string>();
        if (File.Exists(csv))
            foreach (string line in File.ReadLines(csv).Skip(1))
            {
                string[] c = line.Split(',');
                done.Add($"{c[0]}|{c[1]}|{c[2]}");
            }
        else File.WriteAllText(csv, Header + "\n");
        output.WriteLine($"table: {csv} ({done.Count} fits there already)");

        foreach (Fixture f in Fixtures().Where(f => wanted.Length == 0 || wanted.Contains(f.Name)))
        {
            // A short fit of each arm first, discarded, so no timed fit carries the process's compilation of its paths.
            foreach (string arm in arms)
                Bass.Fit(Arm(f, arm, 99, 300, 300) with { Sampling = Arm(f, arm, 99, 300, 300).Sampling with { Chains = 1 } }, f.Data,
                         new BassRunOptions { MaxDegreeOfParallelism = 1 });
            for (int k = 0; k < seeds; k++)
                foreach (string arm in k % 2 == 0 ? arms : reversed)   // the order alternates by seed
                {
                    int seed = 101 + k;
                    if (done.Contains($"{f.Name}|{arm}|{seed}")) continue;
                    BassConfig config = Arm(f, arm, seed, burn, samples);
                    TimeSpan cpu0 = Process.GetCurrentProcess().TotalProcessorTime;
                    long wall0 = Stopwatch.GetTimestamp();
                    BassFit fit = Bass.Fit(config, f.Data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
                    double wall = (Stopwatch.GetTimestamp() - wall0) / (double)Stopwatch.Frequency;
                    double cpu = (Process.GetCurrentProcess().TotalProcessorTime - cpu0).TotalSeconds;
                    BassDiagnostics d = fit.Diagnostics;
                    BassKindCount births = d.MoveKinds[BassMoveKind.Birth];
                    var functionals = new List<(string Name, ConvergenceDiagnostic Diagnostic)> { ("M", d.Count) };
                    if (d.Tau is { } tau) functionals.Add(("tau", tau));
                    functionals.Add(("log posterior", d.LogPosterior));
                    functionals.Add(("sigma2", d.Sigma2));
                    for (int p = 0; p < d.Fitted.Count; p++) functionals.Add(($"f{p}", d.Fitted[p].Diagnostic));
                    var lines = functionals.Select(x => string.Join(",", f.Name, arm, seed.ToString(CultureInfo.InvariantCulture), 4, burn, samples,
                        cpu.ToString("R", CultureInfo.InvariantCulture), wall.ToString("R", CultureInfo.InvariantCulture), x.Name,
                        x.Diagnostic.RHat.Value.ToString("R", CultureInfo.InvariantCulture), x.Diagnostic.BulkEss.Value.ToString("R", CultureInfo.InvariantCulture),
                        x.Diagnostic.TailEss.Value.ToString("R", CultureInfo.InvariantCulture), births.Attempts, births.Accepted));
                    File.AppendAllText(csv, string.Join("\n", lines) + "\n");
                    output.WriteLine($"{f.Name}, {arm}, seed {seed}: {cpu:F1} CPU s ({wall:F1} wall); births {births.Accepted}/{births.Attempts}; " +
                                     string.Join(", ", functionals.Select(x => $"{x.Name} R̂ {x.Diagnostic.RHat.Value:F3} ESS {x.Diagnostic.BulkEss.Value:F0}")));
                }
        }
    }

    /// <summary>
    /// MCMC's climb from empty on the SMC trap's fixture (KHAOS's Ishigami, the raw ridge; <c>BASS_KERNEL_MEASURE=1</c>): a chain at b = 1
    /// from M = 0, under the fit's own birth and under the residual birth at ω_R = ¼, eight seeds each, until its collapsed log posterior
    /// first reaches the band's floor, the 5% quantile of a long reference chain's last 30,000 sweeps (60,000 from empty, the fit's own
    /// birth). The reading: the sweeps and CPU seconds to reach it, and M there; a chain that has not within 20,000 sweeps is reported so.
    /// </summary>
    [Fact]
    public void ResidualBirth_ClimbFromEmpty()
    {
        if (!Enabled) return;
        BassSmcMutationQualificationTests.Fixture khaos = BassSmcMutationQualificationTests.Khaos();
        BassModel baseModel = BassModel.Create(khaos.Config, khaos.Data);
        BassModel residualModel = BassModel.Create(khaos.Config with { Sampling = khaos.Config.Sampling with { Births = new BirthMixture { Residual = Residual } } }, khaos.Data);
        double logTau = Math.Log(baseModel.Target.HeldTau!.Value);
        var reference = new BassChain(baseModel, [], logTau, new Xoshiro256PlusPlus(7001), new Xoshiro256PlusPlus(7002));
        var tail = new List<double>();
        for (int t = 0; t < 60_000; t++)
        {
            BassState s = reference.Step();
            if (t >= 30_000) tail.Add(BassMixing.LogPosterior(baseModel, s));
        }
        double floor = BassMixing.Quantile([.. tail], 0.05);
        output.WriteLine($"the band's floor: ln π {floor:F2} (the reference's last 30,000 sweeps: median {BassMixing.Quantile([.. tail], 0.5):F2}, M {reference.Current.Count})");
        foreach ((string arm, BassModel model) in new[] { ("base", baseModel), ("residual", residualModel) })
        {
            var reached = new List<(int Sweeps, double Seconds, int M)>();
            for (int seed = 0; seed < 8; seed++)
            {
                var chain = new BassChain(model, [], logTau, new Xoshiro256PlusPlus(8000 + 2 * seed), new Xoshiro256PlusPlus(8001 + 2 * seed));
                TimeSpan cpu0 = Process.GetCurrentProcess().TotalProcessorTime;
                int sweep = 0;
                for (; sweep < 20_000; sweep++)
                    if (BassMixing.LogPosterior(model, chain.Step()) >= floor) break;
                double seconds = (Process.GetCurrentProcess().TotalProcessorTime - cpu0).TotalSeconds;
                reached.Add((sweep + 1, seconds, chain.Current.Count));
                BassKindCount births = chain.Structure.Kinds[BassMoveKind.Birth];
                output.WriteLine($"{arm}, seed {seed}: {(sweep < 20_000 ? $"reached in {sweep + 1} sweeps" : "not within 20,000 sweeps")}, {seconds:F2} CPU s, " +
                                 $"M {chain.Current.Count}; births {births.Accepted}/{births.Attempts}");
            }
            output.WriteLine($"{arm}: median sweeps {BassMixing.Quantile([.. reached.Select(r => (double)r.Sweeps)], 0.5):F0}, " +
                             $"median CPU s {BassMixing.Quantile([.. reached.Select(r => r.Seconds)], 0.5):F2}, reached {reached.Count(r => r.Sweeps <= 20_000)}/8");
        }
    }

    /// <summary>
    /// The cost row (<c>BASS_KERNEL_COST=1</c>; the clock is <see cref="BassCost.Tick"/>'s, the cycle counter with
    /// <c>BASS_COST_CLOCK=cycles</c>, indicative): a birth's kernel work, its draw and its ratio, and a death's ratio, under the
    /// fit's own kernel and under the mixture with a history learned over 500 sweeps, at a burned chain's state. The design's
    /// work, the factor and the marginal, is the same under both and is left out; the measurement's CPU seconds read the whole sweep.
    /// Points: the cost report's base (the rich surface), the same at p 50, and KHAOS's Ishigami.
    /// </summary>
    [Fact]
    public void HistoryBirth_CostRow()
    {
        if (Environment.GetEnvironmentVariable("BASS_KERNEL_COST") is not "1") return;
        BassSmcMutationQualificationTests.Fixture khaos = BassSmcMutationQualificationTests.Khaos();
        BassModel ishigami = BassModel.Create(khaos.Config, khaos.Data);
        BassChain ishigamiChain = BassMixing.Chain(ishigami, 0);
        for (int t = 0; t < 2_000; t++) ishigamiChain.Step();
        (BassModel Model, BassChain Chain, int _) atBase = BassCost.Burned(BassCost.Base), atWide = BassCost.Burned(BassCost.Base with { P = 50 });
        (string Name, BassModel Model, BassChain Chain)[] points =
        [
            ("rich, the cost report's base", atBase.Model, atBase.Chain),
            ("rich at p 50", atWide.Model, atWide.Chain),
            ("KHAOS Ishigami, raw ridge", ishigami, ishigamiChain),
        ];
        double sink = 0.0;
        foreach ((string name, BassModel model, BassChain chain) in points)
        {
            var history = new BassInclusionHistory(model.Prior, History);
            history.Learn();
            for (int t = 0; t < 500; t++) history.Record(chain.Step().Terms);
            history.Freeze();
            var mixture = new BassBirthMixture(model.Kernel, History, history);
            BassTermCounts counts = BassMixing.Counts(model, chain.Current);
            BassTerm[] terms = [.. chain.Current.Terms];
            var rng = new Xoshiro256PlusPlus(17);
            double Time(Action step, int reps)
            {
                for (int k = 0; k < reps / 10; k++) step();
                long t0 = BassCost.Tick();
                for (int k = 0; k < reps; k++) step();
                return (BassCost.Tick() - t0) / BassCost.TicksPerSecond / reps * 1e6;
            }
            var at = new BassBirthContext(terms, 1.0);
            BassBirthContext[] reduced = [.. Enumerable.Range(0, terms.Length).Select(k => new BassBirthContext([.. terms.Where((_, i) => i != k)], 1.0))];
            double Birth(Func<IBassProposalKernel> kernelOf, int reps) => Time(() =>
            {
                IBassProposalKernel kernel = kernelOf();
                if (kernel.DrawBirth(rng, counts, at) is { } born) sink += kernel.LogBirthRatio(born, counts, at);
            }, reps);
            double Death(Func<IBassProposalKernel> kernelOf, int reps) => Time(() =>
            {
                IBassProposalKernel kernel = kernelOf();
                for (int k = 0; k < terms.Length; k++) sink += kernel.LogDeathRatio(terms[k], counts, reduced[k]);
            }, Math.Max(1, reps / Math.Max(1, terms.Length))) / Math.Max(1, terms.Length);
            var residualMixture = new BassBirthMixture(model.Kernel, null, null, new BassResidualBirth(model), Residual);
            IBassProposalKernel Cold() => new BassBirthMixture(model.Kernel, null, null, new BassResidualBirth(model), Residual);
            double bb = Birth(() => model.Kernel, 20_000), bm = Birth(() => mixture, 20_000), db = Death(() => model.Kernel, 20_000), dm = Death(() => mixture, 20_000);
            double br = Birth(() => residualMixture, 20_000), dr = Death(() => residualMixture, 20_000), bc = Birth(Cold, 200), dc = Death(Cold, 200 * Math.Max(1, terms.Length));
            output.WriteLine($"{name}: M {terms.Length}; birth {bb:F2} µs base, {bm:F2} µs with the history ({bm / bb:F2}×), {br:F2} µs with the residual warm " +
                             $"({br / bb:F1}×), {bc:F1} µs cold ({bc / bb:F0}×); death {db:F2} µs base, {dm:F2} µs with the history ({dm / db:F2}×), " +
                             $"{dr:F2} µs with the residual warm ({dr / db:F1}×), {dc:F1} µs cold ({dc / db:F0}×) [{(BassCost.Cycles ? "cycles" : "wall")}]");
        }
        Assert.True(double.IsFinite(sink));
    }
}
