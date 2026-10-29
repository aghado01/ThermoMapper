using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime;
using System.Text;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Repo.TestHarness;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// F-M2's Phase 2, the pre-registered pass (the local-move pass brief): seven arms against the baseline on four fixtures,
/// five paired master seeds each, each fit read by its least rank-normalized bulk ESS per second and the arms judged by the
/// rule the user decided on 2026-10-29.
/// <list type="bullet">
/// <item><b>Arms</b> (<see cref="Arms"/>): S10's configuration with the move schedule, the local steps, the fixed-M family
/// and the tempered split rule set as the brief's table. Tempered arms count burn-in and samples in scans.</item>
/// <item><b>Fixtures</b> (<see cref="Fixtures"/>): S10's three and Friedman at n = 2,000, master seeds s₀ + 1,000k for
/// k &lt; 5, the same five for every arm.</item>
/// <item><b>Readings:</b> the bulk ESS and R̂ of M, σ², the collapsed log posterior at b = 1 and the fitted values at the
/// default monitored points, from the shared estimator (<see cref="ChainDiagnostics.Diagnose(double[][])"/>); the fit's
/// least ESS over them, a constant quantity left out; that over the whole <see cref="Bass.Fit"/>'s seconds. Context: the
/// acceptance by move kind, each chain's median M, σ's coverage and the tempered ladders.</item>
/// <item><b>The rule</b> (<see cref="Rule"/>): r_f = exp(mean ln e_arm − mean ln e_A0) over the paired seeds, the score
/// the geometric mean of r_f over the fixtures; an arm clears at a score of at least 1.5 with no r_f below 0.8, and the
/// highest clearing score wins.</item>
/// </list>
/// Every fit runs alone on one thread (MaxDegreeOfParallelism 1) after an untimed warm-up of every arm, the arms shuffled
/// within each (fixture, seed) by a fixed seed; A0 reruns at the end for each fixture's first seed, a timing control.
/// BASS_PASS_SCALE divides the lengths and BASS_PASS_SEEDS keeps the first seeds, for the dry run. The per-fit CSV and the
/// verdict's tables land under <c>artifacts/test-runs/BassMovePass/</c>. The pass and the log posterior's check are
/// qualification runs (<see cref="BassQualification.Enabled"/>); the arms', the data's and the minimum's checks run always.
/// </summary>
public sealed class BassMovePassTests
{
    /// <summary>Master seeds per fixture, their stride, and the seed the arms' shuffle derives from.</summary>
    internal const int Seeds = 5, SeedStride = 1000, OrderSeed = 20261029;

    /// <summary>The rule's thresholds: the score an arm must reach, and the least r_f it may hold.</summary>
    internal const double Clears = 1.5, Floor = 0.8;

    // S10's Friedman data at n = 500, hashed before Friedman took its sample size. The hash reads the default config's JSON
    // too: 6b52ce789967ee0c until N14 added Sampling.StoppingScope to it, the data unchanged (the old JSON still gives it).
    private const string S10Friedman = "e38c2c6c05ff8dbc";

#if DEBUG
    private const string Configuration = "Debug";
#else
    private const string Configuration = "Release";
#endif

    private static readonly string[] QuantityNames = ["M", "sigma2", "logpost", "f1", "f2", "f3"];

    private readonly ITestOutputHelper _out;
    public BassMovePassTests(ITestOutputHelper output) => _out = output;

    /// <summary>An arm: S10's sampling with the brief's changes to the schedule, the local steps, the family and tempering.</summary>
    internal sealed record Arm(string Name, string Label, bool PriorCancelling, bool Steps, bool Family, bool Tempered)
    {
        public BassSampling On(BassSampling sampling) => sampling with
        {
            Moves = PriorCancelling ? new PriorCancelling() : new EqualOverAvailable(),
            LocalSteps = Steps ? new LocalSteps() : null,
            FactorMoves = Family ? new FactorMoves() : null,
            Tempering = Tempered ? new TemperingSpec { Split = new SplitRule() } : null,
        };
    }

    /// <summary>The brief's arms, A0 the baseline.</summary>
    internal static readonly Arm[] Arms =
    [
        new("A0", "baseline", false, false, false, false),
        new("A1", "local steps", false, true, false, false),
        new("A2", "fixed-M family", false, false, true, false),
        new("A3", "both", false, true, true, false),
        new("A4", "prior-cancelling", true, false, false, false),
        new("A5", "prior-cancelling with both", true, true, true, false),
        new("A6", "tempered", false, false, false, true),
        new("A7", "tempered with both", false, true, true, true),
    ];

    /// <summary>A fixture: its data, σ's truth, and s₀.</summary>
    internal sealed record Fixture(string Name, string Label, BassData Data, double Sigma, int Seed);

    /// <summary>S10's three fixtures and Friedman at n = 2,000, from its generator and data seed.</summary>
    internal static Fixture[] Fixtures() =>
    [
        new("friedman", "Friedman", BassQualificationRecoveryTests.Friedman().Data, 1.0, 53),
        new("mixed", "mixed", BassQualificationRecoveryTests.Mixed().Data, 0.5, 57),
        new("functional", "functional", BassQualificationRecoveryTests.Functional().Data, 0.2, 61),
        new("friedman2000", "Friedman, n = 2,000", BassQualificationRecoveryTests.Friedman(2000).Data, 1.0, 65),
    ];

    /// <summary>Arm <paramref name="arm"/> on <paramref name="fixture"/> at master seed s₀ + 1,000<paramref name="k"/>, its
    /// lengths divided by <paramref name="scale"/>.</summary>
    internal static BassConfig Config(Fixture fixture, int k, Arm arm, int scale = 1)
    {
        BassConfig s10 = BassQualificationRecoveryTests.Config(fixture.Seed + SeedStride * k);
        BassSampling sampling = arm.On(s10.Sampling);
        return s10 with { Sampling = sampling with { BurnIn = sampling.BurnIn / scale, Samples = sampling.Samples / scale } };
    }

    [Fact]
    public void Friedman_AtN500_IsS10s()
    {
        Assert.Equal(S10Friedman, DataHash(BassQualificationRecoveryTests.Friedman()));
        Assert.Equal(S10Friedman, DataHash(BassQualificationRecoveryTests.Friedman(500)));

        // n = 2,000 continues the generator's stream, so its first 500 cases are S10's.
        (BassData small, double[] smallTruth) = BassQualificationRecoveryTests.Friedman();
        (BassData large, double[] largeTruth) = BassQualificationRecoveryTests.Friedman(2000);
        Assert.Equal(2000, large.Response.Cases);
        for (int v = 0; v < small.Predictors.Count; v++)
            Assert.Equal(Bits(((ContinuousPredictor)small.Predictors[v]).Values.Span), Bits(((ContinuousPredictor)large.Predictors[v]).Values.Span[..500]));
        Assert.Equal(Bits(((ScalarResponse)small.Response).Y.Span), Bits(((ScalarResponse)large.Response).Y.Span[..500]));
        Assert.Equal(Bits(smallTruth), Bits(largeTruth.AsSpan(0, 500)));
    }

    private static string DataHash((BassData Data, double[] Truth) d)
        => new RunFingerprint().Add(Bass.Fingerprint(new BassConfig(), d.Data)).Add(d.Truth).ToString();

    private static long[] Bits(ReadOnlySpan<double> values)
    {
        var bits = new long[values.Length];
        for (int i = 0; i < bits.Length; i++) bits[i] = BitConverter.DoubleToInt64Bits(values[i]);
        return bits;
    }

    /// <summary>
    /// The arms are what the table says: each sets the schedule, the local steps, the family and tempering as its row, its
    /// resolved echo (the step windows, the factor shares, the split threshold) follows, its lengths and seed are S10's,
    /// and its fingerprint differs from every other arm's, while with its four fields set back it is S10's own.
    /// </summary>
    [Fact]
    public void Arms_AreWhatTheySay()
    {
        foreach (Fixture f in Fixtures())
        {
            BassConfig s10 = BassQualificationRecoveryTests.Config(f.Seed);
            string baseline = Bass.Fingerprint(s10, f.Data);
            BassInputs inputs = BassModel.Create(s10, f.Data).Inputs;
            var prints = new HashSet<string>();
            foreach (Arm arm in Arms)
            {
                BassConfig config = Config(f, 0, arm);
                BassSampling s = config.Sampling;
                Assert.Equal((4, f.Seed, 5000, 2000, 1), (s.Chains, s.MasterSeed, s.BurnIn, s.Samples, s.Thin));
                Assert.Equal(arm.PriorCancelling ? new PriorCancelling(0.4) : (MoveSchedule)new EqualOverAvailable(), s.Moves);
                Assert.Equal(arm.Steps ? new LocalSteps(0.5, 0.02) : null, s.LocalSteps);
                Assert.Equal(arm.Family ? new FactorMoves(1, 1, 1) : null, s.FactorMoves);
                Assert.Equal(arm.Tempered ? new TemperingSpec { Split = new SplitRule(0.9) } : null, s.Tempering);
                Assert.True(s.Tempering is null or { Levels: null, Adaptive: true, ReferenceRung: false });

                BassResolved resolved = BassResolved.Resolve(config, f.Data);
                if (arm.Steps)
                {
                    // w = max(1, ⌈0.02·n_v⌉) over each hinge variable's knots, taken in decimal.
                    static int Window(int knots) => (int)Math.Max(1m, Math.Ceiling(0.02m * knots));
                    Assert.Equal(Enumerable.Range(0, inputs.Cases.ContinuousCount).Select(v => Window(inputs.ContinuousKnots(v).Length)),
                                 resolved.StepWindows!.Continuous);
                    Assert.Equal(Enumerable.Range(0, inputs.Outputs?.Count ?? 0).Select(v => Window(inputs.FunctionalKnots(v).Length)),
                                 resolved.StepWindows.Functional);
                }
                else Assert.Null(resolved.StepWindows);
                if (arm.Family)
                {
                    BassFactorShares shares = resolved.FactorShares!;
                    Assert.Equal(new[] { 1 / 3.0, 1 / 6.0, 1 / 6.0, 1 / 3.0 }, new[] { shares.Change, shares.Prune, shares.Grow, shares.VariableChange },
                                 new ToleranceComparer(1e-15));
                }
                else Assert.Null(resolved.FactorShares);
                Assert.Equal(arm.Tempered ? 0.9 : (double?)null, resolved.SplitThreshold);

                BassConfig reset = config with
                {
                    Sampling = s with { Moves = new EqualOverAvailable(), LocalSteps = null, FactorMoves = null, Tempering = null },
                };
                Assert.Equal(baseline, Bass.Fingerprint(reset, f.Data));
                Assert.True(prints.Add(Bass.Fingerprint(config, f.Data)), $"{f.Name}, {arm.Name}: the fingerprint repeats another arm's.");
            }
            _out.WriteLine($"{f.Name}: eight distinct fingerprints, S10's {baseline} under each arm's fields set back");
        }
    }

    private sealed class ToleranceComparer(double tolerance) : IEqualityComparer<double>
    {
        public bool Equals(double a, double b) => Math.Abs(a - b) <= tolerance;
        public int GetHashCode(double value) => 0;
    }

    /// <summary>
    /// The log posterior read from a retained draw (<see cref="LogPosterior"/>) is the chain's own log density at that
    /// state: each chain replayed as the run takes it, untempered or its ladder burned and scanned, and at three of its
    /// draws the state's structure and τ match the draw's and the chain's b = 1 density matches the value from the draw's
    /// structure and τ, to rounding. A0, A5 and A7 on each fixture, two chains of 255 and 40.
    /// </summary>
    [Fact]
    public void LogPosterior_IsTheChainsOwn()
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int burn = 255, samples = 40;
        foreach (Fixture f in Fixtures())
            foreach (Arm arm in new[] { Arms[0], Arms[5], Arms[7] })
            {
                BassConfig config = Config(f, 0, arm);
                config = config with { Sampling = config.Sampling with { Chains = 2, BurnIn = burn, Samples = samples } };
                BassFit fit = Bass.Fit(config, f.Data);
                BassModel model = BassModel.Create(config, f.Data);
                var fits = new Dictionary<string, BassStructureFit>();
                double worst = 0.0;
                for (int c = 0; c < config.Sampling.Chains; c++)
                {
                    Func<BassState> step;
                    Func<double> density;
                    if (arm.Tempered)
                    {
                        BassLadder ladder = BassMixing.Ladder(model, c, config.Sampling.Tempering!);
                        BassMixing.Burn(ladder, burn);
                        step = ladder.Step;
                        density = () => ladder.Moves.CurrentLogDensity;
                    }
                    else
                    {
                        BassChain chain = BassMixing.Chain(model, c);
                        for (int s = 0; s < burn; s++) chain.Step();
                        step = chain.Step;
                        density = () => chain.Moves.CurrentLogDensity;
                    }
                    for (int s = 0; s < samples; s++)
                    {
                        BassState state = step();
                        if (s is not (0 or samples / 2 or samples - 1)) continue;
                        BassDraw draw = fit.Draws[c][s];
                        Assert.Equal(BassMixing.Key(state.Terms), BassMixing.Key(draw.Structure.TermIds.Select(id => fit.Terms[id])));
                        Assert.Equal(BitConverter.DoubleToInt64Bits(state.Tau), BitConverter.DoubleToInt64Bits(draw.Tau));
                        double own = density(), read = LogPosterior(model, fit.Terms, draw, fits);
                        Assert.Equal(BitConverter.DoubleToInt64Bits(own), BitConverter.DoubleToInt64Bits(Temperature.LogPosterior(model.Target.Evaluate(state))));
                        Assert.True(Math.Abs(read - own) <= 1e-9 * (1.0 + Math.Abs(own)), $"{f.Name}, {arm.Name}, chain {c}, draw {s}: {read:R} against the chain's {own:R}");
                        worst = Math.Max(worst, Math.Abs(read - own));
                    }
                }
                _out.WriteLine($"{f.Name}, {arm.Name}: the draws' log posterior within {worst:G3} of the chains' own");
            }
    }

    /// <summary>The collapsed log π(structure, u | y) at b = 1 of every retained draw, chain by chain, each structure's fit
    /// taken once.</summary>
    internal static double[][] LogPosteriors(BassModel model, BassFit fit)
    {
        var fits = new Dictionary<string, BassStructureFit>();
        return [.. fit.Draws.Select(chain => chain.Select(draw => LogPosterior(model, fit.Terms, draw, fits)).ToArray())];
    }

    /// <summary>
    /// The collapsed log π(structure, u | y) at b = 1 of <paramref name="draw"/>, from its structure and τ: the structure
    /// loaded whole on a fresh design (its least-squares statistics, which do not depend on τ, kept in
    /// <paramref name="fits"/> by its ids), then the target's prior and log marginal at u = ln τ.
    /// </summary>
    internal static double LogPosterior(BassModel model, BassTermTable terms, BassDraw draw, Dictionary<string, BassStructureFit> fits)
    {
        string key = string.Join(",", draw.Structure.TermIds);
        if (!fits.TryGetValue(key, out BassStructureFit fit))
        {
            BassState state = BassMixing.Loaded(model, [.. draw.Structure.TermIds.Select(id => terms[id])], 0.0, out BassRefusal? refusal)
                              ?? throw new InvalidOperationException($"A retained structure lies outside F ({refusal}).");
            fits[key] = fit = state.Fit;
        }
        return Temperature.LogPosterior(model.Target.Evaluate(fit, Math.Log(draw.Tau)));
    }

    /// <summary>The pass: every arm on every fixture and seed, the timing control, the rule and the verdict's tables.</summary>
    [Fact]
    public void Pass()
    {
        if (!BassQualification.Enabled) return;   // heavy
        int scale = Setting("BASS_PASS_SCALE", 1), seeds = Setting("BASS_PASS_SEEDS", Seeds);
        Fixture[] fixtures = Fixtures();
        ArtifactRun run = HarnessArtifacts.Create("test-runs", "BassMovePass", "Maths.Regression.Tests", new Dictionary<string, object?>
        {
            ["Scale"] = scale, ["Seeds"] = seeds, ["OrderSeed"] = OrderSeed, ["Configuration"] = Configuration,
            ["ServerGC"] = GCSettings.IsServerGC, ["Processors"] = Environment.ProcessorCount,
        });
        string csvPath = Path.Combine(run.RunDirectory, "local-move-pass-results.csv");
        using var csv = new StreamWriter(csvPath) { NewLine = "\n" };
        csv.WriteLine(Reading.Header);
        var clock = Stopwatch.StartNew();

        // The warm-up: every arm on every fixture, short and discarded, so no timed fit pays a first compilation.
        foreach (Fixture f in fixtures)
            foreach (Arm arm in Arms) Read(f, 0, arm, Math.Max(20, scale), "warm-up", -1, -1);
        _out.WriteLine($"warm-up {clock.Elapsed.TotalSeconds:F0} s ({Configuration})");

        var pass = new List<Reading>();
        int[] orders = SeedTree.Derive(OrderSeed, seeds * fixtures.Length);
        int block = 0;
        for (int k = 0; k < seeds; k++)
            foreach (Fixture f in fixtures)
            {
                Arm[] order = Shuffled(Arms, orders[block]);
                var held = new List<Reading>();
                for (int position = 0; position < order.Length; position++)
                {
                    Reading r = Read(f, k, order[position], scale, "pass", block, position);
                    csv.WriteLine(r.Row());
                    csv.Flush();
                    held.Add(r);
                }
                // The pairing: the arms of one (fixture, seed) start from the same exact prior draws.
                Assert.Single(held.Select(r => r.StartHash).Distinct());
                pass.AddRange(held);
                _out.WriteLine($"{f.Name}, k = {k}: {string.Join(", ", held.Select(r => $"{r.Arm} {r.Seconds:F1} s"))}");
                block++;
            }

        // The timing control: A0 again for each fixture's first seed, at the end; its draws the same bit for bit.
        var controls = new List<(Reading First, Reading Again)>();
        foreach (Fixture f in fixtures)
        {
            Reading again = Read(f, 0, Arms[0], scale, "control", block++, 0);
            csv.WriteLine(again.Row());
            csv.Flush();
            Reading first = pass.Single(r => r.Fixture == f.Name && r.K == 0 && r.Arm == Arms[0].Name);
            Assert.Equal(first.DrawsHash, again.DrawsHash);
            controls.Add((first, again));
        }

        Verdict verdict = Rule(pass, [.. fixtures.Select(f => f.Name)], seeds);
        string report = Report(verdict, pass, controls, fixtures, scale, seeds, clock.Elapsed);
        _out.WriteLine(report);
        _out.WriteLine($"CSV {csvPath}");
        _out.WriteLine($"verdict {run.WriteText("local-move-pass-verdict.md", report)}");
    }

    private static int Setting(string name, int fallback)
        => int.TryParse(Environment.GetEnvironmentVariable(name), NumberStyles.Integer, CultureInfo.InvariantCulture, out int v) && v > 0 ? v : fallback;

    // The arms in an order drawn from `seed` (Fisher–Yates).
    private static Arm[] Shuffled(Arm[] arms, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        Arm[] order = [.. arms];
        for (int i = order.Length - 1; i > 0; i--)
        {
            int j = (int)(rng.NextDouble() * (i + 1));
            (order[i], order[j]) = (order[j], order[i]);
        }
        return order;
    }

    /// <summary>One fit's readings, a row of the CSV.</summary>
    internal sealed record Reading(
        string Role, string Fixture, int K, int MasterSeed, string Arm, int Block, int Position,
        double Seconds, double CpuSeconds, double BurnSeconds, double SweepsSeconds,
        (string Name, ConvergenceDiagnostic Diagnostic)[] Quantities, double MinEss, string MinQuantity, string Constant,
        double SigmaLower, double SigmaUpper, bool SigmaCovers, int[] MedianM, BassKindCount[] Kinds,
        int[]? Levels, int[]? SteepPairs, double RoundTrips, double IdealRoundTrips, double MaxPairRejection, double Barrier,
        string StartHash, string DrawsHash, string Fingerprint)
    {
        /// <summary>The fit's least ESS per second of its whole run.</summary>
        public double MinEssPerSecond => MinEss / Seconds;

        /// <summary>The sampling phase's share of the run: from the last burn-in sweep to the last sweep.</summary>
        public double SamplingShare => (SweepsSeconds - BurnSeconds) / Seconds;

        public static string Header => string.Join(",", (string[])
            ["role", "fixture", "k", "master_seed", "arm", "block", "position", "seconds", "cpu_seconds", "burn_seconds", "sweeps_seconds", "sampling_share",
             .. QuantityNames.SelectMany(q => new[] { $"ess_{q}", $"rhat_{q}" }),
             "min_ess", "min_quantity", "constant", "min_ess_per_s", "sigma_lower", "sigma_upper", "sigma_covers", "median_M",
             .. Enum.GetValues<BassMoveKind>().SelectMany(k => new[] { "attempts", "accepted", "null", "refused" }.Select(c => $"{k.ToString().ToLowerInvariant()}_{c}")),
             "levels", "steep_pairs", "round_trips_per_scan", "ideal_round_trips_per_scan", "max_pair_rejection", "barrier",
             "start_hash", "draws_hash", "fingerprint"]);

        public string Row() => string.Join(",", (string[])
            [Role, Fixture, I(K), I(MasterSeed), Arm, I(Block), I(Position), F(Seconds), F(CpuSeconds), F(BurnSeconds), F(SweepsSeconds), F(SamplingShare),
             .. Quantities.SelectMany(q => new[] { q.Diagnostic.BulkEss.ToString("R", CultureInfo.InvariantCulture), q.Diagnostic.RHat.ToString("R", CultureInfo.InvariantCulture) }),
             F(MinEss), MinQuantity, Constant, F(MinEssPerSecond), F(SigmaLower), F(SigmaUpper), SigmaCovers ? "1" : "0", string.Join(";", MedianM),
             .. Kinds.SelectMany(k => new[] { k.Attempts, k.Accepted, k.Null, k.Refused }.Select(c => c.ToString(CultureInfo.InvariantCulture))),
             Levels is null ? "" : string.Join(";", Levels), SteepPairs is null ? "" : string.Join(";", SteepPairs),
             F(RoundTrips), F(IdealRoundTrips), F(MaxPairRejection), F(Barrier), StartHash, DrawsHash, Fingerprint]);
    }

    private static string F(double value) => double.IsNaN(value) ? "" : value.ToString("R", CultureInfo.InvariantCulture);
    private static string I(int value) => value.ToString(CultureInfo.InvariantCulture);

    // One fit timed on one thread, then read: everything after Bass.Fit is untimed.
    private static Reading Read(Fixture f, int k, Arm arm, int scale, string role, int block, int position)
    {
        BassConfig config = Config(f, k, arm, scale);
        BassSampling s = config.Sampling;
        long burn = (long)s.Chains * s.BurnIn, planned = burn + (long)s.Chains * s.Samples * s.Thin;
        var phases = new PhaseClock(burn, planned);
        GC.Collect();
        GC.WaitForPendingFinalizers();
        GC.Collect();
        TimeSpan cpu = Process.GetCurrentProcess().TotalProcessorTime;
        phases.Start();
        BassFit fit = Bass.Fit(config, f.Data, new BassRunOptions { MaxDegreeOfParallelism = 1, Progress = phases, ProgressInterval = checked((int)burn) });
        double seconds = phases.Stop();
        double cpuSeconds = (Process.GetCurrentProcess().TotalProcessorTime - cpu).TotalSeconds;

        BassModel model = BassModel.Create(config, f.Data);
        BassDiagnostics d = fit.Diagnostics;
        double[][] m = Trace(fit, draw => draw.Structure.TermIds.Length);
        // The shared estimator: the fit's own readings of M and σ² are Diagnose's on the retained draws.
        Assert.Equal(d.Count, ChainDiagnostics.Diagnose(m));
        Assert.Equal(d.Sigma2, ChainDiagnostics.Diagnose(Trace(fit, draw => draw.Sigma2)));
        (string Name, ConvergenceDiagnostic Diagnostic)[] quantities =
        [
            ("M", d.Count), ("sigma2", d.Sigma2), ("logpost", ChainDiagnostics.Diagnose(LogPosteriors(model, fit))),
            .. d.Fitted.Select((v, p) => ($"f{p + 1}", v.Diagnostic)),
        ];
        Assert.Equal(QuantityNames, quantities.Select(q => q.Name));
        (double minEss, string minQuantity, string constant) = Minimum(quantities, m);

        double[] sd = [.. fit.Draws.SelectMany(c => c).Select(draw => Math.Sqrt(draw.Sigma2)).Order()];
        double lower = PosteriorInterval.Quantile(sd, 0.025), upper = PosteriorInterval.Quantile(sd, 0.975);
        int[] medianM = [.. m.Select(chain => (int)chain.Order().ElementAt(chain.Length / 2))];
        BassTemperingReadout? t = d.Tempering;
        int[]? steep = t?.PairRejection.Select(r => r.Count(x => x > 0.9)).ToArray();
        return new Reading(role, f.Name, k, s.MasterSeed, arm.Name, block, position, seconds, cpuSeconds, phases.BurnEnd, phases.SweepsEnd,
                           quantities, minEss, minQuantity, constant, lower, upper, lower <= f.Sigma && f.Sigma <= upper, medianM,
                           [.. d.MoveKinds.Entries], t?.Levels, steep,
                           t?.RoundTripsPerScan ?? double.NaN, t?.IdealRoundTripRate ?? double.NaN, t?.MaxPairRejection ?? double.NaN,
                           t?.Barrier ?? double.NaN, StartHash(d), BassBitIdentityTests.Hash(fit), fit.Identity.Fingerprint);
    }

    private static double[][] Trace(BassFit fit, Func<BassDraw, double> value) => [.. fit.Draws.Select(chain => chain.Select(value).ToArray())];

    /// <summary>
    /// The fit's least bulk ESS over its quantities, and which set it. A quantity constant over every draw has no ESS and
    /// is left out, named in the third value. M frozen (each split chain constant, not all at one value) has no
    /// rank-normalized estimate, so it reads the shared estimator's ESS of its mean, about the split chains' count over
    /// two: set before any run. σ², the log posterior and the fitted values carry each draw's own τ or coefficients and
    /// cannot freeze; a reading of theirs without an estimate stops the pass.
    /// </summary>
    internal static (double Ess, string Quantity, string Constant) Minimum(IReadOnlyList<(string Name, ConvergenceDiagnostic Diagnostic)> quantities,
                                                                          double[][] m)
    {
        double least = double.PositiveInfinity;
        string which = "";
        var constant = new List<string>();
        foreach ((string name, ConvergenceDiagnostic diagnostic) in quantities)
        {
            EssDiagnostic ess = diagnostic.BulkEss;
            string label = name;
            if (ess.State == DiagnosticState.Constant)
            {
                constant.Add(name);
                continue;
            }
            if (ess.State == DiagnosticState.Frozen && name == "M")
            {
                ess = ChainDiagnostics.MeanEss(m);
                label = "M (frozen)";
            }
            if (ess.State != DiagnosticState.Estimated) throw new InvalidOperationException($"{name}'s bulk ESS is {ess}.");
            if (ess.Value < least) (least, which) = (ess.Value, label);
        }
        return (least, which, string.Join(";", constant));
    }

    /// <summary>The minimum's edge cases on synthetic chains: a constant quantity is left out and named; a frozen M (each
    /// split chain constant, at different values) reads the ESS of its mean; any other reading without an estimate stops.</summary>
    [Fact]
    public void Minimum_LeavesAConstantOut_AndReadsAFrozenMByItsMean()
    {
        var rng = new Xoshiro256PlusPlus(seed: 7);
        double[][] Chains(Func<int, int, double> value) => [.. Enumerable.Range(0, 4).Select(c => Enumerable.Range(0, 100).Select(s => value(c, s)).ToArray())];
        double[][] noise = Chains((_, _) => rng.NextDouble()), constant = Chains((_, _) => 5.0), frozen = Chains((c, s) => c + (s < 50 ? 0 : 10));
        ConvergenceDiagnostic moving = ChainDiagnostics.Diagnose(noise), still = ChainDiagnostics.Diagnose(constant), stuck = ChainDiagnostics.Diagnose(frozen);
        Assert.Equal(DiagnosticState.Constant, still.BulkEss.State);
        Assert.Equal(DiagnosticState.Frozen, stuck.BulkEss.State);

        Assert.Equal((moving.BulkEss.Value, "sigma2", "M"), Minimum([("M", still), ("sigma2", moving)], constant));
        EssDiagnostic mean = ChainDiagnostics.MeanEss(frozen);
        Assert.Equal(DiagnosticState.Estimated, mean.State);
        Assert.True(mean.Value < moving.BulkEss.Value);
        Assert.Equal((mean.Value, "M (frozen)", ""), Minimum([("M", stuck), ("sigma2", moving)], frozen));
        Assert.Throws<InvalidOperationException>(() => Minimum([("M", moving), ("sigma2", stuck)], noise));
    }

    // Each chain's start, an exact prior draw: its terms, τ and attempts.
    private static string StartHash(BassDiagnostics d)
    {
        var hash = new RunFingerprint();
        foreach (BassChainAccount chain in d.Chains)
        {
            hash.Add(chain.Start.Terms.Count).Add(chain.Start.Tau).Add((ulong)chain.Start.Attempts);
            foreach (BassTerm term in chain.Start.Terms) hash.Add(term.ToString());
        }
        return hash.ToString();
    }

    // The run's phases from its progress: the time at the last burn-in sweep over every chain (the ensemble burns every
    // chain before it samples), and at the last sweep, reported after the ensemble's own diagnostics.
    private sealed class PhaseClock(long burn, long planned) : IProgress<BassProgress>
    {
        private readonly Stopwatch _watch = new();

        public double BurnEnd { get; private set; } = double.NaN;

        public double SweepsEnd { get; private set; } = double.NaN;

        public void Start() => _watch.Start();

        public double Stop()
        {
            _watch.Stop();
            return _watch.Elapsed.TotalSeconds;
        }

        public void Report(BassProgress value)
        {
            double now = _watch.Elapsed.TotalSeconds;
            if (value.Sweeps == burn) BurnEnd = now;
            if (value.Sweeps == planned) SweepsEnd = now;
        }
    }

    /// <summary>One arm's standing under the rule: ln r_f per fixture, the per-seed log ratios' standard deviation and
    /// the standard error of ln r_f, the score, whether it clears, and the score with each fixture left out.</summary>
    internal sealed record Standing(string Arm, double[] LogRatio, double[] Spread, double[] StandardError, double Score, bool Clears,
                                    double[] LeaveOneOut);

    /// <summary>The standings, and the winning arm; null when none clears.</summary>
    internal sealed record Verdict(Standing[] Standings, string? Winner);

    /// <summary>
    /// The rule, as decided: per fixture f, r_f = exp(mean over seeds of ln(min-ESS/s of the arm) − mean over seeds of
    /// ln(min-ESS/s of A0)), the seeds paired; the score is the geometric mean of r_f over the fixtures. An arm clears when
    /// its score is at least 1.5 and no r_f lies below 0.8; among the arms that clear, the highest score wins.
    /// </summary>
    internal static Verdict Rule(IReadOnlyList<Reading> pass, string[] fixtures, int seeds)
    {
        double LnE(string arm, string fixture, int k) => Math.Log(pass.Single(r => r.Arm == arm && r.Fixture == fixture && r.K == k).MinEssPerSecond);
        string a0 = Arms[0].Name;
        var standings = new List<Standing>();
        foreach (Arm arm in Arms.Skip(1))
        {
            int n = fixtures.Length;
            double[] logRatio = new double[n], spread = new double[n], error = new double[n];
            for (int i = 0; i < n; i++)
            {
                string f = fixtures[i];
                logRatio[i] = Enumerable.Range(0, seeds).Average(k => LnE(arm.Name, f, k)) - Enumerable.Range(0, seeds).Average(k => LnE(a0, f, k));
                double[] paired = [.. Enumerable.Range(0, seeds).Select(k => LnE(arm.Name, f, k) - LnE(a0, f, k))];
                spread[i] = StandardDeviation(paired);
                error[i] = spread[i] / Math.Sqrt(seeds);
            }
            double score = Math.Exp(logRatio.Average());
            bool clears = score >= Clears && logRatio.All(l => Math.Exp(l) >= Floor);
            double[] leaveOneOut = [.. Enumerable.Range(0, n).Select(i => Math.Exp(logRatio.Where((_, j) => j != i).Average()))];
            standings.Add(new Standing(arm.Name, logRatio, spread, error, score, clears, leaveOneOut));
        }
        return new Verdict([.. standings], standings.Where(s => s.Clears).MaxBy(s => s.Score)?.Arm);
    }

    private static double StandardDeviation(double[] values)
    {
        if (values.Length < 2) return double.NaN;
        double mean = values.Average();
        return Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / (values.Length - 1));
    }

    // The verdict first, then the context tables, as markdown.
    private static string Report(Verdict verdict, List<Reading> pass, List<(Reading First, Reading Again)> controls, Fixture[] fixtures,
                                 int scale, int seeds, TimeSpan elapsed)
    {
        var md = new StringBuilder();
        string Label(string arm) => $"{arm}, {Arms.Single(a => a.Name == arm).Label}";
        IEnumerable<Reading> Cell(string arm, string fixture) => pass.Where(r => r.Arm == arm && r.Fixture == fixture).OrderBy(r => r.K);
        void Row(params object[] cells) => md.Append("| ").AppendJoin(" | ", cells).Append(" |\n");
        void Head(params string[] cells)
        {
            Row(cells);
            Row([.. cells.Select(_ => "---")]);
        }

        md.Append($"Lengths 1/{scale} of S10's, {seeds} seed(s) per fixture, {Configuration}, {Environment.ProcessorCount} logical processors, " +
                  $"server GC {GCSettings.IsServerGC}; {elapsed.TotalMinutes:F1} minutes in all.\n\n");
        md.Append("### The verdict\n\n");
        Head(["arm", "score", .. fixtures.Select(f => $"r_f {f.Label}"), "least r_f", "clears"]);
        foreach (Standing s in verdict.Standings)
            Row([Label(s.Arm), $"{s.Score:F3}", .. s.LogRatio.Select(l => $"{Math.Exp(l):F3}"), $"{s.LogRatio.Min(Math.Exp):F3}", s.Clears ? "yes" : "no"]);
        md.Append(verdict.Winner is { } w
            ? $"\nThe winner is {Label(w)}, score {verdict.Standings.Single(s => s.Arm == w).Score:F3}.\n\n"
            : "\nNo arm clears: every option stays an option, and the defaults stay D01's.\n\n");

        md.Append("### Context for the verdict\n\nThe per-seed log ratios' standard deviation, and the standard error of ln r_f:\n\n");
        Head(["arm", .. fixtures.Select(f => $"{f.Label}: sd, SE")]);
        foreach (Standing s in verdict.Standings)
            Row([Label(s.Arm), .. s.Spread.Select((sd, i) => $"{sd:F3}, {s.StandardError[i]:F3}")]);
        md.Append("\nEach arm's score with one fixture left out:\n\n");
        Head(["arm", .. fixtures.Select(f => $"without {f.Label}")]);
        foreach (Standing s in verdict.Standings) Row([Label(s.Arm), .. s.LeaveOneOut.Select(x => $"{x:F3}")]);

        md.Append("\n### Min-ESS/s by arm and fixture\n\nThe geometric mean over seeds, each seed's value, the quantity that set each seed's minimum, and any constant quantity:\n\n");
        Head("arm", "fixture", "min-ESS/s", "per seed", "set by", "constant");
        foreach (Arm arm in Arms)
            foreach (Fixture f in fixtures)
            {
                Reading[] c = [.. Cell(arm.Name, f.Name)];
                double geo = Math.Exp(c.Average(r => Math.Log(r.MinEssPerSecond)));
                string setBy = string.Join(", ", c.GroupBy(r => r.MinQuantity).OrderByDescending(g => g.Count()).Select(g => $"{g.Key} ×{g.Count()}"));
                string constant = string.Join(", ", c.Where(r => r.Constant.Length > 0).Select(r => $"k{r.K}: {r.Constant}"));
                Row(Label(arm.Name), f.Label, $"{geo:F2}", string.Join(", ", c.Select(r => $"{r.MinEssPerSecond:F2}")), setBy, constant.Length > 0 ? constant : "none");
            }

        md.Append("\n### R̂ by quantity\n\nThe median over seeds, and the largest in brackets (∞ frozen, – constant):\n\n");
        Head(["arm", "fixture", .. QuantityNames]);
        foreach (Arm arm in Arms)
            foreach (Fixture f in fixtures)
            {
                Reading[] c = [.. Cell(arm.Name, f.Name)];
                Row([Label(arm.Name), f.Label, .. QuantityNames.Select((q, i) =>
                {
                    double[] r = [.. c.Select(x => x.Quantities[i].Diagnostic.RHat).Where(h => h.State != DiagnosticState.Constant)
                                      .Select(h => h.State == DiagnosticState.Estimated ? h.Value : double.PositiveInfinity).Order()];
                    return r.Length == 0 ? "–" : $"{Median(r):F2} [{r[^1]:F2}]";
                })]);
            }

        md.Append("\n### Acceptance by move kind\n\nPooled over seeds: accepted over attempted, and over admitted (neither null nor refused) in brackets:\n\n");
        BassMoveKind[] kinds = Enum.GetValues<BassMoveKind>();
        Head(["arm", "fixture", .. kinds.Select(k => k.ToString())]);
        foreach (Arm arm in Arms)
            foreach (Fixture f in fixtures)
            {
                Reading[] c = [.. Cell(arm.Name, f.Name)];
                Row([Label(arm.Name), f.Label, .. kinds.Select((k, i) =>
                {
                    long attempts = c.Sum(r => r.Kinds[i].Attempts), accepted = c.Sum(r => r.Kinds[i].Accepted);
                    long admitted = attempts - c.Sum(r => r.Kinds[i].Null + r.Kinds[i].Refused);
                    return attempts == 0 ? "" : $"{(double)accepted / attempts:F4} ({(admitted > 0 ? (double)accepted / admitted : double.NaN):F4})";
                })]);
            }

        md.Append("\n### σ's coverage and each chain's median M\n\n");
        Head("arm", "fixture", "σ's 95% interval covers the truth", "chains' median M, by seed");
        foreach (Arm arm in Arms)
            foreach (Fixture f in fixtures)
            {
                Reading[] c = [.. Cell(arm.Name, f.Name)];
                Row(Label(arm.Name), f.Label, $"{c.Count(r => r.SigmaCovers)} of {c.Length}", string.Join("; ", c.Select(r => string.Join(",", r.MedianM))));
            }

        md.Append("\n### The tempered arms' ladders\n\nOver seeds: each chain's settled rungs, its steep pairs (rejection above 0.9), the round trips per scan, the ideal rate, and observed over ideal:\n\n");
        Head("arm", "fixture", "rungs per chain [range]", "steep pairs per chain (mean)", "round trips/scan (mean)", "ideal (mean)", "observed/ideal (mean)", "largest pair rejection");
        foreach (Arm arm in Arms.Where(a => a.Tempered))
            foreach (Fixture f in fixtures)
            {
                Reading[] c = [.. Cell(arm.Name, f.Name)];
                int[] rungs = [.. c.SelectMany(r => r.Levels!)];
                Row(Label(arm.Name), f.Label, $"{rungs.Average():F1} [{rungs.Min()}, {rungs.Max()}]", $"{c.SelectMany(r => r.SteepPairs!).Average():F2}",
                    $"{c.Average(r => r.RoundTrips):G3}", $"{c.Average(r => r.IdealRoundTrips):G3}", $"{c.Average(r => r.RoundTrips / r.IdealRoundTrips):G3}",
                    $"{c.Max(r => r.MaxPairRejection):F3}");
            }

        md.Append("\n### Times\n\nThe whole fit's seconds (mean over seeds), the sampling phase's share, and CPU over wall time:\n\n");
        Head("arm", "fixture", "seconds", "sampling share", "CPU/wall");
        foreach (Arm arm in Arms)
            foreach (Fixture f in fixtures)
            {
                Reading[] c = [.. Cell(arm.Name, f.Name)];
                Row(Label(arm.Name), f.Label, $"{c.Average(r => r.Seconds):F2}", $"{c.Average(r => r.SamplingShare):F3}", $"{c.Average(r => r.CpuSeconds / r.Seconds):F3}");
            }
        md.Append($"\nThe pass's fits: {pass.Sum(r => r.Seconds) / 60:F1} minutes timed.\n\nThe timing control, A0 at each fixture's first seed run again at the end:\n\n");
        Head("fixture", "first (s)", "again (s)", "again/first", "draws identical");
        foreach ((Reading first, Reading again) in controls)
            Row(fixtures.Single(f => f.Name == first.Fixture).Label, $"{first.Seconds:F3}", $"{again.Seconds:F3}", $"{again.Seconds / first.Seconds:F3}",
                first.DrawsHash == again.DrawsHash ? "yes" : "no");

        Reading[] low = [.. pass.Where(r => r.MinEss < 10 || r.Constant.Length > 0 || r.MinQuantity.Contains("frozen"))];
        md.Append($"\nFits with a minimum ESS below 10, a constant quantity or a frozen M: {(low.Length == 0 ? "none" : "")}\n");
        foreach (Reading r in low) md.Append($"- {r.Arm}, {r.Fixture}, k = {r.K}: min ESS {r.MinEss:F2} ({r.MinQuantity}), constant {(r.Constant.Length > 0 ? r.Constant : "none")}\n");
        return md.ToString();
    }

    private static double Median(double[] sorted) => sorted.Length % 2 == 1 ? sorted[sorted.Length / 2] : 0.5 * (sorted[sorted.Length / 2 - 1] + sorted[sorted.Length / 2]);
}
