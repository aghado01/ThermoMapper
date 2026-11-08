using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// A BASS ladder whose rungs step at once (<see cref="ParallelTempering{TState}.MaxDegreeOfParallelism"/>) runs as the
/// in-order ladder does, bit for bit.
/// <list type="bullet">
/// <item><b>The stress run.</b> Ladders of 24 rungs, and a self-sizing one, on 16 workers against 1, through the burn-in's
/// rounds and resizes, the freeze and the traced phase: every scan's cold state, then every replica's state, density,
/// moves, blocks, refusals and trace, and the pairs' statistics. The cases cover what a rung's step reaches: the integrated,
/// linear and powered forms, the functional and categorical bases, the ridge's spectra, the g-form's τ block, and the
/// noise laws' blocks with interweaving, whose states carry their own working responses.</item>
/// <item><b>The pinned fixtures.</b> The tempered fits the pinned tier holds (the bit-identity baseline's self-sizing ladder,
/// the fixed-level S8 runs) and the flat histogram's pilots, the rungs in order (the default) and at once
/// (<see cref="BassRunOptions.RungWorkers"/> −1), at one bound.</item>
/// <item><b>Resume.</b> A ladder snapshotted mid-burn-in and a fit checkpointed mid-run continue bit for bit at another
/// degree: the degree is not restart state.</item>
/// <item><b>The split.</b> The run's bound goes to the chains first, and what they leave to each one's rungs when asked.</item>
/// </list>
/// The stress run at scale (<see cref="Stress_FortyEightRungsOnFortyEightWorkers_RunAsInOrder"/>) and the wall-time row
/// (<see cref="WallTime_OneAndFourChainsOfTwelveRungs"/>) run behind <c>BASS_QUALIFICATION</c> = 1.
/// </summary>
public sealed class BassRungParallelismTests
{
    private readonly ITestOutputHelper _out;
    public BassRungParallelismTests(ITestOutputHelper output) => _out = output;

    private const int AtOnce = 16;

    public static TheoryData<string> LadderCases => new()
    {
        "reference", "linear", "powered", "fixed g", "functional", "categorical", "ridge", "self-sizing", "t noise", "skewed noise, interweaving",
    };

    [Theory]
    [MemberData(nameof(LadderCases))]
    public void Ladder_RungsAtOnce_RunAsInOrder_BitForBit(string option)
    {
        (BassModel model, TemperingSpec spec) = LadderCase(option);
        List<string> inOrder = Run(model, spec, degree: 1), atOnce = Run(model, spec, AtOnce);
        AssertSameRecord(inOrder, atOnce);
        _out.WriteLine($"{option}: {inOrder.Count} entries agree; the last rung count {atOnce.Count(e => e.StartsWith("rung", StringComparison.Ordinal))}");
    }

    /// <summary>The stress run at scale, behind <c>BASS_QUALIFICATION</c> = 1: 48 rungs on 48 workers, more than the machine's
    /// cores, through 1,023 scans of burn-in and 2,000 traced, on the cases whose rungs reach the most shared state.</summary>
    [Theory]
    [InlineData("reference")]
    [InlineData("ridge")]
    [InlineData("functional")]
    [InlineData("skewed noise, interweaving")]
    public void Stress_FortyEightRungsOnFortyEightWorkers_RunAsInOrder(string option)
    {
        if (!BassQualification.Enabled) return;
        (BassModel model, TemperingSpec spec) = LadderCase(option, levels: 48);
        var clock = Stopwatch.StartNew();
        List<string> inOrder = Run(model, spec, degree: 1, burn: 1_023, traced: 2_000);
        double serial = clock.Elapsed.TotalSeconds;
        clock.Restart();
        List<string> atOnce = Run(model, spec, degree: 48, burn: 1_023, traced: 2_000);
        AssertSameRecord(inOrder, atOnce);
        _out.WriteLine($"{option}: {inOrder.Count} entries agree over 3,023 scans of 48 rungs; {serial:F1} s in order, {clock.Elapsed.TotalSeconds:F1} s at once");
    }

    [Theory]
    [InlineData("reference")]
    [InlineData("linear")]
    [InlineData("powered")]
    [InlineData("fixed ladder")]
    [InlineData("fixed g")]
    [InlineData("functional")]
    [InlineData("categorical")]
    public void FixedLevels_AsTheirPinsRunThem_TheSameAtOnceAndInOrder(string option)
    {
        (BassConfig config, BassData data) = BassTemperingTests.RunCase(option, samples: 120);
        int chains = config.Sampling.Chains;
        string inOrder = BassTemperingSizingTests.Pin(Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = chains * 6 }));
        string atOnce = BassTemperingSizingTests.Pin(Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = chains * 6, RungWorkers = -1 }));
        Assert.Equal(inOrder, atOnce);
    }

    [Fact]
    public void TheSelfSizingBaseline_AsItsPinRunsIt_TheSameAtOnceAndInOrder()
    {
        (BassConfig config, BassData data) = BassBitIdentityTests.Fixture(BassMixing.Friedman, "equal", tempered: true);
        int chains = config.Sampling.Chains;
        string inOrder = BassBitIdentityTests.Hash(Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = chains * 9 }));
        string atOnce = BassBitIdentityTests.Hash(Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = chains * 9, RungWorkers = -1 }));
        _out.WriteLine($"Friedman tempered: {inOrder} in order, {atOnce} at once");
        Assert.Equal(inOrder, atOnce);
    }

    [Fact]
    public void TheFlatHistogramsPilots_TheSameAtOnceAndInOrder()
    {
        BassModel model = BassFlatHistogramTests.FixtureA(heldTau: false);
        BassFlatHistogramPilot spec = BassFlatHistogramTests.FixtureAPilot;
        foreach (int seed in SeedTree.Derive(7702, 3).Prepend(7701))
        {
            var inOrder = new BassRunOptions { MaxDegreeOfParallelism = spec.Chains * 8 };
            var atOnce = new BassRunOptions { MaxDegreeOfParallelism = spec.Chains * 8, RungWorkers = -1 };
            (TemperingPilot a, BassTemperingReadout x) = BassFlatHistogramFit.Pilot(model, spec, seed, inOrder, new BassSweepCount(inOrder, 0));
            (TemperingPilot b, BassTemperingReadout y) = BassFlatHistogramFit.Pilot(model, spec, seed, atOnce, new BassSweepCount(atOnce, 0));
            Assert.Equal(a.InverseTemperatures, b.InverseTemperatures);
            Assert.Equal((a.Draws, a.ColdRung), (b.Draws, b.ColdRung));
            Assert.Equal(a.Inefficiencies, b.Inefficiencies);
            for (int r = 0; r < a.Rungs; r++) Assert.True(a.Levels(r).SequenceEqual(b.Levels(r)), $"seed {seed}, rung {r}'s levels differ");
            Assert.Equal(x.RoundTripsPerScan, y.RoundTripsPerScan);
            for (int c = 0; c < x.Ladders.Length; c++) Assert.Equal(x.PairRejection[c], y.PairRejection[c]);
        }
    }

    [Fact]
    public void ALadder_ResumedMidBurnIn_AtAnotherDegree_ContinuesBitForBit()
    {
        (BassModel model, TemperingSpec spec) = LadderCase("self-sizing");
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(5)).State;
        var whole = new BassLadder(model, spec, start.Terms, start.LogTau, seed: 6, start.StartLogRate);
        whole.Tempering.MaxDegreeOfParallelism = 8;
        int done = 0, round = 1;
        for (; round <= 16; round *= 2)
        {
            for (int t = 0; t < round; t++) whole.Step();
            done += round;
            whole.Resize();
        }
        BassLadder restored = BassLadder.Restore(model, spec, whole.Snapshot());
        Assert.Equal(1, restored.Tempering.MaxDegreeOfParallelism);
        restored.Tempering.MaxDegreeOfParallelism = 3;
        for (; round <= 64; round *= 2)
        {
            for (int t = 0; t < round; t++) Assert.Equal(Describe(whole.Step()), Describe(restored.Step()));
            whole.Resize();
            restored.Resize();
        }
        whole.EndBurnIn(traceRungs: true);
        restored.EndBurnIn(traceRungs: true);
        for (int t = 0; t < 100; t++) Assert.Equal(Describe(whole.Step()), Describe(restored.Step()));
        AssertSameRecord(Final(whole), Final(restored));
    }

    [Fact]
    public void AFit_CheckpointedInOrder_ResumesAtOnce_BitForBit()
    {
        (BassConfig whole, BassData data) = BassTemperingTests.RunCase("self-sizing", samples: 300);
        BassConfig half = whole with { Sampling = whole.Sampling with { Samples = 150 } };
        int chains = whole.Sampling.Chains;
        BassFit full = Bass.Fit(whole, data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
        BassFit first = Bass.Fit(half, data, new BassRunOptions { MaxDegreeOfParallelism = chains, Checkpoint = true });
        BassFit resumed = Bass.Resume(first.Checkpoint!, half, data, new BassRunOptions { MaxDegreeOfParallelism = chains * AtOnce, RungWorkers = -1 });
        BassRunTests.AssertSameFit(full, resumed);
        BassTemperingResumeTests.AssertSameTempering(full, resumed);
    }

    [Fact]
    public void TheBound_GoesToTheChainsFirst_AndWhatTheyLeaveToTheirRungs_WhenAsked()
    {
        Assert.Equal((4, 4), new BassRunOptions { MaxDegreeOfParallelism = 18 }.Split(4));
        Assert.Equal((2, 9), new BassRunOptions { MaxDegreeOfParallelism = 18 }.Split(2));
        Assert.Equal((3, 1), new BassRunOptions { MaxDegreeOfParallelism = 3 }.Split(4));
        Assert.Equal((1, 1), new BassRunOptions { MaxDegreeOfParallelism = 1 }.Split(1));
        (int units, int within) = new BassRunOptions().Split(1);
        Assert.Equal((1, Environment.ProcessorCount), (units, within));
        Assert.Throws<ArgumentOutOfRangeException>(() => new BassRunOptions().Split(0));
        Assert.Throws<InvalidOperationException>(() => new BassRunOptions { MaxDegreeOfParallelism = 0 }.Split(2));

        // The rungs step in order unless asked; asked, they take what the chains leave, or at most the number asked.
        Assert.Equal(1, new BassRunOptions { MaxDegreeOfParallelism = 18 }.RungWorkersFor(2));
        Assert.Equal(9, new BassRunOptions { MaxDegreeOfParallelism = 18, RungWorkers = -1 }.RungWorkersFor(2));
        Assert.Equal(4, new BassRunOptions { MaxDegreeOfParallelism = 18, RungWorkers = 4 }.RungWorkersFor(2));
        Assert.Equal(1, new BassRunOptions { MaxDegreeOfParallelism = 4, RungWorkers = 12 }.RungWorkersFor(4));
        Assert.Equal(Environment.ProcessorCount, new BassRunOptions { RungWorkers = -1 }.RungWorkersFor(1));
        Assert.Throws<InvalidOperationException>(() => new BassRunOptions { RungWorkers = 0 }.RungWorkersFor(2));
        (BassConfig config, BassData data) = BassTemperingTests.RunCase("reference", samples: 10);
        Assert.Contains("RungWorkers", Assert.Throws<ArgumentOutOfRangeException>(() => Bass.Fit(config, data, new BassRunOptions { RungWorkers = -2 })).Message);
    }

    // ---- The wall-time row ----

    /// <summary>
    /// Wall time of a Friedman fit (p = 10; n = 500, 1,000 sweeps of burn-in and 1,000 draws, and n = 4,000, 300 and 300) on a
    /// fixed ladder of 12 rungs with the reference rung, 1 chain and 4, the rungs in order and on 2 to 12 workers
    /// (<see cref="BassRunOptions.RungWorkers"/> −1 within each bound), three runs each in turn, the median kept, with each arm's CPU over its wall
    /// and its GC pauses, and the load the runs ran under. Rows go to <c>BASS_COST_OUT</c> as <c>rung-parallelism-wall.csv</c>.
    /// Indicative: the machine is shared.
    /// </summary>
    [Fact]
    public void WallTime_OneAndFourChainsOfTwelveRungs()
    {
        if (!BassQualification.Enabled) return;
        var rows = new BassCost.Rows("rung-parallelism", "wall");
        foreach ((string q, string v) in BassCost.Machine()) rows.Add("machine", q, v);
        var arms = new (string Name, int Chains, int Bound, int Cases, int Sweeps)[]
        {
            ("1 chain, rungs in order", 1, 1, 500, 1_000), ("1 chain, 2 rung workers", 1, 2, 500, 1_000), ("1 chain, 4 rung workers", 1, 4, 500, 1_000),
            ("1 chain, 6 rung workers", 1, 6, 500, 1_000), ("1 chain, 12 rung workers", 1, 12, 500, 1_000),
            ("4 chains, rungs in order", 4, 4, 500, 1_000), ("4 chains, 2 rung workers each", 4, 8, 500, 1_000),
            ("4 chains, 4 rung workers each", 4, 16, 500, 1_000),
            ("n 4,000: 1 chain, rungs in order", 1, 1, 4_000, 300), ("n 4,000: 1 chain, 4 rung workers", 1, 4, 4_000, 300),
            ("n 4,000: 1 chain, 12 rung workers", 1, 12, 4_000, 300), ("n 4,000: 4 chains, rungs in order", 4, 4, 4_000, 300),
            ("n 4,000: 4 chains, 4 rung workers each", 4, 16, 4_000, 300),
        };
        var data = new Dictionary<int, BassData> { [500] = WallData, [4_000] = BassQualificationRecoveryTests.Friedman(n: 4_000).Data };
        var seconds = arms.ToDictionary(a => a.Name, _ => new List<double>());
        var cpu = arms.ToDictionary(a => a.Name, _ => 0.0);
        var gc = arms.ToDictionary(a => a.Name, _ => 0.0);
        static BassRunOptions Options(int bound) => new() { MaxDegreeOfParallelism = bound, RungWorkers = -1 };
        foreach ((string name, int chains, int bound, int cases, _) in arms) Bass.Fit(WallConfig(chains, 50), data[cases], Options(bound));
        var load = new BassCost.Load();
        Process self = Process.GetCurrentProcess();
        for (int rep = 0; rep < 3; rep++)
            foreach ((string name, int chains, int bound, int cases, int sweeps) in arms)
            {
                self.Refresh();
                TimeSpan cpu0 = self.TotalProcessorTime, pause0 = GC.GetTotalPauseDuration();
                var clock = Stopwatch.StartNew();
                Bass.Fit(WallConfig(chains, sweeps), data[cases], Options(bound));
                seconds[name].Add(clock.Elapsed.TotalSeconds);
                self.Refresh();
                cpu[name] += (self.TotalProcessorTime - cpu0).TotalSeconds;
                gc[name] += (GC.GetTotalPauseDuration() - pause0).TotalSeconds;
            }
        (string[] atEnd, double cpuOverWall, double wall, double otherCores) = load.End();
        foreach ((string name, int chains, int bound, int cases, int sweeps) in arms)
        {
            double median = seconds[name].OrderBy(s => s).ElementAt(1), total = seconds[name].Sum();
            rows.Add(name, "chains", chains);
            rows.Add(name, "bound", bound);
            rows.Add(name, "cases", cases);
            rows.Add(name, "burn-in and draws, each", sweeps);
            rows.Add(name, "median seconds", median);
            rows.Add(name, "runs", string.Join(" ", seconds[name].Select(s => s.ToString("F3", CultureInfo.InvariantCulture))));
            rows.Add(name, "cpu over wall", cpu[name] / total);
            rows.Add(name, "gc pause fraction", gc[name] / total);
            _out.WriteLine($"{name}: median {median:F3} s ({string.Join(", ", seconds[name].Select(s => s.ToString("F3", CultureInfo.InvariantCulture)))}); " +
                           $"cpu/wall {cpu[name] / total:F2}, gc pause {gc[name] / total:P1}");
        }
        rows.Add("load", "other test hosts at start", string.Join(" ", load.AtStart));
        rows.Add("load", "other test hosts at end", string.Join(" ", atEnd));
        rows.Add("load", "this process's cpu over wall", cpuOverWall);
        rows.Add("load", "other processes' busy cores", otherCores);
        rows.Add("load", "seconds", wall);
        _out.WriteLine($"other processes kept {otherCores:F1} cores busy; this process {cpuOverWall:F1}; rows in {rows.Flush()}");
    }

    private static readonly BassData WallData = BassQualificationRecoveryTests.Friedman().Data;

    private static BassConfig WallConfig(int chains, int sweeps)
    {
        BassConfig config = BassMixing.Get(BassMixing.Friedman, new TemperingSpec { Levels = 12, ReferenceRung = true, Adaptive = false }).Config;
        return config with { Sampling = config.Sampling with { Chains = chains, BurnIn = sweeps, Samples = sweeps } };
    }

    // ---- The stress run ----

    private static readonly BassNoise Skewed = BassNoise.NormalWald(new BassNoisePrior(1.5, 0.5), new BassNoisePrior(0.0, 0.7));

    // The case's model and ladder: `levels` fixed rungs, or a self-sizing ladder within 32.
    private static (BassModel Model, TemperingSpec Spec) LadderCase(string option, int levels = 24)
    {
        var reference = new TemperingSpec { Levels = levels, HotEnd = 0.005, ReferenceRung = true };
        (BassConfig config, BassData data) = option switch
        {
            "t noise" => (BassStudentTTests.Config(BassNoise.StudentT(4)) with { Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 6 } },
                          BassQuantileTests.Data()),
            "skewed noise, interweaving" => (BassStudentTTests.Config(Skewed with { Interweave = true }) with
            {
                Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 6 },
                Sigma2 = Sigma2Prior.InverseGamma(3, 0.5),
            }, BassQuantileTests.Data()),
            "ridge" => BassTemperingTests.RunCase("reference"),
            "self-sizing" => BassTemperingTests.RunCase("self-sizing"),
            _ => BassTemperingTests.RunCase(option),
        };
        TemperingSpec spec = option switch
        {
            "linear" => reference with { Form = TemperedForm.Linear },
            "powered" => new TemperingSpec { Levels = levels, HotEnd = 0.005, Path = TemperingPath.Powered },
            "self-sizing" => new TemperingSpec { ReferenceRung = true, MaxLevels = 32 },
            _ => reference,
        };
        if (option == "ridge") config = config with { Coefficients = BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma(1.5)) };
        config = config with { Sampling = config.Sampling with { Tempering = spec } };
        return (BassModel.Create(config, data), spec);
    }

    // The ladder's run as BassRun drives it: burn-in rounds of 1, 2, 4, … scans, each followed by a resize or an adaptation,
    // then the frozen ladder traced. Every scan's cold state, then the final record.
    private static List<string> Run(BassModel model, TemperingSpec spec, int degree, int burn = 255, int traced = 200)
    {
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(5)).State;
        var ladder = new BassLadder(model, spec, start.Terms, start.LogTau, seed: 6, start.StartLogRate);
        ladder.Tempering.MaxDegreeOfParallelism = degree;
        var record = new List<string>();
        void Scan() => record.Add($"{Describe(ladder.Step())} swaps {ladder.Tempering.SwapAccepts}");
        int done = 0;
        for (int round = 1; done + round <= burn; round *= 2)
        {
            for (int t = 0; t < round; t++) Scan();
            done += round;
            if (ladder.SizesItself && done + 2L * round <= burn) ladder.Resize();
            else ladder.Adapt();
            Assert.Equal(degree, ladder.Tempering.MaxDegreeOfParallelism);
        }
        ladder.EndBurnIn(traceRungs: true);
        for (int t = 0; t < traced; t++) Scan();
        record.AddRange(Final(ladder));
        return record;
    }

    private static List<string> Final(BassLadder ladder)
    {
        ParallelTempering<BassState> pt = ladder.Tempering;
        var record = new List<string>();
        for (int i = 0; i < pt.Replicas.Count; i++)
        {
            ReversibleJumpChain<BassState> r = pt.Replicas[i];
            record.Add($"rung {i}: {Describe(r.Current)} {Bits(r.CurrentLogDensity)} {Bits(r.Temperature.InverseTemperature)} " +
                       string.Join(" ", r.MoveStats().Select(m => $"{m.Key} {m.Attempts}/{m.Accepted}")) + " | " +
                       string.Join(" ", pt.BlockStats(i).Select(m => $"{m.Key} {m.Attempts}/{m.Accepted}")) + " | " +
                       string.Join(" ", ladder.Structures[i].Refusals.Entries.Select(e => $"{e.Category} {e.Count}")) +
                       $" empty births {ladder.Structures[i].EmptyBirths}");
            if (pt.RungTraces is { } traces)
                record.Add($"trace {i}: " + string.Join(" ", traces[i].Draws(RungFunctional.Slope).Concat(traces[i].Draws(RungFunctional.LogRatio)).Select(Bits)));
        }
        foreach (SwapPairStat p in pt.PairStats) record.Add($"pair {p.Attempts}/{p.Accepts} {Bits(p.MeanRejection)}");
        record.Add($"scans {pt.Scans}, swaps {pt.SwapAttempts}/{pt.SwapAccepts}, round trips {pt.RoundTrips}, restarts {pt.AnnealedRestarts}, " +
                   $"reference attempts {ladder.Reference?.Attempts ?? 0}");
        return record;
    }

    // A state by its bits: its terms, u = ln τ, v, and its fit's sums of squares.
    private static string Describe(BassState s)
        => $"{string.Join(";", s.Terms)} {Bits(s.LogTau)} {Bits(s.LogRate)} {Bits(s.Fit.ResidualSquares)} {Bits(s.Fit.FittedSquares)}";

    private static string Bits(double x) => BitConverter.DoubleToInt64Bits(x).ToString("x16");

    private static void AssertSameRecord(List<string> expected, List<string> actual)
    {
        Assert.Equal(expected.Count, actual.Count);
        for (int i = 0; i < expected.Count; i++)
            if (expected[i] != actual[i]) Assert.Fail($"Entry {i} differs: in order '{expected[i]}', at once '{actual[i]}'.");
    }
}
