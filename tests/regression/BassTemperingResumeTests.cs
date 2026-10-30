using System;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Evidence;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// A tempered run's restart state (stage 4's S8): resuming a checkpoint continues the run bit for bit, its ladders, their
/// readout and its evidence included; a ladder restored from its snapshot continues bit for bit before burn-in ends and
/// after; and the result does not depend on the parallelism bound.
/// </summary>
public sealed class BassTemperingResumeTests
{
    private readonly ITestOutputHelper _out;
    public BassTemperingResumeTests(ITestOutputHelper output) => _out = output;

    public static TheoryData<string> Options => new()
    {
        "reference", "linear", "powered", "fixed ladder", "fixed g", "functional", "categorical", "flat + jeffreys", "thinned", "two rungs",
        "three rungs, reference", "self-sizing", "self-sizing, linear", "self-sizing, powered", "self-sizing, fixed g", "self-sizing, functional",
        "self-sizing, categorical", "self-sizing, flat",
    };

    [Theory]
    [MemberData(nameof(Options))]
    public void Resume_IsBitIdentical(string option)
    {
        (BassConfig whole, BassData data) = BassTemperingTests.RunCase(option, samples: 400);
        BassConfig half = whole with { Sampling = whole.Sampling with { Samples = 200 } };
        BassFit full = Bass.Fit(whole, data);
        BassFit first = Bass.Fit(half, data, new BassRunOptions { Checkpoint = true });
        BassFit resumed = Bass.Resume(first.Checkpoint!, half, data);
        Assert.Equal(400, resumed.Diagnostics.SamplesUsed);
        BassRunTests.AssertSameFit(full, resumed);
        AssertSameTempering(full, resumed);
        Assert.Equal(option is "reference" or "linear" or "fixed ladder" or "thinned" or "three rungs, reference"
                     || option.StartsWith("self-sizing", StringComparison.Ordinal) && option is not ("self-sizing, powered" or "self-sizing, flat"),
                     full.Evidence is not null);
        _out.WriteLine($"{option}: {full.Terms.Count} terms, round trips {full.Diagnostics.Tempering!.RoundTripsPerScan:F4}/scan, " +
                       $"evidence {full.Evidence?.LogEvidence.ToString("F4") ?? "null"}");
    }

    [Fact]
    public void Resume_Twice_AndRefusesAnotherLadder()
    {
        (BassConfig whole, BassData data) = BassTemperingTests.RunCase("linear", samples: 300);
        BassConfig third = whole with { Sampling = whole.Sampling with { Samples = 100 } };
        BassFit first = Bass.Fit(third, data, new BassRunOptions { Checkpoint = true });
        BassFit second = Bass.Resume(first.Checkpoint!, third, data, new BassRunOptions { Checkpoint = true });
        BassFit full = Bass.Fit(whole, data);
        BassFit resumed = Bass.Resume(second.Checkpoint!, third, data);
        BassRunTests.AssertSameFit(full, resumed);
        AssertSameTempering(full, resumed);
        Assert.NotNull(resumed.Evidence!.PathSampling);

        // A tempered checkpoint does not resume an untempered run, nor one on another ladder.
        Assert.Contains("fingerprint", Assert.Throws<ArgumentException>(() =>
            Bass.Resume(first.Checkpoint!, third with { Sampling = third.Sampling with { Tempering = null } }, data)).Message);
        Assert.Throws<ArgumentException>(() =>
            Bass.Resume(first.Checkpoint!, third with { Sampling = third.Sampling with { Tempering = third.Sampling.Tempering! with { Levels = 6 } } }, data));
    }

    [Fact]
    public void Fit_DoesNotDependOnTheParallelismBound()
    {
        (BassConfig config, BassData data) = BassTemperingTests.RunCase("reference", samples: 150);
        BassFit reference = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
        foreach (int degree in new[] { 2, 3, -1 })
        {
            BassFit other = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = degree });
            BassRunTests.AssertSameFit(reference, other);
            AssertSameTempering(reference, other);
        }
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Ladder_RestoredFromItsSnapshot_ContinuesBitForBit(bool reference)
    {
        BassModel model = BassModel.Create(BassTemperingTests.Config(null, maxBasis: 8), BassSweepTests.Data(large: true));
        var spec = new TemperingSpec { Levels = 5, HotEnd = 0.01, ReferenceRung = reference };
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(5)).State;
        var ladder = new BassLadder(model, spec, start.Terms, start.LogTau, seed: 6);
        for (int round = 1; round <= 64; round *= 2)
        {
            for (int t = 0; t < round; t++) ladder.Step();
            ladder.Adapt();
        }

        // Mid burn-in: the restored ladder adapts as the original does.
        BassLadder restored = BassLadder.Restore(model, spec, ladder.Snapshot());
        Assert.False(restored.Frozen);
        for (int t = 0; t < 128; t++)
        {
            ladder.Step();
            restored.Step();
        }
        Assert.Equal(ladder.Adapt(), restored.Adapt());
        AssertSameLadder(ladder, restored);

        // After burn-in, with the rungs traced: the swaps, the hand-overs and the traces continue as one.
        for (int t = 0; t < 100; t++) ladder.Step();
        ladder.EndBurnIn(traceRungs: true);
        for (int t = 0; t < 250; t++) ladder.Step();
        restored = BassLadder.Restore(model, spec, ladder.Snapshot());
        Assert.True(restored.Frozen);
        for (int t = 0; t < 400; t++)
        {
            ladder.Step();
            restored.Step();
        }
        AssertSameLadder(ladder, restored);
        Assert.True(ladder.Tempering.SwapAccepts > 200);
    }

    [Theory]
    [InlineData(false, false)]
    [InlineData(true, false)]
    [InlineData(true, true)]
    public void SelfSizingLadder_RestoredAcrossItsResizes_ContinuesBitForBit(bool reference, bool fixedG)
    {
        BassConfig config = BassTemperingTests.Config(null, maxBasis: 8);
        if (fixedG) config = config with { Coefficients = BassCoefficientPrior.G(TauPrior.Held()) };
        BassModel model = BassModel.Create(config, BassSweepTests.Data(large: true));
        var spec = new TemperingSpec { ReferenceRung = reference, MaxLevels = 20 };
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(5)).State;
        var ladder = new BassLadder(model, spec, start.Terms, start.LogTau, seed: 6);
        var counts = new System.Collections.Generic.List<int> { ladder.Tempering.Replicas.Count };
        for (int round = 1; round <= 32; round *= 2)
        {
            for (int t = 0; t < round; t++) ladder.Step();
            ladder.Resize();
            counts.Add(ladder.Tempering.Replicas.Count);
        }

        // Before a sizing round, and after it: the restored ladder sizes, places and equalizes as the original does.
        BassLadder restored = BassLadder.Restore(model, spec, ladder.Snapshot());
        void Both(int scans)
        {
            for (int t = 0; t < scans; t++)
            {
                ladder.Step();
                restored.Step();
            }
        }
        Both(64);
        Assert.Equal(ladder.Resize(), restored.Resize());
        counts.Add(ladder.Tempering.Replicas.Count);
        AssertSameLadder(ladder, restored);
        restored = BassLadder.Restore(model, spec, ladder.Snapshot());
        Both(128);
        Assert.Equal(ladder.Adapt(), restored.Adapt());
        AssertSameLadder(ladder, restored);

        // After burn-in, with the rungs traced.
        Both(50);
        ladder.EndBurnIn(traceRungs: true);
        restored.EndBurnIn(traceRungs: true);
        Both(100);
        restored = BassLadder.Restore(model, spec, ladder.Snapshot());
        Both(300);
        AssertSameLadder(ladder, restored);
        _out.WriteLine($"reference {reference}, fixed g {fixedG}: counts [{string.Join(", ", counts)}]");
        Assert.True(counts.Distinct().Count() > 1, "the ladder never resized");
    }

    private static void AssertSameLadder(BassLadder a, BassLadder b)
    {
        ParallelTempering<BassState> x = a.Tempering, y = b.Tempering;
        Assert.Equal(x.InverseTemperatures, y.InverseTemperatures);
        Assert.Equal((x.Scans, x.SwapAttempts, x.SwapAccepts, x.RoundTrips, x.AnnealedRestarts),
                     (y.Scans, y.SwapAttempts, y.SwapAccepts, y.RoundTrips, y.AnnealedRestarts));
        Assert.Equal(x.PairStats, y.PairStats);
        Assert.Equal(x.Snapshot().GrowthRng, y.Snapshot().GrowthRng);
        for (int i = 0; i < x.Replicas.Count; i++)
        {
            BassState s = x.Replicas[i].Current, t = y.Replicas[i].Current;
            Assert.Equal(s.Terms, t.Terms);
            Assert.Equal(s.LogTau, t.LogTau);
            Assert.Equal(s.Fit, t.Fit);
            Assert.Equal(x.Replicas[i].MoveStats(), y.Replicas[i].MoveStats());
            Assert.Equal(x.BlockStats(i), y.BlockStats(i));
            BassTemperingTests.AssertSameDesign(a.Structures[i].DesignAt(s), b.Structures[i].DesignAt(t));
            Assert.Equal(a.Structures[i].Refusals.Entries, b.Structures[i].Refusals.Entries);
        }
        Assert.Equal(a.Reference?.Attempts, b.Reference?.Attempts);
        Assert.Equal(x.RungTraces is null, y.RungTraces is null);
        if (x.RungTraces is null) return;
        for (int i = 0; i < x.RungTraces.Count; i++)
            for (int f = 0; f < 2; f++)
                Assert.Equal(x.RungTraces[i].Draws(f), y.RungTraces![i].Draws(f));
    }

    /// <summary>The tempering readouts and the evidences of two fits agree bit for bit.</summary>
    internal static void AssertSameTempering(BassFit a, BassFit b)
    {
        BassTemperingReadout? x = a.Diagnostics.Tempering, y = b.Diagnostics.Tempering;
        Assert.Equal(x is null, y is null);
        if (x is not null)
        {
            Assert.Equal((x.Path, x.Form, x.Schedule, x.ReferenceRung), (y!.Path, y.Form, y.Schedule, y.ReferenceRung));
            Assert.Equal((x.RoundTripsPerScan, x.AnnealedRestartsPerScan, x.Barrier, x.MaxPairRejection, x.IdealRoundTripRate, x.AtMaxLevels, x.ReferenceAttemptsPerDraw),
                         (y.RoundTripsPerScan, y.AnnealedRestartsPerScan, y.Barrier, y.MaxPairRejection, y.IdealRoundTripRate, y.AtMaxLevels, y.ReferenceAttemptsPerDraw));
            Assert.Equal(x.Levels, y.Levels);
            Assert.Equal(x.Ladders.Length, y.Ladders.Length);
            for (int c = 0; c < x.Ladders.Length; c++)
            {
                Assert.Equal(x.Ladders[c], y.Ladders[c]);
                Assert.Equal(x.PairRejection[c], y.PairRejection[c]);
                Assert.Equal(x.PairAcceptance[c], y.PairAcceptance[c]);
            }
        }
        Assert.Equal(a.Evidence is null, b.Evidence is null);
        if (a.Evidence is not { } e) return;
        AssertSameEstimate(e.SteppingStone, b.Evidence!.SteppingStone);
        Assert.Equal(e.PathSampling is null, b.Evidence.PathSampling is null);
        if (e.PathSampling is { } ps) AssertSameEstimate(ps, b.Evidence.PathSampling!);
        Assert.Equal(e.Priors, b.Evidence.Priors);
        Assert.Equal((e.MaxPairRejection, e.AtMaxLevels), (b.Evidence.MaxPairRejection, b.Evidence.AtMaxLevels));
    }

    private static void AssertSameEstimate(BassEvidenceEstimate a, BassEvidenceEstimate b)
    {
        Assert.Equal((a.Method, a.LogEvidence, a.MonteCarloError, a.WithinChainError, a.BetweenChainError, a.DiscretizationError),
                     (b.Method, b.LogEvidence, b.MonteCarloError, b.WithinChainError, b.BetweenChainError, b.DiscretizationError));
        Assert.Equal(a.Chains.Count, b.Chains.Count);
        for (int c = 0; c < a.Chains.Count; c++)
        {
            LadderEstimate x = a.Chains[c], y = b.Chains[c];
            Assert.Equal((x.Method, x.LogRatio, x.MonteCarloError, x.DiscretizationError), (y.Method, y.LogRatio, y.MonteCarloError, y.DiscretizationError));
            Assert.Equal(x.Rungs, y.Rungs);
        }
    }
}
