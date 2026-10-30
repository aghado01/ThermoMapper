using System;
using System.IO;
using System.Linq;
using System.Numerics;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Observables;
using Xunit;
using PottsObs = Graphs.Models.Potts.Observables;

namespace VizCore.Tests;

/// <summary>
/// The exact cluster-size histogram (S7): <see cref="AccumulationSpec.ExactClusterSizeHistogram"/>, off by default,
/// counts the FK clusters at each exact size 1..N beside the log-binned bins. Checked against what is already kept —
/// Σ_s s·n_s = N × draws, Σ_s s²·n_s = Σ|c|²'s running sum, and its rebinning is the log-binned histogram — and reduced
/// by <see cref="ClusterSizeSums"/> to S₂, F_n and ⟨Σ c ln c⟩.
/// </summary>
public sealed class SwExactClusterSizeHistogramTests
{
    private static readonly AccumulationSpec Exact = new() { ExactClusterSizeHistogram = true };

    // ── Against what is already kept ─────────────────────────────────────────

    [Theory]
    [MemberData(nameof(SwAccumulatorProbeTests.CaseNames), MemberType = typeof(SwAccumulatorProbeTests))]
    public void Histogram_AgreesWithTheSumsAndTheLogBins(string name)
    {
        SwAccumulatorProbeTests.ProbeCase probe = SwAccumulatorProbeTests.Cases.Single(c => c.Name == name);
        Accumulator acc = SwAccumulatorProbeTests.Run(probe, probe.Accumulation with { ExactClusterSizeHistogram = true });
        long[] counts = acc.ExactClusterSizeHistogram!;
        int n = acc.Spins.Length;

        Assert.Equal(n + 1, counts.Length);
        Assert.Equal(0, counts[0]);

        long mass = 0, squares = 0;
        var rebinned = new long[acc.ClusterSizeHistogram.Length];
        for (int s = 1; s <= n; s++)
        {
            mass += s * counts[s];
            squares += (long)s * s * counts[s];
            rebinned[BitOperations.Log2((uint)s)] += counts[s];
        }

        Assert.Equal((long)n * acc.DrawCount, mass);
        Assert.Equal(acc.RunningSumSqClusterSizes, squares);
        Assert.Equal(acc.SumSqClusterSizesPerDraw.Sum(), squares);
        Assert.Equal(acc.ClusterSizeHistogram, rebinned);
    }

    [Fact]
    public void Histogram_IsAbsentByDefault()
        => Assert.Null(Run(AccumulationSpec.None).ExactClusterSizeHistogram);

    [Fact]
    public void Histogram_ResetsAfterBurnIn()
    {
        // Counted over the measured draws alone: the mass is N per measured draw, burn-in excluded.
        Accumulator acc = Run(Exact);
        long mass = acc.ExactClusterSizeHistogram!.Select((count, s) => s * count).Sum();

        Assert.Equal(16L * SwAccumulatorProbeTests.Budget.Cycles, mass);
    }

    // ── SPCX v11 and Restore ─────────────────────────────────────────────────

    [Fact]
    public void Histogram_SurvivesSerialization_PresentAndAbsent()
    {
        Accumulator on = Run(Exact);
        Assert.Equal(on.ExactClusterSizeHistogram, SwDrawRecordTests.RoundTrip(on).ExactClusterSizeHistogram);
        Assert.Null(SwDrawRecordTests.RoundTrip(Run(AccumulationSpec.None)).ExactClusterSizeHistogram);
    }

    [Fact]
    public void Restore_CarriesTheHistogram_AndRefusesAMismatch()
    {
        var original = Engine(Exact, seed: 23);
        original.BurnIn(20);
        original.Run(60);
        Accumulator cp = original.GetCheckpoint();
        var resumed = Engine(Exact, seed: 2);
        resumed.Restore(cp);
        original.Run(30);
        resumed.Run(30);

        Assert.Equal(original.GetCheckpoint().ExactClusterSizeHistogram, resumed.GetCheckpoint().ExactClusterSizeHistogram);

        Assert.Throws<InvalidOperationException>(() => Engine(AccumulationSpec.None, 1).Restore(cp));
        Assert.Throws<InvalidOperationException>(() => Engine(Exact, 1).Restore(cp with { ExactClusterSizeHistogram = null }));
        Assert.Throws<InvalidOperationException>(() =>
            Engine(Exact, 1).Restore(cp with { ExactClusterSizeHistogram = cp.ExactClusterSizeHistogram![..^1] }));
    }

    // ── The reductions ───────────────────────────────────────────────────────

    // Two draws on N = 5: clusters {3, 1, 1}, then {5}. So n₁ = 2, n₃ = 1, n₅ = 1.
    private static readonly long[] TwoDraws = [0, 2, 0, 1, 0, 1];

    [Fact]
    public void Mean_IsThePerDrawSumOverClusters()
        => Assert.Equal(2.0, ClusterSizeSums.Mean(TwoDraws, 2, _ => 1.0)); // (3 + 1) clusters over 2 draws

    [Fact]
    public void SecondMoment_IsS2()
        => Assert.InRange(Math.Abs(ClusterSizeSums.SecondMoment(TwoDraws, 2, 5) - (11.0 + 25.0) / 2 / 25), 0.0, 1e-15);

    [Theory]
    [InlineData(1.0, 3.0)] // draws: (2 + 0 + 0), (4); mean 3, over n = 1
    [InlineData(2.0, 1.0)] // draws: (1), (3); mean 2, over n = 2
    [InlineData(0.5, 8.0)] // draws: (2.5 + 0.5 + 0.5), (4.5); mean 4, over n = 1/2
    public void FlooredCount_IsFn(double floor, double expected)
        => Assert.InRange(Math.Abs(ClusterSizeSums.FlooredCount(TwoDraws, 2, floor) - expected), 0.0, 1e-15);

    [Fact]
    public void SizeLogSize_IsTheMeanOfSumCLnC()
        => Assert.InRange(
            Math.Abs(ClusterSizeSums.SizeLogSize(TwoDraws, 2) - (3 * Math.Log(3) + 5 * Math.Log(5)) / 2), 0.0, 1e-14);

    [Fact]
    public void Reductions_ReturnZeroWithoutDraws_AndRefuseBadArguments()
    {
        Assert.Equal(0.0, ClusterSizeSums.Mean(TwoDraws, 0, s => s));
        Assert.Throws<ArgumentOutOfRangeException>(() => ClusterSizeSums.FlooredCount(TwoDraws, 2, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => ClusterSizeSums.SecondMoment(TwoDraws, 2, 0));
    }

    [Fact]
    public void SecondMoment_OfARun_IsTheFkSusceptibilityOverN()
    {
        Accumulator acc = Run(Exact);
        double chi = PottsObs.Susceptibility.Reduce(acc.RunningSumSqClusterSizes, acc.DrawCount, acc.Spins.Length);

        Assert.InRange(
            Math.Abs(ClusterSizeSums.SecondMoment(acc.ExactClusterSizeHistogram, acc.DrawCount, acc.Spins.Length) - chi / 16),
            0.0, 1e-12);
    }

    // ── Helpers ──────────────────────────────────────────────────────────────

    private static Accumulator Run(AccumulationSpec accumulation)
        => SwRunner.Run(new SwRunSpec
        {
            Graph = SwAccumulatorProbeTests.Random16(), Temperature = 0.6, Q = 20,
            Accumulation = accumulation, Seed = SwAccumulatorProbeTests.Seed, Budget = SwAccumulatorProbeTests.Budget,
        }).Accumulator;

    private static SwendsenWang<NoCurrencies> Engine(AccumulationSpec accumulation, int seed)
        => new(SwAccumulatorProbeTests.Random16(), 0.6, 20, accumulation, seed);
}
