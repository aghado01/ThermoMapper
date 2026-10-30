using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// S5, the kept states: every k-th draw's FK partition, canonically labelled, and its spins. Each check runs the kept
/// states against what the run already accumulates, so the partitions are the ones the sums and counts describe.
/// </summary>
public sealed class SwKeptStatesTests
{
    private static readonly AccumulationSpec EveryDraw =
        new() { CoMembership = true, ExactClusterSizeHistogram = true, PartitionEvery = 1, SpinsEvery = 1 };

    [Fact]
    public void Off_ByDefault()
    {
        Accumulator acc = Run(AccumulationSpec.None, 0.6);

        Assert.Null(acc.KeptPartitions);
        Assert.Null(acc.KeptSpins);
        Assert.Equal(0, acc.PartitionEvery);
        Assert.Equal(0, acc.SpinsEvery);
    }

    [Theory]
    [InlineData(0.1)]
    [InlineData(0.6)]
    [InlineData(2.0)]
    public void Partitions_ReproduceTheAccumulatedStatistics(double temperature)
    {
        // Σ|c|², c_max, the exact size histogram and the co-membership counts, each rebuilt from the kept partitions.
        Accumulator acc = Run(EveryDraw, temperature);
        var graph = SwAccumulatorProbeTests.Random16();
        int n = acc.Spins.Length;
        Assert.Equal(acc.DrawCount * n, acc.KeptPartitions!.Length);

        var coMembership = new int[acc.CoMembershipCount!.Length];
        var exact = new long[n + 1];
        for (int d = 0; d < acc.DrawCount; d++)
        {
            int[] labels = acc.KeptPartitions[(d * n)..((d + 1) * n)];
            var size = new int[n];
            foreach (int label in labels) size[label]++;
            long sumSq = 0;
            int max = 0;
            foreach (int s in size.Where(s => s > 0))
            {
                sumSq += (long)s * s;
                max = Math.Max(max, s);
                exact[s]++;
            }
            Assert.Equal(acc.SumSqClusterSizesPerDraw[d], sumSq);
            Assert.Equal(acc.MaxClusterSizePerDraw[d], max);

            for (int i = 0; i < n; i++)
                for (int e = graph.RowPointers[i]; e < graph.RowPointers[i + 1]; e++)
                {
                    int j = graph.Targets[e];
                    if (j > i && labels[i] == labels[j]) coMembership[e]++;
                }
        }

        Assert.Equal(acc.CoMembershipCount, coMembership);
        Assert.Equal(acc.ExactClusterSizeHistogram, exact);
    }

    [Fact]
    public void Labels_AreEachClustersSmallestNode()
    {
        Accumulator acc = Run(EveryDraw, 0.6);
        int n = acc.Spins.Length;

        for (int d = 0; d < acc.DrawCount; d++)
        {
            int[] labels = acc.KeptPartitions![(d * n)..((d + 1) * n)];
            var smallest = new Dictionary<int, int>();
            for (int i = 0; i < n; i++)
                smallest.TryAdd(labels[i], i);
            foreach ((int label, int first) in smallest)
                Assert.Equal(label, first);
        }
    }

    [Fact]
    public void Spins_ColourEachClusterUniformly_AndGiveTheRecordsM()
    {
        Accumulator acc = Run(EveryDraw, 0.6);
        int n = acc.Spins.Length, q = acc.Q;

        for (int d = 0; d < acc.DrawCount; d++)
        {
            byte[] spins = acc.KeptSpins![(d * n)..((d + 1) * n)];
            int[] labels = acc.KeptPartitions![(d * n)..((d + 1) * n)];
            for (int i = 0; i < n; i++)
                Assert.Equal(spins[labels[i]], spins[i]);

            int nMax = spins.GroupBy(s => s).Max(g => g.Count());
            Assert.Equal(acc.MagPerDraw[d], ((double)nMax * q / n - 1.0) / (q - 1.0));
        }

        // The last kept spins are the run's final spins.
        Assert.Equal(acc.Spins.Select(s => (byte)s).ToArray(), acc.KeptSpins![^n..]);
    }

    [Fact]
    public void Thinning_KeepsEveryKthDraw_OfTheSameStream()
    {
        Accumulator all = Run(EveryDraw, 0.6);
        Accumulator thin = Run(EveryDraw with { PartitionEvery = 3, SpinsEvery = 7 }, 0.6);
        int n = all.Spins.Length;

        Assert.Equal(3, thin.PartitionEvery);
        Assert.Equal(7, thin.SpinsEvery);
        Assert.Equal((all.DrawCount + 2) / 3 * n, thin.KeptPartitions!.Length);
        Assert.Equal((all.DrawCount + 6) / 7 * n, thin.KeptSpins!.Length);
        for (int j = 0; j * n < thin.KeptPartitions.Length; j++)
            Assert.Equal(all.KeptPartitions![(3 * j * n)..((3 * j + 1) * n)], thin.KeptPartitions[(j * n)..((j + 1) * n)]);
        for (int j = 0; j * n < thin.KeptSpins.Length; j++)
            Assert.Equal(all.KeptSpins![(7 * j * n)..((7 * j + 1) * n)], thin.KeptSpins[(j * n)..((j + 1) * n)]);
    }

    [Fact]
    public void Serialization_RoundTrips_PresentAndAbsent()
    {
        Accumulator on = Run(EveryDraw with { PartitionEvery = 3, SpinsEvery = 2 }, 0.6);
        Accumulator round = SwDrawRecordTests.RoundTrip(on);

        Assert.Equal(on.KeptPartitions, round.KeptPartitions);
        Assert.Equal(on.KeptSpins, round.KeptSpins);
        Assert.Equal(3, round.PartitionEvery);
        Assert.Equal(2, round.SpinsEvery);

        Accumulator off = SwDrawRecordTests.RoundTrip(Run(AccumulationSpec.None, 0.6));
        Assert.Null(off.KeptPartitions);
        Assert.Null(off.KeptSpins);
        Assert.Equal(0, off.PartitionEvery);
        Assert.Equal(0, off.SpinsEvery);
    }

    [Fact]
    public void Restore_AppendsTheKeptStates_AndRefusesAMismatch()
    {
        var spec = new AccumulationSpec { PartitionEvery = 3, SpinsEvery = 2 };
        var original = Engine(spec, seed: 23);
        original.BurnIn(20);
        original.Run(61);
        Accumulator cp = original.GetCheckpoint();
        var resumed = Engine(spec, seed: 2);
        resumed.Restore(cp);
        original.Run(30);
        resumed.Run(30);

        Assert.Equal(original.GetCheckpoint().KeptPartitions, resumed.GetCheckpoint().KeptPartitions);
        Assert.Equal(original.GetCheckpoint().KeptSpins, resumed.GetCheckpoint().KeptSpins);

        Assert.Throws<InvalidOperationException>(() => Engine(AccumulationSpec.None, 1).Restore(cp));
        Assert.Throws<InvalidOperationException>(() => Engine(spec with { PartitionEvery = 4 }, 1).Restore(cp));
        Assert.Throws<InvalidOperationException>(() => Engine(spec, 1).Restore(cp with { KeptSpins = null }));
        Assert.Throws<InvalidOperationException>(() =>
            Engine(spec, 1).Restore(cp with { KeptPartitions = cp.KeptPartitions![..^1] }));
    }

    [Fact]
    public void Construction_Refuses_NegativeIntervals_AndSpinsPastAByte()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Engine(new AccumulationSpec { PartitionEvery = -1 }, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new SwendsenWang<NoCurrencies>(SwAccumulatorProbeTests.Random16(), 0.6, 300, new AccumulationSpec { SpinsEvery = 1 }, 1));
    }

    // ── Helpers ──────────────────────────────────────────────────────────────

    private static Accumulator Run(AccumulationSpec accumulation, double temperature)
        => SwRunner.Run(new SwRunSpec
        {
            Graph = SwAccumulatorProbeTests.Random16(), Temperature = temperature, Q = 20,
            Accumulation = accumulation, Seed = SwAccumulatorProbeTests.Seed, Budget = SwAccumulatorProbeTests.Budget,
        }).Accumulator;

    private static SwendsenWang<NoCurrencies> Engine(AccumulationSpec accumulation, int seed)
        => new(SwAccumulatorProbeTests.Random16(), 0.6, 20, accumulation, seed);
}
