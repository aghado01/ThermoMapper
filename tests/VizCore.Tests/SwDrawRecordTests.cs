using System;
using System.IO;
using System.Linq;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The per-draw record (S1 with S3): E, m, Σ|c|² and c_max for each measured draw, kept on the
/// <see cref="Accumulator"/>. The running sums are its views, so each series summed in order must equal its sum bit for
/// bit; the record starts after burn-in, travels through SPCX v9, and a restored run appends to it.
/// </summary>
public sealed class SwDrawRecordTests
{
    // ── The series against the sums ──────────────────────────────────────────

    [Theory]
    [MemberData(nameof(SwAccumulatorProbeTests.CaseNames), MemberType = typeof(SwAccumulatorProbeTests))]
    public void Series_SummedInOrder_EqualTheRunningSums(string name)
    {
        AssertSeriesAreTheSums(SwAccumulatorProbeTests.Run(SwAccumulatorProbeTests.Cases.Single(c => c.Name == name)));
    }

    /// <summary>Each series has one entry per draw and folds, in draw order, to its running sum bit for bit.</summary>
    internal static void AssertSeriesAreTheSums(Accumulator acc)
    {
        Assert.Equal(acc.DrawCount, acc.EnergyPerDraw.Length);
        Assert.Equal(acc.DrawCount, acc.MagPerDraw.Length);
        Assert.Equal(acc.DrawCount, acc.SumSqClusterSizesPerDraw.Length);
        Assert.Equal(acc.DrawCount, acc.MaxClusterSizePerDraw.Length);

        double e = 0.0, e2 = 0.0, m = 0.0, m2 = 0.0, sumSq = 0.0, sumSqExcl = 0.0;
        for (int t = 0; t < acc.DrawCount; t++)
        {
            double energy = acc.EnergyPerDraw[t], mag = acc.MagPerDraw[t];
            long maxSize = acc.MaxClusterSizePerDraw[t];
            e += energy;
            e2 += energy * energy;
            m += mag;
            m2 += mag * mag;
            sumSq += acc.SumSqClusterSizesPerDraw[t];
            sumSqExcl += acc.SumSqClusterSizesPerDraw[t] - maxSize * maxSize;
        }

        AssertBits(acc.RunningSumEnergy, e);
        AssertBits(acc.RunningSumEnergySq, e2);
        AssertBits(acc.RunningSumMag, m);
        AssertBits(acc.RunningSumMagSq, m2);
        AssertBits(acc.RunningSumSqClusterSizes, sumSq);
        AssertBits(acc.RunningSumSqClusterSizesExcl, sumSqExcl);
    }

    // ── Burn-in and what each value describes ────────────────────────────────

    [Fact]
    public void Record_StartsAfterBurnIn_AsTheTailOfAnUnburnedRun()
    {
        Accumulator burned = Run(SwAccumulatorProbeTests.Random16(), 0.6, 20, new RunBudget(50, 400));
        Accumulator unburned = Run(SwAccumulatorProbeTests.Random16(), 0.6, 20, new RunBudget(0, 450));

        Assert.Equal(unburned.EnergyPerDraw[50..], burned.EnergyPerDraw);
        Assert.Equal(unburned.MagPerDraw[50..], burned.MagPerDraw);
        Assert.Equal(unburned.SumSqClusterSizesPerDraw[50..], burned.SumSqClusterSizesPerDraw);
        Assert.Equal(unburned.MaxClusterSizePerDraw[50..], burned.MaxClusterSizePerDraw);
    }

    [Fact]
    public void Magnetization_DescribesTheSpinsLeavingTheDraw()
    {
        const int q = 20;
        Accumulator acc = Run(SwAccumulatorProbeTests.Random16(), 0.6, q, new RunBudget(50, 400));

        int[] buckets = new int[q];
        foreach (int s in acc.Spins) buckets[s]++;
        double m = ((double)buckets.Max() * q / acc.Spins.Length - 1.0) / (q - 1.0);

        AssertBits(m, acc.MagPerDraw[^1]);
    }

    [Fact]
    public void Energy_DescribesTheSpinsEnteringTheDraw()
    {
        // One stream, one draw longer: the shorter run's final spins enter the longer run's last draw.
        CsrGraph graph = SwAccumulatorProbeTests.Random16();
        Accumulator shorter = Run(graph, 0.6, 20, new RunBudget(50, 400));
        Accumulator longer = Run(graph, 0.6, 20, new RunBudget(50, 401));

        Assert.Equal(shorter.EnergyPerDraw, longer.EnergyPerDraw[..400]);
        AssertBits(Energy(graph, shorter.Spins), longer.EnergyPerDraw[400]);
    }

    // ── SPCX v9 ──────────────────────────────────────────────────────────────

    [Fact]
    public void Record_SurvivesSerialization()
    {
        Accumulator acc = SwAccumulatorProbeTests.Run(
            SwAccumulatorProbeTests.Cases.Single(c => c.Name == "random16/T0.6/full"));

        Accumulator round = RoundTrip(acc);

        Assert.Equal(acc.EnergyPerDraw, round.EnergyPerDraw);
        Assert.Equal(acc.MagPerDraw, round.MagPerDraw);
        Assert.Equal(acc.SumSqClusterSizesPerDraw, round.SumSqClusterSizesPerDraw);
        Assert.Equal(acc.MaxClusterSizePerDraw, round.MaxClusterSizePerDraw);
        Assert.Equal(SwAccumulatorProbeTests.Fingerprint(acc), SwAccumulatorProbeTests.Fingerprint(round));
    }

    [Fact]
    public void EmptyRecord_SurvivesSerialization()
    {
        // A hand-built frame carries sums alone; its record is empty, never null.
        var frame = new Accumulator
        {
            Temperature = 1.0, Q = 4, DrawCount = 3,
            Spins = new int[2], ClusterSizeHistogram = new long[2],
            RngState0 = 1, RngState1 = 2, RngState2 = 3, RngState3 = 4,
            RunningSumSqClusterSizes = 0, RunningSumSqClusterSizesExcl = 0,
            RunningSumEnergy = 0, RunningSumEnergySq = 0,
            RunningSumMag = 0, RunningSumMagSq = 0,
        };

        Accumulator round = RoundTrip(frame);

        Assert.Empty(round.EnergyPerDraw);
        Assert.Empty(round.MagPerDraw);
        Assert.Empty(round.SumSqClusterSizesPerDraw);
        Assert.Empty(round.MaxClusterSizePerDraw);
    }

    [Theory]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    public void Serializer_RefusesAnOlderFile_NamingItsVersion(int version)
    {
        using var ms = new MemoryStream();
        using (var w = new BinaryWriter(ms, System.Text.Encoding.UTF8, leaveOpen: true))
        {
            w.Write(0x58435053u); // "SPCX"
            w.Write(version);
        }
        ms.Position = 0;

        var ex = Assert.Throws<InvalidDataException>(() => AccumulatorSerializer.Instance.ReadFrom(ms));
        Assert.Contains($"version {version}", ex.Message);
    }

    [Fact]
    public void Serializer_RefusesUnequalSeriesLengths()
    {
        Accumulator acc = Run(SwAccumulatorProbeTests.Triangle(), 1.0, 4, new RunBudget(0, 10));

        using var ms = new MemoryStream();
        Assert.Throws<ArgumentException>(() =>
            AccumulatorSerializer.Instance.WriteTo(acc with { MagPerDraw = acc.MagPerDraw[..9] }, ms));
    }

    // ── The engine: read surface, growth, Restore ────────────────────────────

    [Fact]
    public void ReadSurface_IsTheRecord()
    {
        var engine = Engine(seed: 5);
        engine.BurnIn(20);
        engine.Run(60);
        Accumulator cp = engine.GetCheckpoint();

        Assert.Equal(cp.EnergyPerDraw, engine.EnergyPerDraw.ToArray());
        Assert.Equal(cp.MagPerDraw, engine.MagPerDraw.ToArray());
        Assert.Equal(cp.SumSqClusterSizesPerDraw, engine.SumSqClusterSizesPerDraw.ToArray());
        Assert.Equal(cp.MaxClusterSizePerDraw, engine.MaxClusterSizePerDraw.ToArray());
    }

    [Fact]
    public void BareDraws_GrowTheRecordPastItsCapacity()
    {
        // Run reserves its count up front; bare Draws double the capacity instead. Same stream, same record.
        var reserved = Engine(seed: 9);
        reserved.Run(37);
        var grown = Engine(seed: 9);
        for (int i = 0; i < 37; i++) grown.Draw();

        Accumulator a = reserved.GetCheckpoint(), b = grown.GetCheckpoint();
        Assert.Equal(a.EnergyPerDraw, b.EnergyPerDraw);
        Assert.Equal(a.MagPerDraw, b.MagPerDraw);
        Assert.Equal(a.SumSqClusterSizesPerDraw, b.SumSqClusterSizesPerDraw);
        Assert.Equal(a.MaxClusterSizePerDraw, b.MaxClusterSizePerDraw);
        AssertSeriesAreTheSums(b);
    }

    [Fact]
    public void Restore_AppendsToTheRecord()
    {
        var original = Engine(seed: 11);
        original.BurnIn(30);
        original.Run(120);
        var resumed = Engine(seed: 99);
        resumed.Restore(original.GetCheckpoint());

        original.Run(80);
        resumed.Run(80);

        Accumulator a = original.GetCheckpoint(), b = resumed.GetCheckpoint();
        Assert.Equal(200, b.EnergyPerDraw.Length);
        Assert.Equal(a.EnergyPerDraw, b.EnergyPerDraw);
        Assert.Equal(a.MagPerDraw, b.MagPerDraw);
        Assert.Equal(a.SumSqClusterSizesPerDraw, b.SumSqClusterSizesPerDraw);
        Assert.Equal(a.MaxClusterSizePerDraw, b.MaxClusterSizePerDraw);
        Assert.Equal(SwAccumulatorProbeTests.Fingerprint(a), SwAccumulatorProbeTests.Fingerprint(b));
        AssertSeriesAreTheSums(b);
    }

    [Fact]
    public void Restore_RefusesARecordThatMissesTheDrawCount()
    {
        var engine = Engine(seed: 3);
        engine.Run(40);
        Accumulator cp = engine.GetCheckpoint();

        Assert.Throws<InvalidOperationException>(() =>
            Engine(seed: 3).Restore(cp with { MaxClusterSizePerDraw = cp.MaxClusterSizePerDraw[..39] }));
        Assert.Throws<InvalidOperationException>(() =>
            Engine(seed: 3).Restore(cp with { DrawCount = 41 }));
    }

    // ── Helpers ──────────────────────────────────────────────────────────────

    private static Accumulator Run(CsrGraph graph, double temperature, int q, RunBudget budget)
        => SwRunner.Run(new SwRunSpec
        {
            Graph = graph, Temperature = temperature, Q = q,
            Seed = SwAccumulatorProbeTests.Seed, Budget = budget,
        }).Accumulator;

    private static SwendsenWang<NoCurrencies> Engine(int seed)
        => new(SwAccumulatorProbeTests.Random16(), 0.6, 20, AccumulationSpec.None, seed);

    /// <summary>H = −Σ J δ, walked in the cycle's own edge order so the sum is bit-comparable.</summary>
    private static double Energy(CsrGraph graph, int[] spins)
    {
        double h = 0.0;
        for (int i = 0; i < graph.NodeCount; i++)
            for (int e = graph.RowPointers[i]; e < graph.RowPointers[i + 1]; e++)
            {
                int j = graph.Targets[e];
                if (j > i && spins[i] == spins[j]) h -= graph.Weights[e];
            }
        return h;
    }

    internal static Accumulator RoundTrip(Accumulator acc)
    {
        using var ms = new MemoryStream();
        AccumulatorSerializer.Instance.WriteTo(acc, ms);
        ms.Position = 0;
        return AccumulatorSerializer.Instance.ReadFrom(ms);
    }

    internal static void AssertBits(double expected, double actual)
        => Assert.Equal(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(actual));
}
