using System;
using System.IO;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The burn-in trace (S4): <see cref="AccumulationSpec.BurnInTrace"/>, off by default, keeps the per-draw record of the
/// draws the reset after burn-in discards. Its length is the burn-in count, it is the head of the same stream run
/// without burn-in, the measured record is unchanged by turning it on, and it travels through SPCX v10.
/// </summary>
public sealed class SwBurnInTraceTests
{
    private static readonly AccumulationSpec Trace = new() { BurnInTrace = true };
    private static readonly RunBudget Budget = new(BurnIn: 50, Cycles: 400);

    [Fact]
    public void Trace_IsAbsentByDefault()
    {
        Accumulator acc = Run(AccumulationSpec.None, Budget);

        Assert.Null(acc.BurnInEnergyPerDraw);
        Assert.Null(acc.BurnInMagPerDraw);
        Assert.Null(acc.BurnInSumSqClusterSizesPerDraw);
        Assert.Null(acc.BurnInMaxClusterSizePerDraw);
    }

    [Fact]
    public void Trace_LengthIsTheBurnInCount()
    {
        Accumulator acc = Run(Trace, Budget);

        Assert.Equal(50, acc.BurnInEnergyPerDraw!.Length);
        Assert.Equal(50, acc.BurnInMagPerDraw!.Length);
        Assert.Equal(50, acc.BurnInSumSqClusterSizesPerDraw!.Length);
        Assert.Equal(50, acc.BurnInMaxClusterSizePerDraw!.Length);
        Assert.Equal(400, acc.EnergyPerDraw.Length);
    }

    [Fact]
    public void Trace_IsEmptyWithoutBurnIn()
    {
        Accumulator acc = Run(Trace, new RunBudget(0, 10));

        Assert.Empty(acc.BurnInEnergyPerDraw!);
        Assert.Empty(acc.BurnInMaxClusterSizePerDraw!);
    }

    [Fact]
    public void Trace_IsTheHeadOfTheUnburnedStream()
    {
        Accumulator traced = Run(Trace, Budget);
        Accumulator unburned = Run(AccumulationSpec.None, new RunBudget(0, 450));

        Assert.Equal(unburned.EnergyPerDraw[..50], traced.BurnInEnergyPerDraw);
        Assert.Equal(unburned.MagPerDraw[..50], traced.BurnInMagPerDraw);
        Assert.Equal(unburned.SumSqClusterSizesPerDraw[..50], traced.BurnInSumSqClusterSizesPerDraw);
        Assert.Equal(unburned.MaxClusterSizePerDraw[..50], traced.BurnInMaxClusterSizePerDraw);
    }

    [Fact]
    public void Trace_LeavesTheMeasuredRecordUnchanged()
    {
        Accumulator off = Run(AccumulationSpec.None, Budget);
        Accumulator on = Run(Trace, Budget);

        Assert.Equal(off.EnergyPerDraw, on.EnergyPerDraw);
        Assert.Equal(off.MagPerDraw, on.MagPerDraw);
        Assert.Equal(off.SumSqClusterSizesPerDraw, on.SumSqClusterSizesPerDraw);
        Assert.Equal(off.MaxClusterSizePerDraw, on.MaxClusterSizePerDraw);
        Assert.Equal(SwAccumulatorProbeTests.Fingerprint(off), SwAccumulatorProbeTests.Fingerprint(on));
    }

    // ── SPCX v10 ─────────────────────────────────────────────────────────────

    [Fact]
    public void Trace_SurvivesSerialization_PresentAndAbsent()
    {
        Accumulator on = Run(Trace, Budget);
        Accumulator roundOn = SwDrawRecordTests.RoundTrip(on);

        Assert.Equal(on.BurnInEnergyPerDraw, roundOn.BurnInEnergyPerDraw);
        Assert.Equal(on.BurnInMagPerDraw, roundOn.BurnInMagPerDraw);
        Assert.Equal(on.BurnInSumSqClusterSizesPerDraw, roundOn.BurnInSumSqClusterSizesPerDraw);
        Assert.Equal(on.BurnInMaxClusterSizePerDraw, roundOn.BurnInMaxClusterSizePerDraw);
        Assert.Equal(on.EnergyPerDraw, roundOn.EnergyPerDraw);

        Accumulator roundOff = SwDrawRecordTests.RoundTrip(Run(AccumulationSpec.None, Budget));
        Assert.Null(roundOff.BurnInEnergyPerDraw);
        Assert.Null(roundOff.BurnInMaxClusterSizePerDraw);
    }

    [Fact]
    public void Serializer_RefusesAPartialTrace()
    {
        Accumulator on = Run(Trace, Budget);

        using var ms = new MemoryStream();
        Assert.Throws<ArgumentException>(() =>
            AccumulatorSerializer.Instance.WriteTo(on with { BurnInMagPerDraw = null }, ms));
    }

    // ── The engine: repeated burn-ins, Restore ───────────────────────────────

    [Fact]
    public void EveryDiscardedDraw_JoinsTheTrace()
    {
        // Bare draws before a burn-in are discarded by its reset too, and a second burn-in appends.
        var traced = Engine(Trace, seed: 13);
        for (int i = 0; i < 5; i++) traced.Draw();
        traced.BurnIn(20);
        traced.BurnIn(30);
        var plain = Engine(AccumulationSpec.None, seed: 13);
        plain.Run(55);

        Accumulator t = traced.GetCheckpoint(), p = plain.GetCheckpoint();
        Assert.Equal(0, t.DrawCount);
        Assert.Equal(p.EnergyPerDraw, t.BurnInEnergyPerDraw);
        Assert.Equal(p.MagPerDraw, t.BurnInMagPerDraw);
        Assert.Equal(p.SumSqClusterSizesPerDraw, t.BurnInSumSqClusterSizesPerDraw);
        Assert.Equal(p.MaxClusterSizePerDraw, t.BurnInMaxClusterSizePerDraw);
    }

    [Fact]
    public void Restore_CarriesTheTrace()
    {
        var original = Engine(Trace, seed: 17);
        original.BurnIn(25);
        original.Run(40);
        Accumulator cp = original.GetCheckpoint();
        var resumed = Engine(Trace, seed: 1);
        resumed.Restore(cp);

        Accumulator round = resumed.GetCheckpoint();
        Assert.Equal(cp.BurnInEnergyPerDraw, round.BurnInEnergyPerDraw);
        Assert.Equal(cp.BurnInMagPerDraw, round.BurnInMagPerDraw);
        Assert.Equal(cp.BurnInSumSqClusterSizesPerDraw, round.BurnInSumSqClusterSizesPerDraw);
        Assert.Equal(cp.BurnInMaxClusterSizePerDraw, round.BurnInMaxClusterSizePerDraw);
    }

    [Fact]
    public void Restore_RefusesATraceThatMissesTheCapability()
    {
        var traced = Engine(Trace, seed: 19);
        traced.BurnIn(10);
        Accumulator withTrace = traced.GetCheckpoint();
        var plain = Engine(AccumulationSpec.None, seed: 19);
        plain.BurnIn(10);
        Accumulator withoutTrace = plain.GetCheckpoint();

        Assert.Throws<InvalidOperationException>(() => Engine(AccumulationSpec.None, 1).Restore(withTrace));
        Assert.Throws<InvalidOperationException>(() => Engine(Trace, 1).Restore(withoutTrace));
        Assert.Throws<InvalidOperationException>(() =>
            Engine(Trace, 1).Restore(withTrace with { BurnInSumSqClusterSizesPerDraw = withTrace.BurnInSumSqClusterSizesPerDraw![..9] }));
    }

    // ── Helpers ──────────────────────────────────────────────────────────────

    private static Accumulator Run(AccumulationSpec accumulation, RunBudget budget)
        => SwRunner.Run(new SwRunSpec
        {
            Graph = SwAccumulatorProbeTests.Random16(), Temperature = 0.6, Q = 20,
            Accumulation = accumulation, Seed = SwAccumulatorProbeTests.Seed, Budget = budget,
        }).Accumulator;

    private static SwendsenWang<NoCurrencies> Engine(AccumulationSpec accumulation, int seed)
        => new(SwAccumulatorProbeTests.Random16(), 0.6, 20, accumulation, seed);
}
