using System;
using System.Text.Json;
using Maths.Samplers.FlatHistogram;
using Xunit;
using static Maths.Samplers.Tests.FlatHistogramToys;

namespace Maths.Samplers.Tests;

/// <summary>
/// The flat-histogram run's configuration and resume (N17 Phase 1, step 4). The spec is values only: it validates, it
/// round-trips through JSON unchanged, and it says when it needs a pilot. A run taken in pieces, each resuming the last
/// one's checkpoint, mid-learning, at the learning's end and mid-production, ends bit for bit where the run taken whole
/// ends; a checkpoint of another run, or of a finished one, is refused.
/// </summary>
public sealed class FlatHistogramSpecTests
{
    private const double Depth = 10.0;

    private static FlatHistogramSpec Explicit => new()
    {
        Chains = 3,
        LearningIterations = 2_000,
        ProductionDraws = 1_000,
        Window = new FlatHistogramWindow(-30.0, -0.02),
        LowSlope = 0.05,
        Nodes = 31,
        GainScale = 120.0,
    };

    [Fact]
    public void Spec_RoundTripsThroughJson()
    {
        FlatHistogramSpec[] specs =
        [
            Explicit,
            new() { LearningIterations = 10, ProductionDraws = 5, Start = FlatHistogramStart.Wham, GainStart = 77, WindowQuantile = 0.05, WindowMargin = 0.5 },
        ];
        foreach (FlatHistogramSpec spec in specs)
        {
            string json = JsonSerializer.Serialize(spec);
            Assert.DoesNotContain(nameof(FlatHistogramSpec.NeedsPilot), json);   // computed, not configured
            FlatHistogramSpec back = JsonSerializer.Deserialize<FlatHistogramSpec>(json)!;
            Assert.Equal(spec, back);
            Assert.Equal(json, JsonSerializer.Serialize(back));
        }
    }

    [Fact]
    public void Spec_Validates_AndSaysWhenItNeedsAPilot()
    {
        Explicit.Validate();
        Assert.False(Explicit.NeedsPilot);
        Assert.True((Explicit with { Start = FlatHistogramStart.Wham }).NeedsPilot);
        Assert.True((Explicit with { Window = null }).NeedsPilot);
        Assert.True((Explicit with { LowSlope = null }).NeedsPilot);
        Assert.True((Explicit with { Nodes = null }).NeedsPilot);

        FlatHistogramSpec[] bad =
        [
            Explicit with { Chains = 0 },
            Explicit with { LearningIterations = -1 },
            Explicit with { ProductionDraws = 0 },
            Explicit with { Start = (FlatHistogramStart)7 },
            Explicit with { Window = new FlatHistogramWindow(0.0, 0.0) },
            Explicit with { Window = new FlatHistogramWindow(double.NaN, 1.0) },
            Explicit with { LowSlope = 1.5 },
            Explicit with { LowSlope = -0.1 },
            Explicit with { Nodes = 1 },
            Explicit with { GainScale = 0.5 },
            Explicit with { GainScale = double.PositiveInfinity },
            Explicit with { GainStart = -1 },
            Explicit with { WindowQuantile = 0.5 },
            Explicit with { WindowQuantile = 0.0 },
            Explicit with { WindowMargin = -1.0 },
        ];
        foreach (FlatHistogramSpec spec in bad)
        {
            Assert.Throws<ArgumentException>(spec.Validate);
            Assert.Throws<ArgumentException>(() => MulticanonicalRun.Run(spec, DoubleWell(Depth), new WellModel(Depth, 0.25), 1));
        }
        Assert.Throws<ArgumentException>(() => MulticanonicalRun.Run(Explicit with { Nodes = null }, DoubleWell(Depth), new WellModel(Depth, 0.25), 1));
    }

    [Fact]
    public void ARunInPieces_EndsBitForBit_WhereTheWholeRunEnds()
    {
        FlatHistogramSpec spec = Explicit;
        var model = new WellModel(Depth, 0.25);
        MulticanonicalResult<double> whole = MulticanonicalRun.Run(spec, DoubleWell(Depth), model, 23);

        // Mid-learning, the learning's end exactly, mid-production, the rest.
        MulticanonicalResult<double> a = MulticanonicalRun.Run(spec, DoubleWell(Depth), model, 23, iterations: 700, checkpoint: true);
        Assert.False(a.Complete);
        Assert.Null(a.Production);
        Assert.Equal(700, a.Learning.Iterations);
        MulticanonicalResult<double> b = MulticanonicalRun.Run(spec, DoubleWell(Depth), model, 23, resumeFrom: a.Checkpoint, iterations: 1_300, checkpoint: true);
        Assert.Equal(2_000, b.Learning.Iterations);
        Assert.Null(b.Production);
        MulticanonicalResult<double> c = MulticanonicalRun.Run(spec, DoubleWell(Depth), model, 23, resumeFrom: b.Checkpoint, iterations: 400, checkpoint: true);
        Assert.Equal(400, c.Draws!.Draws);
        Assert.False(c.Complete);
        MulticanonicalResult<double> d = MulticanonicalRun.Run(spec, DoubleWell(Depth), model, 23, resumeFrom: c.Checkpoint, checkpoint: true);
        Assert.True(d.Complete);

        Assert.Equal(whole.Weights.LogWeights, d.Weights.LogWeights);
        Assert.Equal(whole.Learning, d.Learning with { Weights = whole.Learning.Weights });
        Assert.Equal(whole.WeightTrace.Count, d.WeightTrace.Count);
        for (int i = 0; i < whole.WeightTrace.Count; i++)
        {
            Assert.Equal(whole.WeightTrace[i].Iteration, d.WeightTrace[i].Iteration);
            Assert.Equal(whole.WeightTrace[i].Weights.LogWeights, d.WeightTrace[i].Weights.LogWeights);
        }
        AssertSameDiagnostics(whole.LearningDiagnostics, d.LearningDiagnostics);
        AssertSameDiagnostics(whole.ProductionDiagnostics!, d.ProductionDiagnostics!);
        for (int ch = 0; ch < spec.Chains; ch++)
        {
            Assert.Equal(whole.Draws!.Levels(ch).ToArray(), d.Draws!.Levels(ch).ToArray());
            for (int f = 0; f < 3; f++) Assert.Equal(whole.Draws.Functional(ch, f), d.Draws.Functional(ch, f));
        }
        Assert.Equal(whole.Production!.FunctionalMean, d.Production!.FunctionalMean);
        Assert.Equal(whole.Production.FunctionalVariance, d.Production.FunctionalVariance);
        Assert.Equal(whole.Production.AcceptanceRate, d.Production.AcceptanceRate);
        Assert.Equal(whole.Production.Verdict.Outcome, d.Production.Verdict.Outcome);
    }

    [Fact]
    public void Resume_RefusesAnotherRunsCheckpoint_AndAFinishedOne()
    {
        FlatHistogramSpec spec = Explicit;
        var model = new WellModel(Depth, 0.25);
        MulticanonicalResult<double> part = MulticanonicalRun.Run(spec, DoubleWell(Depth), model, 23, iterations: 100, checkpoint: true);
        Assert.Throws<ArgumentException>(() => MulticanonicalRun.Run(spec, DoubleWell(Depth), model, 24, resumeFrom: part.Checkpoint));
        Assert.Throws<ArgumentException>(() => MulticanonicalRun.Run(spec with { GainScale = 200.0 }, DoubleWell(Depth), model, 23, resumeFrom: part.Checkpoint));
        Assert.Throws<ArgumentException>(() => MulticanonicalRun.Run(spec, DoubleWell(Depth), new WellModel(Depth, 0.3), 23, resumeFrom: part.Checkpoint));
        Assert.Throws<ArgumentOutOfRangeException>(() => MulticanonicalRun.Run(spec, DoubleWell(Depth), model, 23, iterations: 0));

        MulticanonicalResult<double> done = MulticanonicalRun.Run(spec, DoubleWell(Depth), model, 23, checkpoint: true);
        Assert.Throws<InvalidOperationException>(() => MulticanonicalRun.Run(spec, DoubleWell(Depth), model, 23, resumeFrom: done.Checkpoint));
    }

    private static void AssertSameDiagnostics(FlatHistogramDiagnostics x, FlatHistogramDiagnostics y)
    {
        Assert.Equal(x.Occupancy, y.Occupancy);
        Assert.Equal((x.InWindow, x.BelowWindow, x.AboveWindow), (y.InWindow, y.BelowWindow, y.AboveWindow));
        Assert.Equal(x.RoundTrips, y.RoundTrips);
        Assert.Equal(x.LastEnd, y.LastEnd);
    }
}
