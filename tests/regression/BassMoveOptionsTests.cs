using System;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// F-M2's Phase 1 options as configuration (the local-move options brief, Step 1): the local steps and the fixed-M
/// factor moves on <see cref="BassSampling"/>, the split rule on <see cref="TemperingSpec"/>. Each is off by default and
/// its record's defaults are the decided values; each validates by name, resolves on the data and round-trips through
/// JSON.
/// </summary>
public sealed class BassMoveOptionsTests
{
    [Fact]
    public void Options_AreOffByDefault_AndTheirRecordsCarryTheDecidedValues()
    {
        var config = new BassConfig();
        Assert.Null(config.Sampling.LocalSteps);
        Assert.Null(config.Sampling.FactorMoves);
        Assert.Null(new TemperingSpec().Split);
        Assert.Equal((0.5, 0.02), (new LocalSteps().Probability, new LocalSteps().WindowFraction));
        Assert.Equal((1.0, 1.0, 1.0), (new FactorMoves().Change, new FactorMoves().PruneOrGrow, new FactorMoves().VariableChange));
        Assert.Equal(0.9, new SplitRule().Threshold);

        BassResolved resolved = BassResolved.Resolve(config, BassSweepTests.Data());
        Assert.Null(resolved.StepWindows);
        Assert.Null(resolved.FactorShares);
        Assert.Null(resolved.SplitThreshold);
        Assert.DoesNotContain(resolved.AutoPicks, p => p.Setting == "LocalSteps windows");
    }

    [Theory]
    [InlineData(1, 1)]
    [InlineData(3, 1)]
    [InlineData(20, 1)]
    [InlineData(50, 1)]
    [InlineData(51, 2)]
    [InlineData(60, 2)]
    [InlineData(150, 3)]   // 0.02·150 is 3.0000000000000004 in binary; decimal keeps it 3
    [InlineData(500, 10)]
    public void Window_IsTheCeilingOfTheFraction_AtLeastOne(int values, int expected)
        => Assert.Equal(expected, new LocalSteps().Window(values));

    [Fact]
    public void Options_RefuseOutOfRangeFields_ByName()
    {
        string Reason(BassSampling sampling) => Assert.Throws<ArgumentException>(new BassConfig { Sampling = sampling }.Validate).Message;
        Assert.StartsWith("Sampling.LocalSteps.Probability", Reason(new BassSampling { LocalSteps = new LocalSteps(0) }));
        Assert.StartsWith("Sampling.LocalSteps.Probability", Reason(new BassSampling { LocalSteps = new LocalSteps(1.5) }));
        Assert.StartsWith("Sampling.LocalSteps.WindowFraction", Reason(new BassSampling { LocalSteps = new LocalSteps(0.5, 0) }));
        Assert.StartsWith("Sampling.LocalSteps.WindowFraction", Reason(new BassSampling { LocalSteps = new LocalSteps(0.5, double.NaN) }));
        Assert.StartsWith("Sampling.FactorMoves.Change", Reason(new BassSampling { FactorMoves = new FactorMoves(Change: -1) }));
        Assert.StartsWith("Sampling.FactorMoves.PruneOrGrow", Reason(new BassSampling { FactorMoves = new FactorMoves(PruneOrGrow: double.PositiveInfinity) }));
        Assert.StartsWith("Sampling.FactorMoves.VariableChange", Reason(new BassSampling { FactorMoves = new FactorMoves(VariableChange: double.NaN) }));
        Assert.StartsWith("Sampling.FactorMoves", Reason(new BassSampling { FactorMoves = new FactorMoves(1, 0, 0) }));
        Assert.StartsWith("Sampling.FactorMoves.Change",
                          Reason(new BassSampling { FactorMoves = new FactorMoves(Change: 0), LocalSteps = new LocalSteps() }));
        Assert.StartsWith("Sampling.Tempering.Split.Threshold", Reason(new BassSampling { Tempering = new TemperingSpec { Split = new SplitRule(1) } }));
        Assert.StartsWith("Sampling.Tempering.Split.Threshold", Reason(new BassSampling { Tempering = new TemperingSpec { Split = new SplitRule(0) } }));
        Assert.StartsWith("Sampling.Tempering.Split", Reason(new BassSampling { Tempering = new TemperingSpec { Levels = 6, Split = new SplitRule() } }));

        // A family without the change keeps the prune or grow and the variable change, and needs no local steps.
        new BassConfig { Sampling = new BassSampling { FactorMoves = new FactorMoves(Change: 0) } }.Validate();
    }

    [Fact]
    public void Options_ResolveOnTheData()
    {
        BassConfig config = BassSweepTests.Config() with
        {
            Sampling = new BassSampling
            {
                LocalSteps = new LocalSteps(0.5, 0.5), FactorMoves = new FactorMoves(2, 1, 1),
                Tempering = new TemperingSpec { Split = new SplitRule(0.8) },
            },
        };
        BassResolved resolved = BassResolved.Resolve(config, BassSweepTests.Data());   // three knots per variable
        Assert.Equal([2, 2], resolved.StepWindows!.Continuous);
        Assert.Empty(resolved.StepWindows.Functional);
        Assert.Equal(new BassFactorShares(0.5, 0.125, 0.125, 0.25), resolved.FactorShares);
        Assert.Equal(0.8, resolved.SplitThreshold);
        Assert.Contains(resolved.AutoPicks, p => p.Setting == "LocalSteps windows" && p.Value == "continuous [2, 2]");

        // The functional fixture's output variable has 20 locations; the windows follow each variable's own knots.
        BassMixing.Fixture functional = BassMixing.Get(BassMixing.Functional);
        BassResolved f = BassResolved.Resolve(functional.Config with { Sampling = functional.Config.Sampling with { LocalSteps = new LocalSteps() } },
                                              functional.Data);
        var model = BassModel.Create(functional.Config, functional.Data);
        Assert.Equal(Enumerable.Range(0, 2).Select(v => Math.Max(1, (int)Math.Ceiling(0.02m * model.Inputs.ContinuousKnots(v).Length))),
                     f.StepWindows!.Continuous);
        Assert.Equal([1], f.StepWindows.Functional);
        Assert.Null(f.FactorShares);
        Assert.Null(f.SplitThreshold);
    }

    [Fact]
    public void Options_RoundTripThroughJson()
    {
        var config = new BassConfig
        {
            Sampling = new BassSampling
            {
                LocalSteps = new LocalSteps(0.25, 0.05), FactorMoves = new FactorMoves(1, 2, 0.5),
                Tempering = new TemperingSpec { Split = new SplitRule(0.85) },
            },
        };
        string json = config.ToJson();
        Assert.Equal(config, BassConfig.FromJson(json));
        Assert.Equal(json, BassConfig.FromJson(json).ToJson());
        Assert.Contains("\"LocalSteps\":{\"Probability\":0.25,\"WindowFraction\":0.05}", json);
        Assert.Contains("\"FactorMoves\":{\"Change\":1,\"PruneOrGrow\":2,\"VariableChange\":0.5}", json);
        Assert.Contains("\"Split\":{\"Threshold\":0.85}", json);

        // Off reads back off; a record left empty takes its decided values.
        Assert.Equal(new BassConfig(), BassConfig.FromJson(new BassConfig().ToJson()));
        BassConfig partial = BassConfig.FromJson("{\"Sampling\":{\"LocalSteps\":{},\"FactorMoves\":{}}}");
        Assert.Equal(new LocalSteps(), partial.Sampling.LocalSteps);
        Assert.Equal(new FactorMoves(), partial.Sampling.FactorMoves);

        // The fingerprint moves with each option.
        BassData data = BassSweepTests.Data();
        string baseline = Bass.Fingerprint(new BassConfig(), data);
        Assert.NotEqual(baseline, Bass.Fingerprint(new BassConfig { Sampling = new BassSampling { LocalSteps = new LocalSteps() } }, data));
        Assert.NotEqual(baseline, Bass.Fingerprint(new BassConfig { Sampling = new BassSampling { FactorMoves = new FactorMoves() } }, data));
    }
}
