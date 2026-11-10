using System;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Xunit;

namespace Maths.Regression.Tests;

public sealed class RandomEffectsConfigurationTests
{
    private static CategoricalPredictor First => new("first", [0,0,1,1,2,2], ["a","b","c"]);
    private static CategoricalPredictor Second => new("second", [0,1,0,1,0,1], ["d","e"]);

    [Fact]
    public void NamedGroupings_AreCopiedAndValidatedBesideThePredictors()
    {
        BassData plain = BassSweepTests.Data();
        var first = new CategoricalPredictor("first", [.. Enumerable.Range(0, plain.Response.Cases).Select(i => i % 3)], ["a","b","c"]);
        var second = new CategoricalPredictor("second", [.. Enumerable.Range(0, plain.Response.Cases).Select(i => i % 2)], ["d","e"]);
        CategoricalPredictor[] source = [first, second];
        var data = new BassData(plain.Predictors, plain.Response, source);
        source[0] = second;
        Assert.Equal(["first", "second"], data.Groupings.Select(g => g.Name));
        Assert.Empty(plain.Groupings);
        Assert.Contains("used twice", Assert.Throws<ArgumentException>(() => new BassData(plain.Predictors, plain.Response, [first, first])).Message);
        var duplicate = new CategoricalPredictor(plain.Predictors[0].Name, first.Codes.ToArray(), first.Levels);
        Assert.Contains("used twice", Assert.Throws<ArgumentException>(() => new BassData(plain.Predictors, plain.Response, [duplicate])).Message);
        Assert.Contains("cases", Assert.Throws<ArgumentException>(() => new BassData(plain.Predictors, plain.Response, [First])).Message);
    }

    [Fact]
    public void NamedBlocksAndDraws_ResolveByNameWithDistinctScales()
    {
        BassRandomIntercept[] specs = [new("second", 1.7), new("first", 0.4, ScalePrior: new BassGroupScalePrior(0.8))];
        RandomEffectsBlock block = RandomEffectsBlock.Of([First, Second], specs);
        Assert.Equal(["second", "first"], block.Groupings.Select(g => g.Name));
        Assert.Equal([1.7, 0.4], block.Scales.ToArray());
        var groups = new BassGroupDraw[] { new("first", 0.9, new double[3]), new("second", 2.3, new double[2]) };
        RandomEffectsBlock at = block.At(groups);
        Assert.Equal([2.3, 0.9], at.Scales.ToArray());
        Assert.Equal([1.7, 0.4], block.Scales.ToArray());
        Assert.Equal([2.3, 0.25], at.At(1, 0.25).Scales.ToArray());
        Assert.Same(at, at.At([2.3, 0.9]));
        Assert.Throws<ArgumentException>(() => block.At([groups[0], groups[0]]));
        Assert.Throws<ArgumentException>(() => block.At([new BassGroupDraw("other", 1, new double[3]), groups[1]]));
        Assert.Throws<ArgumentException>(() => block.At([groups[0] with { Effects = new double[2] }, groups[1]]));
        Assert.Throws<ArgumentOutOfRangeException>(() => at.At(1, 0));
        var config = new BassConfig { RandomIntercepts = specs };
        string json = config.ToJson();
        Assert.Contains("\"RandomIntercepts\"", json);
        Assert.Equal(json, BassConfig.FromJson(json).ToJson());
        Assert.Throws<System.Text.Json.JsonException>(() => BassConfig.FromJson("{\"RandomIntercept\":{\"Scale\":1.3}}"));
    }

    [Fact]
    public void Declarations_RequireKnownNamesAndAConsistentParentRelation()
    {
        Assert.Throws<ArgumentException>(() => RandomEffectsBlock.Of([First], [new BassRandomIntercept("missing", 1)]));
        Assert.Throws<ArgumentException>(() => RandomEffectsBlock.Of([First], [new BassRandomIntercept("first", 1, RandomInterceptConvention.WithinParent, ParentGrouping: "missing")]));
        BassRandomIntercept nested = new("first", 1, RandomInterceptConvention.WithinParent, ParentGrouping: "second");
        Assert.Contains("more than one", Assert.Throws<ArgumentException>(() => RandomEffectsBlock.Of([First, Second], [nested])).Message);
        Assert.Throws<ArgumentException>(() => (new BassConfig { RandomIntercepts = [nested with { ParentGrouping = null }] }).Validate());
        Assert.Throws<ArgumentException>(() => (new BassConfig { RandomIntercepts = [nested with { Convention = RandomInterceptConvention.ObservationCentred }] }).Validate());
        Assert.Throws<ArgumentException>(() => (new BassConfig { RandomIntercepts = [new("first", 1), new("first", 2)] }).Validate());
        Assert.Throws<ArgumentException>(() => (new BassConfig { RandomIntercepts =
            [nested, new("second", 1, RandomInterceptConvention.WithinParent, ParentGrouping: "first")] }).Validate());
    }

    [Fact]
    [Trait("Tier", "Pinned")]
    public void UngroupedFingerprints_KeepThePreMigrationIdentity()
    {
        BassData data = BassSweepTests.Data();
        Assert.Equal("77a50a29e3922982", Bass.Fingerprint(new BassConfig(), data));
        Assert.Equal("5663b1f7fddc270e", Bass.Fingerprint(BassSweepTests.Config(maxBasis: 3), data));
    }
}
