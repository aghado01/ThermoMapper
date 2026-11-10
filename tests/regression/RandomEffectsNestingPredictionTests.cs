using System;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Xunit;

namespace Maths.Regression.Tests;

public sealed class RandomEffectsNestingPredictionTests
{
    private static (BassFit Fit, BassData Data) Fit(bool parentBlock, bool nested = true)
    {
        const int n = 24;
        double[] x = [.. Enumerable.Range(0, n).Select(i => ((7*i)%n + 0.25)/n)];
        double[] y = [.. x.Select((v,i) => 0.7 + 1.3*v + (i < 12 ? -0.4 : 0.4) + 0.2*Math.Sin(i))];
        // Declaration order differs from encounter order, with unused declared levels in both columns.
        var parent = new CategoricalPredictor("parent", [.. Enumerable.Range(0,n).Select(i => i < 12 ? 2 : 1)], ["unused", "B", "A"]);
        var child = new CategoricalPredictor("child", [.. Enumerable.Range(0,n).Select(i => new[] {4,2,3,1}[i/6])], ["unused", "B1", "A1", "B0", "A0"]);
        var data = new BassData([new ContinuousPredictor("x", x, new InputDomain(0,1))], new ScalarResponse(y), [child, parent]);
        BassRandomIntercept inner = new("child", 0.6, nested ? RandomInterceptConvention.WithinParent : RandomInterceptConvention.ObservationCentred,
            ParentGrouping: nested ? "parent" : null);
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxBasis = 1, MaxInteraction = 1, Legendre = new LegendreBasis { MaxDegree = 2 } },
            Coefficients = BassCoefficientPrior.RawRidge(TauPrior.Held(0.7)),
            Intercept = GPriorCentre.Zero, Sigma2 = Sigma2Prior.InverseGamma(3,2),
            RandomIntercepts = parentBlock ? [inner, new("parent", 1.1)] : [inner],
            Sampling = new BassSampling { Chains = 1, BurnIn = 40, Samples = 40, MasterSeed = 210610 },
        };
        return (Bass.Fit(config, data), data);
    }

    private static PredictionCases Cases(string[]? parents, string[]? children, int count = 1)
    {
        PredictionColumn[] x = [new ContinuousColumn("x", [.. Enumerable.Repeat(0.4, count)])];
        return new PredictionCases([.. x,
            .. parents is null ? Array.Empty<PredictionColumn>() : new PredictionColumn[] {new CategoricalColumn("parent", parents)},
            .. children is null ? Array.Empty<PredictionColumn>() : new PredictionColumn[] {new CategoricalColumn("child", children)}]);
    }

    private static BassPredictOptions Options(bool meanSurface, UnseenLevel unseen = UnseenLevel.Refuse)
        => new() { MeanSurface = meanSurface, Predictive = true, Seed = 210611, UnseenLevels = unseen };

    [Theory]
    [InlineData(true, false, UnseenLevel.Refuse)]
    [InlineData(true, true, UnseenLevel.Refuse)]
    [InlineData(false, false, UnseenLevel.Refuse)]
    [InlineData(false, true, UnseenLevel.Refuse)]
    [InlineData(true, false, UnseenLevel.OutsideEverySubset)]
    [InlineData(true, true, UnseenLevel.OutsideEverySubset)]
    [InlineData(false, false, UnseenLevel.OutsideEverySubset)]
    [InlineData(false, true, UnseenLevel.OutsideEverySubset)]
    public void ExplicitObservedContradiction_IsRefused(bool parentBlock, bool meanSurface, UnseenLevel unseen)
    {
        var (fit, data) = Fit(parentBlock);
        ArgumentException error = Assert.Throws<ArgumentException>(() => fit.Predict(data, Cases(["B"], ["A0"]), Options(meanSurface, unseen)));
        Assert.Contains("child", error.Message);
        Assert.Contains("A0", error.Message);
        Assert.Contains("parent", error.Message);
        Assert.Contains("case 0", error.Message);
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void EveryCase_IsChecked(bool parentBlock)
    {
        var (fit, data) = Fit(parentBlock);
        ArgumentException error = Assert.Throws<ArgumentException>(() => fit.Predict(data,
            Cases(["A", "B"], ["A0", "A0"], count: 2), Options(true)));
        Assert.Contains("case 1", error.Message);
    }

    [Theory]
    [InlineData(true, false)]
    [InlineData(true, true)]
    [InlineData(false, false)]
    [InlineData(false, true)]
    public void ValidPairs_AndOmittedOrPopulationGroups_KeepTheirReadings(bool parentBlock, bool meanSurface)
    {
        var (fit, data) = Fit(parentBlock);
        BassPredictOptions options = Options(meanSurface);
        BassPrediction population = fit.Predict(data, Cases(null, null), options);
        BassPrediction valid = fit.Predict(data, Cases(["A"], ["A0"]), options);
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        for (int d = 0; d < valid.DrawCount; d++)
        {
            double effect = draws[d].Groups!.Single(g => g.Grouping == "child").Effects[3];
            if (parentBlock) effect += draws[d].Groups!.Single(g => g.Grouping == "parent").Effects[1];
            Assert.InRange(Math.Abs(valid.Latent[d,0,0] - population.Latent[d,0,0] - effect), 0, 1e-12);
        }
        if (meanSurface) Assert.True(double.IsFinite(valid.MeanSurface![0,0]));
        // Missing parents are still population reads, without inferred or newly imposed effects.
        BassPrediction omittedParent = fit.Predict(data, Cases(null, ["A0"]), options);
        BassPrediction omittedChild = fit.Predict(data, Cases(["B"], null), options);
        BassPredictOptions outside = Options(meanSurface, UnseenLevel.OutsideEverySubset);
        foreach (string parent in new[] { "new-parent", "unused" })
            Same(omittedParent, fit.Predict(data, Cases([parent], ["A0"]), outside));
        foreach (string child in new[] { "new-child", "unused" })
            Same(omittedChild, fit.Predict(data, Cases(["B"], [child]), outside));
        Same(omittedParent, fit.Predict(data, Cases(null, ["A0"]), outside));
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void CrossedGroups_AllowAnyObservedCombination(bool meanSurface)
    {
        var (fit, data) = Fit(parentBlock: true, nested: false);
        BassPrediction prediction = fit.Predict(data, Cases(["B"], ["A0"]), Options(meanSurface));
        Assert.Equal(fit.Draws[0].Count, prediction.DrawCount);
        if (meanSurface) Assert.True(double.IsFinite(prediction.MeanSurface![0,0]));
    }

    private static void Same(BassPrediction expected, BassPrediction actual)
    {
        Assert.Equal(expected.Latent, actual.Latent);
        Assert.Equal(expected.Predictive, actual.Predictive);
        Assert.Equal(expected.MeanSurface, actual.MeanSurface);
    }
}
