using System;
using System.Diagnostics;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class RandomEffectsMultipleQualificationTests(ITestOutputHelper output)
{
    private const double Sigma = 0.35;
    private sealed record Generated(BassData Data, double[][] Effects, double Offset, double[] Surface,
                                    PredictionCases Population, PredictionCases Grouped, double[] TestSurface, double[] TestMean);

    private static double Truth(double a, double b) => 1.1 + 1.5 * (a - 0.5) + 0.9 * (b - 0.5) * (b - 0.5) + 0.8 * (a - 0.5) * (b - 0.5);

    // The first effects are two-point, the second shifted exponential. Neither follows the fitted Gaussian law.
    // Child means belong to the parent; both are transformed into the declared conventions before generating y.
    private static Generated Generate(string kind, int seed)
    {
        const int n = 320, held = 400;
        int firstLevels = kind == "nested" ? 4 : 6, secondLevels = kind == "nested" ? 12 : 7;
        var rng = new Xoshiro256PlusPlus(seed);
        var first = new int[n]; var second = new int[n];
        for (int i = 0; i < n; i++)
        {
            second[i] = i < 3 * secondLevels ? i % secondLevels : Math.Min(secondLevels - 1, (int)(secondLevels * Math.Pow(rng.NextDouble(), 1.4)));
            first[i] = kind == "nested" ? second[i] / 3 : i < 3 * firstLevels ? i % firstLevels : Math.Min(firstLevels - 1, (int)(firstLevels * Math.Pow(rng.NextDouble(), 1.6)));
        }
        double[][] effects = [[.. Enumerable.Range(0, firstLevels).Select(_ => rng.NextDouble() < 0.5 ? -0.8 : 0.8)],
                              [.. Enumerable.Range(0, secondLevels).Select(_ => -0.8 * Math.Log(1 - rng.NextDouble()) - 0.8)]];
        double offset = 0;
        if (kind == "nested")
            for (int p = 0; p < firstLevels; p++)
            {
                int parent = p;
                int[] rows = [.. Enumerable.Range(0, n).Where(i => first[i] == parent)];
                double mean = rows.Average(i => effects[1][second[i]]);
                effects[0][p] += mean;
                for (int g = 3 * p; g < 3 * p + 3; g++) effects[1][g] -= mean;
            }
        for (int k = 0; k < (kind == "nested" ? 1 : 2); k++)
        {
            int[] codes = k == 0 ? first : second;
            double mean = codes.Average(g => effects[k][g]);
            offset += mean;
            for (int g = 0; g < effects[k].Length; g++) effects[k][g] -= mean;
        }
        double[] x0 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())], x1 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        double[] surface = [.. x0.Zip(x1, Truth)];
        double[] y = [.. surface.Select((f, i) => f + offset + effects[0][first[i]] + effects[1][second[i]] + Sigma * Normal.SampleStandard(rng))];
        string[] names0 = [.. Enumerable.Range(0, firstLevels).Select(g => $"A{g}")], names1 = [.. Enumerable.Range(0, secondLevels).Select(g => $"B{g}")];
        var data = new BassData([new ContinuousPredictor("x0", x0, new InputDomain(0, 1)), new ContinuousPredictor("x1", x1, new InputDomain(0, 1))],
            new ScalarResponse(y), [new CategoricalPredictor("first", first, names0), new CategoricalPredictor("second", second, names1)]);
        double[] t0 = [.. Enumerable.Range(0, held).Select(_ => rng.NextDouble())], t1 = [.. Enumerable.Range(0, held).Select(_ => rng.NextDouble())];
        int[] g1 = [.. Enumerable.Range(0, held).Select(_ => (int)(rng.NextDouble() * secondLevels))];
        int[] g0 = kind == "nested" ? [.. g1.Select(g => g / 3)] : [.. Enumerable.Range(0, held).Select(_ => (int)(rng.NextDouble() * firstLevels))];
        var population = new PredictionCases([new ContinuousColumn("x0", t0), new ContinuousColumn("x1", t1)]);
        var grouped = new PredictionCases([.. population.Columns, new CategoricalColumn("first", [.. g0.Select(g => names0[g])]),
            new CategoricalColumn("second", [.. g1.Select(g => names1[g])])]);
        double[] testSurface = [.. t0.Zip(t1, Truth).Select(f => f + offset)];
        return new Generated(data, effects, offset, surface, population, grouped, testSurface,
            [.. testSurface.Select((f, i) => f + effects[0][g0[i]] + effects[1][g1[i]])]);
    }

    private static (double Low, double High) Interval(double[] x)
    {
        Array.Sort(x);
        return (PosteriorInterval.Quantile(x, 0.05), PosteriorInterval.Quantile(x, 0.95));
    }

    private static double Error(BassFit fit, BassData data, PredictionCases cases, double[] truth)
    {
        BassPrediction prediction = fit.Predict(data, cases, new BassPredictOptions { Predictive = false });
        double squares = 0;
        for (int i = 0; i < truth.Length; i++)
        {
            double mean = Enumerable.Range(0, prediction.DrawCount).Sum(d => prediction.Weights[d] * prediction.Latent[d, i, 0]);
            squares += (mean - truth[i]) * (mean - truth[i]);
        }
        return Math.Sqrt(squares / truth.Length);
    }

    [Theory]
    [InlineData("crossed")]
    [InlineData("nested")]
    public void NonGaussianEffects_AreRecovered_WithBothScalesUnknown(string kind)
    {
        if (!BassQualification.Enabled) return;
        const int datasets = 4;
        int covered = 0, intervals = 0, scalesCovered = 0, sigmaCovered = 0;
        double z2 = 0;
        for (int seed = 0; seed < datasets; seed++)
        {
            Generated g = Generate(kind, 210520 + seed);
            BassConfig plain = KhaosReachQualificationTests.Config(null, g.Data.Response.Cases, 210530 + seed, burnIn: 1500, samples: 700);
            plain = plain with { Structure = plain.Structure with { MaxBasis = 8, Legendre = new LegendreBasis { MaxDegree = 3 } } };
            BassConfig config = plain with { RandomIntercepts = [new("first", 1, ScalePrior: new BassGroupScalePrior()),
                new("second", 1, kind == "nested" ? RandomInterceptConvention.WithinParent : RandomInterceptConvention.ObservationCentred,
                    new BassGroupScalePrior(), kind == "nested" ? "first" : null)] };
            BassFit fit = Bass.Fit(config, g.Data, new BassRunOptions { MaxDegreeOfParallelism = 2 });
            BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
            BassModel model = BassModel.Create(config, g.Data);
            for (int k = 0; k < 2; k++)
            {
                for (int level = 0; level < g.Effects[k].Length; level++)
                {
                    int group = k, index = level;
                    double[] values = [.. draws.Select(d => d.Groups![group].Effects[index])];
                    double mean = values.Average(), variance = values.Sum(v => (v - mean) * (v - mean)) / (values.Length - 1);
                    z2 += (mean - g.Effects[k][level]) * (mean - g.Effects[k][level]) / variance;
                    var (lo, hi) = Interval(values);
                    if (lo <= g.Effects[k][level] && g.Effects[k][level] <= hi) covered++;
                    intervals++;
                }
                int block = k;
                var (low, high) = Interval([.. draws.Select(d => d.Groups![block].Scale)]);
                // Under misspecification the reported scale target is the realized contrast RMS / σ, not a claimed Gaussian population parameter.
                double target = Math.Sqrt(g.Effects[k].Sum(u => u * u) / model.RandomEffects!.Layout(k).Dimension) / Sigma;
                if (low <= target && target <= high) scalesCovered++;
                output.WriteLine($"{kind}, dataset {seed}, group {k}: 90% s interval [{low:F3}, {high:F3}], realized contrast RMS/σ {target:F3}");
            }
            var (sigmaLow, sigmaHigh) = Interval([.. draws.Select(d => Math.Sqrt(d.Sigma2))]);
            if (sigmaLow <= Sigma && Sigma <= sigmaHigh) sigmaCovered++;
            BassPrediction population = fit.Predict(g.Data, PredictionCases.Of(g.Data.Predictors), new BassPredictOptions { Predictive = false });
            double[] levels = [.. Enumerable.Range(0, population.DrawCount).Select(d => Enumerable.Range(0, g.Data.Response.Cases).Average(i => population.Latent[d, i, 0]))];
            double meanLevel = levels.Average(), levelSd = Math.Sqrt(levels.Sum(v => (v - meanLevel) * (v - meanLevel)) / (levels.Length - 1));
            double levelTarget = g.Surface.Average() + g.Offset;
            BassData bare = new(g.Data.Predictors, g.Data.Response);
            BassFit without = Bass.Fit(plain, bare, new BassRunOptions { MaxDegreeOfParallelism = 2 });
            double grouped = Error(fit, g.Data, g.Grouped, g.TestMean), ungrouped = Error(without, bare, g.Population, g.TestMean);
            double surface = Error(fit, g.Data, g.Population, g.TestSurface), alternative = Error(without, bare, g.Population, g.TestSurface);
            output.WriteLine($"{kind}, dataset {seed}: level {meanLevel:F4} ± {levelSd:F4}, target {levelTarget:F4}; 90% σ [{sigmaLow:F3}, {sigmaHigh:F3}]; " +
                $"held-out RMSE {grouped:F4} (without blocks {ungrouped:F4}), population surface {surface:F4} (without {alternative:F4})");
            Assert.InRange(Math.Abs(meanLevel - levelTarget), 0, 2 * levelSd);
            Assert.InRange(grouped, 0, 0.5 * ungrouped);
            Assert.InRange(surface, 0, 1.1 * alternative);
        }
        output.WriteLine($"{kind}: 90% effects covered {covered}/{intervals}, standardized error mean square {z2/intervals:F3}; scale intervals {scalesCovered}/{2*datasets}, σ intervals {sigmaCovered}/{datasets}");
        Assert.InRange((double)covered / intervals, 0.75, 1);
        Assert.InRange(z2 / intervals, 0.15, 2.5);
        Assert.InRange(scalesCovered, 6, 8);
        Assert.InRange(sigmaCovered, 2, 4);
    }

    [Fact]
    public void Cost_AtTheRetainedRank()
    {
        if (!BassQualification.Enabled) return;
        const int cases = 512, updates = 2000, columns = 4000;
        foreach (int levels in new[] { 2, 4, 8, 16 })
        {
            var first = new CategoricalPredictor("first", [.. Enumerable.Range(0, cases).Select(i => i % levels)], [.. Enumerable.Range(0, levels).Select(i => $"A{i}")]);
            var second = new CategoricalPredictor("second", [.. Enumerable.Range(0, cases).Select(i => i / levels % levels)], [.. Enumerable.Range(0, levels).Select(i => $"B{i}")]);
            RandomEffectsBlock block = RandomEffectsBlock.Of([first, second], [new("first", 0.8), new("second", 1.3)]);
            var reading = block;
            for (int j = 0; j < 100; j++) reading = block.At(0, 0.5 + j % 17 * 0.03);
            var watch = Stopwatch.StartNew();
            for (int j = 0; j < updates; j++) reading = block.At(0, 0.5 + j % 17 * 0.03);
            double update = watch.Elapsed.TotalMicroseconds / updates;
            double[] column = [.. Enumerable.Range(0, cases).Select(i => Math.Sin(i))];
            watch.Restart();
            for (int j = 0; j < columns; j++) { column[j % cases] = 0.1; reading.Whiten(column); }
            output.WriteLine($"n={cases}, q={2*(levels-1)}, r={block.Rank}: scale update {update:F2} µs; column {watch.Elapsed.TotalMicroseconds/columns:F2} µs");
            Assert.Equal(2 * (levels - 1), block.Rank);
            Assert.True(double.IsFinite(reading.LogDeterminant));
        }
    }
}
