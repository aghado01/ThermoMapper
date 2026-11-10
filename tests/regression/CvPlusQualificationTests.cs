using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;
using Unit = Maths.Regression.Tests.CrossConformalTests.Unit;

namespace Maths.Regression.Tests;

public sealed class CvPlusQualificationTests(ITestOutputHelper output)
{
    private static Unit[] Draw(int n, Xoshiro256PlusPlus rng) => [.. Enumerable.Range(0, n).Select(_ =>
    {
        double x = 4 * rng.NextDouble() - 2;
        return new Unit(x, ConformalBayesTests.Truth(x, 0.5, rng));
    })];

    private static (double Mean, double Se) Summary(IReadOnlyList<double> values)
    {
        double mean = values.Average();
        return (mean, Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / (values.Count - 1) / values.Count));
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void HeldOutCoverageWidthAndApplicableBounds(bool orderReading)
    {
        if (!BassQualification.Enabled) return;
        const int n = 40, tests = 5, reps = 500;
        const double alpha = 0.1;
        int[] folds = [2, 5, 10, 40];
        var coverage = folds.Select(_ => new List<double>()).ToArray();
        var width = new double[folds.Length];
        var bounds = new CvPlusBounds[folds.Length];
        var splitCoverage = new List<double>();
        double splitWidth = 0;
        Func<double, double> First(Unit[] u) { double y = u[0].Y; return _ => y; }
        Func<Unit[], Func<double, double>> fit = orderReading ? First : CrossConformalQualificationTests.Line;
        CvPlusAssumptions assumptions = orderReading ? CvPlusAssumptions.MergedRanks : CvPlusAssumptions.Both;
        foreach (int seed in SeedTree.Derive(20261119, reps))
        {
            var rng = new Xoshiro256PlusPlus(seed);
            Unit[] train = Draw(n, rng), test = Draw(tests, rng);
            double[] xs = test.Select(u => u.X).ToArray();
            for (int f = 0; f < folds.Length; f++)
            {
                int k = folds[f];
                int[] labels = [.. Enumerable.Range(0, n).Select(i => i / (n / k))];
                FoldReadings readings = CvPlusTests.Readings(train, labels, k, xs, fit);
                CvPlusIntervals result = CrossConformal.CvPlus(readings, assumptions, alpha);
                bounds[f] = result.Bounds;
                int covered = 0;
                for (int t = 0; t < tests; t++)
                {
                    if (result.Intervals[t].Contains(test[t].Y)) covered++;
                    width[f] += result.Intervals[t].Width;
                }
                coverage[f].Add(covered / (double)tests);
                if (k == n)
                {
                    var jack = CrossConformal.JackknifePlus(CrossConformalTests.ResidualReadings(train, xs, fit), alpha);
                    for (int t = 0; t < tests; t++)
                    {
                        Assert.Equal(jack.Intervals[t].Lower, result.Intervals[t].Lower);
                        Assert.Equal(jack.Intervals[t].Upper, result.Intervals[t].Upper);
                    }
                }
            }
            Func<double, double> splitFit = fit(train[..20]);
            Unit[] cal = train[20..];
            var split = SplitConformal.Residuals(new(cal.Select(u => splitFit(u.X)).ToArray()),
                                                 cal.Select(u => u.Y).ToArray(), new(xs.Select(splitFit).ToArray()), alpha);
            int splitHits = 0;
            for (int t = 0; t < tests; t++)
            {
                if (split.Intervals[t].Contains(test[t].Y)) splitHits++;
                splitWidth += split.Intervals[t].Width;
            }
            splitCoverage.Add(splitHits / (double)tests);
        }
        output.WriteLine($"{(orderReading ? "first-unit order-reading fit, merged assumptions only" : "symmetric line, both assumption sets")}: {reps * tests} held-out units, n={n}, alpha={alpha}");
        for (int f = 0; f < folds.Length; f++)
        {
            var (mean, se) = Summary(coverage[f]);
            var b = bounds[f];
            Assert.True(mean >= b.Coverage!.Value - 4 * se);
            Assert.Equal(!orderReading, b.Tournament.Applies);
            double rankBound = Math.Max(b.Merged.Value, b.OneFold.Value);
            Assert.Equal(orderReading ? rankBound : Math.Max(rankBound, b.Tournament.Value), b.Coverage.Value);
            output.WriteLine($"K={folds[f]}, m={b.FoldSize}, level={b.Level}: coverage {mean:F5} (replication SE {se:F5}), width {width[f] / (reps * tests):F4}; merged {b.Merged.Value:F5}, one fold {b.OneFold.Value:F5}, tournament {b.Tournament.Value:F5} (applies {b.Tournament.Applies}), selected {b.Coverage:F5}");
        }
        var splitSummary = Summary(splitCoverage);
        output.WriteLine($"split comparator at equal data budget, 20 training + 20 calibration: coverage {splitSummary.Mean:F5} (SE {splitSummary.Se:F5}), width {splitWidth / (reps * tests):F4}, lower bound {1 - 2.0 / 21:F5}");
        Assert.True(splitSummary.Mean >= 1 - 2.0 / 21 - 4 * splitSummary.Se);
    }

    private static double Best(Action action)
    {
        action();
        double best = double.PositiveInfinity;
        for (int i = 0; i < 3; i++)
        {
            long start = BassCost.Tick();
            action();
            best = Math.Min(best, (BassCost.Tick() - start) / BassCost.TicksPerSecond * 1000);
        }
        return best;
    }

    [Theory]
    [InlineData(5)]
    [InlineData(10)]
    [InlineData(1000)]
    public void Cost_FoldFitsAndClosedEndpoints(int folds)
    {
        if (!BassQualification.Enabled) return;
        const int n = 1000, tests = 200;
        var rng = new Xoshiro256PlusPlus(20261120);
        Unit[] train = Draw(n, rng), test = Draw(tests, rng);
        double[] xs = test.Select(u => u.X).ToArray();
        int[] labels = [.. Enumerable.Range(0, n).Select(i => i / (n / folds))];
        FoldReadings readings = null!;
        double fitting = Best(() => readings = CvPlusTests.Readings(train, labels, folds, xs, CrossConformalQualificationTests.Line));
        double endpoints = Best(() => CrossConformal.CvPlus(readings, CvPlusAssumptions.Both));
        CvPlusBounds bounds = CrossConformal.CvPlus(readings, CvPlusAssumptions.Both).Bounds;
        var splitCalibration = new ResidualReadings(new double[n]);
        var splitTests = new ResidualReadings(new double[tests]);
        double split = Best(() => SplitConformal.Residuals(splitCalibration, readings.Response, splitTests));
        output.WriteLine($"n={n}, K={folds}, m={n / folds}, {tests} cases, {(BassCost.Cycles ? "cycles" : "wall")}, best of three: fold fits and predictions {fitting:F3} ms; CV+ endpoints {endpoints:F3} ms ({1000 * endpoints / tests:F3} µs/case); split residual calibration/endpoints {split:F3} ms");
        output.WriteLine($"level {bounds.Level}, merged {bounds.Merged.Value:F6}, one fold {bounds.OneFold.Value:F6}, tournament {bounds.Tournament.Value:F6}, both assumption sets supplied so selected {bounds.Coverage:F6}; simple linear fitting, all input construction outside the endpoint clock");
    }
}
