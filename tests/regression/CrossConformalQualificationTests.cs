using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;
using Unit = Maths.Regression.Tests.CrossConformalTests.Unit;

namespace Maths.Regression.Tests;

public sealed class CrossConformalQualificationTests(ITestOutputHelper output)
{
    // Sort the inputs so arithmetic and the fit are invariant under relabeling. The model omits the truth's curvature and skew.
    internal static Func<double, double> Line(Unit[] input)
    {
        Unit[] u = [.. input.OrderBy(v => v.X).ThenBy(v => v.Y)];
        double x = u.Average(v => v.X), y = u.Average(v => v.Y);
        double xx = u.Sum(v => (v.X - x) * (v.X - x));
        double slope = xx == 0 ? 0 : u.Sum(v => (v.X - x) * (v.Y - y)) / xx;
        return v => y + slope * (v - x);
    }

    private static Unit[] Draw(int n, Xoshiro256PlusPlus rng) => [.. Enumerable.Range(0, n).Select(_ =>
    {
        double x = 4 * rng.NextDouble() - 2;
        return new Unit(x, ConformalBayesTests.Truth(x, 0.5, rng));
    })];

    private static (double Mean, double Se) Summary(IReadOnlyList<double> values)
    {
        double mean = values.Average();
        double variance = values.Sum(v => (v - mean) * (v - mean)) / (values.Count - 1);
        return (mean, Math.Sqrt(variance / values.Count));
    }

    [Fact]
    public void HeldOutCoverageAndWidth_JackknifePairwiseAndSplit()
    {
        if (!BassQualification.Enabled) return;
        const int n = 39, reps = 500, tests = 5;
        const double alpha = 0.1;
        string[] names = ["jackknife+", "pairwise residual", "split residual"];
        var hits = names.Select(_ => new List<double>()).ToArray();
        var widths = new double[3];
        double[] grid = [.. Enumerable.Range(0, 1001).Select(g => -10 + 0.025 * g)];
        int clipped = 0;
        foreach (int seed in SeedTree.Derive(20261113, reps))
        {
            var rng = new Xoshiro256PlusPlus(seed);
            Unit[] train = Draw(n, rng), test = Draw(tests, rng);
            double[] xs = [.. test.Select(u => u.X)];
            LeaveOneOutReadings readings = CrossConformalTests.ResidualReadings(train, xs, Line);
            var jack = CrossConformal.JackknifePlus(readings, alpha);
            double[] own = [.. Enumerable.Range(0, n).Select(i => -Math.Abs(train[i].Y - readings.HeldOut[i]))];
            var atTruth = new double[n, tests];
            for (int i = 0; i < n; i++)
                for (int t = 0; t < tests; t++) atTruth[i, t] = -Math.Abs(test[t].Y - readings.Predictions[i, t]);
            var pairTruth = CrossConformal.Pairwise(own, atTruth, [.. test.Select(u => new[] { u.Y })], alpha);
            var atGrid = new double[n, grid.Length * tests];
            for (int i = 0; i < n; i++)
                for (int t = 0; t < tests; t++)
                    for (int g = 0; g < grid.Length; g++) atGrid[i, t * grid.Length + g] = -Math.Abs(grid[g] - readings.Predictions[i, t]);
            var pair = CrossConformal.Pairwise(own, atGrid, [.. Enumerable.Repeat(grid, tests)], alpha);
            Func<double, double> splitFit = Line(train[..20]);
            Unit[] cal = train[20..];
            var split = SplitConformal.Residuals(new(cal.Select(u => splitFit(u.X)).ToArray()),
                                                 cal.Select(u => u.Y).ToArray(), new(xs.Select(splitFit).ToArray()), alpha);
            var counts = new double[3];
            for (int t = 0; t < tests; t++)
            {
                if (jack.Intervals[t].Contains(test[t].Y)) counts[0]++;
                if (pairTruth.Sets[t].Accepts(0)) counts[1]++;
                if (split.Intervals[t].Contains(test[t].Y)) counts[2]++;
                Assert.True(!pairTruth.Sets[t].Accepts(0) || jack.Intervals[t].Contains(test[t].Y));
                if (pair.Sets[t].OpenBelow || pair.Sets[t].OpenAbove) clipped++;
                widths[0] += jack.Intervals[t].Width;
                widths[1] += pair.Sets[t].Width;
                widths[2] += split.Intervals[t].Width;
            }
            for (int i = 0; i < 3; i++) hits[i].Add(counts[i] / tests);
        }
        double proved = 1 - 7.0 / 40;
        for (int i = 0; i < 3; i++)
        {
            var (mean, se) = Summary(hits[i]);
            double bound = i == 2 ? 0.9 : proved;
            Assert.True(mean >= bound - 4 * se);
            Assert.InRange(mean, 0.9 - 4 * se, 0.9 + 4 * se + 0.025);
            output.WriteLine($"{names[i]}: {reps * tests} held-out units, coverage {mean:F5} (replication SE {se:F5}), mean width {widths[i] / (reps * tests):F4}; proved lower bound {bound:F3}");
        }
        Assert.Equal(0, clipped);
        output.WriteLine("pairwise width uses grid spacing 0.025; equal data budgets, jackknife+ uses 39 fits, split uses 20 training and 19 calibration units");
    }

    [Fact]
    public void AsymmetricCqr_EachTailUnderWrongQuantileFits()
    {
        if (!BassQualification.Enabled) return;
        const int reps = 500, tests = 5, calibration = 39;
        var lower = new List<double>();
        var upper = new List<double>();
        var covered = new List<double>();
        var symmetricUpper = new List<double>();
        double width = 0, symmetricWidth = 0;
        foreach (int seed in SeedTree.Derive(20261114, reps))
        {
            var rng = new Xoshiro256PlusPlus(seed);
            Unit[] train = Draw(30, rng), cal = Draw(calibration, rng), test = Draw(tests, rng);
            Func<double, double> f = Line(train);
            QuantileReadings Read(Unit[] u) => new(u.Select(v => f(v.X) - 3).ToArray(), u.Select(v => f(v.X) - 2).ToArray());
            double[] response = [.. cal.Select(u => u.Y)];
            var asym = SplitConformal.AsymmetricQuantiles(Read(cal), response, Read(test));
            var symmetric = SplitConformal.Quantiles(Read(cal), response, Read(test));
            double lo = 0, hi = 0, hit = 0, symHi = 0;
            for (int t = 0; t < tests; t++)
            {
                var d = asym.Intervals[t];
                if (test[t].Y < d.Lower) lo++;
                if (test[t].Y > d.Upper) hi++;
                if (d.Contains(test[t].Y)) hit++;
                if (test[t].Y > symmetric.Intervals[t].Upper) symHi++;
                width += d.Width;
                symmetricWidth += symmetric.Intervals[t].Width;
            }
            lower.Add(lo / tests); upper.Add(hi / tests); covered.Add(hit / tests); symmetricUpper.Add(symHi / tests);
        }
        var low = Summary(lower); var high = Summary(upper); var cov = Summary(covered); var sym = Summary(symmetricUpper);
        Assert.InRange(low.Mean, 0.05 - 4 * low.Se, 0.05 + 4 * low.Se);
        Assert.InRange(high.Mean, 0.05 - 4 * high.Se, 0.05 + 4 * high.Se);
        Assert.True(cov.Mean >= 0.9 - 4 * cov.Se);
        output.WriteLine($"{reps * tests} held-out units, asymmetric CQR lower misses {low.Mean:F5} (SE {low.Se:F5}), upper {high.Mean:F5} (SE {high.Se:F5}), coverage {cov.Mean:F5} (SE {cov.Se:F5}), mean width {width / (reps * tests):F4}");
        output.WriteLine($"bounds 0.05 per tail and 0.90 two-sided; symmetric CQR upper misses {sym.Mean:F5} (SE {sym.Se:F5}), mean width {symmetricWidth / (reps * tests):F4}");
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
    [InlineData(1000)]
    [InlineData(4000)]
    public void DensityCost_ThousandFitsAgainstTheSplit(int draws)
    {
        if (!BassQualification.Enabled) return;
        const int n = 1000, g = 200;
        double[] grid = [.. Enumerable.Range(0, g).Select(i => -4 + 8.0 * i / (g - 1))];
        var own = new double[draws, 1];
        var candidates = new double[draws, g];
        var cal = new double[draws, n];
        for (int s = 0; s < draws; s++)
        {
            own[s, 0] = -2 - (s % 7) / 10.0;
            for (int j = 0; j < g; j++) candidates[s, j] = -0.5 * grid[j] * grid[j] - (s % 7) / 10.0;
            for (int i = 0; i < n; i++) cal[s, i] = own[s, 0] - (i % 11) / 5.0;
        }
        var held = new CandidateLogDensity(own, [[0.0]], []);
        var test = new CandidateLogDensity(candidates, [grid], []);
        var calibrationSeam = new CandidateLogDensity(cal, [.. Enumerable.Range(0, n).Select(i => new[] { (double)i })], []);
        double cross = Best(() => CrossConformal.Density(n, _ => (held, test)));
        double split = Best(() => SplitConformal.Density(calibrationSeam, test));
        var rng = new Xoshiro256PlusPlus(20261115);
        Unit[] data = Draw(n, rng);
        double refits = Best(() => CrossConformalTests.ResidualReadings(data, [0.2], Line));
        LeaveOneOutReadings readings = CrossConformalTests.ResidualReadings(data, [0.2], Line);
        double ends = Best(() => CrossConformal.JackknifePlus(readings));
        output.WriteLine($"n={n}, S={draws}, G={g}, {(BassCost.Cycles ? "cycles" : "wall")}, best of three: n-fold predictive reductions/ranks {cross:F2} ms ({1000 * cross / g:F2} µs/candidate), split {split:F2} ms");
        output.WriteLine($"linear leave-one-out fitting and predictions {refits:F2} ms for {n} fits; jackknife+ two order statistics {ends:F4} ms; density matrices built before clock, one fold retained at a time");
    }
}
