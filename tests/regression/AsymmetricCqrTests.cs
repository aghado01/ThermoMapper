using System;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class AsymmetricCqrTests(ITestOutputHelper output)
{
    [Fact]
    public void CqrAsymCoverage_ExactTailsOverEveryOrdering()
    {
        // cqrAsym_coverage_half, cqrLo_isUpperSet, cqrHi_isLowerSet. Fits deliberately cross and miss the truth.
        double[] y = [-2, -0.5, 1, 3, 6];
        double[] lo = [3, -2, 4, -1, 1], hi = [-4, 4, 0, 8, 5.5];
        int cases = 0;
        foreach (double alpha in new[] { 0.1, 0.4, 0.8 })
        {
            int lowerMiss = 0, upperMiss = 0, covered = 0, total = 0;
            foreach (int[] p in CrossConformalTests.Permutations(y.Length))
            {
                int[] c = p[..^1];
                int t = p[^1];
                var result = SplitConformal.AsymmetricQuantiles(new(c.Select(i => lo[i]).ToArray(), c.Select(i => hi[i]).ToArray()),
                                                                c.Select(i => y[i]).ToArray(), new([lo[t]], [hi[t]]), alpha);
                var d = result.Intervals[0];
                if (y[t] < d.Lower) lowerMiss++;
                if (y[t] > d.Upper) upperMiss++;
                if (d.Contains(y[t])) covered++;
                Assert.Equal(1 - (d.LowerLevel + d.UpperLevel) / (double)y.Length, d.Coverage);
                total++;
                cases++;
            }
            int k = (int)Math.Floor(alpha / 2 * y.Length);
            Assert.Equal(total * k, lowerMiss * y.Length);
            Assert.Equal(total * k, upperMiss * y.Length);
            Assert.True(covered * y.Length >= total * (y.Length - 2 * k));
            Assert.True(covered / (double)total >= 1 - alpha);
        }
        output.WriteLine($"{cases} orderings with distinct scores: each tail misses exactly k/m and total coverage meets 1-2k/m");
    }

    [Fact]
    public void CqrAsymOrdConnected_ClosedEndsStrataTiesAndFarRanges()
    {
        var rng = new Xoshiro256PlusPlus(20261112);
        int probes = 0, empty = 0;
        for (int trial = 0; trial < 1000; trial++)
        {
            int n = rng.NextInt(40);
            double scale = Math.ScaleB(1, rng.NextInt(1001) - 500);
            double[] Values(int count) => [.. Enumerable.Range(0, count).Select(_ => scale * (rng.NextInt(13) - 6))];
            double[] y = Values(n), lo = Values(n), hi = Values(n);
            int[] strata = [.. Enumerable.Range(0, n).Select(_ => rng.NextInt(3))];
            var cal = new QuantileReadings(lo, hi, strata);
            var test = new QuantileReadings(Values(4), Values(4), [0, 1, 2, 3]);
            double alpha = 0.05 + 0.85 * rng.NextDouble();
            var result = SplitConformal.AsymmetricQuantiles(cal, y, test, alpha);
            foreach (var d in result.Intervals)
            {
                if (d.Empty) { empty++; Assert.Equal(0, d.Width); }
                int t = d.Case, m = 1 + strata.Count(h => h == d.Stratum), k = (int)Math.Floor(alpha / 2 * m);
                Assert.Equal(k, d.LowerLevel);
                Assert.Equal(k, d.UpperLevel);
                foreach (double v in new[] { -100 * scale, -scale, 0, scale, 100 * scale, d.Lower, d.Upper }.Where(double.IsFinite))
                {
                    int lowerAtMost = Enumerable.Range(0, n).Count(i => strata[i] == d.Stratum && y[i] - lo[i] <= v - test.Lower[t]);
                    int upperAtMost = Enumerable.Range(0, n).Count(i => strata[i] == d.Stratum && hi[i] - y[i] <= test.Upper[t] - v);
                    Assert.Equal(lowerAtMost + 1 > k && upperAtMost + 1 > k, d.Contains(v));
                    probes++;
                }
            }
        }
        Assert.True(empty > 0);
        output.WriteLine($"{probes} count comparisons, {empty} empty intervals, ties, unseen strata, and dyadic scales across 1e-150 to 1e150");
    }

    [Fact]
    public void AcceptsZero_SmallCalibrationAndFiniteRounding()
    {
        foreach (int n in new[] { 0, 1, 18, 19, 29 })
        {
            double[] y = [.. Enumerable.Range(0, n).Select(i => (double)i)];
            var result = SplitConformal.AsymmetricQuantiles(new(new double[n], new double[n]), y, new([0.0], [0.0]));
            var d = result.Intervals[0];
            if (n < 19)
            {
                Assert.True(d.WholeSpace);
                Assert.Equal(double.NegativeInfinity, d.Lower);
                Assert.Equal(double.PositiveInfinity, d.Upper);
                Assert.Equal(1, d.Coverage);
            }
            else
            {
                Assert.Equal(1, d.LowerLevel);
                Assert.Equal(0, d.Lower);
                Assert.Equal(n - 1, d.Upper);
                Assert.True(d.Contains(d.Lower) && d.Contains(d.Upper));
                Assert.Equal(1 - 2.0 / (n + 1), d.Coverage);
            }
        }
        Assert.Throws<ArgumentException>(() => SplitConformal.AsymmetricQuantiles(new([0], [0, 1]), [0], new([0], [0])));
        Assert.Throws<ArgumentOutOfRangeException>(() => SplitConformal.AsymmetricQuantiles(new([0], [0]), [0], new([0], [0]), 0));
        Assert.True(SplitConformal.AsymmetricQuantiles(new([0], [0]), [0], new([0], [0]), double.Epsilon).Intervals[0].WholeSpace);
        Assert.Throws<ArgumentException>(() => SplitConformal.AsymmetricQuantiles(
            new([-double.MaxValue], [double.MaxValue]), [double.MaxValue], new([0], [0])));
    }
}
