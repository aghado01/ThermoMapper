using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class ConformalExactEndsTests(ITestOutputHelper output)
{
    private const double Tolerance = 1e-8;

    private static double[][] Responses(IEnumerable<double> y) => [.. y.Select(v => new[] { v })];

    private static CandidateLogDensity Gaussian(double[][] grids, int[]? strata = null)
    {
        var cdf = new double[1, grids.Sum(g => g.Length)];
        var logDensity = new double[1, cdf.Length];
        int j = 0;
        foreach (var grid in grids)
            foreach (double y in grid)
            {
                cdf[0, j] = Normal.StandardLogCdf(y);
                logDensity[0, j++] = -0.5 * y * y;
            }
        return new(logDensity, grids, [], strata) { LogCdf = cdf };
    }

    [Fact]
    public void FiniteRanksTiesAndUnboundedTails_AgreeWithHandOrderedResponses()
    {
        double[] ys = [11.25, -12.5, 0.37, 0.37, 9.5];
        var fit = SplitConformal.CalibrateDistribution(Gaussian(Responses(ys)));
        foreach (Func<double, double>? inverse in new Func<double, double>?[] { null, Normal.StandardQuantile })
        {
            var d = fit.AtCase(y => Gaussian([[y]]), predictiveQuantile: inverse);
            double[] expected = [.. ys.Order()];
            for (int j = 1; j < d.Units; j++)
                Assert.InRange(d.Quantile((double)j / d.Units), expected[j - 1] - Tolerance, expected[j - 1] + Tolerance);
            // Q̄ can just reach the end of a plateau. A one-ulp step in p moves the quantile to the next PIT.
            Assert.Equal(double.NegativeInfinity, d.Quantile(0.5 / d.Units));
            Assert.InRange(d.Quantile(3.5 / d.Units), 0.37 - Tolerance, 0.37 + Tolerance);
            Assert.InRange(d.Quantile(Math.BitIncrement(3.5 / d.Units)), 9.5 - Tolerance, 9.5 + Tolerance);
            Assert.Equal(double.PositiveInfinity, d.Quantile(Math.BitIncrement(5.5 / d.Units)));
            Assert.Equal((1.0 + 0.3 * 3) / 6, d.At(0.37, 0.3));
            var ends = d.Interval(0.4);
            Assert.InRange(ends.Lower, ys.Min() - Tolerance, ys.Min());
            Assert.InRange(ends.Upper, ys.Max(), ys.Max() + Tolerance);
        }
        var absent = fit.AtCase(y => Gaussian([[y]], [3]));
        Assert.Equal(new GridInterval(double.NegativeInfinity, double.PositiveInfinity), absent.Interval(0.9));
        Assert.Equal(double.NegativeInfinity, absent.Quantile(0.5));
        Assert.Equal(double.PositiveInfinity, absent.Quantile(Math.BitIncrement(0.5)));
        Assert.Equal(0.8, absent.At(-17, 0.8));
    }

    [Fact]
    public void CalibrationTie_OnAFlatPartOfTheTestCdf_UsesBothEndsOfThePlateau()
    {
        // Calibration PIT 0.5 three times. The test law puts equal uniform mass on [-2,-1] and [1,2].
        var fit = SplitConformal.CalibrateDistribution(Gaussian(Responses([0, 0, 0])));
        CandidateLogDensity At(double y)
        {
            double f = y < -1 ? 0.5 * Math.Clamp(y + 2, 0, 1) : 0.5 + 0.5 * Math.Clamp(y - 1, 0, 1);
            return new(new double[1, 1], [[y]], []) { LogCdf = new double[,] { { Math.Log(f) } } };
        }
        var d = fit.AtCase(At);
        Assert.Equal(0.5, d.At(0.0, 0.5));
        Assert.InRange(d.Quantile(0.5), -1 - Tolerance, -1 + Tolerance);
        Assert.InRange(d.Quantile(0.51), 1 - Tolerance, 1 + Tolerance);
        var ends = d.Interval(0.5);
        Assert.InRange(ends.Lower, -1 - Tolerance, -1);
        Assert.InRange(ends.Upper, 1, 1 + Tolerance);
    }

    [Fact]
    public void WeightedStrata_ClosedQuantilesAndEndsAgainstThePublicFineGrid()
    {
        var rng = new Xoshiro256PlusPlus(20261109);
        const double step = 0.002;
        double[][] grid = [[.. Enumerable.Range(0, 6001).Select(g => -6.0 + step * g)]];
        double worst = 0;
        for (int h = 0; h < 3; h++)
        {
            double[] beta = [-0.8, 2.2, 4.0];
            double[] weights = [1, 3, 0];
            double[] xs = [.. Enumerable.Range(0, 60).Select(_ => 2 * rng.NextDouble() - 1)];
            double[] ys = [.. xs.Select(x => 0.5 * x + 2 * rng.NextDouble() - 1)];
            CandidateLogDensity Seam(double[] x, double[][] grids, int[] strata)
            {
                var plain = ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, x, grids);
                return new(plain.Values, grids, plain.Latents, strata) { LogCdf = plain.LogCdf, Weights = weights };
            }
            var cal = Seam(xs, Responses(ys), [.. Enumerable.Range(0, xs.Length).Select(i => i % 3)]);
            CandidateLogDensity Test(double[][] g) => Seam([0.17], g, [h]);
            var d = SplitConformal.CalibrateDistribution(cal).AtCase(y => Test([[y]]));
            var gridded = SplitConformal.Distribution(cal, Test(grid)).Distributions.Single();
            foreach (double p in new[] { 0.05, 0.1, 0.5, 0.9, 0.95 })
            {
                double q = d.Quantile(p), gap = gridded.Quantile(p) - q;
                Assert.InRange(gap, -Tolerance, step + Tolerance);
                worst = Math.Max(worst, Math.Abs(gap));
            }
            GridInterval closed = d.Interval(0.2), sampled = gridded.Intervals(0.2).Single();
            Assert.InRange(sampled.Lower - closed.Lower, 0, step + Tolerance);
            Assert.InRange(closed.Upper - sampled.Upper, 0, step + Tolerance);
            worst = Math.Max(worst, Math.Max(sampled.Lower - closed.Lower, closed.Upper - sampled.Upper));
            for (int g = 0; g < grid[0].Length; g += 137)
                Assert.Equal(gridded.At(g, 0.31), d.At(grid[0][g], 0.31));
        }
        output.WriteLine($"weighted two-component predictive, three strata, grid step {step}: largest quantile or end difference {worst:G8}");
    }

    [Fact]
    public void BassPredictiveInverseAndCdfRoot_OnTheSameFit()
    {
        BassData data = BassStudentTTests.Data();
        BassFit fit = Bass.Fit(BassStudentTTests.Config(null, burnIn: 50, samples: 20, chains: 1), data);
        PredictionCases cases = PredictionCases.Of([new ContinuousPredictor("x0", [0.1, 0.4, 0.7, 0.9]),
                                                   new ContinuousPredictor("x1", [0.8, 0.3, 0.5, 0.2])]);
        CandidateLogDensity Seam(double[][] grids) => fit.CandidateLogDensity(data, cases, grids, cdf: true);
        var cal = SplitConformal.CalibrateDistribution(Seam(Responses([-0.7, 0.3, -0.2, 0.9])));
        BassPrediction pred = fit.Predict(data, cases, new BassPredictOptions { Predictive = false });
        const int t = 2;
        CandidateLogDensity Test(double y) => Seam(Responses([0, 0, y, 0]));
        int inversions = 0;
        double Inverse(double p) { inversions++; return pred.Quantile(p)[t, 0]; }
        var native = cal.AtCase(Test, t, Inverse);
        var root = cal.AtCase(Test, t);
        var grid = SplitConformal.Distribution(Seam(Responses([-0.7, 0.3, -0.2, 0.9])), Test(0.08)).Distributions[t];
        Assert.Equal(grid.At(0, 0.3), native.At(0.08, 0.3));
        foreach (double p in new[] { 0.2, 0.5, 0.8 })
            Assert.InRange(Math.Abs(native.Quantile(p) - root.Quantile(p)), 0, Tolerance);
        Assert.Equal(3, inversions);
    }

    [Fact]
    public void FormsAndAlignment_CheckedOnConstructionAndEveryNewProbe()
    {
        var cal = Gaussian(Responses([-1.0, 0.0, 1.0])) with { Forms = ConformalForms.Split };
        var fullOnly = cal with { Forms = ConformalForms.Full, FormsReason = "test reason" };
        Assert.Contains("test reason", Assert.Throws<ArgumentException>(() => SplitConformal.CalibrateDistribution(fullOnly)).Message);
        var fit = SplitConformal.CalibrateDistribution(cal);
        CandidateLogDensity Probe(double y) => Gaussian([[y]]) with { Forms = ConformalForms.Split };
        var d = fit.AtCase(y => y == 0.0 ? Probe(y) : Probe(y) with { Forms = ConformalForms.Full, FormsReason = "probe reason" });
        Assert.Contains("probe reason", Assert.Throws<ArgumentException>(() => d.Quantile(0.25)).Message);
        Assert.Throws<ArgumentException>(() => fit.AtCase(y => Probe(y) with { Weights = [1] }));
        Assert.Throws<ArgumentException>(() => fit.AtCase(y => Probe(y + 1)));
        Assert.Throws<ArgumentException>(() => fit.AtCase(y => Probe(y) with { LogCdf = null }));
        Assert.Throws<ArgumentException>(() => fit.AtCase(y => Gaussian([[y]], [y == 0 ? 0 : 1])).At(1, 0.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => fit.AtCase(Probe, tolerance: double.NaN));
        double[] weights = [1];
        var weighted = SplitConformal.CalibrateDistribution(cal with { Weights = weights });
        weights[0] = 2;
        Assert.Throws<ArgumentException>(() => weighted.AtCase(y => Probe(y) with { Weights = weights }));
    }

    [Fact]
    public void UnattainableTolerance_AndSaturatedCalibrationPit()
    {
        var fit = SplitConformal.CalibrateDistribution(Gaussian(Responses([0.37, 0.37, 0.37])));
        var tiny = fit.AtCase(y => Gaussian([[y]]), tolerance: double.Epsilon);
        Assert.Contains("adjacent doubles", Assert.Throws<InvalidOperationException>(() => tiny.Interval(0.5)).Message);
        // A calibration PIT of 0 always passes the lower comparison. The upper comparison includes the test law's
        // entire F = 0 plateau, here the negative half-line of an exponential predictive.
        var zero = Gaussian(Responses([0, 0, 0])) with { LogCdf = new double[,] { { double.NegativeInfinity, double.NegativeInfinity, double.NegativeInfinity } } };
        var d = SplitConformal.CalibrateDistribution(zero).AtCase(y => Gaussian([[y]]) with
        {
            LogCdf = new double[,] { { y <= 0 ? double.NegativeInfinity : Math.Log(1 - Math.Exp(-y)) } }
        });
        var ends = d.Interval(0.5);
        Assert.Equal(double.NegativeInfinity, ends.Lower);
        Assert.InRange(ends.Upper, 0, Tolerance);
    }

    [Theory]
    [InlineData(false, false)]
    [InlineData(true, false)]
    [InlineData(true, true)]
    public void Bisection_OnEachPublicDensityRule_ReturnsOppositeWitnesses(bool full, bool inflate)
    {
        // All draws have the same two-peaked log density -(|y|-2)², up to an additive constant.
        // At the three calibration responses |y|-2 is ±0.4, so there are four crossings and the full rule's errors vanish.
        CandidateLogDensity Seam(double[][] grids)
        {
            var values = new double[4, grids.Sum(g => g.Length)];
            int col = 0;
            foreach (var g in grids)
                foreach (double y in g)
                {
                    for (int s = 0; s < 4; s++) values[s, col] = -Math.Pow(Math.Abs(y) - 2, 2);
                    col++;
                }
            return new(values, grids, []);
        }
        var cal = Seam(Responses([-2.4, 1.6, 2.4]));
        var training = new PointwiseLogDensity(cal.Values, UnitPartition.EachObservation(3), cal.Latents);
        var options = new ConformalOptions { Alpha = 0.5, Inflate = inflate };
        ConformalSet Set(double[] grid) => full ? ConformalBayes.FullStep(training, Seam([grid]), options).Sets.Single()
                                               : SplitConformal.Density(cal, Seam([grid]), options.Alpha).Sets.Single();
        ConformalSet coarse = Set([-3, -2, 0, 2, 3]);
        Assert.Equal([new GridInterval(-2, -2), new GridInterval(2, 2)], coarse.Intervals);
        var bounds = coarse.RefineBoundaries(y => Set([y]), Tolerance);
        Assert.Equal(4, bounds.Count);
        double[] ends = [-2.4, -1.6, 1.6, 2.4];
        for (int i = 0; i < ends.Length; i++)
        {
            var (a, r) = bounds[i];
            Assert.InRange(Math.Abs(a - r), 0, Tolerance);
            Assert.InRange(ends[i], Math.Min(a, r) - 1e-14, Math.Max(a, r) + 1e-14);
            foreach ((double y, bool accepted) in new[] { (a, true), (r, false) })
            {
                ConformalSet at = Set([y]);
                double pi = inflate ? at.Margin!.InflatedRanks[0] / (double)at.Units : at.PValues![0];
                Assert.Equal(accepted, pi > at.Level / (double)at.Units);
            }
        }
        Assert.Throws<ArgumentException>(() => coarse.RefineBoundaries(y => Set([y + 1]), Tolerance));
        var open = Set([-2, 0, 2]);
        Assert.True(open.OpenBelow && open.OpenAbove);
        Assert.Equal(2, open.RefineBoundaries(y => Set([y]), Tolerance).Count);
        Assert.Empty(Set([0, 3]).RefineBoundaries(_ => throw new InvalidOperationException(), Tolerance));
        output.WriteLine($"{(full ? "full" : "split")} {(inflate ? "inflated" : "plain")}: four π crossings; widest bracket {bounds.Max(b => Math.Abs(b.Accepted - b.Rejected)):G6}");
    }

    [Fact]
    public void HeldOutPitAndClosedDeciles_AreUniformOnSkewedTruth()
    {
        if (!BassQualification.Enabled) return;
        const int m = 20, reps = 500, held = 8;
        var pits = new List<double>();
        var deciles = new int[10];
        int covered = 0;
        double width = 0;
        foreach (int seed in SeedTree.Derive(20261109, reps))
        {
            var rng = new Xoshiro256PlusPlus(seed);
            (double[] X, double[] Y) Draw(int count)
            {
                double[] x = [.. Enumerable.Range(0, count).Select(_ => 4 * rng.NextDouble() - 2)];
                return (x, [.. x.Select(v => ConformalBayesTests.Truth(v, 0.5, rng))]);
            }
            var train = Draw(30);
            var cal = Draw(m - 1);
            var test = Draw(held);
            double slope = ConformalBayesTests.ExactDraws(train.X, train.Y, 1, rng)[0];
            CandidateLogDensity Seam(double[] x, double[][] grids)
                => ConformalPredictiveDistributionTests.Seam([slope], 0.7, null, x, grids);
            var fit = SplitConformal.CalibrateDistribution(Seam(cal.X, Responses(cal.Y)));
            for (int t = 0; t < held; t++)
            {
                double x = test.X[t], y = test.Y[t];
                var d = fit.AtCase(v => Seam([x], [[v]]), predictiveQuantile: p => slope * x + 0.7 * Normal.StandardQuantile(p));
                pits.Add(d.At(y, rng.NextDouble()));
                int bin = 0;
                while (bin < 9 && y > d.Quantile((bin + 1) / 10.0)) bin++;
                deciles[bin]++;
                var ends = d.Interval(0.1);
                if (ends.Lower <= y && y <= ends.Upper) covered++;
                width += ends.Upper - ends.Lower;
            }
        }
        var (share, se, ks) = ConformalPredictiveDistributionTests.Uniformity(pits, m, reps);
        double[] bins = [.. deciles.Select(c => (double)c / pits.Count)];
        Assert.All(share.Concat(bins), v => Assert.InRange(v, 0.1 - 4 * se, 0.1 + 4 * se));
        double coverage = covered / (double)pits.Count;
        Assert.InRange(coverage, 0.9 - 4 * se, 0.9 + 4 * se);
        output.WriteLine($"{reps} independent fits, {pits.Count} held-out PITs, m {m}: bins {string.Join(" ", share.Select(v => v.ToString("F4")))}; KS {ks:F4}");
        output.WriteLine($"closed decile bins {string.Join(" ", bins.Select(v => v.ToString("F4")))}; bin SE {se:F4}; coverage {covered}/{pits.Count} = {coverage:F4}; mean width {width / pits.Count:F4}");
    }
}
