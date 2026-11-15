using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// A split set's own miscoverage bound (E16, slice 18's follow-up; Vovk 2012, Proposition 2a): given its calibration, the split interval's
/// miscoverage is Beta(k, m − k), exactly with continuous scores, so its β-quantile is exceeded with chance 1 − β; with tied scores the law
/// is conservative; per stratum it reads the stratum's m and k; the asymmetric interval's union bound holds; the quantile is the binomial
/// tail's root; the full step, the cross-fitted forms and the weighted rule report none.
/// </summary>
public sealed class TrainingConditionalTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);
    private static double Phi(double z) => Math.Exp(Normal.StandardLogCdf(z));

    [Fact]
    public void TheSplitIntervalsMiscoverage_IsBetaKMMinusK_GivenItsCalibration()
    {
        // 99 calibration residuals of N(0, 1) about the true centre, α 0.1: m 100, k 10. Each replication's interval [−c, c] misses a new unit
        // with chance 2Φ(−c) given its calibration; over 4,000 replications that miscoverage passes Beta(10, 90)'s β-quantile with chance 1 − β.
        const int n = 99, replications = 4000;
        var rng = new Xoshiro256PlusPlus(2026_1139);
        double[] betas = [0.5, 0.9, 0.99];
        var above = new int[betas.Length];
        double mean = 0.0;
        double[]? bounds = null;
        var centre = new ResidualReadings(new double[n]);
        var at = new ResidualReadings([0.0]);
        for (int r = 0; r < replications; r++)
        {
            double[] y = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))];
            ConformalInterval interval = SplitConformal.Residuals(centre, y, at, 0.1).Intervals[0];
            double miss = Phi(interval.Lower) + Phi(-interval.Upper);
            mean += miss / replications;
            bounds ??= [.. betas.Select(b => interval.TrainingConditionalMiscoverage(b)!.Value)];
            for (int j = 0; j < betas.Length; j++) if (miss > bounds[j]) above[j]++;
        }
        double sd = Math.Sqrt(10.0 * 90.0 / (100.0 * 100.0 * 101.0));
        output.WriteLine($"mean miscoverage {mean:F5} (k/m 0.1, se {sd / Math.Sqrt(replications):F5}); bounds "
                         + string.Join(", ", betas.Select((b, j) => $"β {b}: {bounds![j]:F4}, passed {above[j] / (double)replications:F4} (1 − β {1 - b:F2})")));
        Assert.True(Math.Abs(mean - 0.1) <= 3.5 * sd / Math.Sqrt(replications), $"mean {mean:F5}");
        for (int j = 0; j < betas.Length; j++)
        {
            double share = above[j] / (double)replications, se = Math.Sqrt(betas[j] * (1 - betas[j]) / replications);
            Assert.True(Math.Abs(share - (1 - betas[j])) <= 3.5 * se, $"β {betas[j]}: passed {share:F4}");
        }
    }

    [Fact]
    public void TiedScores_AreConservative_AndEachStratumReadsItsOwnMAndK()
    {
        // Responses rounded to integers, Y = round(1.5Z), so the residuals tie; strata of 49 and 199 calibration units, m 50 and 200, k 5 and
        // 20. Given its calibration the interval [−c, c] misses with chance P(|Y| > c) = 2Φ(−(c + ½)/1.5).
        const int replications = 4000;
        const double beta = 0.9;
        var rng = new Xoshiro256PlusPlus(2026_1140);
        int[] strata = [.. Enumerable.Repeat(0, 49), .. Enumerable.Repeat(1, 199)];
        var centre = new ResidualReadings(new double[strata.Length], null, strata);
        var at = new ResidualReadings([0.0, 0.0], null, [0, 1]);
        var above = new int[2];
        double[]? bounds = null;
        for (int r = 0; r < replications; r++)
        {
            double[] y = [.. strata.Select(_ => Math.Round(1.5 * Gaussian(rng)))];
            ConformalIntervals read = SplitConformal.Residuals(centre, y, at, 0.1);
            bounds ??= [.. read.Intervals.Select(i => i.TrainingConditionalMiscoverage(beta)!.Value)];
            for (int h = 0; h < 2; h++)
            {
                ConformalInterval interval = read.Intervals[h];
                if (r == 0)
                {
                    Assert.Equal((h == 0 ? 50 : 200, h == 0 ? 5 : 20), (interval.Units, interval.Level));
                    Assert.Equal(interval.TrainingConditionalMiscoverage(beta), read.Calibration.TrainingConditionalMiscoverage(h, beta));
                }
                double miss = 2.0 * Phi(-(interval.Upper + 0.5) / 1.5);
                if (miss > bounds[h]) above[h]++;
            }
        }
        double se = Math.Sqrt(beta * (1 - beta) / replications);
        output.WriteLine($"tied scores at β {beta}: bounds {bounds![0]:F4} and {bounds[1]:F4}, passed {above[0] / (double)replications:F4} and "
                         + $"{above[1] / (double)replications:F4} (1 − β {1 - beta:F2}, se {se:F4})");
        for (int h = 0; h < 2; h++) Assert.True(above[h] / (double)replications <= 1 - beta + 3.5 * se, $"stratum {h}");
    }

    [Fact]
    public void TheAsymmetricInterval_TakesTheTailsUnionBound()
    {
        // CQR's two tails at α/2 on quantile readings of 0: the interval runs from the 5th smallest to the 5th largest of 99 responses N(0, 1),
        // missing with chance Φ(L) + Φ(−U) given its calibration; the bound sums each tail's Beta(5, 95) quantile at (1 + β)/2.
        const int n = 99, replications = 4000;
        const double beta = 0.9;
        var rng = new Xoshiro256PlusPlus(2026_1141);
        var calibration = new QuantileReadings(new double[n], new double[n]);
        var test = new QuantileReadings([0.0], [0.0]);
        int above = 0;
        double bound = double.NaN;
        for (int r = 0; r < replications; r++)
        {
            double[] y = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))];
            AsymmetricConformalInterval interval = SplitConformal.AsymmetricQuantiles(calibration, y, test, 0.1).Intervals[0];
            bound = interval.TrainingConditionalMiscoverage(beta)!.Value;
            if (Phi(interval.Lower) + Phi(-interval.Upper) > bound) above++;
        }
        double each = TrainingConditional.Miscoverage(100, 5, (1 + beta) / 2), se = Math.Sqrt(beta * (1 - beta) / replications);
        output.WriteLine($"asymmetric at β {beta}: bound {bound:F4} (each tail {each:F4}), passed {above / (double)replications:F4} (1 − β {1 - beta:F2})");
        Assert.Equal(2 * each, bound, 15);
        Assert.True(above / (double)replications <= 1 - beta + 3.5 * se);
    }

    [Fact]
    public void TheEquiTailedSet_IsBetaTwoKMMinusTwoK_ItsTailsOneScore()
    {
        // The split predictive distribution of the right fit, y = 0.5x + N(0, 1): 99 calibration units, α 0.2, k 10 per tail. Both tails read
        // the PIT, uniform at any x, so given the calibration the set misses a new unit with chance u₍ₖ₎ + 1 − u₍ₘ₋ₖ₎, Beta(20, 80); read
        // here at x = 0.4 from the set's grid-free ends.
        const int n = 99, replications = 2000;
        const double alpha = 0.2, beta = 0.9;
        var rng = new Xoshiro256PlusPlus(2026_1143);
        double[] slope = [.. Enumerable.Repeat(0.5, 10)];
        int above = 0;
        double bound = double.NaN;
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
            SplitPitCalibration fit = SplitConformal.CalibrateDistribution(
                ConformalPredictiveDistributionTests.Seam(slope, 1.0, null, x, [.. y.Select(v => new[] { v })]));
            SplitPredictiveDistribution at = fit.AtCase(v => ConformalPredictiveDistributionTests.Seam(slope, 1.0, null, [0.4], [[v]]));
            GridInterval set = at.Interval(alpha);
            bound = at.TrainingConditionalMiscoverage(alpha, beta)!.Value;
            if (Phi(set.Lower - 0.2) + Phi(0.2 - set.Upper) > bound) above++;
        }
        double se = Math.Sqrt(beta * (1 - beta) / replications);
        output.WriteLine($"equi-tailed at α {alpha}, β {beta}: bound {bound:F4} (Beta(20, 80)), passed {above / (double)replications:F4} (1 − β {1 - beta:F2}, se {se:F4})");
        Assert.Equal(TrainingConditional.Miscoverage(100, 20, beta), bound);
        Assert.True(Math.Abs(above / (double)replications - (1 - beta)) <= 3.5 * se);

        // The grid form reads the same; the full step and a weighed case read none.
        double[] grid = [.. Enumerable.Range(0, 41).Select(g => -4.0 + 0.2 * g)];
        double[] cx = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))], cy = [.. cx.Select(v => 0.5 * v + Gaussian(rng))];
        CandidateLogDensity calibration = ConformalPredictiveDistributionTests.Seam(slope, 1.0, null, cx, [.. cy.Select(v => new[] { v })]);
        CandidateLogDensity test = ConformalPredictiveDistributionTests.Seam(slope, 1.0, null, [0.4], [grid]);
        Assert.Equal(bound, SplitConformal.Distribution(calibration, test).Distributions[0].TrainingConditionalMiscoverage(alpha, beta));
        PointwiseLogDensity training = ConformalPredictiveDistributionTests.Training(slope, 1.0, null, cx, cy);
        Assert.Null(ConformalBayes.Distribution(training, test).Distributions[0].TrainingConditionalMiscoverage(alpha, beta));
        SplitPitCalibration pits = SplitConformal.CalibrateDistribution(calibration);
        Assert.Null(pits.OnGrid(test, ShiftWeights.ByCase([.. cx.Select(v => 1.0 + v * v)], [2.0], [grid])).Distributions[0].TrainingConditionalMiscoverage(alpha, beta));
    }

    [Fact]
    public void TheQuantile_IsTheBinomialTailsRoot_AndOnlySplitSetsReportIt()
    {
        // Beta(k, m − k)'s CDF at q is P(Bin(m − 1, q) ≥ k), summed here term by term.
        static double Tail(int trials, int k, double q)
        {
            double total = 0.0;
            for (int j = k; j <= trials; j++)
            {
                double logChoose = 0.0;
                for (int i = 1; i <= j; i++) logChoose += Math.Log((trials - j + i) / (double)i);
                total += Math.Exp(logChoose + j * Math.Log(q) + (trials - j) * Math.Log(1.0 - q));
            }
            return total;
        }
        foreach ((int m, int k) in new[] { (100, 10), (1000, 1), (50, 49), (31, 3) })
            foreach (double beta in new[] { 0.5, 0.9, 0.999 })
            {
                double q = TrainingConditional.Miscoverage(m, k, beta);
                Assert.Equal(beta, Tail(m - 1, k, q), 9);
            }
        Assert.Equal(0.0, TrainingConditional.Miscoverage(9, 0, 0.9));
        Assert.Throws<ArgumentOutOfRangeException>(() => TrainingConditional.Miscoverage(0, 0, 0.9));
        Assert.Throws<ArgumentOutOfRangeException>(() => TrainingConditional.Miscoverage(10, 10, 0.9));
        Assert.Throws<ArgumentOutOfRangeException>(() => TrainingConditional.Miscoverage(10, -1, 0.9));
        foreach (double bad in new[] { 0.0, 1.0, double.NaN })
            Assert.Throws<ArgumentOutOfRangeException>(() => TrainingConditional.Miscoverage(10, 1, bad));

        // The split density form and the group form report it at their m and k, a small stratum 0; the weighted rule, the full step: none.
        var rng = new Xoshiro256PlusPlus(2026_1142);
        double[] x = [.. Enumerable.Range(0, 30).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
        double[] slope = [0.45, 0.5, 0.55];
        double[] grid = [.. Enumerable.Range(0, 21).Select(g => -3.0 + 0.3 * g)];
        CandidateLogDensity calibration = ConformalPredictiveDistributionTests.Seam(slope, 1.0, null, x, [.. y.Select(v => new[] { v })]);
        CandidateLogDensity candidates = ConformalPredictiveDistributionTests.Seam(slope, 1.0, null, [0.2], [grid]);
        ConformalSet split = SplitConformal.Density(calibration, candidates).Sets[0];
        Assert.Equal((31, 3), (split.Units, split.Level));
        Assert.Equal(TrainingConditional.Miscoverage(31, 3, 0.9), split.TrainingConditionalMiscoverage(0.9));
        CandidateLogDensity few = ConformalPredictiveDistributionTests.Seam(slope, 1.0, null, x[..5], [.. y[..5].Select(v => new[] { v })]);
        Assert.Equal(0.0, SplitConformal.Density(few, candidates).Sets[0].TrainingConditionalMiscoverage(0.9));
        ConformalSet group = GroupConformal.Density(calibration, [.. Enumerable.Range(0, 30).Select(i => i / 2)], candidates).Sets[0];
        Assert.Equal(TrainingConditional.Miscoverage(16, 1, 0.9), group.TrainingConditionalMiscoverage(0.9));
        ConformalSet weighted = SplitConformal.Density(calibration, candidates, 0.1, ShiftWeights.ByCase([.. x.Select(v => 1.0 + v * v)], [2.0], [grid])).Sets[0];
        Assert.NotNull(weighted.WeightedAlpha);
        Assert.Null(weighted.TrainingConditionalMiscoverage(0.9));
        PointwiseLogDensity training = ConformalPredictiveDistributionTests.Training(slope, 1.0, null, x, y);
        Assert.Null(ConformalBayes.FullStep(training, candidates).Sets[0].TrainingConditionalMiscoverage(0.9));
    }
}
