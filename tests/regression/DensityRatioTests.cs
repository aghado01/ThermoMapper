using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The ratio estimator (the conformal layer's slice 14 follow-up): a logistic regression telling the two samples apart recovers ln w where it
/// is linear in the features, against the laws' own ratio; its weights give the weighted rule coverage within the total variation they
/// lose (<c>sum_weight_refused_estimated_le</c>); separated samples stay finite under the ridge; what it cannot fit is refused.
/// </summary>
public sealed class DensityRatioTests(ITestOutputHelper output)
{
    private static double[,] Column(double[] v) => To(v.Select(x => new[] { x }).ToArray());

    private static double[,] To(double[][] rows)
    {
        var m = new double[rows.Length, rows[0].Length];
        for (int i = 0; i < rows.Length; i++)
            for (int j = 0; j < rows[0].Length; j++) m[i, j] = rows[i][j];
        return m;
    }

    [Fact]
    public void AGaussianMeanShift_ItsLogRatioRecovered_AnIdleFeatureIgnored()
    {
        // Source N(0, 1), target N(1, 1): ln w(x) = x − ½. A second feature, N(0, 1) in both, carries nothing.
        var rng = new Xoshiro256PlusPlus(41);
        const int n = 4000;
        double[][] source = [.. Enumerable.Range(0, n).Select(_ => new[] { Normal.SampleStandard(rng), Normal.SampleStandard(rng) })];
        double[][] target = [.. Enumerable.Range(0, n).Select(_ => new[] { 1.0 + Normal.SampleStandard(rng), Normal.SampleStandard(rng) })];
        DensityRatio ratio = DensityRatio.Logistic(To(source), To(target));
        foreach (double x in new[] { -1.0, 0.0, 1.0, 2.0 })
        {
            double fitted = ratio.LogRatio([x, 0.0]);
            output.WriteLine($"x {x}: ln w {fitted:F4} against {x - 0.5:F4}");
            Assert.InRange(fitted, x - 0.5 - 0.15, x - 0.5 + 0.15);
            Assert.InRange(ratio.LogRatio([x, 1.0]) - fitted, -0.1, 0.1);
        }
        Assert.InRange(ratio.Iterations, 1, 20);
    }

    [Fact]
    public void APowerShift_WithLogFeatures_ItsLogRatioRecovered()
    {
        // Source uniform on (0, 1), target with density 2x: ln w(x) = ln 2 + ln x, linear in ln x.
        var rng = new Xoshiro256PlusPlus(42);
        const int n = 4000;
        double[] source = [.. Enumerable.Range(0, n).Select(_ => Math.Log(1.0 - rng.NextDouble()))];
        double[] target = [.. Enumerable.Range(0, n).Select(_ => 0.5 * Math.Log(1.0 - rng.NextDouble()))];
        DensityRatio ratio = DensityRatio.Logistic(Column(source), Column(target));
        foreach (double x in new[] { 0.1, 0.5, 0.9 })
        {
            double fitted = ratio.LogRatio([Math.Log(x)]), truth = Math.Log(2.0 * x);
            output.WriteLine($"x {x}: ln w {fitted:F4} against {truth:F4}");
            Assert.InRange(fitted, truth - 0.12, truth + 0.12);
        }
    }

    [Fact]
    public void TheEstimatedWeights_CoverWithinTheTotalVariationTheyLose()
    {
        // Slice 14's covariate shift (training x uniform, the test's density 2x, noise growing with x, a one-scale fit), the ratio estimated
        // in each replication from 400 unlabeled inputs of each law, on ln x: coverage at least 1 − α less the mean TV of the data set's
        // normalized true and estimated weights.
        const int n = 199, replications = 2000;
        const double alpha = 0.1;
        var rng = new Xoshiro256PlusPlus(2026_1117);
        int plain = 0, estimated = 0;
        double tv = 0.0;
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            double[] y = [.. x.Select(v => Math.Sin(2.0 * Math.PI * v) + (0.1 + 0.6 * v) * Normal.SampleStandard(rng))];
            double xt = Math.Sqrt(rng.NextDouble()), yt = Math.Sin(2.0 * Math.PI * xt) + (0.1 + 0.6 * xt) * Normal.SampleStandard(rng);
            double[] sourceSample = [.. Enumerable.Range(0, 400).Select(_ => Math.Log(1.0 - rng.NextDouble()))];
            double[] targetSample = [.. Enumerable.Range(0, 400).Select(_ => 0.5 * Math.Log(1.0 - rng.NextDouble()))];
            DensityRatio ratio = DensityRatio.Logistic(Column(sourceSample), Column(targetSample));
            ShiftWeights weights = ratio.Weights(Column([.. x.Select(v => Math.Log(v))]), Column([Math.Log(xt)]));
            double[] w = [.. x.Select(v => 2.0 * v)];
            double wt = 2.0 * xt, st = w.Sum() + wt, se = weights.Units.Sum() + weights.Candidates[0];
            double gap = Math.Max(0.0, wt / st - weights.Candidates[0] / se);
            for (int i = 0; i < n; i++) gap += Math.Max(0.0, w[i] / st - weights.Units[i] / se);
            tv += gap;
            var cal = new double[1, n];
            for (int i = 0; i < n; i++) cal[0, i] = Normal.MeanSd(Math.Sin(2.0 * Math.PI * x[i]), 0.35).LogPdf(y[i]);
            var calibration = new CandidateLogDensity(cal, [.. y.Select(v => new[] { v })], []);
            var test = new CandidateLogDensity(new double[,] { { Normal.MeanSd(Math.Sin(2.0 * Math.PI * xt), 0.35).LogPdf(yt) } }, [[yt]], []);
            if (SplitConformal.Density(calibration, test, alpha).Sets[0].Accepts(0)) plain++;
            if (SplitConformal.Density(calibration, test, alpha, weights).Sets[0].Accepts(0)) estimated++;
        }
        double p = plain / (double)replications, q = estimated / (double)replications, meanTv = tv / replications;
        double se1 = Math.Sqrt(alpha * (1 - alpha) / replications);
        output.WriteLine($"unweighted {p:F4}, the estimated ratio {q:F4} (mean TV {meanTv:F4}, se {se1:F4})");
        Assert.True(q >= 1 - alpha - meanTv - 3.5 * se1);
        Assert.True(q > p);
    }

    [Fact]
    public void ASeededSweep_ShiftsSizesAndScales_ConvergesInvariantAndRecovers()
    {
        // 60 trials of a Gaussian mean shift δ in [−3, 3] with unequal samples of 50 to 3,000 points, the feature in units scaled by 10⁻³ to
        // 10³: the fit converges, its ratio reads alike in any units of the feature, and recovers ln w(x) = δx − δ²/2 within six of its
        // standard errors at the source's quartiles and median.
        var rng = new Xoshiro256PlusPlus(2026_1118);
        double worst = 0.0;
        for (int trial = 0; trial < 60; trial++)
        {
            double delta = -3.0 + 6.0 * rng.NextDouble(), units = Math.Pow(10.0, -3.0 + 6.0 * rng.NextDouble());
            int n0 = 50 + (int)(2950 * rng.NextDouble()), n1 = 50 + (int)(2950 * rng.NextDouble());
            double[] source = [.. Enumerable.Range(0, n0).Select(_ => Normal.SampleStandard(rng))];
            double[] target = [.. Enumerable.Range(0, n1).Select(_ => delta + Normal.SampleStandard(rng))];
            DensityRatio plain = DensityRatio.Logistic(Column(source), Column(target));
            DensityRatio scaled = DensityRatio.Logistic(Column([.. source.Select(v => v * units)]), Column([.. target.Select(v => v * units)]));
            Assert.InRange(plain.Iterations, 1, 60);
            double sd = Math.Sqrt(1.0 / n0 + 1.0 / n1) * (1.0 + Math.Abs(delta)) * Math.Exp(0.25 * delta * delta);
            foreach (double x in new[] { -0.674, 0.0, 0.674 })
            {
                double fitted = plain.LogRatio([x]);
                Assert.Equal(fitted, scaled.LogRatio([x * units]), 8);
                double error = Math.Abs(fitted - (delta * x - 0.5 * delta * delta));
                worst = Math.Max(worst, error / sd);
                Assert.True(error <= 6.0 * sd, $"trial {trial}: δ {delta:F2}, n {n0}/{n1}, x {x}: error {error:F4} against {sd:F4}");
            }
        }
        output.WriteLine($"the largest error in its standard-error units: {worst:F2}");
    }

    [Fact]
    public void SeparatedSamples_StayFinite_AndWhatCannotBeFit_IsRefused()
    {
        DensityRatio ratio = DensityRatio.Logistic(Column([-3.0, -2.0, -1.5, -1.0]), Column([1.0, 1.5, 2.0, 3.0]));
        Assert.True(double.IsFinite(ratio.LogRatio([2.0])) && ratio.LogRatio([2.0]) > ratio.LogRatio([-2.0]));
        Assert.Throws<ArgumentException>(() => DensityRatio.Logistic(Column([1.0, 1.0]), Column([1.0, 1.0])));
        Assert.Throws<ArgumentException>(() => DensityRatio.Logistic(new double[0, 1], Column([1.0])));
        Assert.Throws<ArgumentException>(() => DensityRatio.Logistic(Column([1.0, 2.0]), To([[1.0, 2.0]])));
        Assert.Throws<ArgumentOutOfRangeException>(() => DensityRatio.Logistic(Column([1.0, 2.0]), Column([1.5, 2.5]), ridge: 0.0));
        Assert.Throws<ArgumentException>(() => ratio.LogRatio([1.0, 2.0]));
    }
}
