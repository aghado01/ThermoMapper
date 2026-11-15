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
/// lose (<c>sum_weight_refused_estimated_le</c>); separated samples stay finite under the ridge; what it cannot fit is refused. Its spline
/// form (τ through richer features) recovers a log ratio no linear fit holds, covers under a shift the linear fit misses, and reads alike in
/// any units of its inputs.
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

    // The total variation between the data set's normalized true and estimated weights, the units' and the one case's.
    private static double TotalVariation(double[] units, double candidate, ShiftWeights estimated)
    {
        double st = units.Sum() + candidate, se = estimated.Units.Sum() + estimated.Candidates[0];
        double gap = Math.Max(0.0, candidate / st - estimated.Candidates[0] / se);
        for (int i = 0; i < units.Length; i++) gap += Math.Max(0.0, units[i] / st - estimated.Units[i] / se);
        return gap;
    }

    [Fact]
    public void AVarianceShift_TheSplineRatioRecoversItsLogRatio_WhereALinearOneCannot()
    {
        // Source N(0, 1), target N(0, 1.6²): ln w(x) = −ln 1.6 + (1 − 1/1.6²)x²/2, even in x, which log odds linear in x read as flat.
        var rng = new Xoshiro256PlusPlus(43);
        const int n = 4000;
        double[] source = [.. Enumerable.Range(0, n).Select(_ => Normal.SampleStandard(rng))];
        double[] target = [.. Enumerable.Range(0, n).Select(_ => 1.6 * Normal.SampleStandard(rng))];
        DensityRatio spline = DensityRatio.Spline(Column(source), Column(target)), linear = DensityRatio.Logistic(Column(source), Column(target));
        output.WriteLine($"the spline's ridge {spline.Ridge}, {spline.Iterations} Newton steps");
        double linearWorst = 0.0;
        foreach (double x in new[] { -2.0, -1.0, 0.0, 1.0, 2.0 })
        {
            double truth = -Math.Log(1.6) + (1.0 - 1.0 / 2.56) * x * x / 2.0, fitted = spline.LogRatio([x]), flat = linear.LogRatio([x]);
            output.WriteLine($"x {x}: ln w {truth:F4}; the spline's {fitted:F4}, the linear fit's {flat:F4}");
            Assert.InRange(fitted, truth - 0.2, truth + 0.2);
            linearWorst = Math.Max(linearWorst, Math.Abs(flat - truth));
        }
        Assert.True(linearWorst > 0.4, $"the linear fit's largest error {linearWorst:F4}");
        Assert.Equal(1, spline.Inputs);
        Assert.Contains(spline.Ridge, new[] { 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0 });
    }

    [Fact]
    public void ANonlinearShift_TheSplineWeightsCover_WhereTheLinearOnesDoNot()
    {
        // Training x uniform on (0, 1), the test's x triangular, density 4 min(x, 1 − x): ln w is even about ½, so log odds linear in x read
        // it as flat. The noise grows toward the middle, where the test's inputs gather, so the unweighted rule undercovers. Each ratio is
        // fitted once from 2,000 unlabeled inputs of each law and held fixed: the spline's weights cover within the total variation they
        // lose (sum_weight_refused_estimated_le), the linear fit's undercover with the unweighted rule.
        const int n = 199, replications = 4000;
        const double alpha = 0.1;
        var rng = new Xoshiro256PlusPlus(2026_1120);
        static double Sd(double v) => 0.1 + 1.6 * Math.Min(v, 1.0 - v);
        static double W(double v) => 4.0 * Math.Min(v, 1.0 - v);
        double Triangular() => 0.5 * (rng.NextDouble() + rng.NextDouble());
        DensityRatio spline = DensityRatio.Spline(Column([.. Enumerable.Range(0, 2000).Select(_ => rng.NextDouble())]),
                                                  Column([.. Enumerable.Range(0, 2000).Select(_ => Triangular())]));
        DensityRatio linear = DensityRatio.Logistic(Column([.. Enumerable.Range(0, 2000).Select(_ => rng.NextDouble())]),
                                                    Column([.. Enumerable.Range(0, 2000).Select(_ => Triangular())]));
        int plain = 0, bySpline = 0, byLinear = 0;
        double tvSpline = 0.0, tvLinear = 0.0;
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            double[] y = [.. x.Select(v => Math.Sin(2.0 * Math.PI * v) + Sd(v) * Normal.SampleStandard(rng))];
            double xt = Triangular(), yt = Math.Sin(2.0 * Math.PI * xt) + Sd(xt) * Normal.SampleStandard(rng);
            ShiftWeights splineWeights = spline.Weights(Column(x), Column([xt])), linearWeights = linear.Weights(Column(x), Column([xt]));
            double[] w = [.. x.Select(W)];
            tvSpline += TotalVariation(w, W(xt), splineWeights);
            tvLinear += TotalVariation(w, W(xt), linearWeights);
            var cal = new double[1, n];
            for (int i = 0; i < n; i++) cal[0, i] = Normal.MeanSd(Math.Sin(2.0 * Math.PI * x[i]), 0.5).LogPdf(y[i]);
            var calibration = new CandidateLogDensity(cal, [.. y.Select(v => new[] { v })], []);
            var test = new CandidateLogDensity(new double[,] { { Normal.MeanSd(Math.Sin(2.0 * Math.PI * xt), 0.5).LogPdf(yt) } }, [[yt]], []);
            if (SplitConformal.Density(calibration, test, alpha).Sets[0].Accepts(0)) plain++;
            if (SplitConformal.Density(calibration, test, alpha, splineWeights).Sets[0].Accepts(0)) bySpline++;
            if (SplitConformal.Density(calibration, test, alpha, linearWeights).Sets[0].Accepts(0)) byLinear++;
        }
        double p = plain / (double)replications, qs = bySpline / (double)replications, ql = byLinear / (double)replications;
        double se = Math.Sqrt(alpha * (1 - alpha) / replications);
        output.WriteLine($"unweighted {p:F4}; the spline's weights {qs:F4} (mean TV {tvSpline / replications:F4}, ridge {spline.Ridge}); "
                         + $"the linear fit's {ql:F4} (mean TV {tvLinear / replications:F4}); se {se:F4}");
        Assert.True(qs >= 1 - alpha - tvSpline / replications - 3.5 * se);
        Assert.True(p < 1 - alpha - 3.5 * se && ql < 1 - alpha - 3.5 * se);
        Assert.True(qs > ql + 3.5 * se && qs > p + 3.5 * se);
    }

    [Fact]
    public void TheSplineRatio_ReadsAlikeInAnyUnits_HoldsItsEdges_AndRefusesWhatItCannotFit()
    {
        // Two inputs, a Gaussian shifted and widened and a uniform squared; the second sample in units 10³ and 10⁻³ with offsets.
        var rng = new Xoshiro256PlusPlus(44);
        double[][] source = [.. Enumerable.Range(0, 600).Select(_ => new[] { Normal.SampleStandard(rng), rng.NextDouble() })];
        double[][] target = [.. Enumerable.Range(0, 500).Select(_ => new[] { 0.3 + 1.5 * Normal.SampleStandard(rng), Math.Pow(rng.NextDouble(), 2.0) })];
        static double[] Rescaled(double[] row) => [1e3 * row[0] + 7.0, 1e-3 * row[1] - 2.0];
        DensityRatio ratio = DensityRatio.Spline(To(source), To(target), knots: 4);
        DensityRatio scaled = DensityRatio.Spline(To([.. source.Select(Rescaled)]), To([.. target.Select(Rescaled)]), knots: 4);
        Assert.Equal(2, ratio.Inputs);
        Assert.Equal(ratio.Ridge, scaled.Ridge);
        foreach (double[] row in new[] { new[] { -1.5, 0.1 }, new[] { 0.0, 0.5 }, new[] { 0.7, 0.9 }, new[] { 2.5, 0.3 } })
            Assert.Equal(ratio.LogRatio(row), scaled.LogRatio(Rescaled(row)), 8);
        double top = source.Concat(target).Max(r => r[0]), bottom = source.Concat(target).Min(r => r[1]);
        Assert.Equal(ratio.LogRatio([top, 0.5]), ratio.LogRatio([top + 100.0, 0.5]), 12);
        Assert.Equal(ratio.LogRatio([0.0, bottom]), ratio.LogRatio([0.0, bottom - 3.0]), 12);
        Assert.Throws<ArgumentException>(() => ratio.LogRatio([1.0]));

        // A binary input leaves one column that varies, the rest constant and left out: it fits, the flag reading as a step.
        double[][] flagged0 = [.. source.Select((r, i) => new[] { r[0], i % 4 == 0 ? 1.0 : 0.0 })];
        double[][] flagged1 = [.. target.Select((r, i) => new[] { r[0], i % 2 == 0 ? 1.0 : 0.0 })];
        DensityRatio flagged = DensityRatio.Spline(To(flagged0), To(flagged1));
        double step = flagged.LogRatio([0.0, 1.0]) - flagged.LogRatio([0.0, 0.0]);
        output.WriteLine($"the flag's step {step:F4} against ln 3 = {Math.Log(3.0):F4}");
        Assert.InRange(step, Math.Log(3.0) - 0.4, Math.Log(3.0) + 0.4);

        double[,] s = Column([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]), t = Column([1.5, 2.5, 3.5, 4.5, 5.5, 6.5]);
        Assert.Throws<ArgumentOutOfRangeException>(() => DensityRatio.Spline(s, t, knots: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => DensityRatio.Spline(s, t, folds: 1));
        Assert.Throws<ArgumentException>(() => DensityRatio.Spline(s, t, ridges: []));
        Assert.Throws<ArgumentOutOfRangeException>(() => DensityRatio.Spline(s, t, ridges: [1.0, 0.0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => DensityRatio.Spline(s, t, ridges: [double.NaN]));
        Assert.Throws<ArgumentException>(() => DensityRatio.Spline(Column([1.0, 2.0, 3.0]), t));
        Assert.Throws<ArgumentException>(() => DensityRatio.Spline(Column([1.0, 1.0, 1.0, 1.0, 1.0]), Column([1.0, 1.0, 1.0, 1.0, 1.0])));
        Assert.Throws<ArgumentException>(() => DensityRatio.Spline(s, To([[1.0, 2.0], [2.0, 3.0], [3.0, 4.0], [4.0, 5.0], [5.0, 6.0]])));
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
