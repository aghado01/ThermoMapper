using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The weighted rule under shift (E16, the conformal layer's slice 14): refused units carry weight at most α (<c>sum_weight_refused_le</c>);
/// equal weights give CB-1's sets bit for bit in every form (<c>weighted_eq_cb1</c>); under a known covariate shift and a response shift
/// the true ratio covers where the unweighted sets fall short (<c>coverage_weighted</c>), and an estimated ratio loses at most its total
/// variation (<c>sum_weight_refused_estimated_le</c>). With the qualification switch (<c>BASS_QUALIFICATION</c> = 1): the full step under
/// covariate shift on exact draws, and a drifting stream under weights by position (<c>refused_le_swap_tv</c>).
/// </summary>
public sealed class WeightedConformalTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    [Fact]
    public void TheWeightedLemma_RefusedUnitsCarryWeightAtMostAlpha()
    {
        // 3,000 random data sets with ties and some zero weights: each unit in turn is the candidate among the rest, and the refused units'
        // normalized weight is at most α; a score is accepted exactly when it reaches the threshold at its weight.
        var rng = new Xoshiro256PlusPlus(14);
        double worst = 0.0;
        for (int r = 0; r < 3000; r++)
        {
            int m = 2 + (int)(rng.NextDouble() * 12);
            double alpha = 0.05 + 0.4 * rng.NextDouble();
            double[] s = [.. Enumerable.Range(0, m).Select(_ => Math.Round(3.0 * rng.NextDouble(), 1))];
            double[] w = [.. Enumerable.Range(0, m).Select(_ => rng.NextDouble() < 0.15 ? 0.0 : rng.NextDouble())];
            double total = w.Sum();
            if (!(total > 0.0)) continue;
            double refused = 0.0;
            for (int u = 0; u < m; u++)
            {
                int[] rest = [.. Enumerable.Range(0, m).Where(j => j != u)];
                var fit = new WeightedCalibration([.. rest.Select(j => s[j])], [.. rest.Select(j => w[j])], null, alpha);
                double p = fit.PValue(0, s[u], w[u]);
                bool accepted = p > alpha;
                Assert.Equal(accepted, s[u] >= fit.Threshold(0, w[u]));
                if (!accepted) refused += w[u] / total;
            }
            worst = Math.Max(worst, refused - alpha);
            Assert.True(refused <= alpha + 1e-12, $"refused weight {refused} at α {alpha}");
        }
        output.WriteLine($"the refused weight's largest excess over α: {worst:F4}");
    }

    // ---- Equal weights: CB-1's sets, bit for bit ----

    private static (CandidateLogDensity Calibration, CandidateLogDensity Candidates, double[] Response, double[] Mu) Seams(int n, int cases, double[] grid, Xoshiro256PlusPlus rng)
    {
        double[] mu = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))];
        double[] y = [.. mu.Select(v => v + 0.5 * Gaussian(rng))];
        var cal = new double[3, n];
        var cand = new double[3, cases * grid.Length];
        double[] tmu = [.. Enumerable.Range(0, cases).Select(_ => Gaussian(rng))];
        for (int d = 0; d < 3; d++)
        {
            for (int i = 0; i < n; i++) cal[d, i] = Normal.MeanSd(mu[i] + 0.1 * d, 0.6).LogPdf(y[i]);
            for (int t = 0; t < cases; t++)
                for (int g = 0; g < grid.Length; g++) cand[d, t * grid.Length + g] = Normal.MeanSd(tmu[t] + 0.1 * d, 0.6).LogPdf(grid[g]);
        }
        return (new CandidateLogDensity(cal, [.. y.Select(v => new[] { v })], []),
                new CandidateLogDensity(cand, [.. Enumerable.Repeat(grid, cases)], []), y, mu);
    }

    [Fact]
    public void EqualWeights_GiveCb1sSets_BitForBit_InEveryForm()
    {
        var rng = new Xoshiro256PlusPlus(3);
        double[] grid = [.. Enumerable.Range(0, 41).Select(g => -3.0 + 0.15 * g)];
        (CandidateLogDensity cal, CandidateLogDensity cand, double[] y, double[] mu) = Seams(39, 4, grid, rng);
        var equal = ShiftWeights.ByCase([.. Enumerable.Repeat(0.37, 39)], [0.37, 0.37, 0.37, 0.37], cand.Grids);
        ConformalPrediction a = SplitConformal.Density(cal, cand), b = SplitConformal.Density(cal, cand, weights: equal);
        for (int t = 0; t < 4; t++)
        {
            Assert.Null(b.Sets[t].WeightedAlpha);
            Assert.Equal(a.Sets[t].PValues, b.Sets[t].PValues);
        }
        // One case weighed otherwise: only it takes the weighted rule.
        var one = ShiftWeights.ByCase([.. Enumerable.Repeat(0.37, 39)], [0.37, 0.9, 0.37, 0.37], cand.Grids);
        ConformalPrediction c = SplitConformal.Density(cal, cand, weights: one);
        Assert.Equal(a.Sets[0].PValues, c.Sets[0].PValues);
        Assert.Equal(0.1, c.Sets[1].WeightedAlpha);

        var quantiles = new QuantileReadings([.. mu.Select(v => v - 0.6)], [.. mu.Select(v => v + 0.6)]);
        var test = new QuantileReadings([-0.6, 0.0, 0.4], [0.6, 1.2, 1.4]);
        var flat = new ShiftWeights([.. Enumerable.Repeat(2.0, 39)], [2.0, 2.0, 2.0]);
        Assert.Equal(SplitConformal.Quantiles(quantiles, y, test).Intervals, SplitConformal.Quantiles(quantiles, y, test, weights: flat).Intervals);
        Assert.Equal(SplitConformal.AsymmetricQuantiles(quantiles, y, test).Intervals, SplitConformal.AsymmetricQuantiles(quantiles, y, test, weights: flat).Intervals);
        var residuals = new ResidualReadings(mu, [.. mu.Select(_ => 0.7)]);
        var at = new ResidualReadings([0.0, 0.3, -1.0], [0.7, 0.7, 0.7]);
        Assert.Equal(SplitConformal.Residuals(residuals, y, at).Intervals, SplitConformal.Residuals(residuals, y, at, weights: flat).Intervals);

        // The full step on exact draws.
        double[] x = [.. Enumerable.Range(0, 29).Select(_ => Gaussian(rng))];
        double[] yy = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
        double[] beta = ConformalBayesTests.ExactDraws(x, yy, 200, rng);
        (PointwiseLogDensity training, CandidateLogDensity candidates) = Conjugate(x, yy, beta, [0.3, -1.2], grid);
        ConformalPrediction full = ConformalBayes.FullStep(training, candidates);
        ConformalPrediction weighed = ConformalBayes.FullStep(training, candidates, new ConformalOptions
        {
            Weights = ShiftWeights.ByCase([.. Enumerable.Repeat(1.5, 29)], [1.5, 1.5], candidates.Grids),
        });
        for (int t = 0; t < 2; t++)
        {
            Assert.Null(weighed.Sets[t].WeightedAlpha);
            Assert.Equal(full.Sets[t].PValues, weighed.Sets[t].PValues);
            Assert.Equal(full.Sets[t].SplitHalfMass, weighed.Sets[t].SplitHalfMass);
        }
    }

    // The conjugate model's terms, y = βx + N(0, 1): the training units' and the cases' on their grids.
    private static (PointwiseLogDensity Training, CandidateLogDensity Candidates) Conjugate(double[] x, double[] y, double[] beta, double[] cases, double[] grid)
    {
        int s = beta.Length, n = x.Length;
        var tv = new double[s, n];
        var cv = new double[s, cases.Length * grid.Length];
        for (int d = 0; d < s; d++)
        {
            for (int i = 0; i < n; i++) tv[d, i] = ConformalBayesTests.LogNormal(y[i], beta[d] * x[i]);
            for (int t = 0; t < cases.Length; t++)
                for (int g = 0; g < grid.Length; g++) cv[d, t * grid.Length + g] = ConformalBayesTests.LogNormal(grid[g], beta[d] * cases[t]);
        }
        return (new PointwiseLogDensity(tv, UnitPartition.EachObservation(n), []),
                new CandidateLogDensity(cv, [.. Enumerable.Repeat(grid, cases.Length)], []));
    }

    // ---- Coverage under shift ----

    // Training inputs uniform on (0, 1); the test's from the density 2x, so w(x) = 2x. The noise grows with x, and the fit, the right curve
    // with one scale, is too narrow where the test lands.
    private static double Curve(double x) => Math.Sin(2.0 * Math.PI * x);
    private static double Noise(double x) => 0.1 + 0.6 * x;

    [Fact]
    public void CovariateShift_TheTrueRatioCovers_AnEstimatedOneLosesAtMostItsTotalVariation()
    {
        const int n = 199, replications = 4000;
        const double alpha = 0.1;
        var rng = new Xoshiro256PlusPlus(2026_1113);
        int[] plain = new int[3], trueRatio = new int[3], estimated = new int[3];
        double tv = 0.0;
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            double[] y = [.. x.Select(v => Curve(v) + Noise(v) * Gaussian(rng))];
            double xt = Math.Sqrt(rng.NextDouble()), yt = Curve(xt) + Noise(xt) * Gaussian(rng);
            double[] w = [.. x.Select(v => 2.0 * v)], we = [.. x.Select(v => 0.5 + v)];
            double wt = 2.0 * xt, wet = 0.5 + xt;
            // The data set's normalized true and estimated weights, test included: their total variation.
            double st = w.Sum() + wt, se = we.Sum() + wet, gap = Math.Max(0.0, wt / st - wet / se);
            for (int i = 0; i < n; i++) gap += Math.Max(0.0, w[i] / st - we[i] / se);
            tv += gap;

            // The density score, one draw: N(y; sin 2πx, 0.35²).
            var cal = new double[1, n];
            for (int i = 0; i < n; i++) cal[0, i] = Normal.MeanSd(Curve(x[i]), 0.35).LogPdf(y[i]);
            var calibration = new CandidateLogDensity(cal, [.. y.Select(v => new[] { v })], []);
            var test = new CandidateLogDensity(new double[,] { { Normal.MeanSd(Curve(xt), 0.35).LogPdf(yt) } }, [[yt]], []);
            // CQR's quantiles and the residual's centre, the same curve.
            var q = new QuantileReadings([.. x.Select(v => Curve(v) - 0.45)], [.. x.Select(v => Curve(v) + 0.45)]);
            var qt = new QuantileReadings([Curve(xt) - 0.45], [Curve(xt) + 0.45]);
            var c = new ResidualReadings([.. x.Select(Curve)]);
            var ct = new ResidualReadings([Curve(xt)]);
            bool[] Read(ShiftWeights? weights) =>
            [
                SplitConformal.Density(calibration, test, alpha, weights).Sets[0].Accepts(0),
                SplitConformal.Quantiles(q, y, qt, alpha, weights).Intervals[0].Contains(yt),
                SplitConformal.Residuals(c, y, ct, alpha, weights).Intervals[0].Contains(yt),
            ];
            bool[] a = Read(null), b = Read(new ShiftWeights(w, [wt])), e = Read(new ShiftWeights(we, [wet]));
            for (int f = 0; f < 3; f++)
            {
                if (a[f]) plain[f]++;
                if (b[f]) trueRatio[f]++;
                if (e[f]) estimated[f]++;
            }
        }
        double se1 = Math.Sqrt(alpha * (1 - alpha) / replications), meanTv = tv / replications;
        string[] forms = ["density", "CQR", "residual"];
        for (int f = 0; f < 3; f++)
        {
            double p = plain[f] / (double)replications, t = trueRatio[f] / (double)replications, e = estimated[f] / (double)replications;
            output.WriteLine($"{forms[f]}: unweighted {p:F4}, the true ratio {t:F4}, an estimated one {e:F4} (mean TV {meanTv:F4}); se {se1:F4}");
            Assert.True(t >= 1 - alpha - 3.5 * se1, $"{forms[f]}: the true ratio covers {t:F4}");
            Assert.True(e >= 1 - alpha - meanTv - 3.5 * se1, $"{forms[f]}: the estimated ratio covers {e:F4}");
            Assert.True(p < 1 - alpha - 3.0 * se1, $"{forms[f]}: the unweighted sets should fall short under the shift, {p:F4}");
        }
    }

    [Fact]
    public void ResponseShift_TheTrueRatioCovers()
    {
        // Training y ~ N(0, 1) and test y ~ N(1, 1), the inputs' law given y unchanged, x | y ~ N(0.7y, 0.5²): w(y) = e^(y − 1/2). The fit is
        // the training law's own predictive, y | x ~ N(0.946x, 0.338). Each candidate's weight moves with its value.
        const int n = 199, replications = 4000;
        const double alpha = 0.1, slope = 0.7 / 0.74, spread = 0.338;
        var rng = new Xoshiro256PlusPlus(2026_1114);
        int plain = 0, weighted = 0;
        for (int r = 0; r < replications; r++)
        {
            double[] y = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))];
            double[] x = [.. y.Select(v => 0.7 * v + 0.5 * Gaussian(rng))];
            double yt = 1.0 + Gaussian(rng), xt = 0.7 * yt + 0.5 * Gaussian(rng);
            var cal = new double[1, n];
            for (int i = 0; i < n; i++) cal[0, i] = Normal.MeanSd(slope * x[i], Math.Sqrt(spread)).LogPdf(y[i]);
            var calibration = new CandidateLogDensity(cal, [.. y.Select(v => new[] { v })], []);
            var test = new CandidateLogDensity(new double[,] { { Normal.MeanSd(slope * xt, Math.Sqrt(spread)).LogPdf(yt) } }, [[yt]], []);
            var weights = new ShiftWeights([.. y.Select(v => Math.Exp(v - 0.5))], [Math.Exp(yt - 0.5)]);
            if (SplitConformal.Density(calibration, test, alpha).Sets[0].Accepts(0)) plain++;
            if (SplitConformal.Density(calibration, test, alpha, weights).Sets[0].Accepts(0)) weighted++;
        }
        double p = plain / (double)replications, q = weighted / (double)replications, se = Math.Sqrt(alpha * (1 - alpha) / replications);
        output.WriteLine($"response shift: unweighted {p:F4}, the true ratio {q:F4} (se {se:F4})");
        Assert.True(q >= 1 - alpha - 3.5 * se);
    }

    [Fact]
    public void WhatTheWeightsDoNotServe_IsRefused()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new ShiftWeights([1.0, -0.1], [1.0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ShiftWeights([1.0, double.NaN], [1.0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ShiftWeights([1.0], [double.PositiveInfinity]));
        var rng = new Xoshiro256PlusPlus(8);
        double[] grid = [-1.0, 0.0, 1.0];
        (CandidateLogDensity cal, CandidateLogDensity cand, _, _) = Seams(19, 2, grid, rng);
        Assert.Throws<ArgumentException>(() => SplitConformal.Density(cal, cand, weights: new ShiftWeights([.. Enumerable.Repeat(1.0, 19)], [1.0, 2.0])));
        Assert.Throws<ArgumentException>(() => SplitConformal.Density(cal, cand, weights: new ShiftWeights([1.0, 2.0], new double[6])));
        double[] x = [.. Enumerable.Range(0, 19).Select(_ => Gaussian(rng))];
        double[] y = [.. x.Select(v => v + Gaussian(rng))];
        (PointwiseLogDensity training, CandidateLogDensity candidates) = Conjugate(x, y, ConformalBayesTests.ExactDraws(x, y, 50, rng), [0.5], grid);
        var weights = ShiftWeights.ByCase([.. x.Select(v => 1.0 + v * v)], [2.0], candidates.Grids);
        foreach (ConformalOptions bad in new[]
        {
            new ConformalOptions { Weights = weights, MarginBounds = MarginBounds.Kept },
            new ConformalOptions { Weights = weights, FallbackShare = 0.1 },
            new ConformalOptions { Weights = weights, Score = FullStepScore.LeaveOneOut },
        })
            Assert.Contains("weighted rule", Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(training, candidates, bad)).Message);
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(training, candidates, new ConformalOptions { Weights = new ShiftWeights([1.0], [1.0, 1.0, 1.0]) }));
        // Unequal weights: the full step's weighted rule, its p-values the weight at or below over the total.
        ConformalSet set = ConformalBayes.FullStep(training, candidates, new ConformalOptions { Weights = weights }).Sets[0];
        Assert.Equal(0.1, set.WeightedAlpha);
        Assert.Equal(0.9, set.Coverage, 12);
        // Its margin record holds the weighted inflated rule's readings (slice 14's follow-up 3); the band was not asked for, so the coverage
        // under error is not read.
        Assert.NotNull(set.Margin!.WeightedInflated);
        Assert.Null(set.Margin.WeightedFragile);
        Assert.True(double.IsNaN(set.CoverageUnderError));
    }

    // ---- Behind the qualification switch ----

    [Fact]
    public void Measurement_TheFullStepUnderCovariateShift_OnExactDraws()
    {
        // n 49, exact draws of the conjugate model y = βx + N(0, 1) under the conformal tests' heteroscedastic truth, training x ~ N(0, 1) and
        // test x ~ N(1, 1), so w(x) = e^(x − 1/2); 2,000 replications at the truth.
        if (!BassQualification.Enabled) return;
        const int n = 49, replications = 2000;
        var rng = new Xoshiro256PlusPlus(2026_1115);
        int plain = 0, weighted = 0;
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))];
            double[] y = [.. x.Select(v => ConformalBayesTests.Truth(v, 1.0, rng))];
            double xt = 1.0 + Gaussian(rng), yt = ConformalBayesTests.Truth(xt, 1.0, rng);
            (PointwiseLogDensity training, CandidateLogDensity candidates) = Conjugate(x, y, ConformalBayesTests.ExactDraws(x, y, 400, rng), [xt], [yt]);
            if (ConformalBayes.FullStep(training, candidates).Sets[0].Accepts(0)) plain++;
            var weights = new ShiftWeights([.. x.Select(v => Math.Exp(v - 0.5))], [Math.Exp(xt - 0.5)]);
            if (ConformalBayes.FullStep(training, candidates, new ConformalOptions { Weights = weights }).Sets[0].Accepts(0)) weighted++;
        }
        double p = plain / (double)replications, q = weighted / (double)replications, se = Math.Sqrt(0.09 / replications);
        output.WriteLine($"the full step under covariate shift: unweighted {p:F4}, the true ratio {q:F4} (se {se:F4}; nominal 0.9000)");
        Assert.True(q >= 0.9 - 3.5 * se);
    }

    [Fact]
    public void Measurement_ADriftingStream_UnderWeightsByPosition()
    {
        // 199 calibration units in time order, then the test unit; the response's level drifts upward by 0.006 a step. A fixed fit, the
        // curve without the drift, scored by the residual. Weights ρ^(n − i) by position (Barber et al.'s), the test's 1, against equal
        // weights: a reading, the bound's TV term unknown.
        if (!BassQualification.Enabled) return;
        const int n = 199, replications = 4000;
        var rng = new Xoshiro256PlusPlus(2026_1116);
        var counts = new Dictionary<double, int> { [1.0] = 0, [0.99] = 0, [0.97] = 0 };
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, n + 1).Select(_ => rng.NextDouble())];
            double[] y = [.. x.Select((v, i) => Curve(v) + 0.006 * i + 0.3 * Gaussian(rng))];
            var c = new ResidualReadings([.. x[..n].Select(Curve)]);
            var ct = new ResidualReadings([Curve(x[n])]);
            foreach (double rho in counts.Keys.ToArray())
            {
                var weights = new ShiftWeights([.. Enumerable.Range(0, n).Select(i => Math.Pow(rho, n - i))], [1.0]);
                if (SplitConformal.Residuals(c, y[..n], ct, 0.1, weights).Intervals[0].Contains(y[n])) counts[rho]++;
            }
        }
        output.WriteLine("a drifting stream: " + string.Join(", ", counts.Select(kv => $"ρ {kv.Key}: coverage {kv.Value / (double)replications:F4}")));
    }
}
