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
/// The weighted predictive distribution (E16, slice 14's follow-up 1): with the tied mass randomized its PIT is exactly uniform when the
/// test unit lands on each unit with chance proportional to its weight (<c>sum_mul_wPitCdf</c>), checked by enumeration with ties and zero
/// weights; equal weights give the counts' distribution bit for bit; under covariate and response shift with the true ratio the split
/// form's PIT is uniform where the unweighted one is not (<c>pit_weighted</c>), and an estimated ratio moves it by at most its total
/// variation (<c>wPit_estimated</c>). With the qualification switch (<c>BASS_QUALIFICATION</c> = 1): the full step under covariate shift
/// on exact draws.
/// </summary>
public sealed class WeightedDistributionTests(ITestOutputHelper output)
{
    private static readonly LatentReport[] Theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];

    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    // A one-draw seam whose split score at each column is the given CDF value; each case one column.
    private static CandidateLogDensity Scores(double[] cdf)
    {
        var values = new double[1, cdf.Length];
        var logCdf = new double[1, cdf.Length];
        for (int c = 0; c < cdf.Length; c++) logCdf[0, c] = Math.Log(cdf[c]);
        return new CandidateLogDensity(values, [.. cdf.Select(_ => new[] { 0.0 })], Theta) { LogCdf = logCdf };
    }

    [Fact]
    public void TheWeightedPit_IsExactlyUniform_ByEnumeration()
    {
        // 2,000 data sets of 2 to 9 units, scores on a coarse lattice so many tie, weights with zeros. The test unit is unit u with chance
        // w_u/Σw, and its PIT T = Q⁻ + V·(Q⁺ − Q⁻) with V uniform: P(T ≤ x) = Σ_u q_u·clamp((x − Q⁻_u)/(Q⁺_u − Q⁻_u)) is x at every x, and
        // the ties counted for the candidate (V = 1) super-uniform.
        var rng = new Xoshiro256PlusPlus(1401);
        double worst = 0.0, worstSuper = double.NegativeInfinity;
        int used = 0;
        for (int r = 0; r < 2000; r++)
        {
            int n = 2 + rng.NextInt(8);
            double[] score = [.. Enumerable.Range(0, n).Select(_ => 0.1 * (1 + rng.NextInt(5)))];
            double[] w = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble() < 0.2 ? 0.0 : rng.NextDouble())];
            double total = w.Sum();
            if (!(total > 0.0)) continue;
            used++;
            var lower = new double[n];
            var upper = new double[n];
            for (int u = 0; u < n; u++)
            {
                int[] rest = [.. Enumerable.Range(0, n).Where(i => i != u)];
                ConformalDistribution d = SplitConformal.Distribution(Scores([.. rest.Select(i => score[i])]), Scores([score[u]]),
                                                                      new ShiftWeights([.. rest.Select(i => w[i])], [w[u]])).Distributions[0];
                lower[u] = d.Lower(0);
                upper[u] = d.Upper(0);
            }
            for (int k = 0; k <= 20; k++)
            {
                double x = k / 20.0, p = 0.0, plain = 0.0;
                for (int u = 0; u < n; u++)
                {
                    double q = w[u] / total, spread = upper[u] - lower[u];
                    p += q * (spread > 0.0 ? Math.Clamp((x - lower[u]) / spread, 0.0, 1.0) : (x >= lower[u] ? 1.0 : 0.0));
                    if (upper[u] <= x) plain += q;
                }
                worst = Math.Max(worst, Math.Abs(p - x));
                worstSuper = Math.Max(worstSuper, plain - x);
                Assert.True(Math.Abs(p - x) < 1e-12, $"data set {r}: P(T ≤ {x}) = {p}");
                Assert.True(plain <= x + 1e-12, $"data set {r}: P(Q⁺ ≤ {x}) = {plain}");
            }
        }
        output.WriteLine($"{used} data sets: |P(T ≤ x) − x| at most {worst:E2}; P(Q⁺ ≤ x) − x at most {worstSuper:E2}");
    }

    [Fact]
    public void EqualWeights_GiveTheCountsDistribution_BitForBit()
    {
        var rng = new Xoshiro256PlusPlus(1402);
        double[] beta = [.. Enumerable.Range(0, 40).Select(_ => 0.5 + 0.1 * Gaussian(rng))];
        double[] x = [.. Enumerable.Range(0, 30).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
        double[] grid = [.. Enumerable.Range(0, 31).Select(g => -3.0 + 0.2 * g)];
        CandidateLogDensity calibration = ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, x, [.. y.Select(v => new[] { v })]);
        CandidateLogDensity candidates = ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, [0.3, -1.0], [grid, grid]);
        ConformalDistributions plain = SplitConformal.Distribution(calibration, candidates);
        var flat = ShiftWeights.ByCase([.. x.Select(_ => 0.4)], [0.4, 0.4], candidates.Grids);
        var perCase = ShiftWeights.ByCase([Enumerable.Repeat(1.7, 30).ToArray(), Enumerable.Repeat(0.2, 30).ToArray()], [1.7, 0.2], candidates.Grids);
        foreach (ShiftWeights weights in new[] { flat, perCase })
        {
            ConformalDistributions weighed = SplitConformal.Distribution(calibration, candidates, weights);
            for (int t = 0; t < 2; t++)
            {
                Assert.Null(weighed.Distributions[t].Weighted);
                Assert.Equal(plain.Distributions[t].Below, weighed.Distributions[t].Below);
                Assert.Equal(plain.Distributions[t].Tied, weighed.Distributions[t].Tied);
            }
        }

        PointwiseLogDensity training = ConformalPredictiveDistributionTests.Training(beta, 1.0, null, x, y);
        ConformalDistributions full = ConformalBayes.Distribution(training, candidates, new ConformalOptions());
        ConformalDistributions fullWeighed = ConformalBayes.Distribution(training, candidates, new ConformalOptions { Weights = flat });
        for (int t = 0; t < 2; t++)
        {
            Assert.Null(fullWeighed.Distributions[t].Weighted);
            Assert.Equal(full.Distributions[t].Below, fullWeighed.Distributions[t].Below);
            Assert.Equal(full.Distributions[t].Tied, fullWeighed.Distributions[t].Tied);
        }
    }

    [Fact]
    public void PerCaseWeights_GiveEachCaseItsOwnCall_InBothForms()
    {
        var rng = new Xoshiro256PlusPlus(1403);
        double[] beta = [.. Enumerable.Range(0, 40).Select(_ => 0.5 + 0.1 * Gaussian(rng))];
        double[] x = [.. Enumerable.Range(0, 25).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
        double[] grid = [.. Enumerable.Range(0, 21).Select(g => -2.0 + 0.2 * g)], at = [0.3, -1.0, 1.2];
        CandidateLogDensity calibration = ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, x, [.. y.Select(v => new[] { v })]);
        CandidateLogDensity candidates = ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, at, [grid, grid, grid]);
        PointwiseLogDensity training = ConformalPredictiveDistributionTests.Training(beta, 1.0, null, x, y);
        double[][] byCase = [.. at.Select(_ => x.Select(_ => 0.1 + rng.NextDouble()).ToArray())];
        double[] own = [.. at.Select(_ => 0.1 + rng.NextDouble())];
        var weights = ShiftWeights.ByCase(byCase, own, candidates.Grids);
        ConformalDistributions split = SplitConformal.Distribution(calibration, candidates, weights);
        ConformalDistributions full = ConformalBayes.Distribution(training, candidates, new ConformalOptions { Weights = weights });
        for (int t = 0; t < at.Length; t++)
        {
            CandidateLogDensity one = ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, [at[t]], [grid]);
            var alone = ShiftWeights.ByCase(byCase[t], [own[t]], [grid]);
            ConformalDistribution s = SplitConformal.Distribution(calibration, one, alone).Distributions[0];
            ConformalDistribution f = ConformalBayes.Distribution(training, one, new ConformalOptions { Weights = alone }).Distributions[0];
            Assert.NotNull(split.Distributions[t].Weighted);
            Assert.Equal(s.Weighted!.Lower, split.Distributions[t].Weighted!.Lower);
            Assert.Equal(s.Weighted.Upper, split.Distributions[t].Weighted!.Upper);
            Assert.Equal(f.Weighted!.Lower, full.Distributions[t].Weighted!.Lower);
            Assert.Equal(f.Weighted.Upper, full.Distributions[t].Weighted!.Upper);
            Assert.True(split.Distributions[t].Monotone, "the split form under weights fixed along a case's grid rises");
            Assert.Equal(0.9, split.Distributions[t].Coverage(0.1), 12);
        }
    }

    // ---- The split form under shift ----

    private static double Curve(double x) => Math.Sin(2.0 * Math.PI * x);
    private static double Noise(double x) => 0.1 + 0.6 * x;

    // One draw's seam: ln N(y; mean, sd²) and ln Φ((y − mean)/sd) at each case's grid.
    private static CandidateLogDensity Normal1(double[] mean, double sd, double[][] grids)
    {
        int columns = grids.Sum(g => g.Length);
        var values = new double[1, columns];
        var logCdf = new double[1, columns];
        for (int t = 0, c = 0; t < mean.Length; t++)
            foreach (double v in grids[t])
            {
                values[0, c] = Normal.MeanSd(mean[t], sd).LogPdf(v);
                logCdf[0, c++] = Normal.StandardLogCdf((v - mean[t]) / sd);
            }
        return new CandidateLogDensity(values, grids, Theta) { LogCdf = logCdf };
    }

    // The PITs' largest deviation from uniform at the deciles, in standard errors.
    private static double Deviation(List<double> pits)
    {
        double worst = 0.0;
        for (int k = 1; k <= 9; k++)
        {
            double x = k / 10.0, p = pits.Count(v => v <= x) / (double)pits.Count, se = Math.Sqrt(x * (1 - x) / pits.Count);
            worst = Math.Max(worst, Math.Abs(p - x) / se);
        }
        return worst;
    }

    [Fact]
    public void CovariateShift_TheWeightedPitIsUniform_AnEstimatedRatioWithinItsTotalVariation()
    {
        // Slice 14's covariate shift: training x uniform, test x of density 2x, so w(x) = 2x; the fit N(sin 2πx, 0.35²) is too narrow where
        // the test lands. 4,000 replications of 199 calibration units and a test case, the PIT read at the true y with τ uniform.
        const int n = 199, replications = 4000;
        var rng = new Xoshiro256PlusPlus(2026_1121);
        List<double> plain = [], weighted = [], estimated = [];
        int covered = 0;
        double tv = 0.0;
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            double[] y = [.. x.Select(v => Curve(v) + Noise(v) * Gaussian(rng))];
            double xt = Math.Sqrt(rng.NextDouble()), yt = Curve(xt) + Noise(xt) * Gaussian(rng);
            CandidateLogDensity calibration = Normal1([.. x.Select(Curve)], 0.35, [.. y.Select(v => new[] { v })]);
            CandidateLogDensity test = Normal1([Curve(xt)], 0.35, [[yt]]);
            double[] w = [.. x.Select(v => 2.0 * v)], we = [.. x.Select(v => 0.5 + v)];
            double st = w.Sum() + 2.0 * xt, se = we.Sum() + 0.5 + xt, gap = Math.Max(0.0, 2.0 * xt / st - (0.5 + xt) / se);
            for (int i = 0; i < n; i++) gap += Math.Max(0.0, w[i] / st - we[i] / se);
            tv += gap;
            double tau = rng.NextDouble();
            plain.Add(SplitConformal.Distribution(calibration, test).Distributions[0].At(0, tau));
            ConformalDistribution d = SplitConformal.Distribution(calibration, test, new ShiftWeights(w, [2.0 * xt])).Distributions[0];
            weighted.Add(d.At(0, tau));
            if (d.Accepts(0, 0.1)) covered++;
            estimated.Add(SplitConformal.Distribution(calibration, test, new ShiftWeights(we, [0.5 + xt])).Distributions[0].At(0, tau));
        }
        double meanTv = tv / replications, coverage = covered / (double)replications, seCover = Math.Sqrt(0.09 / replications);
        double shift = Enumerable.Range(1, 9).Max(k => Math.Abs(estimated.Count(v => v <= k / 10.0) / (double)replications - k / 10.0));
        output.WriteLine($"PIT's largest deviation at the deciles, in se: unweighted {Deviation(plain):F2}, the true ratio {Deviation(weighted):F2}, "
                         + $"an estimated one {Deviation(estimated):F2} (largest |P − x| {shift:F4}, mean TV {meanTv:F4}); the weighted set covers {coverage:F4}");
        Assert.True(Deviation(weighted) < 4.0);
        Assert.True(Deviation(plain) > 5.0);
        Assert.True(shift <= meanTv + 4.0 * Math.Sqrt(0.25 / replications));
        Assert.True(coverage >= 0.9 - 3.5 * seCover);
    }

    [Fact]
    public void ResponseShift_TheWeightedPitIsUniform()
    {
        // Slice 14's response shift: training y ~ N(0, 1), test y ~ N(1, 1), x | y ~ N(0.7y, 0.5²) unchanged, so w(y) = e^(y − 1/2), each
        // candidate's weight at its own value; the fit is the training law's predictive y | x ~ N(0.946x, 0.338).
        const int n = 199, replications = 4000;
        const double slope = 0.7 / 0.74, sd = 0.5814;
        var rng = new Xoshiro256PlusPlus(2026_1122);
        List<double> plain = [], weighted = [];
        for (int r = 0; r < replications; r++)
        {
            double[] y = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))];
            double[] x = [.. y.Select(v => 0.7 * v + 0.5 * Gaussian(rng))];
            double yt = 1.0 + Gaussian(rng), xt = 0.7 * yt + 0.5 * Gaussian(rng);
            CandidateLogDensity calibration = Normal1([.. x.Select(v => slope * v)], sd, [.. y.Select(v => new[] { v })]);
            CandidateLogDensity test = Normal1([slope * xt], sd, [[yt]]);
            double tau = rng.NextDouble();
            plain.Add(SplitConformal.Distribution(calibration, test).Distributions[0].At(0, tau));
            weighted.Add(SplitConformal.Distribution(calibration, test, new ShiftWeights([.. y.Select(v => Math.Exp(v - 0.5))], [Math.Exp(yt - 0.5)]))
                                       .Distributions[0].At(0, tau));
        }
        output.WriteLine($"response shift, the PIT's largest deviation at the deciles in se: unweighted {Deviation(plain):F2}, the true ratio {Deviation(weighted):F2}");
        Assert.True(Deviation(weighted) < 4.0);
        Assert.True(Deviation(plain) > 5.0);
    }

    [Fact]
    public void WhatTheWeightsDoNotServe_IsRefused()
    {
        var rng = new Xoshiro256PlusPlus(1404);
        double[] beta = [.. Enumerable.Range(0, 20).Select(_ => 0.5 + 0.1 * Gaussian(rng))];
        double[] x = [.. Enumerable.Range(0, 12).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
        double[] grid = [-1.0, 0.0, 1.0];
        CandidateLogDensity calibration = ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, x, [.. y.Select(v => new[] { v })]);
        CandidateLogDensity candidates = ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, [0.3, 0.6], [grid, grid]);
        PointwiseLogDensity training = ConformalPredictiveDistributionTests.Training(beta, 1.0, null, x, y);
        double[] units = [.. x.Select(v => 1.0 + v * v)];
        var good = ShiftWeights.ByCase(units, [2.0, 0.5], candidates.Grids);
        ShiftWeights[] bad =
        [
            ShiftWeights.ByCase(units[..11], [2.0, 0.5], candidates.Grids),
            new ShiftWeights(units, [2.0, 0.5]),
            ShiftWeights.PerCase([units], new double[6]),
        ];
        foreach (ShiftWeights weights in bad)
        {
            Assert.Throws<ArgumentException>(() => SplitConformal.Distribution(calibration, candidates, weights));
            Assert.Throws<ArgumentException>(() => ConformalBayes.Distribution(training, candidates, new ConformalOptions { Weights = weights }));
        }
        foreach (ConformalOptions options in new[]
        {
            new ConformalOptions { Weights = good, MarginBounds = MarginBounds.Carried },
            new ConformalOptions { Weights = good, Inflate = true },
        })
            Assert.Contains("weighted ranks", Assert.Throws<ArgumentException>(() => ConformalBayes.Distribution(training, candidates, options)).Message);
        // Unequal weights: each case reads its weighted ranks, the error readouts unread.
        ConformalDistribution d = ConformalBayes.Distribution(training, candidates, new ConformalOptions { Weights = good }).Distributions[1];
        Assert.NotNull(d.Weighted);
        Assert.Null(d.Margin);
        Assert.Equal(0.9, d.Coverage(0.1), 12);
        Assert.False(d.WholeSpace(0.1));
    }

    // ---- Behind the qualification switch ----

    [Fact]
    public void Measurement_TheFullStepsPitUnderCovariateShift_OnExactDraws()
    {
        // n 49, exact draws of the conjugate model y = βx + N(0, 1) under the conformal tests' heteroscedastic truth, training x ~ N(0, 1)
        // and test x ~ N(1, 1), so w(x) = e^(x − 1/2); 2,000 replications at the truth, the PIT read at the true y with τ uniform.
        if (!BassQualification.Enabled) return;
        const int n = 49, replications = 2000;
        var rng = new Xoshiro256PlusPlus(2026_1123);
        List<double> plain = [], weighted = [];
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))];
            double[] y = [.. x.Select(v => ConformalBayesTests.Truth(v, 1.0, rng))];
            double xt = 1.0 + Gaussian(rng), yt = ConformalBayesTests.Truth(xt, 1.0, rng);
            double[] beta = ConformalBayesTests.ExactDraws(x, y, 400, rng);
            PointwiseLogDensity training = ConformalPredictiveDistributionTests.Training(beta, 1.0, null, x, y);
            CandidateLogDensity test = ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, [xt], [[yt]]);
            double tau = rng.NextDouble();
            plain.Add(ConformalBayes.Distribution(training, test).Distributions[0].At(0, tau));
            var weights = new ShiftWeights([.. x.Select(v => Math.Exp(v - 0.5))], [Math.Exp(xt - 0.5)]);
            weighted.Add(ConformalBayes.Distribution(training, test, new ConformalOptions { Weights = weights }).Distributions[0].At(0, tau));
        }
        output.WriteLine($"the full step's PIT under covariate shift, largest deviation at the deciles in se: unweighted {Deviation(plain):F2}, "
                         + $"the true ratio {Deviation(weighted):F2}");
        Assert.True(Deviation(weighted) < 4.0);
    }
}
