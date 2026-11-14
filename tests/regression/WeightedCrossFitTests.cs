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
/// The cross-fitted forms under shift weights (E16, slice 14's follow-up 2): with likelihood ratios travelling with the units, the weighted
/// jackknife+ and the n-fold cross-conformal set miss at most 2α of the weight over every placement of the test unit
/// (<c>jackknifePlus_likelihoodRatio</c>, <c>crossFit_likelihoodRatio</c>), checked by enumeration with a symmetric fit; equal weights give
/// the counts' forms bit for bit; under covariate shift the weighted jackknife+ and CV+ cover where the unweighted fall short; and CV+'s
/// weighted bounds read as the Lean states them (<c>cvPlus_oneFold_likelihoodRatio</c>, <c>cvPlus_merged_weighted</c>,
/// <c>coverage_across_likelihoodRatio</c>).
/// </summary>
public sealed class WeightedCrossFitTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    // The leave-one-out readings of the mean of the training responses, a symmetric fit: row i the mean without unit i, at the unit and at
    // the test case.
    private static LeaveOneOutReadings MeanReadings(double[] y)
    {
        int n = y.Length;
        double sum = y.Sum();
        double[] heldOut = [.. y.Select(v => (sum - v) / (n - 1))];
        var predictions = new double[n, 1];
        for (int i = 0; i < n; i++) predictions[i, 0] = heldOut[i];
        return new LeaveOneOutReadings(y, heldOut, predictions);
    }

    [Fact]
    public void LikelihoodRatioWeights_MissAtMostTwoAlpha_ByEnumeration()
    {
        // 1,500 bags of 3 to 8 units, responses on a coarse lattice so residuals tie, weights with zeros. The test unit is unit u with chance
        // w_u/Σw, the fit the mean of the rest's responses, its leave-one-out fits the means without each training unit; over the placements
        // the weighted jackknife+ and the weighted pairwise set (scores −|residual|) each miss at most 2α of the weight.
        var rng = new Xoshiro256PlusPlus(1501);
        double worstJack = 0.0, worstPair = 0.0;
        int used = 0;
        for (int r = 0; r < 1500; r++)
        {
            int size = 3 + rng.NextInt(6);
            double alpha = 0.05 + 0.4 * rng.NextDouble();
            double[] y = [.. Enumerable.Range(0, size).Select(_ => 0.5 * rng.NextInt(7))];
            double[] w = [.. Enumerable.Range(0, size).Select(_ => rng.NextDouble() < 0.15 ? 0.0 : rng.NextDouble())];
            double total = w.Sum();
            if (!(total > 0.0)) continue;
            used++;
            double jackMiss = 0.0, pairMiss = 0.0;
            for (int u = 0; u < size; u++)
            {
                int[] rest = [.. Enumerable.Range(0, size).Where(i => i != u)];
                double[] ty = [.. rest.Select(i => y[i])];
                LeaveOneOutReadings readings = MeanReadings(ty);
                var weights = new ShiftWeights([.. rest.Select(i => w[i])], [w[u]]);
                if (!CrossConformal.JackknifePlus(readings, alpha, weights).Intervals[0].Contains(y[u])) jackMiss += w[u] / total;
                double[] heldOut = [.. ty.Select((v, i) => -Math.Abs(v - readings.HeldOut[i]))];
                var candidate = new double[ty.Length, 1];
                for (int i = 0; i < ty.Length; i++) candidate[i, 0] = -Math.Abs(y[u] - readings.Predictions[i, 0]);
                if (!CrossConformal.Pairwise(heldOut, candidate, [[y[u]]], alpha, weights).Sets[0].Accepts(0)) pairMiss += w[u] / total;
            }
            worstJack = Math.Max(worstJack, jackMiss / (2.0 * alpha));
            worstPair = Math.Max(worstPair, pairMiss / (2.0 * alpha));
            Assert.True(jackMiss <= 2.0 * alpha + 1e-12, $"bag {r}: the jackknife+ misses {jackMiss} at α {alpha}");
            Assert.True(pairMiss <= 2.0 * alpha + 1e-12, $"bag {r}: the pairwise set misses {pairMiss} at α {alpha}");
        }
        output.WriteLine($"{used} bags: the largest share of 2α missed, jackknife+ {worstJack:F3}, pairwise {worstPair:F3}");
    }

    [Fact]
    public void EqualWeights_GiveTheCountsForms_BitForBit()
    {
        var rng = new Xoshiro256PlusPlus(1502);
        double[] y = [.. Enumerable.Range(0, 30).Select(_ => Gaussian(rng))];
        var predictions = new double[30, 3];
        for (int i = 0; i < 30; i++)
            for (int t = 0; t < 3; t++) predictions[i, t] = 0.1 * t + 0.05 * Gaussian(rng);
        var readings = new LeaveOneOutReadings(y, [.. y.Select(_ => 0.05 * Gaussian(rng))], predictions);
        var flat = new ShiftWeights([.. y.Select(_ => 0.7)], [0.7, 0.7, 0.7]);
        Assert.Equal(CrossConformal.JackknifePlus(readings).Intervals, CrossConformal.JackknifePlus(readings, weights: flat).Intervals);

        double[] heldOut = [.. y.Select(_ => Gaussian(rng))];
        double[] grid = [-1.0, 0.0, 1.0];
        var scores = new double[30, 6];
        for (int i = 0; i < 30; i++)
            for (int c = 0; c < 6; c++) scores[i, c] = Gaussian(rng);
        var columns = new ShiftWeights([.. y.Select(_ => 2.5)], [.. Enumerable.Repeat(2.5, 6)]);
        ConformalPrediction a = CrossConformal.Pairwise(heldOut, scores, [grid, grid]), b = CrossConformal.Pairwise(heldOut, scores, [grid, grid], weights: columns);
        ConformalPrediction c2 = CrossConformal.Tournament(scores, scores, [grid, grid]), d = CrossConformal.Tournament(scores, scores, [grid, grid], weights: columns);
        for (int t = 0; t < 2; t++)
        {
            Assert.Null(b.Sets[t].WeightedAlpha);
            Assert.Equal(a.Sets[t].PValues, b.Sets[t].PValues);
            Assert.Null(d.Sets[t].WeightedAlpha);
            Assert.Equal(c2.Sets[t].PValues, d.Sets[t].PValues);
        }

        int[] foldOf = [.. Enumerable.Range(0, 30).Select(i => i % 5)];
        var foldPredictions = new double[5, 3];
        for (int k = 0; k < 5; k++)
            for (int t = 0; t < 3; t++) foldPredictions[k, t] = 0.1 * t + 0.05 * Gaussian(rng);
        var folds = new FoldReadings(y, readings.HeldOut, foldOf, foldPredictions);
        CvPlusIntervals plain = CrossConformal.CvPlus(folds, CvPlusAssumptions.Both);
        CvPlusIntervals equal = CrossConformal.CvPlus(folds, CvPlusAssumptions.Both, weights: new CvPlusWeights(flat, CvPlusWeighting.LikelihoodRatio));
        for (int t = 0; t < 3; t++)
        {
            Assert.Null(equal.Intervals[t].WeightedAlpha);
            Assert.Equal(plain.Intervals[t].Lower, equal.Intervals[t].Lower);
            Assert.Equal(plain.Intervals[t].Upper, equal.Intervals[t].Upper);
            Assert.Same(equal.Bounds, equal.Intervals[t].Bounds);
        }

        // Per-case weights: each case what it gets alone.
        double[][] byCase = [.. Enumerable.Range(0, 3).Select(_ => y.Select(_ => 0.1 + rng.NextDouble()).ToArray())];
        double[] own = [0.4, 1.3, 0.8];
        JackknifePlusIntervals all = CrossConformal.JackknifePlus(readings, weights: ShiftWeights.PerCase(byCase, own));
        for (int t = 0; t < 3; t++)
        {
            var column = new double[30, 1];
            for (int i = 0; i < 30; i++) column[i, 0] = predictions[i, t];
            JackknifePlusInterval alone = CrossConformal.JackknifePlus(readings with { Predictions = column }, weights: new ShiftWeights(byCase[t], [own[t]])).Intervals[0];
            Assert.Equal(alone with { Case = t }, all.Intervals[t]);
            Assert.Equal(0.8, all.Intervals[t].Coverage, 12);
        }
    }

    // ---- Covariate shift ----

    private static double Curve(double x) => Math.Sin(2.0 * Math.PI * x);
    private static double Noise(double x) => 0.1 + 0.6 * x;

    // Least squares y ≈ a + bx on the units in `keep`, read at `at`.
    private static double Line(double[] x, double[] y, IEnumerable<int> keep, double at)
    {
        double n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
        foreach (int i in keep)
        {
            n++; sx += x[i]; sy += y[i]; sxx += x[i] * x[i]; sxy += x[i] * y[i];
        }
        double b = (n * sxy - sx * sy) / (n * sxx - sx * sx), a = (sy - b * sx) / n;
        return a + b * at;
    }

    [Fact]
    public void CovariateShift_TheWeightedJackknifePlusAndCvPlusCover()
    {
        // Training x uniform, test x of density 2x, so w(x) = 2x (JAWS's setting); y = sin 2πx + (0.1 + 0.6x)ε, fitted by a line, which is
        // wrong, by leave-one-out and by 10 folds. 2,000 replications of 99 units and a test case at α 0.1.
        const int n = 100, replications = 2000, folds = 10;
        const double alpha = 0.1;
        var rng = new Xoshiro256PlusPlus(2026_1124);
        int[] covered = new int[4];
        double oneFold = 0.0, merged = 0.0;
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            double[] y = [.. x.Select(v => Curve(v) + Noise(v) * Gaussian(rng))];
            double xt = Math.Sqrt(rng.NextDouble()), yt = Curve(xt) + Noise(xt) * Gaussian(rng);
            var predictions = new double[n, 1];
            double[] heldOut = new double[n];
            for (int i = 0; i < n; i++)
            {
                IEnumerable<int> rest = Enumerable.Range(0, n).Where(j => j != i);
                heldOut[i] = Line(x, y, rest, x[i]);
                predictions[i, 0] = Line(x, y, rest, xt);
            }
            var readings = new LeaveOneOutReadings(y, heldOut, predictions);
            var weights = new ShiftWeights([.. x.Select(v => 2.0 * v)], [2.0 * xt]);
            if (CrossConformal.JackknifePlus(readings, alpha).Intervals[0].Contains(yt)) covered[0]++;
            if (CrossConformal.JackknifePlus(readings, alpha, weights).Intervals[0].Contains(yt)) covered[1]++;

            int[] foldOf = [.. Enumerable.Range(0, n).Select(i => i % folds)];
            var foldPredictions = new double[folds, 1];
            double[] foldHeldOut = new double[n];
            for (int k = 0; k < folds; k++)
            {
                int[] rest = [.. Enumerable.Range(0, n).Where(j => foldOf[j] != k)];
                foldPredictions[k, 0] = Line(x, y, rest, xt);
                for (int i = 0; i < n; i++) if (foldOf[i] == k) foldHeldOut[i] = Line(x, y, rest, x[i]);
            }
            var foldReadings = new FoldReadings(y, foldHeldOut, foldOf, foldPredictions);
            if (CrossConformal.CvPlus(foldReadings, CvPlusAssumptions.Both, alpha).Intervals[0].Contains(yt)) covered[2]++;
            CvPlusInterval weighed = CrossConformal.CvPlus(foldReadings, CvPlusAssumptions.Both, alpha,
                                                           new CvPlusWeights(weights, CvPlusWeighting.LikelihoodRatio, RatioBound: 2.0)).Intervals[0];
            if (weighed.Contains(yt)) covered[3]++;
            oneFold += weighed.Bounds.OneFold.Value;
            merged += weighed.Bounds.Merged.Value;
        }
        double se = Math.Sqrt(alpha * (1 - alpha) / replications);
        double[] rate = [.. covered.Select(c => c / (double)replications)];
        output.WriteLine($"covariate shift: jackknife+ unweighted {rate[0]:F4}, weighted {rate[1]:F4}; CV+ unweighted {rate[2]:F4}, weighted {rate[3]:F4} "
                         + $"(its bounds on average: one fold {oneFold / replications:F4}, merged {merged / replications:F4}); se {se:F4}");
        Assert.True(rate[1] >= 1.0 - 2.0 * alpha - 3.5 * se);
        Assert.True(rate[1] > rate[0] + 3.0 * se, "the weights should help the jackknife+ under the shift");
        Assert.True(rate[3] > rate[2] + 3.0 * se, "the weights should help CV+ under the shift");
    }

    [Fact]
    public void CvPlusWeightedBounds_ReadAsTheLeanStatesThem()
    {
        // Three folds of two units; unit weights 1, 2 | 1, 1 | 3, 1, so W_k = 3, 2, 4 and W = 9.
        double[] y = [0.1, -0.2, 0.4, 0.0, -0.3, 0.2];
        int[] foldOf = [0, 0, 1, 1, 2, 2];
        var folds = new FoldReadings(y, [0.0, 0.1, -0.1, 0.05, 0.0, 0.1], foldOf, new double[,] { { 0.0 }, { 0.1 }, { -0.05 } });
        double[] units = [1.0, 2.0, 1.0, 1.0, 3.0, 1.0];
        const double alpha = 0.2;

        // Position weights, the test's the largest (3): c = (α(W + w_t) + (K − 1)w_t)/(W + K·w_t) = (0.2·12 + 6)/18.
        CvPlusInterval position = CrossConformal.CvPlus(folds, CvPlusAssumptions.Both, alpha,
                                                        new CvPlusWeights(new ShiftWeights(units, [3.0]), CvPlusWeighting.Position)).Intervals[0];
        double c = (0.2 * 12.0 + 6.0) / 18.0;
        Assert.Equal(Math.Max(0.0, 1.0 - 2.0 * c), position.Bounds.Merged.Value, 12);
        Assert.True(position.Bounds.Merged.Applies);
        Assert.False(position.Bounds.OneFold.Applies);
        Assert.Equal(0.2, position.WeightedAlpha);
        CvPlusInterval lighter = CrossConformal.CvPlus(folds, CvPlusAssumptions.Both, alpha,
                                                       new CvPlusWeights(new ShiftWeights(units, [2.0]), CvPlusWeighting.Position)).Intervals[0];
        Assert.False(lighter.Bounds.Merged.Applies, "the test position's weight must be the largest");
        Assert.Null(lighter.Coverage);

        // Likelihood ratios, the test's 2: one fold's share α·W_tot/(W_max + w_t) = 0.2·11/6; across folds 2α + max(4, 2·w̄)/W_tot.
        CvPlusInterval ratio = CrossConformal.CvPlus(folds, CvPlusAssumptions.Both, alpha,
                                                     new CvPlusWeights(new ShiftWeights(units, [2.0]), CvPlusWeighting.LikelihoodRatio, RatioBound: 3.0)).Intervals[0];
        Assert.Equal(1.0 - 0.2 * 11.0 / 6.0, ratio.Bounds.OneFold.Value, 12);
        Assert.Equal(Math.Max(0.0, 1.0 - 0.4 - 6.0 / 11.0), ratio.Bounds.Tournament.Value, 12);
        // Merged (cvPlus_merged_likelihoodRatio): W = 11 with the test's 2, w_max = 3, so 2c(3) = 2·8.2/17; refined at γ = 2, the largest
        // weight at most αW = 2.2: L = 6.2, N = 15, e = Σw²/W = 21/11, 2L(W + 2e)/N².
        double simple = Math.Min(1.0, 2.0 * 8.2 / 17.0), refined = 2.0 * 6.2 * (11.0 + 2.0 * 21.0 / 11.0) / 225.0;
        Assert.Equal(1.0 - Math.Min(simple, refined), ratio.Bounds.Merged.Value, 12);
        Assert.True(ratio.Bounds.Merged.Applies);
        Assert.Equal(new[] { 1.0 - 0.2 * 11.0 / 6.0, Math.Max(0.0, 1.0 - 0.4 - 6.0 / 11.0), 1.0 - Math.Min(simple, refined) }.Max(), ratio.Coverage!.Value, 12);
        // Every weight above αW: no unit can be refused, the merged bound 1.
        CvPlusInterval heavy = CrossConformal.CvPlus(folds, CvPlusAssumptions.MergedRanks, 0.1,
                                                     new CvPlusWeights(new ShiftWeights([1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [1.5]), CvPlusWeighting.LikelihoodRatio)).Intervals[0];
        Assert.Equal(1.0, heavy.Bounds.Merged.Value);
        CvPlusInterval unbounded = CrossConformal.CvPlus(folds, CvPlusAssumptions.MergedRanks, alpha,
                                                         new CvPlusWeights(new ShiftWeights(units, [2.0]), CvPlusWeighting.LikelihoodRatio)).Intervals[0];
        Assert.False(unbounded.Bounds.Tournament.Applies, "the across-folds bound reads a bound on the ratio");
        Assert.True(unbounded.Bounds.OneFold.Applies);
    }

    [Fact]
    public void WhatTheWeightsDoNotServe_IsRefused()
    {
        double[] y = [0.1, -0.2, 0.4, 0.0];
        var readings = new LeaveOneOutReadings(y, [0.0, 0.1, -0.1, 0.05], new double[,] { { 0.0 }, { 0.1 }, { -0.05 }, { 0.0 } });
        Assert.Throws<ArgumentException>(() => CrossConformal.JackknifePlus(readings, 0.2, new ShiftWeights([1.0, 1.0, 2.0], [1.0])));
        Assert.Throws<ArgumentException>(() => CrossConformal.JackknifePlus(readings, 0.2, new ShiftWeights([1.0, 1.0, 2.0, 1.0], [1.0, 2.0])));
        Assert.Throws<ArgumentException>(() => CrossConformal.Pairwise(y, new double[4, 2], [[0.0, 1.0]], 0.2, new ShiftWeights([1.0, 1.0, 2.0, 1.0], [1.0])));
        Assert.Throws<ArgumentException>(() => CrossConformal.Tournament(new double[4, 2], new double[4, 2], [[0.0, 1.0]], 0.2,
                                                                         weights: ShiftWeights.PerCase([[1.0, 1.0, 2.0, 1.0], [1.0, 1.0, 1.0, 1.0]], [1.0, 1.0])));
        var folds = new FoldReadings(y, [0.0, 0.1, -0.1, 0.05], [0, 0, 1, 1], new double[,] { { 0.0 }, { 0.1 } });
        var units = new ShiftWeights([1.0, 2.0, 1.0, 1.0], [2.0]);
        Assert.Throws<ArgumentException>(() => CrossConformal.CvPlus(folds, CvPlusAssumptions.Both, 0.2, new CvPlusWeights(units, CvPlusWeighting.Position, RatioBound: 2.0)));
        Assert.Throws<ArgumentOutOfRangeException>(() => CrossConformal.CvPlus(folds, CvPlusAssumptions.Both, 0.2, new CvPlusWeights(units, CvPlusWeighting.LikelihoodRatio, RatioBound: 0.0)));
        Assert.Throws<ArgumentOutOfRangeException>(() => CrossConformal.CvPlus(folds, CvPlusAssumptions.Both, 0.2, new CvPlusWeights(units, (CvPlusWeighting)9)));
        Assert.Throws<ArgumentException>(() => CrossConformal.CvPlus(folds, CvPlusAssumptions.Both, 0.2, new CvPlusWeights(new ShiftWeights([1.0, 2.0], [2.0]), CvPlusWeighting.LikelihoodRatio)));
        // A weighted cross-fitted set reports 1 − 2α.
        ConformalSet weighedSet = CrossConformal.Pairwise(y, new double[4, 2], [[0.0, 1.0]], 0.2, new ShiftWeights([1.0, 1.0, 2.0, 1.0], [1.0, 3.0])).Sets[0];
        Assert.Equal(0.2, weighedSet.WeightedAlpha);
        Assert.Equal(0.6, weighedSet.Coverage, 12);
        // Weights of zero total refuse nothing: the whole line.
        JackknifePlusInterval zero = CrossConformal.JackknifePlus(readings, 0.2, new ShiftWeights([0.0, 0.0, 0.0, 0.0], [0.0])).Intervals[0];
        Assert.True(zero.WholeSpace || zero.WeightedAlpha is null);
    }
}
