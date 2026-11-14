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
/// Posterior conformal prediction given the memberships (E16, the conformal layer's slice 15a; Zhang and Candès, arXiv 2409.19712): the
/// randomized weights' rule is exact given the draw (<c>coverage_pcp</c>), checked by enumerating every draw and every position of the test;
/// per-case weights give each case what it gets called alone, in every form; on a two-cluster mixture with known memberships the posterior
/// weights cover per cluster where the unweighted split does not, the KL weights within their gap (<c>pcp_cluster_gap</c>).
/// </summary>
public sealed class PosteriorConformalTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    // Every composition of m into k parts.
    private static IEnumerable<int[]> Compositions(int m, int k)
    {
        if (k == 1) { yield return [m]; yield break; }
        for (int first = 0; first <= m; first++)
            foreach (int[] rest in Compositions(m - first, k - 1)) yield return [first, .. rest];
    }

    [Fact]
    public void GivenTheDraw_TheRuleIsExact_ByEnumeration()
    {
        // 400 random data sets of 6 points in 3 clusters, m = 3 (every draw L), ties among the scores: for each L, the test's position j
        // drawn with chance ∝ P(L | π(xⱼ)), as it is given L for exchangeable data, and the refused chance is at most α.
        var rng = new Xoshiro256PlusPlus(15);
        const int points = 6, k = 3, m = 3;
        double worst = double.NegativeInfinity;
        int checkedDraws = 0;
        for (int r = 0; r < 400; r++)
        {
            double alpha = 0.1 + 0.4 * rng.NextDouble();
            var pi = new double[points, k];
            for (int i = 0; i < points; i++)
                for (int c = 0; c < k; c++) pi[i, c] = 0.05 + rng.NextDouble();
            double[] y = [.. Enumerable.Range(0, points).Select(_ => Math.Round(2.0 * rng.NextDouble(), 1))];
            foreach (int[] draw in Compositions(m, k))
            {
                double chance = 0.0, refused = 0.0;
                for (int j = 0; j < points; j++)
                {
                    int[] rest = [.. Enumerable.Range(0, points).Where(i => i != j)];
                    var units = new double[rest.Length, k];
                    var test = new double[1, k];
                    double sum = 0.0;
                    for (int c = 0; c < k; c++) sum += pi[j, c];
                    for (int c = 0; c < k; c++) test[0, c] = pi[j, c];
                    for (int u = 0; u < rest.Length; u++)
                        for (int c = 0; c < k; c++) units[u, c] = pi[rest[u], c];
                    double p = Math.Exp(PosteriorConformal.LogKernel([.. Enumerable.Range(0, k).Select(c => pi[j, c] / sum)], draw));
                    PosteriorWeights weights = PosteriorConformal.GivenDraws(units, test, [draw]);
                    bool accepted = SplitConformal.Residuals(new ResidualReadings([.. rest.Select(_ => 0.0)]), [.. rest.Select(i => y[i])],
                                                             new ResidualReadings([0.0]), alpha, weights.Weights).Intervals[0].Contains(y[j]);
                    chance += p;
                    if (!accepted) refused += p;
                }
                worst = Math.Max(worst, refused / chance - alpha);
                Assert.True(refused / chance <= alpha + 1e-12, $"given L = ({string.Join(", ", draw)}), the refused chance {refused / chance} at α {alpha}");
                checkedDraws++;
            }
        }
        output.WriteLine($"{checkedDraws} draws: the refused chance given L, its largest excess over α: {worst:F4}");
    }

    [Fact]
    public void TheDraws_AreSeededPerCase_AndFollowTheMemberships()
    {
        double[,] units = { { 0.9, 0.1 }, { 0.2, 0.8 }, { 0.5, 0.5 } };
        double[,] cases = { { 0.3, 0.7 }, { 1.0, 0.0 }, { 0.6, 0.4 } };
        PosteriorWeights a = PosteriorConformal.Weights(units, cases, 12, seed: 4), b = PosteriorConformal.Weights(units, cases, 12, seed: 4);
        PosteriorWeights first = PosteriorConformal.Weights(units, new double[,] { { 0.3, 0.7 } }, 12, seed: 4);
        Assert.Equal(a.Draws!, b.Draws!);
        Assert.Equal(a.Draws![0], first.Draws![0]);
        Assert.All(a.Draws!, d => Assert.Equal(12, d.Sum()));
        Assert.Equal(new[] { 12, 0 }, a.Draws![1]);
        Assert.Null(a.ClusterGaps);
        for (int t = 0; t < 3; t++)
        {
            IReadOnlyList<double> w = a.Weights.UnitsFor(t);
            double own = a.Weights.Candidates[t], total = w.Sum() + own;
            Assert.Equal(own / total, a.SelfWeights[t], 12);
            Assert.Equal(total * total / (w.Sum(v => v * v) + own * own), a.EffectiveSampleSizes[t], 9);
        }
        // Case 1's memberships put nothing on cluster 2, so its draw is all cluster 1: the units weigh 0.9¹², 0.2¹² and 0.5¹² to its 1.
        Assert.Equal(Math.Pow(0.2 / 0.9, 12), a.Weights.UnitsFor(1)[1] / a.Weights.UnitsFor(1)[0], 12);

        // The draws' mean is mπ.
        const int many = 4000;
        var same = new double[many, 2];
        for (int t = 0; t < many; t++) { same[t, 0] = 0.3; same[t, 1] = 0.7; }
        PosteriorWeights spread = PosteriorConformal.Weights(units, same, 12, seed: 9);
        double mean = spread.Draws!.Average(d => d[0]), se = Math.Sqrt(12 * 0.3 * 0.7 / many);
        output.WriteLine($"the draws' mean count of cluster 1: {mean:F4} (3.6 expected, se {se:F4})");
        Assert.True(Math.Abs(mean - 3.6) < 4.0 * se);

        // The KL weights: no draws, the test's own the largest, a gap per case.
        PosteriorWeights kl = PosteriorConformal.Weights(units, cases, 12, weighting: PosteriorWeighting.KullbackLeibler);
        Assert.Null(kl.Draws);
        Assert.NotNull(kl.ClusterGaps);
        for (int t = 0; t < 3; t++)
        {
            Assert.Equal(1.0, kl.Weights.Candidates[t]);
            Assert.All(kl.Weights.UnitsFor(t), w => Assert.True(w <= 1.0));
        }
        // Case 1, (1, 0): unit 0 at KL ln(1/0.9) weighs 0.9¹², unit 1 0.2¹²; its gap is Σ qᵢ·2|π₁(i) − 1|.
        IReadOnlyList<double> w1 = kl.Weights.UnitsFor(1);
        Assert.Equal(Math.Pow(0.9, 12), w1[0], 12);
        double totalKl = w1.Sum() + 1.0;
        Assert.Equal((w1[0] * 0.2 + w1[1] * 1.6 + w1[2] * 1.0) / totalKl, kl.ClusterGaps![1], 12);
    }

    // ---- Per-case weights: each case gets what it gets alone ----

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

    // Case t's columns of a candidate seam whose cases share one grid.
    private static CandidateLogDensity CaseOf(CandidateLogDensity seam, double[] grid, int t)
    {
        int draws = seam.Values.GetLength(0);
        var values = new double[draws, grid.Length];
        for (int d = 0; d < draws; d++)
            for (int g = 0; g < grid.Length; g++) values[d, g] = seam.Values[d, t * grid.Length + g];
        return new CandidateLogDensity(values, [grid], []);
    }

    [Fact]
    public void PerCaseWeights_GiveEachCaseItsOwnCall_InEveryForm()
    {
        var rng = new Xoshiro256PlusPlus(151);
        double[] grid = [.. Enumerable.Range(0, 41).Select(g => -3.0 + 0.15 * g)];
        const int n = 39, cases = 4;
        (CandidateLogDensity cal, CandidateLogDensity cand, double[] y, double[] mu) = Seams(n, cases, grid, rng);
        double[][] byCase = [.. Enumerable.Range(0, cases).Select(_ => Enumerable.Range(0, n).Select(_ => 0.2 + rng.NextDouble()).ToArray())];
        double[] own = [.. Enumerable.Range(0, cases).Select(_ => 0.2 + rng.NextDouble())];

        ConformalPrediction all = SplitConformal.Density(cal, cand, weights: ShiftWeights.ByCase(byCase, own, cand.Grids));
        for (int t = 0; t < cases; t++)
        {
            ConformalSet alone = SplitConformal.Density(cal, CaseOf(cand, grid, t), weights: ShiftWeights.ByCase(byCase[t], [own[t]], [grid])).Sets[0];
            Assert.Equal(0.1, all.Sets[t].WeightedAlpha);
            Assert.Equal(alone.PValues, all.Sets[t].PValues);
        }

        var quantiles = new QuantileReadings([.. mu.Select(v => v - 0.6)], [.. mu.Select(v => v + 0.6)]);
        double[] lo = [-0.6, 0.0, 0.4, -1.1], hi = [0.6, 1.2, 1.4, 0.2], centre = [0.0, 0.3, -1.0, 0.8];
        var residuals = new ResidualReadings(mu, [.. mu.Select(_ => 0.7)]);
        var perCase = ShiftWeights.PerCase(byCase, own);
        ConformalIntervals q = SplitConformal.Quantiles(quantiles, y, new QuantileReadings(lo, hi), weights: perCase);
        AsymmetricConformalIntervals aq = SplitConformal.AsymmetricQuantiles(quantiles, y, new QuantileReadings(lo, hi), weights: perCase);
        ConformalIntervals r = SplitConformal.Residuals(residuals, y, new ResidualReadings(centre, [.. centre.Select(_ => 0.7)]), weights: perCase);
        for (int t = 0; t < cases; t++)
        {
            var shared = new ShiftWeights(byCase[t], [own[t]]);
            Assert.Equal(SplitConformal.Quantiles(quantiles, y, new QuantileReadings([lo[t]], [hi[t]]), weights: shared).Intervals[0] with { Case = t }, q.Intervals[t]);
            Assert.Equal(SplitConformal.AsymmetricQuantiles(quantiles, y, new QuantileReadings([lo[t]], [hi[t]]), weights: shared).Intervals[0] with { Case = t }, aq.Intervals[t]);
            Assert.Equal(SplitConformal.Residuals(residuals, y, new ResidualReadings([centre[t]], [0.7]), weights: shared).Intervals[0] with { Case = t }, r.Intervals[t]);
            Assert.Equal(0.1, r.Intervals[t].WeightedAlpha);
        }

        // Each case's weights all one value, a different one per case: CB-1's sets, bit for bit.
        double[] level = [0.4, 1.0, 2.5, 7.0];
        var flat = ShiftWeights.ByCase([.. level.Select(v => Enumerable.Repeat(v, n).ToArray())], level, cand.Grids);
        ConformalPrediction plain = SplitConformal.Density(cal, cand), equal = SplitConformal.Density(cal, cand, weights: flat);
        for (int t = 0; t < cases; t++)
        {
            Assert.Null(equal.Sets[t].WeightedAlpha);
            Assert.Equal(plain.Sets[t].PValues, equal.Sets[t].PValues);
        }

        // The full step on exact draws.
        double[] x = [.. Enumerable.Range(0, 29).Select(_ => Gaussian(rng))];
        double[] yy = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
        double[] beta = ConformalBayesTests.ExactDraws(x, yy, 200, rng);
        double[] at = [0.3, -1.2, 0.8];
        (PointwiseLogDensity training, CandidateLogDensity candidates) = Conjugate(x, yy, beta, at, grid);
        double[][] fullByCase = [.. at.Select(_ => x.Select(_ => 0.2 + rng.NextDouble()).ToArray())];
        double[] fullOwn = [.. at.Select(_ => 0.2 + rng.NextDouble())];
        ConformalPrediction full = ConformalBayes.FullStep(training, candidates, new ConformalOptions { Weights = ShiftWeights.ByCase(fullByCase, fullOwn, candidates.Grids) });
        for (int t = 0; t < at.Length; t++)
        {
            (PointwiseLogDensity tr, CandidateLogDensity one) = Conjugate(x, yy, beta, [at[t]], grid);
            ConformalSet alone = ConformalBayes.FullStep(tr, one, new ConformalOptions { Weights = ShiftWeights.ByCase(fullByCase[t], [fullOwn[t]], [grid]) }).Sets[0];
            Assert.Equal(0.1, full.Sets[t].WeightedAlpha);
            Assert.Equal(alone.PValues, full.Sets[t].PValues);
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

    // ---- A two-cluster mixture with known memberships ----

    // x uniform on (0, 1); cluster 1's membership 1/(1 + e^((x − ½)/0.06)), so the clusters meet around ½; the residual's spread 0.2 in
    // cluster 1 and 1.0 in cluster 2. The fit is the true curve, scored by the absolute residual: a mixture model of the scores with the
    // true memberships.
    private static double Curve(double x) => Math.Sin(2.0 * Math.PI * x);
    private static double Membership(double x) => 1.0 / (1.0 + Math.Exp((x - 0.5) / 0.06));

    [Fact]
    public void TwoClusters_ThePosteriorWeightsCoverPerCluster_WhereTheSplitDoesNot()
    {
        const int n = 199, replications = 3000, m = 20;
        const double alpha = 0.1;
        var rng = new Xoshiro256PlusPlus(2026_1117);
        string[] rules = ["unweighted", "Mondrian", "randomized", "KL"];
        var covered = new int[4, 2];
        var widths = new double[4, 2];
        var finite = new int[4, 2];
        var tests = new int[2];
        var gaps = new double[2];
        var ess = new double[2];
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            double[] y = [.. x.Select(v => Curve(v) + (rng.NextDouble() < Membership(v) ? 0.2 : 1.0) * Gaussian(rng))];
            double xt = rng.NextDouble(), yt = Curve(xt) + (rng.NextDouble() < Membership(xt) ? 0.2 : 1.0) * Gaussian(rng);
            int region = Membership(xt) >= 0.5 ? 0 : 1;
            tests[region]++;
            var units = new double[n, 2];
            for (int i = 0; i < n; i++) { units[i, 0] = Membership(x[i]); units[i, 1] = 1.0 - units[i, 0]; }
            var test = new double[,] { { Membership(xt), 1.0 - Membership(xt) } };
            PosteriorWeights randomized = PosteriorConformal.Weights(units, test, m, seed: r);
            PosteriorWeights kl = PosteriorConformal.Weights(units, test, m, weighting: PosteriorWeighting.KullbackLeibler);
            gaps[region] += kl.ClusterGaps![0];
            ess[region] += randomized.EffectiveSampleSizes[0];
            int[] strata = [.. x.Select(v => Membership(v) >= 0.5 ? 0 : 1)];
            ConformalInterval[] read =
            [
                SplitConformal.Residuals(new ResidualReadings([.. x.Select(Curve)]), y, new ResidualReadings([Curve(xt)]), alpha).Intervals[0],
                SplitConformal.Residuals(new ResidualReadings([.. x.Select(Curve)], null, strata), y, new ResidualReadings([Curve(xt)], null, [region]), alpha).Intervals[0],
                SplitConformal.Residuals(new ResidualReadings([.. x.Select(Curve)]), y, new ResidualReadings([Curve(xt)]), alpha, randomized.Weights).Intervals[0],
                SplitConformal.Residuals(new ResidualReadings([.. x.Select(Curve)]), y, new ResidualReadings([Curve(xt)]), alpha, kl.Weights).Intervals[0],
            ];
            for (int f = 0; f < 4; f++)
            {
                if (read[f].Contains(yt)) covered[f, region]++;
                if (!read[f].WholeSpace) { widths[f, region] += read[f].Width; finite[f, region]++; }
            }
        }
        double nominal = 1 - alpha, seAll = Math.Sqrt(alpha * nominal / replications);
        double Rate(int f, int h) => covered[f, h] / (double)tests[h];
        double Marginal(int f) => (covered[f, 0] + covered[f, 1]) / (double)replications;
        for (int f = 0; f < 4; f++)
            output.WriteLine($"{rules[f]}: marginal {Marginal(f):F4}; cluster 1 {Rate(f, 0):F4}, cluster 2 {Rate(f, 1):F4}; mean width {widths[f, 0] / finite[f, 0]:F3} and "
                             + $"{widths[f, 1] / finite[f, 1]:F3}, whole space {tests[0] - finite[f, 0]} and {tests[1] - finite[f, 1]}");
        output.WriteLine($"KL gap by cluster {gaps[0] / tests[0]:F4}, {gaps[1] / tests[1]:F4}; randomized mean ESS {ess[0] / tests[0]:F1}, {ess[1] / tests[1]:F1}; "
                         + $"{tests[0]} and {tests[1]} tests; se {seAll:F4}");

        // The randomized weights are exact marginally; the KL weights within their mean gap, marginally and per cluster (the cluster is the
        // test's dominant membership, a function of its input, under which the bound given the inputs holds).
        Assert.True(Marginal(2) >= nominal - 3.5 * seAll, $"randomized marginal {Marginal(2):F4}");
        Assert.True(Marginal(3) >= nominal - (gaps[0] + gaps[1]) / replications - 3.5 * seAll, $"KL marginal {Marginal(3):F4}");
        for (int h = 0; h < 2; h++)
        {
            double se = Math.Sqrt(alpha * nominal / tests[h]);
            Assert.True(Rate(3, h) >= nominal - gaps[h] / tests[h] - 3.5 * se, $"KL cluster {h + 1}: {Rate(3, h):F4}");
        }
        // The unweighted split falls short in the wide cluster, and the posterior weights narrow the narrow one.
        Assert.True(Rate(0, 1) < nominal - 3.0 * Math.Sqrt(alpha * nominal / tests[1]), $"unweighted cluster 2 {Rate(0, 1):F4}");
        Assert.True(widths[2, 0] / finite[2, 0] < 0.7 * widths[0, 0] / finite[0, 0]);
    }

    [Fact]
    public void WhatTheMembershipsDoNotServe_IsRefused()
    {
        double[,] units = { { 0.5, 0.5 }, { 0.2, 0.8 } };
        double[,] one = { { 0.4, 0.6 } };
        Assert.Throws<ArgumentOutOfRangeException>(() => PosteriorConformal.Weights(new double[,] { { 0.5, -0.1 } }, one, 5));
        Assert.Throws<ArgumentOutOfRangeException>(() => PosteriorConformal.Weights(units, new double[,] { { double.NaN, 1.0 } }, 5));
        Assert.Throws<ArgumentException>(() => PosteriorConformal.Weights(units, new double[,] { { 0.0, 0.0 } }, 5));
        Assert.Throws<ArgumentException>(() => PosteriorConformal.Weights(units, new double[,] { { 0.2, 0.3, 0.5 } }, 5));
        Assert.Throws<ArgumentException>(() => PosteriorConformal.Weights(new double[2, 0], new double[1, 0], 5));
        Assert.Throws<ArgumentOutOfRangeException>(() => PosteriorConformal.Weights(units, one, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => PosteriorConformal.Weights(units, one, 5, weighting: (PosteriorWeighting)7));
        Assert.Throws<ArgumentException>(() => PosteriorConformal.Weights(units, one, 5, grids: [[0.0], [1.0]]));
        Assert.Throws<ArgumentException>(() => PosteriorConformal.GivenDraws(units, new double[,] { { 1.0, 0.0 } }, [[2, 1]]));
        Assert.Throws<ArgumentException>(() => PosteriorConformal.GivenDraws(units, one, [[2, 1], [1, 2]]));
        Assert.Throws<ArgumentException>(() => ShiftWeights.PerCase([[1.0, 2.0], [1.0]], [1.0, 1.0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => ShiftWeights.PerCase([[1.0, -2.0]], [1.0]));
        // Per-case weights for another number of cases, or of units.
        var residuals = new ResidualReadings([0.0, 0.0, 0.0]);
        Assert.Throws<ArgumentException>(() => SplitConformal.Residuals(residuals, [0.1, 0.5, 0.9], new ResidualReadings([0.0, 0.0]),
                                                                        weights: ShiftWeights.PerCase([[1.0, 1.0, 2.0]], [1.0, 1.0])));
        Assert.Throws<ArgumentException>(() => SplitConformal.Residuals(residuals, [0.1, 0.5, 0.9], new ResidualReadings([0.0]),
                                                                        weights: ShiftWeights.PerCase([[1.0, 2.0]], [1.0])));
    }
}
