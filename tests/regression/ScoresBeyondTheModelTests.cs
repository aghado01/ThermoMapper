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
/// Scores beyond the model's own (E16, the conformal layer's slice 17): conformal Bayesian model averaging across noise laws (CBMA), its
/// weights against each law's exact evidence and its sets against each law's own; the residual scores, distance to origin, distance to
/// average and RoBAS-EB's shrinkage, under a right and a shifted predictor; and Wu et al.'s highest-density sets with their threshold by the
/// high-probability rule, covering with confidence β over calibration sets.
/// </summary>
public sealed class ScoresBeyondTheModelTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    // ---- CBMA across noise laws ----

    // y = βx + N(0, σ²), β ~ N(0, 1): exact posterior draws and the exact log evidence.
    private static double[] Draws(double[] x, double[] y, double sigma, int s, Xoshiro256PlusPlus rng)
    {
        double sxx = x.Sum(v => v * v), sxy = x.Zip(y, (a, b) => a * b).Sum(), v = 1.0 / (1.0 + sxx / (sigma * sigma)), m = v * sxy / (sigma * sigma);
        return [.. Enumerable.Range(0, s).Select(_ => m + Math.Sqrt(v) * Gaussian(rng))];
    }

    private static double LogEvidence(double[] x, double[] y, double sigma)
    {
        int n = x.Length;
        double s2 = sigma * sigma, sxx = x.Sum(v => v * v), sxy = x.Zip(y, (a, b) => a * b).Sum(), syy = y.Sum(v => v * v);
        double logDet = n * Math.Log(s2) + Math.Log(1.0 + sxx / s2), quad = (syy - sxy * sxy / (s2 + sxx)) / s2;
        return -0.5 * n * Math.Log(2 * Math.PI) - 0.5 * logDet - 0.5 * quad;
    }

    private static double LogNormal(double y, double mean, double sigma)
        => -0.5 * Math.Log(2 * Math.PI) - Math.Log(sigma) - 0.5 * (y - mean) * (y - mean) / (sigma * sigma);

    // Each unit's log density at its own response under each draw, a split seam's values.
    private static double[,] AtResponses(double[] x, double[] y, double[] beta, double sigma)
    {
        var values = new double[beta.Length, x.Length];
        for (int d = 0; d < beta.Length; d++)
            for (int i = 0; i < x.Length; i++) values[d, i] = LogNormal(y[i], beta[d] * x[i], sigma);
        return values;
    }

    private static (PointwiseLogDensity Training, CandidateLogDensity Candidates) Seams(double[] x, double[] y, double[] beta, double sigma,
                                                                                         double[] cases, double[] grid)
    {
        int s = beta.Length, n = x.Length;
        var tv = new double[s, n];
        var cv = new double[s, cases.Length * grid.Length];
        for (int d = 0; d < s; d++)
        {
            for (int i = 0; i < n; i++) tv[d, i] = LogNormal(y[i], beta[d] * x[i], sigma);
            for (int t = 0; t < cases.Length; t++)
                for (int g = 0; g < grid.Length; g++) cv[d, t * grid.Length + g] = LogNormal(grid[g], beta[d] * cases[t], sigma);
        }
        return (new PointwiseLogDensity(tv, UnitPartition.EachObservation(n), []),
                new CandidateLogDensity(cv, [.. Enumerable.Repeat(grid, cases.Length)], []));
    }

    [Fact]
    public void Cbma_WithOneModel_IsTheFullStep_AndItsWeightsAreTheEvidences()
    {
        var rng = new Xoshiro256PlusPlus(2001);
        double[] x = [.. Enumerable.Range(0, 25).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.7 * v + Gaussian(rng))];
        double[] grid = [.. Enumerable.Range(0, 41).Select(g => -4.0 + 0.2 * g)];
        var one = Seams(x, y, Draws(x, y, 1.0, 4000, rng), 1.0, [0.5, -1.0], grid);
        ConformalPrediction full = ConformalBayes.FullStep(one.Training, one.Candidates);
        AveragedPrediction averaged = ModelAveraging.FullStep([one], [0.0]);
        for (int t = 0; t < 2; t++)
        {
            Assert.Equal(full.Sets[t].PValues, averaged.Prediction.Sets[t].PValues);
            Assert.All(Enumerable.Range(0, grid.Length), g => Assert.Equal(1.0, averaged.ModelWeights[t][g, 0]));
        }

        // Two noise laws, σ 1 and 1.5, equal priors, on 8 units so the evidence leaves both plausible: q_k(y) against the exact
        // p(M_k | Z₁:ₙ, (x, y)) from each law's evidence with the candidate.
        double[] x8 = x[..8], y8 = y[..8];
        var narrow = Seams(x8, y8, Draws(x8, y8, 1.0, 4000, rng), 1.0, [0.5, -1.0], grid);
        var wide = Seams(x8, y8, Draws(x8, y8, 1.5, 4000, rng), 1.5, [0.5, -1.0], grid);
        double[] logWeights = [LogEvidence(x8, y8, 1.0), LogEvidence(x8, y8, 1.5)];
        AveragedPrediction two = ModelAveraging.FullStep([narrow, wide], logWeights);
        double worst = 0.0, lowest = 1.0, highest = 0.0;
        double[] cases = [0.5, -1.0];
        for (int t = 0; t < 2; t++)
            for (int g = 0; g < grid.Length; g += 2)
            {
                double[] xa = [.. x8, cases[t]], ya = [.. y8, grid[g]];
                double a = LogEvidence(xa, ya, 1.0), b = LogEvidence(xa, ya, 1.5), exact = 1.0 / (1.0 + Math.Exp(b - a));
                worst = Math.Max(worst, Math.Abs(two.ModelWeights[t][g, 0] - exact));
                lowest = Math.Min(lowest, exact);
                highest = Math.Max(highest, exact);
            }
        output.WriteLine($"CBMA's weight on σ 1 against the exact posterior model probability (from {lowest:F3} to {highest:F3} over the grid): the largest difference {worst:F4} (4,000 draws)");
        Assert.True(highest - lowest > 0.2, "the check needs weights that move with the candidate");
        Assert.True(worst < 0.02);
    }

    [Fact]
    public void Cbma_Covers_AndTracksTheRightModel()
    {
        // Data from y = 0.7x + N(0, 1); models y = βx and y = βx², each with σ 1, their exact evidence and 400 exact draws; n 30, 500
        // replications, α 0.1. A model wrong only in its noise scale ranks the units as the right one does, its density falling with
        // |residual| alike, so its conformal set barely moves; one wrong in its mean does not.
        const int n = 30, replications = 500;
        var rng = new Xoshiro256PlusPlus(2026_1132);
        double[] grid = [.. Enumerable.Range(0, 81).Select(g => -6.0 + 0.15 * g)];
        int[] covered = new int[3];
        double[] width = new double[3];
        double weight = 0.0;
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.7 * v + Gaussian(rng))];
            double xt = Gaussian(rng), yt = 0.7 * xt + Gaussian(rng);
            double[] x2 = [.. x.Select(v => v * v)];
            var right = Seams(x, y, Draws(x, y, 1.0, 400, rng), 1.0, [xt], grid);
            var wrong = Seams(x2, y, Draws(x2, y, 1.0, 400, rng), 1.0, [xt * xt], grid);
            AveragedPrediction averaged = ModelAveraging.FullStep([right, wrong], [LogEvidence(x, y, 1.0), LogEvidence(x2, y, 1.0)]);
            ConformalSet[] sets =
            [
                averaged.Prediction.Sets[0],
                ConformalBayes.FullStep(right.Training, right.Candidates).Sets[0],
                ConformalBayes.FullStep(wrong.Training, wrong.Candidates).Sets[0],
            ];
            int nearest = Enumerable.Range(0, grid.Length).MinBy(g => Math.Abs(grid[g] - yt));
            for (int f = 0; f < 3; f++)
            {
                if (sets[f].Accepts(nearest)) covered[f]++;
                width[f] += sets[f].Width;
            }
            weight += averaged.ModelWeights[0][nearest, 0];
        }
        double se = Math.Sqrt(0.09 / replications);
        string[] names = ["CBMA", "βx alone", "βx² alone"];
        for (int f = 0; f < 3; f++)
            output.WriteLine($"{names[f]}: coverage {covered[f] / (double)replications:F4}, mean width {width[f] / replications:F3}");
        output.WriteLine($"CBMA's mean weight on βx at the truth: {weight / replications:F4}; se {se:F4}");
        Assert.True(covered[0] / (double)replications >= 0.9 - 3.5 * se);
        Assert.True(width[0] < 0.9 * width[2], "CBMA should track the right model, narrower than the wrong one's set");
        Assert.True(Math.Abs(width[0] - width[1]) < 0.05 * width[1], "CBMA should be near the right model's own set");
    }

    // ---- The residual scores ----

    [Fact]
    public void TheResidualScores_AreTheirDefinitions()
    {
        // 300 random residual vectors, some shifted far: each score against its two-pass definition.
        var rng = new Xoshiro256PlusPlus(2002);
        for (int c = 0; c < 300; c++)
        {
            int m = 2 + rng.NextInt(20);
            double shift = rng.NextDouble() < 0.3 ? 1e4 * Gaussian(rng) : Gaussian(rng);
            double[] r = [.. Enumerable.Range(0, m).Select(_ => shift + Gaussian(rng))];
            var got = new double[m];
            foreach (ResidualScore score in Enum.GetValues<ResidualScore>())
            {
                ResidualShrinkage.Score(r, score, got);
                for (int j = 0; j < m; j++)
                {
                    double[] others = [.. r.Where((_, i) => i != j)];
                    int n = others.Length;
                    double mean = others.Average(), sigma2 = others.Sum(v => (v - mean) * (v - mean)) / n;
                    double upsilon2 = Math.Max(mean * mean - sigma2 / n, 0.0), shrink = upsilon2 + sigma2 / n > 0 ? upsilon2 / (upsilon2 + sigma2 / n) : 0.0;
                    double expected = score switch
                    {
                        ResidualScore.DistanceToOrigin => Math.Abs(r[j]),
                        ResidualScore.DistanceToAverage => Math.Abs(r[j] - mean),
                        _ => Math.Abs(r[j] - shrink * mean),
                    };
                    Assert.True(Math.Abs(got[j] - expected) <= 1e-9 * (1.0 + Math.Abs(expected) + Math.Abs(shift)), $"{score}, case {c}, unit {j}: {got[j]} against {expected}");
                }
            }
        }
    }

    [Fact]
    public void TheResidualScores_Cover_AndShrinkageTracksTheBetterReference()
    {
        // A fixed predictor f(x) = 0 and residuals N(shift, 1): unbiased (shift 0) and biased (shift 3); n 15, α 0.1, 2,000 replications.
        const int n = 15, replications = 2000;
        var rng = new Xoshiro256PlusPlus(2026_1133);
        double[] grid = [.. Enumerable.Range(0, 401).Select(g => -8.0 + 0.05 * g)];
        ResidualScore[] scores = Enum.GetValues<ResidualScore>();
        foreach (double shift in new[] { 0.0, 3.0 })
        {
            int[] covered = new int[3];
            double[] width = new double[3];
            for (int r = 0; r < replications; r++)
            {
                double[] residuals = [.. Enumerable.Range(0, n).Select(_ => shift + Gaussian(rng))];
                double yt = shift + Gaussian(rng);
                int nearest = (int)Math.Round((yt + 8.0) / 0.05);
                for (int f = 0; f < 3; f++)
                {
                    ConformalSet set = ResidualShrinkage.Prediction(residuals, [0.0], [grid], scores[f], 0.1).Sets[0];
                    if (nearest >= 0 && nearest < grid.Length && set.Accepts(nearest)) covered[f]++;
                    width[f] += set.Width;
                }
            }
            double se = Math.Sqrt(0.09 / replications);
            output.WriteLine($"shift {shift}: " + string.Join("; ", Enumerable.Range(0, 3).Select(f => $"{scores[f]} coverage {covered[f] / (double)replications:F4}, width {width[f] / replications:F3}"))
                             + $"; se {se:F4}");
            for (int f = 0; f < 3; f++) Assert.True(covered[f] / (double)replications >= 0.9 - 3.5 * se, $"{scores[f]} at shift {shift}");
            double origin = width[0], average = width[1], shrinkage = width[2];
            if (shift == 0.0) Assert.True(shrinkage < average, "unbiased: shrinkage should beat distance to average");
            else Assert.True(shrinkage < origin, "shifted: shrinkage should beat distance to origin");
        }
    }

    // ---- Wu et al.'s high-density sets by decision risk ----

    [Fact]
    public void HighDensitySets_CoverWithConfidenceBeta()
    {
        // A fit to 50 training units of y = 0.7x + N(0, 1), 200 exact draws; 60 calibration units; the set's own miscoverage read on 1,000
        // fresh units; 300 calibration sets, α 0.1, β 0.9: the miscoverage is at most α in at least β of them.
        const int replications = 300, tests = 1000;
        const double alpha = 0.1, beta = 0.9;
        var rng = new Xoshiro256PlusPlus(2026_1134);
        double[] xTrain = [.. Enumerable.Range(0, 50).Select(_ => Gaussian(rng))], yTrain = [.. xTrain.Select(v => 0.7 * v + Gaussian(rng))];
        double[] draws = Draws(xTrain, yTrain, 1.0, 200, rng);
        double[] xTest = [.. Enumerable.Range(0, tests).Select(_ => Gaussian(rng))], yTest = [.. xTest.Select(v => 0.7 * v + Gaussian(rng))];
        var candidates = new CandidateLogDensity(AtResponses(xTest, yTest, draws, 1.0), [.. yTest.Select(v => new[] { v })], []);
        int held = 0;
        double miscoverage = 0.0;
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, 60).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.7 * v + Gaussian(rng))];
            var calibrationSeam = new CandidateLogDensity(AtResponses(x, y, draws, 1.0), [.. y.Select(v => new[] { v })], []);
            HighDensitySets sets = RiskControlledSets.Density(calibrationSeam, candidates, alpha, beta);
            double missed = Enumerable.Range(0, tests).Count(t => !sets.Accepts[t][0]) / (double)tests;
            miscoverage += missed;
            if (missed <= alpha) held++;
        }
        double rate = held / (double)replications, se = Math.Sqrt(beta * (1 - beta) / replications);
        output.WriteLine($"the set's miscoverage at most α in {rate:F4} of the calibration sets (β {beta}, se {se:F4}); mean miscoverage {miscoverage / replications:F4}");
        Assert.True(rate >= beta - 3.5 * se);
    }

    [Fact]
    public void WhatTheScoresDoNotServe_IsRefused()
    {
        var rng = new Xoshiro256PlusPlus(2003);
        double[] x = [.. Enumerable.Range(0, 10).Select(_ => Gaussian(rng))], y = [.. x.Select(v => v + Gaussian(rng))];
        double[] grid = [-1.0, 0.0, 1.0];
        var a = Seams(x, y, Draws(x, y, 1.0, 50, rng), 1.0, [0.2], grid);
        var b = Seams(x[..9], y[..9], Draws(x[..9], y[..9], 1.0, 50, rng), 1.0, [0.2], grid);
        Assert.Throws<ArgumentException>(() => ModelAveraging.FullStep([], []));
        Assert.Throws<ArgumentException>(() => ModelAveraging.FullStep([a], [0.0, 0.0]));
        Assert.Throws<ArgumentException>(() => ModelAveraging.FullStep([a, b], [0.0, 0.0]));
        Assert.Throws<ArgumentException>(() => ModelAveraging.FullStep([a], [double.NegativeInfinity]));
        Assert.Throws<ArgumentException>(() => ResidualShrinkage.Prediction([0.1, double.NaN], [0.0], [grid]));
        Assert.Throws<ArgumentException>(() => ResidualShrinkage.Prediction([0.1, 0.2], [0.0, 1.0], [grid]));
        Assert.Throws<ArgumentException>(() => ResidualShrinkage.Prediction([0.1, 0.2], [0.0], [[1.0, 0.0]]));
        Assert.Throws<ArgumentOutOfRangeException>(() => ResidualShrinkage.Prediction([0.1, 0.2], [0.0], [grid], (ResidualScore)7));
        var calibration = new CandidateLogDensity(a.Training.Values, [.. y.Select(v => new[] { v })], []);
        Assert.Throws<ArgumentOutOfRangeException>(() => RiskControlledSets.Density(calibration, a.Candidates, 0.1, 1.0));
    }
}
