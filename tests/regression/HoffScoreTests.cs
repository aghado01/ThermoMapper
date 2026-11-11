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

/// <summary>
/// Hoff's leave-one-out score as the full step's opt-in (E16 slice 12a; the user's decision 9 of 2026-11-10), against the Lean's
/// slice 12. Each test names its declaration.
/// <list type="bullet">
/// <item>On the discrete model with exact posterior weights: the score is Hoff's ratio p(z)/p(z₋ᵢ), the accept rule reads only the Bᵢ,
/// and CB-1 covers by enumeration with ties.</item>
/// <item>On the conjugate model with exact draws: the draw identity cᵢ = A/B̂ᵢ, the closed-form leave-one-out predictive as the
/// reference, and equivariance under permuted units.</item>
/// <item>The inflated rule's coverage under the shift form's bounds, and the refusals.</item>
/// </list>
/// With the qualification switch: decision 9's measurement, width and coverage against the in-sample score on the conjugate model and
/// on BASS fits, and the cost row.
/// </summary>
public sealed class HoffScoreTests(ITestOutputHelper output)
{
    private static readonly ConformalOptions Hoff = new() { Score = FullStepScore.LeaveOneOut };

    // ---- The discrete model, exact posterior weights (ConformalBayesTests' table) ----

    // Hoff's score in the ratio form, ln p(z) − ln p(z₋ᵢ), under the uniform prior over θ.
    private static double LogEvidence(IEnumerable<int> units)
    {
        int[] u = [.. units];
        double[] log = [.. Enumerable.Range(0, 3).Select(t => u.Sum(v => ConformalBayesTests.LogP(t, v)) - Math.Log(3))];
        double max = log.Max();
        return max + Math.Log(log.Sum(l => Math.Exp(l - max)));
    }

    [Fact]
    public void OnTheDiscreteModel_TheScoreIsHoffsRatio_AndTheRuleReadsOnlyTheB()
    {
        // hoffScore_eq, hoffScore_none, hoffScore_some, hoff_accepts_iff: with exact posterior weights the draws give Hoff's score itself.
        int checkedSets = 0;
        foreach (int[] units in ConformalBayesTests.Multisets(6, 4))
            for (int j = 0; j < units.Length; j++)
            {
                (PointwiseLogDensity training, CandidateLogDensity candidate) = ConformalBayesTests.Split(units, j);
                ConformalScores scores = ConformalBayes.Scores(training, candidate, 0, FullStepScore.LeaveOneOut);
                double all = LogEvidence(units);
                Assert.Equal(all - LogEvidence(units.Where((_, i) => i != j)), scores.Candidate, 12);
                int[] rest = [.. Enumerable.Range(0, units.Length).Where(i => i != j)];
                for (int a = 0; a < rest.Length; a++)
                    Assert.Equal(all - LogEvidence(units.Where((_, i) => i != rest[a])), scores.Scores[a], 12);
                // The rule from the Bᵢ: Bᵢ ≥ 1 exactly when unit i scores at most the candidate.
                double[] weights = ConformalBayesTests.Posterior(rest.Select(i => units[i]));
                for (int k = 1; k < units.Length; k++)
                {
                    ConformalSet set = ConformalBayes.FullStep(training, candidate, Hoff with { Alpha = (k + 0.5) / units.Length }).Sets[0];
                    int atLeastOne = rest.Count(i => Enumerable.Range(0, 3).Sum(t => weights[t] * Math.Exp(ConformalBayesTests.LogP(t, units[j]) - ConformalBayesTests.LogP(t, units[i]))) >= 1.0 - 1e-12);
                    Assert.Equal(atLeastOne >= k, set.Accepts(0));
                }
                checkedSets++;
            }
        output.WriteLine($"{checkedSets} (data set, candidate) pairs of four units");
    }

    [Fact]
    public void Cb1_ByEnumeration_WithTheLeaveOneOutScore()
    {
        // hoff_coverage, with coverage_orderings: over each data set's orderings the candidate is each unit equally often, so coverage is
        // the share accepted with each in turn as the candidate. At least 1 − k/m whatever the ties, exactly with distinct scores; equal
        // units tie exactly.
        int distinctSets = 0, tiedSets = 0, ambiguous = 0;
        for (int m = 2; m <= 5; m++)
            foreach (int[] units in ConformalBayesTests.Multisets(6, m))
            {
                var below = new int[m];
                var tied = new int[m];
                var scores = new double[m];
                for (int j = 0; j < m; j++)
                {
                    (PointwiseLogDensity training, CandidateLogDensity candidate) = ConformalBayesTests.Split(units, j);
                    ConformalSet set = ConformalBayes.FullStep(training, candidate, Hoff with { Alpha = (m - 0.5) / m }).Sets[0];
                    below[j] = set.Below![0];
                    tied[j] = set.Tied![0];
                    scores[j] = LogEvidence(units) - LogEvidence(units.Where((_, i) => i != j));
                }
                bool distinct = units.Distinct().Count() == m, near = false;
                for (int a = 0; a < m; a++)
                    for (int b = a + 1; b < m; b++)
                        if (units[a] != units[b] && Math.Abs(scores[a] - scores[b]) < 1e-9) near = true;
                if (near) { ambiguous++; continue; }
                if (distinct) distinctSets++; else tiedSets++;
                for (int a = 0; a < m; a++)
                    if (units.Count(u => u == units[a]) > 1) Assert.True(tied[a] >= 1, "equal units must tie exactly");
                for (int k = 1; k < m; k++)
                {
                    int accepted = 0;
                    for (int j = 0; j < m; j++) if (below[j] + tied[j] + 1 > k) accepted++;
                    Assert.True(accepted >= m - k, $"hoff_coverage: {accepted}/{m} at k = {k}, units [{string.Join(",", units)}]");
                    if (distinct) Assert.Equal(m - k, accepted);
                }
            }
        output.WriteLine($"data sets: {distinctSets} distinct, {tiedSets} with ties, {ambiguous} skipped with unequal units within 1e-9");
        Assert.True(distinctSets > 50 && tiedSets > 300 && ambiguous <= 10);
    }

    // ---- The conjugate model, exact draws: β ~ N(0, 10), y = βx + N(0, 1) ----

    private static (double[] X, double[] Y) Data(int n, Xoshiro256PlusPlus rng)
    {
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++) { x[i] = 4 * rng.NextDouble() - 2; y[i] = ConformalBayesTests.Truth(x[i], 0.5, rng); }
        return (x, y);
    }

    private static (PointwiseLogDensity Training, CandidateLogDensity Candidates) Seams(double[] x, double[] y, double[] beta, double xt, double[] grid)
    {
        int s = beta.Length, n = x.Length;
        var ll = new double[s, n];
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++) ll[d, i] = ConformalBayesTests.LogNormal(y[i], beta[d] * x[i]);
        var lc = new double[s, grid.Length];
        for (int d = 0; d < s; d++)
            for (int g = 0; g < grid.Length; g++) lc[d, g] = ConformalBayesTests.LogNormal(grid[g], beta[d] * xt);
        return (new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), ConformalBayesTests.Theta),
                new CandidateLogDensity(lc, [grid], ConformalBayesTests.Theta));
    }

    // ln N(y; Mx, 1 + x²/P), the predictive at x given units with precision P = 0.1 + Σx² and mean M = Σxy/P.
    private static double LogPredictive(double y, double x, IEnumerable<(double X, double Y)> given)
    {
        double precision = 0.1, xy = 0.0;
        foreach ((double a, double b) in given) { precision += a * a; xy += a * b; }
        double mean = xy / precision * x, variance = 1.0 + x * x / precision;
        return -0.5 * Math.Log(2 * Math.PI * variance) - 0.5 * (y - mean) * (y - mean) / variance;
    }

    [Fact]
    public void OnExactDraws_TheDrawIdentityHolds_AndMeetsTheClosedFormLeaveOneOutPredictive()
    {
        // hoffScore_none, hoffScore_some, postUExcept_some: the candidate scores A = Σ W̄ p(y | θ), each unit A/B̂ᵢ with B̂ᵢ = Σ W̄ p(y | θ)/p(zᵢ | θ),
        // and both meet their closed forms within a few of the layer's own errors.
        var rng = new Xoshiro256PlusPlus(seed: 1110);
        (double[] x, double[] y) = Data(19, rng);
        double[] beta = ConformalBayesTests.ExactDraws(x, y, 40_000, rng);
        double xt = 0.7;
        int worstUnit = 0;
        double worstStandardized = 0.0;
        foreach (double yt in new[] { -3.0, 0.4, 1.1, 5.0 })
        {
            (PointwiseLogDensity training, CandidateLogDensity candidates) = Seams(x, y, beta, xt, [yt]);
            ConformalScores scores = ConformalBayes.Scores(training, candidates, 0, FullStepScore.LeaveOneOut);
            double logA = LogSumExp(beta.Select(b => ConformalBayesTests.LogNormal(yt, b * xt))) - Math.Log(beta.Length);
            Assert.Equal(logA, scores.Candidate, 10);
            Assert.Equal(0.0, scores.CandidateError);
            double closedA = LogPredictive(yt, xt, x.Zip(y));
            Assert.True(Math.Abs(scores.Candidate - closedA) < 0.01, $"A at {yt}: {scores.Candidate:R} against {closedA:R}");
            for (int i = 0; i < x.Length; i++)
            {
                double logB = LogSumExp(beta.Select(b => ConformalBayesTests.LogNormal(yt, b * xt) - ConformalBayesTests.LogNormal(y[i], b * x[i]))) - Math.Log(beta.Length);
                Assert.Equal(logA - logB, scores.Scores[i], 9);
                var others = x.Zip(y).Where((_, j) => j != i).Append((xt, yt));
                double closed = LogPredictive(y[i], x[i], others), standardized = Math.Abs(scores.Scores[i] - closed) / scores.Errors[i];
                if (standardized > worstStandardized) { worstStandardized = standardized; worstUnit = i; }
                Assert.True(standardized < 5.0, $"unit {i} at {yt}: {scores.Scores[i]:R} against {closed:R}, error {scores.Errors[i]:E2}");
            }
            Assert.InRange(scores.EffectiveSampleSize, 1.0, beta.Length);
        }
        output.WriteLine($"the largest gap to the closed form, in the layer's own errors: {worstStandardized:F2} (unit {worstUnit})");
    }

    private static double LogSumExp(IEnumerable<double> values)
    {
        double[] v = [.. values];
        double max = v.Max();
        return max + Math.Log(v.Sum(z => Math.Exp(z - max)));
    }

    [Fact]
    public void PermutingTheUnits_PermutesTheScores_BitForBit()
    {
        // hoffScore_equivariant.
        var rng = new Xoshiro256PlusPlus(seed: 7);
        (double[] x, double[] y) = Data(15, rng);
        double[] beta = ConformalBayesTests.ExactDraws(x, y, 600, rng);
        int[] order = [.. Enumerable.Range(0, 15).OrderBy(_ => rng.NextDouble())];
        (PointwiseLogDensity a, CandidateLogDensity c) = Seams(x, y, beta, 0.3, [0.2, 2.5]);
        (PointwiseLogDensity b, CandidateLogDensity c2) = Seams([.. order.Select(i => x[i])], [.. order.Select(i => y[i])], beta, 0.3, [0.2, 2.5]);
        for (int g = 0; g < 2; g++)
        {
            ConformalScores sa = ConformalBayes.Scores(a, c, g, FullStepScore.LeaveOneOut), sb = ConformalBayes.Scores(b, c2, g, FullStepScore.LeaveOneOut);
            for (int j = 0; j < 15; j++) Assert.Equal(sa.Scores[order[j]], sb.Scores[j]);
            Assert.Equal(sa.Candidate, sb.Candidate);
        }
        Assert.Equal(ConformalBayes.FullStep(a, c, Hoff).Sets[0].PValues, ConformalBayes.FullStep(b, c2, Hoff).Sets[0].PValues);
    }

    [Fact]
    public void TheInflatedRule_CoversUnderTheShiftFormsBounds()
    {
        // coverage_inflated_shift: exact draws, few of them (S 40), n 19, α 0.1, k 2. The inflated set covers at least 1 − k/m − ν within its
        // binomial error, under per-position and uniform bounds; the candidate's position carries no bound, Â being common to every unit.
        const int n = 19, s = 40, replications = 4000;
        int[] seeds = SeedTree.Derive(2026_1110, replications);
        int plain = 0, inflated = 0, uniform = 0;
        foreach (int seed in seeds)
        {
            var rng = new Xoshiro256PlusPlus(seed);
            (double[] x, double[] y) = Data(n, rng);
            double xt = 4 * rng.NextDouble() - 2, yt = ConformalBayesTests.Truth(xt, 0.5, rng);
            (PointwiseLogDensity training, CandidateLogDensity candidates) = Seams(x, y, ConformalBayesTests.ExactDraws(x, y, s, rng), xt, [yt]);
            ConformalSet set = ConformalBayes.FullStep(training, candidates, Hoff).Sets[0];
            ConformalSet wide = ConformalBayes.FullStep(training, candidates, Hoff with { InflationBounds = InflationBounds.Uniform }).Sets[0];
            Assert.Equal(0.0, set.Margin!.CandidateErrors[0]);
            Assert.False(set.Margin.Computed);
            Assert.True(double.IsNaN(set.SplitHalfMass));
            if (set.AcceptsPlain(0)) plain++;
            if (set.AcceptsInflated(0)) inflated++;
            if (wide.AcceptsInflated(0)) uniform++;
            Assert.True(!set.AcceptsPlain(0) || set.AcceptsInflated(0), "the inflated set holds the plain one");
        }
        double se = Math.Sqrt(0.09 / replications);
        output.WriteLine($"plain {plain / (double)replications:F4}, inflated per position {inflated / (double)replications:F4}, uniform {uniform / (double)replications:F4} "
                         + $"(1 − k/m 0.9000, ν 0.01, se {se:F4})");
        Assert.True(inflated / (double)replications >= 0.89 - 3.5 * se);
        Assert.True(uniform >= inflated);
    }

    [Fact]
    public void WhatTheOptionDoesNotServe_IsRefusedWithItsReason()
    {
        var rng = new Xoshiro256PlusPlus(seed: 3);
        (double[] x, double[] y) = Data(12, rng);
        (PointwiseLogDensity training, CandidateLogDensity candidates) = Seams(x, y, ConformalBayesTests.ExactDraws(x, y, 100, rng), 0.5, [0.0, 1.0]);
        Assert.Contains("carried bounds", Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(training, candidates, Hoff with { MarginBounds = MarginBounds.Carried })).Message);
        Assert.Contains("fallback", Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(training, candidates, Hoff with { FallbackShare = 0.5 })).Message);
        Assert.Throws<ArgumentOutOfRangeException>(() => ConformalBayes.FullStep(training, candidates, new ConformalOptions { Score = (FullStepScore)5 }));
        // A unit with zero density under a retained draw: its weight is unbounded there.
        var ll = (double[,])training.Values.Clone();
        ll[17, 4] = double.NegativeInfinity;
        var holed = new PointwiseLogDensity(ll, UnitPartition.EachObservation(12), ConformalBayesTests.Theta);
        Assert.Contains("zero density", Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(holed, candidates, Hoff)).Message);
        // Off, nothing moves: the default and an explicit in-sample score agree bit for bit.
        Assert.Equal(ConformalBayes.FullStep(training, candidates).Sets[0].PValues,
                     ConformalBayes.FullStep(training, candidates, new ConformalOptions { Score = FullStepScore.InSample }).Sets[0].PValues);
    }

    [Fact]
    public void TheUnitsKHat_IsAnOptIn_NotReadByDefault()
    {
        // k̂ over the deciding units is read only when asked for: NaN by default, which Untrusted counts as none; asked for, a number or −∞
        // where the inflated rule accepts.
        var rng = new Xoshiro256PlusPlus(seed: 9);
        (double[] x, double[] y) = Data(19, rng);
        double[] grid = [.. Enumerable.Range(0, 81).Select(g => -4.0 + 0.1 * g)];
        (PointwiseLogDensity training, CandidateLogDensity candidates) = Seams(x, y, ConformalBayesTests.ExactDraws(x, y, 60, rng), 0.4, grid);
        ConformalSet off = ConformalBayes.FullStep(training, candidates, Hoff).Sets[0];
        ConformalSet on = ConformalBayes.FullStep(training, candidates, Hoff with { LeaveOneOutParetoK = true }).Sets[0];
        Assert.All(off.Margin!.ParetoK, k => Assert.True(double.IsNaN(k)));
        Assert.Equal(0, off.Margin.Untrusted(off.Level));
        int read = 0;
        for (int g = 0; g < grid.Length; g++)
        {
            Assert.Equal(off.PValues![g], on.PValues![g]);
            if (on.Margin!.InflatedRanks[g] > on.Level) { Assert.False(double.IsNaN(on.Margin.ParetoK[g])); read++; }
        }
        Assert.True(read > 10);
    }

    [Fact]
    public void ACandidateOutsideEveryDrawsSupport_ScoresMinusInfinity()
    {
        var rng = new Xoshiro256PlusPlus(seed: 4);
        (double[] x, double[] y) = Data(12, rng);
        (PointwiseLogDensity training, _) = Seams(x, y, ConformalBayesTests.ExactDraws(x, y, 50, rng), 0.5, [0.0]);
        var lc = new double[50, 2];
        for (int d = 0; d < 50; d++) { lc[d, 0] = double.NegativeInfinity; lc[d, 1] = -1.0; }
        ConformalSet set = ConformalBayes.FullStep(training, new CandidateLogDensity(lc, [[0.0, 1.0]], ConformalBayesTests.Theta), Hoff).Sets[0];
        Assert.Equal(double.NegativeInfinity, set.CandidateScores![0]);
        Assert.Equal(1.0 / 13, set.PValues![0]);
        Assert.Equal(0.0, set.EffectiveSampleSizes![0]);
    }

    [Fact]
    public void Scores_SeededSweep_FarRanges_WeightedDraws_Ties_StrataOfOne_AgainstReference()
    {
        // Each unit's ln cᵢ = ln A − ln Bᵢ and the candidate's ln A, against log-sum-exp references, over a seeded sweep: log terms down to
        // −700 (the inverse likelihoods up to e⁷⁰⁰, past the product's floor), draw weights far in a tail, a candidate equal to a unit
        // (an exact tie), and strata of one.
        var rng = new Xoshiro256PlusPlus(seed: 20261110);
        int ties = 0, inLogs = 0;
        for (int trial = 0; trial < 60; trial++)
        {
            int s = 5 + (int)(80 * rng.NextDouble()), n = 1 + (int)(12 * rng.NextDouble());
            double spread = trial % 3 == 0 ? 700.0 : trial % 3 == 1 ? 40.0 : 3.0;
            var ll = new double[s, n];
            for (int d = 0; d < s; d++)
                for (int i = 0; i < n; i++) ll[d, i] = -spread * rng.NextDouble();
            var lc = new double[s, 3];
            for (int d = 0; d < s; d++)
            {
                lc[d, 0] = -spread * rng.NextDouble();
                lc[d, 1] = ll[d, 0];   // equal to unit 0's column: an exact tie
                lc[d, 2] = d == 0 ? -1.0 : -800.0;   // one draw carries the candidate
            }
            double[]? weights = trial % 2 == 0 ? null : [.. Enumerable.Range(0, s).Select(_ => Math.Exp(-30.0 * rng.NextDouble()))];
            if (weights is not null) { double total = weights.Sum(); for (int d = 0; d < s; d++) weights[d] /= total; }
            int[] strata = [.. Enumerable.Range(0, n).Select(i => trial % 4 == 0 ? i : 0)];
            var training = new PointwiseLogDensity(ll, new UnitPartition([.. Enumerable.Range(0, n)], strata), ConformalBayesTests.Theta) { Weights = weights };
            var candidates = new CandidateLogDensity(lc, [[0.0, 1.0, 2.0]], ConformalBayesTests.Theta, [strata[0]]) { Weights = weights };
            double LogW(int d) => weights is null ? -Math.Log(s) : Math.Log(weights[d]);
            for (int c = 0; c < 3; c++)
            {
                ConformalScores scores = ConformalBayes.Scores(training, candidates, c, FullStepScore.LeaveOneOut);
                double logA = LogSumExp(Enumerable.Range(0, s).Select(d => LogW(d) + lc[d, c]));
                Assert.Equal(logA, scores.Candidate, 9);
                for (int a = 0; a < scores.Units.Length; a++)
                {
                    int i = scores.Units[a];
                    double logB = LogSumExp(Enumerable.Range(0, s).Select(d => LogW(d) + (lc[d, c] - ll[d, i])));
                    double expected = logA - logB, read = scores.Scores[a];
                    Assert.True(Math.Abs(read - expected) <= 1e-9 * Math.Max(1.0, Math.Abs(expected)), $"trial {trial}, candidate {c}, unit {i}: {read:R} against {expected:R}");
                    if (c == 1 && i == 0) { Assert.Equal(scores.Candidate, read); ties++; }
                    if (Math.Abs(expected) > 230) inLogs++;
                }
            }
        }
        output.WriteLine($"{ties} exact ties held; {inLogs} unit scores beyond e^±230");
        Assert.True(ties > 20);
    }

    // ---- Decision 9's measurement, with the qualification switch ----

    [Theory]
    [InlineData("outside", 400)]
    [InlineData("outside", 4000)]
    [InlineData("prior predictive", 400)]
    [InlineData("prior predictive", 4000)]
    public void Measurement_OnTheConjugateModel_WidthAndCoverageAgainstTheInSampleScore(string truth, int s)
    {
        // Both scores on the same draws: n 19, m 20, k 2, exact draws; coverage over 4,000 replications and widths over 1,000 more on a grid
        // of step 0.05, with the paired difference's standard error.
        // The truth is slice 0's, outside the model, or the model's own prior predictive (β ~ N(0, 10), y = βx + ζ), the setting of
        // hoff_expectedSize_le, where Hoff's set is the smallest of those covering as much on every data set given exact posterior weights.
        // Two S part the score's own width from its Monte Carlo error. A reading for decision 9: no default is proposed.
        if (!BassQualification.Enabled) return;
        const int n = 19, replications = 5000, widths = 1000;
        int[] seeds = SeedTree.Derive(truth == "outside" ? 2026_1111 : 2026_1113, replications);
        double[] grid = [.. Enumerable.Range(0, 601).Select(g => -12.0 + 0.05 * g)];
        int coveredIn = 0, coveredOut = 0, differ = 0;
        double widthIn = 0.0, widthOut = 0.0;
        var paired = new List<double>();
        for (int r = 0; r < replications; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            double[] x, y;
            double xt, yt;
            if (truth == "outside")
            {
                (x, y) = Data(n, rng);
                xt = 4 * rng.NextDouble() - 2;
                yt = ConformalBayesTests.Truth(xt, 0.5, rng);
            }
            else
            {
                double beta = Math.Sqrt(10.0) * Normal.SampleStandard(rng);
                x = [.. Enumerable.Range(0, n).Select(_ => 4 * rng.NextDouble() - 2)];
                y = [.. x.Select(v => beta * v + Normal.SampleStandard(rng))];
                xt = 4 * rng.NextDouble() - 2;
                yt = beta * xt + Normal.SampleStandard(rng);
            }
            double[] at = r < widths ? [.. grid.Select(g => g + (truth == "outside" ? 0.0 : Math.Round(yt)))] : [yt];
            (PointwiseLogDensity training, CandidateLogDensity candidates) = Seams(x, y, ConformalBayesTests.ExactDraws(x, y, s, rng), xt, at);
            ConformalSet inSample = ConformalBayes.FullStep(training, candidates).Sets[0], hoff = ConformalBayes.FullStep(training, candidates, Hoff).Sets[0];
            if (r < widths)
            {
                widthIn += inSample.Width;
                widthOut += hoff.Width;
                paired.Add(hoff.Width - inSample.Width);
                for (int g = 0; g < at.Length; g++) if (inSample.Accepts(g) != hoff.Accepts(g)) differ++;
                continue;
            }
            if (inSample.Accepts(0)) coveredIn++;
            if (hoff.Accepts(0)) coveredOut++;
        }
        int count = replications - widths;
        double mean = paired.Average(), pairedSe = Math.Sqrt(paired.Sum(d => (d - mean) * (d - mean)) / (widths - 1) / widths);
        output.WriteLine($"conjugate model, truth {truth}, S {s}: coverage in-sample {coveredIn / (double)count:F4}, leave-one-out {coveredOut / (double)count:F4} "
                         + $"(nominal 0.9000, se {Math.Sqrt(0.09 / count):F4}); mean width in-sample {widthIn / widths:F3}, leave-one-out {widthOut / widths:F3} "
                         + $"({100 * (widthOut / widthIn - 1):+0.00;-0.00}%, the paired difference {mean:+0.0000;-0.0000} with se {pairedSe:F4}); "
                         + $"grid values decided differently {differ} of {widths * grid.Length}");
        Assert.True(Math.Abs(coveredOut / (double)count - 0.9) <= 3.5 * Math.Sqrt(0.09 / count));
    }

    [Theory]
    [InlineData("a healthy fit")]
    [InlineData("a Student t fit")]
    [InlineData("a wrong fit")]
    public void Measurement_OnBassFits_WidthAndCoverageAgainstTheInSampleScore(string kind)
    {
        // Slice 4's fits: n 99 (m 100, k 10), 20 fits, 50 held-out cases each from a truth outside every model, and widths on 20 cases a
        // fit over a grid of step 0.04; both scores on the same fit and seams, with the least ESS and the largest k̂ at the true value.
        if (!BassQualification.Enabled) return;
        const int fits = 20, tests = 50, widths = 20;
        double[] grid = [.. Enumerable.Range(0, 201).Select(g => -4.0 + 0.04 * g)];
        int[] seeds = SeedTree.Derive(2026_1112, fits);
        int coveredIn = 0, coveredOut = 0, lowEss = 0, highK = 0;
        double widthIn = 0.0, widthOut = 0.0, leastEss = double.PositiveInfinity;
        long ticksIn = 0, ticksOut = 0;
        var perFit = new List<(double In, double Out)>();
        foreach (int seed in seeds)
        {
            var rng = new Xoshiro256PlusPlus(seed);
            (double[] x0, double[] x1, double[] y) = BassCandidateTests.Draw(99, rng);
            (double[] t0, double[] t1, double[] yt) = BassCandidateTests.Draw(tests + widths, rng);
            bool blind = kind == "a wrong fit";
            BassPredictor[] Predictors(double[] a, double[] b)
                => blind ? [new ContinuousPredictor("x1", b)] : [new ContinuousPredictor("x0", a), new ContinuousPredictor("x1", b)];
            var data = new BassData(Predictors(x0, x1), new ScalarResponse(y));
            BassFit fit = Bass.Fit(BassConformalDistributionTests.ConfigFor(kind), data);
            double[][] grids = [.. yt[..tests].Select(v => new[] { v }), .. Enumerable.Repeat(grid, widths)];
            PointwiseLogDensity training = fit.PointwiseLogDensity(data);
            CandidateLogDensity candidates = fit.CandidateLogDensity(data, PredictionCases.Of(Predictors(t0, t1)), grids);
            long a0 = System.Diagnostics.Stopwatch.GetTimestamp();
            ConformalPrediction inSample = ConformalBayes.FullStep(training, candidates);
            long a1 = System.Diagnostics.Stopwatch.GetTimestamp();
            ConformalPrediction hoff = ConformalBayes.FullStep(training, candidates, Hoff with { LeaveOneOutParetoK = true });
            long a2 = System.Diagnostics.Stopwatch.GetTimestamp();
            ticksIn += a1 - a0;
            ticksOut += a2 - a1;
            int fitIn = 0, fitOut = 0;
            for (int t = 0; t < tests; t++)
            {
                if (inSample.Sets[t].Accepts(0)) { coveredIn++; fitIn++; }
                if (hoff.Sets[t].Accepts(0)) { coveredOut++; fitOut++; }
                double ess = hoff.Sets[t].EffectiveSampleSizes![0];
                leastEss = Math.Min(leastEss, ess);
                if (ess < 50) lowEss++;
                if (hoff.Sets[t].Margin!.ParetoK[0] > hoff.Sets[t].Margin!.ParetoKThreshold) highK++;
            }
            perFit.Add((fitIn / (double)tests, fitOut / (double)tests));
            for (int t = tests; t < tests + widths; t++)
            {
                widthIn += inSample.Sets[t].Width;
                widthOut += hoff.Sets[t].Width;
            }
        }
        int count = fits * tests;
        double SeOf(Func<(double In, double Out), double> pick)
        {
            double mean = perFit.Average(pick);
            return Math.Sqrt(perFit.Sum(f => (pick(f) - mean) * (pick(f) - mean)) / (fits - 1) / fits);
        }
        output.WriteLine($"{kind}: coverage in-sample {coveredIn / (double)count:F4} (se {SeOf(f => f.In):F4}), leave-one-out {coveredOut / (double)count:F4} (se {SeOf(f => f.Out):F4}); "
                         + $"mean width in-sample {widthIn / (fits * widths):F3}, leave-one-out {widthOut / (fits * widths):F3} ({100 * (widthOut / widthIn - 1):+0.00;-0.00}%)");
        output.WriteLine($"{kind}: the leave-one-out weights at the true value: least ESS {leastEss:F1}, below 50 at {lowEss} of {count}; k̂ above its threshold at {highK} of "
                         + $"the inflated sets that held the truth; the full step's time over the fits {ticksIn * 1e3 / System.Diagnostics.Stopwatch.Frequency:F0} ms in-sample, "
                         + $"{ticksOut * 1e3 / System.Diagnostics.Stopwatch.Frequency:F0} ms leave-one-out");
    }

    // ---- The cost row ----

    [Theory]
    [InlineData(1000, 1000)]
    [InlineData(1000, 4000)]
    public void CostRow(int n, int s)
    {
        // The full step on one grid of 200 at n 1,000, both scores, and the leave-one-out scorer's own parts: its training side (the units'
        // inverse likelihoods, once per stratum) and its readings per grid. Cycle clock, best of five, indicative.
        if (!BassQualification.Enabled) return;
        var rng = new Xoshiro256PlusPlus(n + s);
        var ll = new double[s, n];
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++) ll[d, i] = -1 - 4 * rng.NextDouble();
        var lc = new double[s, 200];
        for (int d = 0; d < s; d++)
            for (int g = 0; g < 200; g++) lc[d, g] = -1 - 6 * rng.NextDouble() * (1 + g / 50.0);
        var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), ConformalBayesTests.Theta);
        var candidates = new CandidateLogDensity(lc, [[.. Enumerable.Range(0, 200).Select(g => (double)g)]], ConformalBayesTests.Theta);
        int[] rows = [.. Enumerable.Range(0, s)], members = [.. Enumerable.Range(0, n)];
        int level = (int)Math.Floor(0.1 * (n + 1));
        LeaveOneOutScorer scorer = null!;
        var reading = new ConformalBayes.Reading(200);
        double side = Best(() => scorer = new LeaveOneOutScorer(ll, rows, null, members, null, level, 3.5));
        double withoutK = Best(() => scorer.Rank(lc, 0, 200, reading));
        var bare = scorer;
        int deciding = 0, accepted = 0;
        for (int g = 0; g < 200; g++)
        {
            if (reading.InflatedRanks[g] <= level) continue;
            accepted++;
            for (int i = 0; i < n; i++)
                if (Math.Abs(bare.UnitScores[i * 200 + g] - reading.Own[g]) <= 3.5 * bare.UnitErrors[i * 200 + g]) deciding++;
        }
        double inSample = Best(() => ConformalBayes.FullStep(training, candidates));
        double hoff = Best(() => ConformalBayes.FullStep(training, candidates, Hoff));
        long k0 = BassCost.Tick();
        ConformalBayes.FullStep(training, candidates, Hoff with { LeaveOneOutParetoK = true });
        double withK = (BassCost.Tick() - k0) / BassCost.TicksPerSecond * 1e3;
        output.WriteLine($"n {n}, S {s}, G 200 ({(BassCost.Cycles ? "cycles" : "wall")} clock, indicative, best of 5): the full step in-sample {inSample:F1} ms, "
                         + $"leave-one-out {hoff:F1} ms ({withK:F0} ms once with k̂ over {deciding} deciding units at {accepted} accepted values); "
                         + $"the leave-one-out training side {side:F1} ms, its readings per grid {withoutK:F1} ms ({1e3 * withoutK / 200:F0} µs per candidate)");
    }

    private static double Best(Action call)
    {
        call();
        double best = double.PositiveInfinity;
        for (int r = 0; r < 5; r++)
        {
            long t0 = BassCost.Tick();
            call();
            best = Math.Min(best, (BassCost.Tick() - t0) / BassCost.TicksPerSecond * 1e3);
        }
        return best;
    }
}
