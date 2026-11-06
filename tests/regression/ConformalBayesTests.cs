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
/// The conformal layer's core (E16, slice 0) against its Lean contract (<c>lean/Enthymemata/ConformalValidity.lean</c>), each
/// test naming its declarations. CB-1 and its per-stratum form by enumeration on a discrete model whose posterior is held
/// exactly, as weighted draws; coverage on a conjugate model with exact draws under noise the model does not share; the small
/// stratum's whole space; CB-4's refusals, with the counterexample they guard; and the scores against a reference in logs over a
/// seeded sweep with far ranges, ties and strata of one.
/// </summary>
public sealed class ConformalBayesTests(ITestOutputHelper output)
{
    private static readonly LatentReport[] Theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];

    // ---- the discrete model: θ ∈ {0, 1, 2} under a uniform prior, a unit (x, y) with x ∈ {0, 1}, y ∈ {0, 1, 2} ----------------

    // p(y | x, θ), [θ, x, y]; every unit's column over θ is distinct, so only equal units tie by construction.
    private static readonly double[,,] Table =
    {
        { { 0.50, 0.30, 0.20 }, { 0.25, 0.45, 0.30 } },
        { { 0.20, 0.20, 0.60 }, { 0.55, 0.15, 0.30 } },
        { { 0.10, 0.65, 0.25 }, { 0.35, 0.30, 0.35 } },
    };

    private static double LogP(int theta, int unit) => Math.Log(Table[theta, unit / 3, unit % 3]);

    // The exact posterior over θ given the units, as the three draws' weights.
    private static double[] Posterior(IEnumerable<int> units)
    {
        var log = new double[3];
        foreach (int u in units)
            for (int t = 0; t < 3; t++) log[t] += LogP(t, u);
        double max = log.Max();
        double[] w = log.Select(l => Math.Exp(l - max)).ToArray();
        double total = w.Sum();
        return w.Select(v => v / total).ToArray();
    }

    // The full step with unit `test` of `units` as the candidate and the rest as training, each unit's stratum `strata[i]`.
    private static (PointwiseLogDensity Training, CandidateLogDensity Candidate) Split(int[] units, int test, int[]? strata = null)
    {
        int[] rest = [.. Enumerable.Range(0, units.Length).Where(i => i != test)];
        double[] weights = Posterior(rest.Select(i => units[i]));
        var ll = new double[3, rest.Length];
        for (int t = 0; t < 3; t++)
            for (int j = 0; j < rest.Length; j++) ll[t, j] = LogP(t, units[rest[j]]);
        var partition = new UnitPartition([.. Enumerable.Range(0, rest.Length)], strata is null ? null : [.. rest.Select(i => strata[i])]);
        var training = new PointwiseLogDensity(ll, partition, Theta) { Weights = weights };
        var lc = new double[3, 1];
        for (int t = 0; t < 3; t++) lc[t, 0] = LogP(t, units[test]);
        var candidate = new CandidateLogDensity(lc, [[units[test] % 3]], Theta, strata is null ? null : [strata[test]]) { Weights = weights };
        return (training, candidate);
    }

    private static IEnumerable<int[]> Multisets(int alphabet, int size, int from = 0)
    {
        if (size == 0)
        {
            yield return [];
            yield break;
        }
        for (int a = from; a < alphabet; a++)
            foreach (int[] rest in Multisets(alphabet, size - 1, a)) yield return [a, .. rest];
    }

    [Fact]
    public void Cb1_ByEnumeration_CoversAtLeastOneMinusKOverM_ExactlyWithDistinctScores()
    {
        // coverage_orderings, coverage_orderings_eq, coverage_orderings_le, coverage_strict_le, coverage_tieBreak, card_rankLe_le.
        // Over the orderings of a data set the test unit is each unit equally often, so coverage is the share of units accepted
        // when each in turn is the candidate; the posterior is exact, so the scores are the symmetric in-sample predictive.
        int distinctSets = 0, tiedSets = 0, straddling = 0, ambiguous = 0;
        for (int m = 2; m <= 5; m++)
            foreach (int[] units in Multisets(6, m))
            {
                var below = new int[m];
                var tied = new int[m];
                var scores = new double[m];
                for (int j = 0; j < m; j++)
                {
                    (PointwiseLogDensity training, CandidateLogDensity candidate) = Split(units, j);
                    ConformalSet set = ConformalBayes.FullStep(training, candidate, new ConformalOptions { Alpha = (m - 0.5) / m }).Sets[0];
                    Assert.Equal(m, set.Units);
                    Assert.Equal(m - 1, set.Level);
                    below[j] = set.Below![0];
                    tied[j] = set.Tied![0];
                    scores[j] = set.CandidateScores![0];
                    Assert.Equal((below[j] + tied[j] + 1) / (double)m, set.PValues![0]);
                }
                // Distinct: every pair of unequal units apart beyond rounding; equal units tie exactly, by the shared arithmetic.
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
                    int accepted = 0, strict = 0, atMostK = 0, straddle = 0;
                    double tieBroken = 0.0;
                    for (int j = 0; j < m; j++)
                    {
                        (PointwiseLogDensity training, CandidateLogDensity candidate) = Split(units, j);
                        ConformalSet set = ConformalBayes.FullStep(training, candidate, new ConformalOptions { Alpha = (k + 0.5) / m }).Sets[0];
                        Assert.Equal(k, set.Level);
                        int rankLe = below[j] + tied[j] + 1, rankLt = below[j];
                        Assert.Equal(rankLe > k, set.Accepts(0));
                        if (set.Accepts(0)) accepted++;
                        if (below[j] >= k) strict++;
                        if (rankLe <= k) atMostK++;
                        if (rankLt < k && k < rankLe) straddle++;
                        // A uniform position among the ties: the label order of coverage_tieBreak.
                        tieBroken += Math.Clamp((rankLe - k) / (double)(tied[j] + 1), 0.0, 1.0);
                    }
                    Assert.True(atMostK <= k, "card_rankLe_le");
                    Assert.True(accepted >= m - k, $"coverage_orderings: {accepted}/{m} at k = {k}, units [{string.Join(",", units)}]");
                    Assert.True(accepted <= m - k + straddle, "coverage_orderings_le");
                    Assert.True(strict <= m - k, "coverage_strict_le");
                    Assert.Equal(m - k, tieBroken, 9);
                    if (distinct) Assert.Equal(m - k, accepted);
                    if (straddle > 0) straddling++;
                }
            }
        output.WriteLine($"data sets: {distinctSets} distinct, {tiedSets} with ties, {ambiguous} with unequal units within 1e-9 (skipped); {straddling} (set, k) with straddling ties");
        Assert.True(distinctSets == 56 && tiedSets > 300 && straddling > 50);
        Assert.True(ambiguous <= 5);
    }

    [Fact]
    public void Cb1_PerStratum_RanksWithinTheCandidatesStratum()
    {
        // coverage_law with C the stratum, and coverage_stratum: the candidate's stratum (x = 0) is exchangeable within itself,
        // the score reads the other stratum's units through the posterior, and coverage is taken over orderings inside it.
        int[] other = [3, 5, 4];   // stratum 1: x = 1
        int checkedSets = 0;
        for (int m0 = 2; m0 <= 5; m0++)
            foreach (int[] ys in Multisets(3, m0))
            {
                int[] units = [.. ys, .. other];
                int[] strata = [.. units.Select(u => u / 3)];
                for (int k = 1; k < m0; k++)
                {
                    int accepted = 0;
                    var candidateScores = new List<double>();
                    for (int j = 0; j < m0; j++)
                    {
                        (PointwiseLogDensity training, CandidateLogDensity candidate) = Split(units, j, strata);
                        ConformalSet set = ConformalBayes.FullStep(training, candidate, new ConformalOptions { Alpha = (k + 0.5) / m0 }).Sets[0];
                        Assert.Equal(m0, set.Units);
                        Assert.Equal(0, set.Stratum);
                        if (set.Accepts(0)) accepted++;
                        candidateScores.Add(set.CandidateScores![0]);
                    }
                    Assert.True(accepted >= m0 - k, $"coverage_stratum: {accepted}/{m0} at k = {k}");
                    bool near = false;
                    for (int a = 0; a < m0; a++)
                        for (int b = a + 1; b < m0; b++)
                            if (ys[a] != ys[b] && Math.Abs(candidateScores[a] - candidateScores[b]) < 1e-9) near = true;
                    if (ys.Distinct().Count() == m0 && !near) Assert.Equal(m0 - k, accepted);
                    checkedSets++;
                }
            }
        Assert.True(checkedSets > 50);
    }

    [Fact]
    public void SmallStratum_IsTheWholeSpace_BelowCeilOneOverAlpha()
    {
        // CB-1's small strata: k = ⌊α·m⌋ is 0 below m = ⌈1/α⌉, and the set is everything, reported so, with nothing computed.
        foreach ((double alpha, int threshold) in new[] { (0.1, 10), (0.15, 7), (0.5, 2) })
            for (int m = 1; m <= threshold + 1; m++)
            {
                int n = m - 1;
                var ll = new double[2, Math.Max(n, 1)];
                for (int i = 0; i < ll.GetLength(1); i++) { ll[0, i] = -1.0 - 0.1 * i; ll[1, i] = -2.0 + 0.05 * i; }
                // Stratum 1 holds the n units the candidate is ranked against; a lone unit sits in stratum 0 when n = 0.
                var partition = new UnitPartition([.. Enumerable.Range(0, ll.GetLength(1))], [.. Enumerable.Repeat(n == 0 ? 0 : 1, ll.GetLength(1))]);
                var training = new PointwiseLogDensity(ll, partition, Theta);
                var candidate = new CandidateLogDensity(new double[,] { { -1.0, -3.0 }, { -1.5, -0.5 } }, [[0.0, 1.0]], Theta, [1]);
                ConformalSet set = ConformalBayes.FullStep(training, candidate, new ConformalOptions { Alpha = alpha }).Sets[0];
                Assert.Equal(m, set.Units);
                Assert.Equal(m < threshold, set.WholeSpace);
                if (set.WholeSpace)
                {
                    Assert.Equal(0, set.Level);
                    Assert.Null(set.PValues);
                    Assert.Empty(set.Intervals);
                    Assert.Equal(double.PositiveInfinity, set.Width);
                    Assert.True(set.Accepts(0) && set.Accepts(1));
                }
                else Assert.True(set.Level >= 1);
            }
    }

    // ---- CB-4 -----------------------------------------------------------------------------------------------------------------

    [Fact]
    public void Cb4_RefusesWhatBreaksTheScoresSymmetry()
    {
        var ll = new double[,] { { -1.0, -2.0 }, { -1.5, -1.0 } };
        var units = UnitPartition.EachObservation(2);
        var lc = new double[,] { { -1.2 }, { -0.9 } };
        var candidate = new CandidateLogDensity(lc, [[0.5]], Theta);
        CandidateLogDensity Candidate(IReadOnlyList<LatentReport> latents) => new(lc, [[0.5]], latents);
        var options = new ConformalOptions { Alpha = 0.4 };

        // The collapsed content: leave-unit-out terms (LeavesUnitsOut).
        LatentReport[] collapsed = [new("coefficients", LatentScope.Global, LatentTreatment.IntegratedGivenOthers)];
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(new PointwiseLogDensity(ll, units, collapsed), Candidate(collapsed), options));
        // An observation-level latent conditioned on, or mixed (mixed_counterexample, accepts_latent_of_mixed).
        foreach (LatentTreatment treatment in new[] { LatentTreatment.Conditioned, LatentTreatment.Mixed })
        {
            LatentReport[] latents = [.. Theta, new("v", LatentScope.Observation, treatment)];
            Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(new PointwiseLogDensity(ll, units, latents), Candidate(latents), options));
        }
        // Integrated, it is read; the candidate must report the same treatment.
        LatentReport[] integrated = [.. Theta, new("v", LatentScope.Observation, LatentTreatment.Integrated)];
        ConformalBayes.FullStep(new PointwiseLogDensity(ll, units, integrated), Candidate(integrated), options);
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(new PointwiseLogDensity(ll, units, integrated), candidate, options));
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(new PointwiseLogDensity(ll, units, integrated),
            Candidate([.. Theta, new("v", LatentScope.Observation, LatentTreatment.Conditioned)]), options));
        // Draws that do not align: their count, their chains, their weights.
        var training = new PointwiseLogDensity(ll, units, Theta);
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(training, new CandidateLogDensity(new double[,] { { -1.0 } }, [[0.5]], Theta), options));
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(training, candidate with { ChainLengths = [1, 1] }, options));
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(training, candidate with { Weights = [0.5, 0.5] }, options));
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(training with { Weights = [0.4, 0.6] }, candidate with { Weights = [0.5, 0.5] }, options));
        ConformalBayes.FullStep(training with { Weights = [0.4, 0.6] }, candidate with { Weights = [0.4, 0.6] }, options);
        Assert.Throws<ArgumentOutOfRangeException>(() => ConformalBayes.FullStep(training, candidate, new ConformalOptions { Alpha = 1.0 }));
    }

    [Fact]
    public void Cb4_TheMixedReading_TheRefusalGuards_CoversZero()
    {
        // mixed_counterexample: one θ, a fair latent bit v, y = v with probability 3/4, two units. The data unit read at its
        // retained latent (two draws of v from p(v | y), weights 3/4 and 1/4) scores 5/8; the candidate read integrated, 1/2. So
        // at k = 1 the candidate is never accepted, coverage 0 against 1/2. Labelled honestly the layer refuses it; the
        // integrated reading (integratedScore_equivariant) ties the two units and covers 1.
        var weights = new[] { 0.75, 0.25 };
        double[,] retained = { { Math.Log(0.75) }, { Math.Log(0.25) } };
        double[,] integrated = { { Math.Log(0.5) }, { Math.Log(0.5) } };
        double[,] lc = { { Math.Log(0.5) }, { Math.Log(0.5) } };
        LatentReport[] conditioned = [new("v", LatentScope.Observation, LatentTreatment.Conditioned)];
        LatentReport[] honest = [new("v", LatentScope.Observation, LatentTreatment.Integrated)];
        var units = UnitPartition.EachObservation(1);
        var options = new ConformalOptions { Alpha = 0.5 };
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(new PointwiseLogDensity(retained, units, conditioned) { Weights = weights },
            new CandidateLogDensity(lc, [[0.0]], conditioned) { Weights = weights }, options));

        ConformalSet mixed = ConformalBayes.FullStep(new PointwiseLogDensity(retained, units, honest) { Weights = weights },
            new CandidateLogDensity(lc, [[0.0]], honest) { Weights = weights }, options).Sets[0];
        ConformalScores scores = ConformalBayes.Scores(new PointwiseLogDensity(retained, units, honest) { Weights = weights },
            new CandidateLogDensity(lc, [[0.0]], honest) { Weights = weights }, 0);
        Assert.Equal(Math.Log(5.0 / 8), scores.Scores[0], 14);
        Assert.Equal(Math.Log(0.5), scores.Candidate, 14);
        Assert.False(mixed.Accepts(0));   // whichever unit is the candidate: coverage 0

        ConformalSet symmetric = ConformalBayes.FullStep(new PointwiseLogDensity(integrated, units, honest) { Weights = weights },
            new CandidateLogDensity(lc, [[0.0]], honest) { Weights = weights }, options).Sets[0];
        Assert.Equal(1, symmetric.Tied![0]);
        Assert.True(symmetric.Accepts(0));
    }

    // ---- coverage from truth outside the model ----------------------------------------------------------------------------------

    // y | x, β ~ N(βx, 1) with β ~ N(0, 10): the posterior of β given (x, y) and S exact draws from it.
    private static double[] ExactDraws(double[] x, double[] y, int s, Xoshiro256PlusPlus rng)
    {
        double precision = 0.1, xy = 0.0;
        for (int i = 0; i < x.Length; i++) { precision += x[i] * x[i]; xy += x[i] * y[i]; }
        double mean = xy / precision, sd = 1.0 / Math.Sqrt(precision);
        var beta = new double[s];
        for (int d = 0; d < s; d++) beta[d] = mean + sd * Normal.SampleStandard(rng);
        return beta;
    }

    private static double LogNormal(double y, double mean) => -0.5 * Math.Log(2 * Math.PI) - 0.5 * (y - mean) * (y - mean);

    // A skewed, heteroscedastic, curved truth the linear Gaussian model does not share: y = sin(1.5x) + (1 + scale·|x|)(E − 1).
    private static double Truth(double x, double scale, Xoshiro256PlusPlus rng)
        => Math.Sin(1.5 * x) + (1 + scale * Math.Abs(x)) * (-Math.Log(1 - rng.NextDouble()) - 1);

    [Fact]
    public void Coverage_ExactDraws_UnderNoiseTheModelDoesNotShare_WithinBinomialError()
    {
        // coverage_given / coverage_law on a continuous space, up to the draws' Monte Carlo error (CB-3): n = 19, m = 20, α = 0.1,
        // k = 2, so 1 − k/m = 0.9 exactly with distinct scores. Coverage reads π at the true y; width, a grid's set.
        const int n = 19, s = 400, replications = 4000, widths = 150;
        int[] seeds = SeedTree.Derive(2026_1106, replications);
        int covered = 0, bayesCovered = 0, intervals = 0, open = 0;
        double width = 0.0, bayesWidth = 0.0;
        double[] grid = [.. Enumerable.Range(0, 1001).Select(g => -25.0 + g / 20.0)];
        for (int r = 0; r < replications; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            var x = new double[n];
            var y = new double[n];
            for (int i = 0; i < n; i++) { x[i] = 4 * rng.NextDouble() - 2; y[i] = Truth(x[i], 0.5, rng); }
            double xt = 4 * rng.NextDouble() - 2, yt = Truth(xt, 0.5, rng);
            double[] beta = ExactDraws(x, y, s, rng);
            var ll = new double[s, n];
            for (int d = 0; d < s; d++)
                for (int i = 0; i < n; i++) ll[d, i] = LogNormal(y[i], beta[d] * x[i]);
            var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), Theta);
            double[] at = r < widths ? grid : [yt];
            var lc = new double[s, at.Length];
            for (int d = 0; d < s; d++)
                for (int g = 0; g < at.Length; g++) lc[d, g] = LogNormal(at[g], beta[d] * xt);
            ConformalSet set = ConformalBayes.FullStep(training, new CandidateLogDensity(lc, [at], Theta)).Sets[0];
            Assert.Equal(2, set.Level);
            if (r >= widths && set.Accepts(0)) covered++;
            // The model's own 90% predictive interval, the mixture of N(βₛx, 1), read from its draws' moments.
            double mean = beta.Average() * xt, sd = Math.Sqrt(1 + beta.Select(b => b * xt).Select(v => (v - mean) * (v - mean)).Sum() / (s - 1));
            if (r >= widths && Math.Abs(yt - mean) <= 1.6449 * sd) bayesCovered++;
            if (r < widths)
            {
                width += set.Width;
                bayesWidth += 2 * 1.6449 * sd;
                intervals += set.Intervals.Count;
                if (set.OpenBelow || set.OpenAbove) open++;
            }
        }
        int count = replications - widths;
        double coverage = covered / (double)count, se = Math.Sqrt(0.9 * 0.1 / count);
        output.WriteLine($"conformal coverage {coverage:F4} (nominal 0.9000, se {se:F4}), mean width {width / widths:F3}, {intervals / (double)widths:F2} intervals per set, {open} open at the grid's ends");
        output.WriteLine($"the model's own 90% interval: coverage {bayesCovered / (double)count:F4}, mean width {bayesWidth / widths:F3}");
        Assert.True(Math.Abs(coverage - 0.9) <= 3.5 * se, $"coverage {coverage:F4}");
        Assert.Equal(0, open);
    }

    [Fact]
    public void Coverage_WithinEachStratum_WhereTheModelsWidthIsWrongInOne()
    {
        // coverage_stratum (Mondrian): x < 0 and x ≥ 0 as strata, 12 training units in each, noise 0.2 and 3 times the base. The
        // homoscedastic model's predictive is too wide in one stratum and too narrow in the other; ranked within the stratum, each
        // covers 1 − k/m = 1 − 2/13 within its binomial error.
        const int perStratum = 12, s = 300, replications = 3000;
        var lowEss = new int[2];
        int[] seeds = SeedTree.Derive(611, replications);
        var covered = new int[2];
        var tries = new int[2];
        var widthSum = new double[2];
        for (int r = 0; r < replications; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            int n = 2 * perStratum;
            var x = new double[n];
            var y = new double[n];
            var strata = new int[n];
            for (int i = 0; i < n; i++)
            {
                strata[i] = i < perStratum ? 0 : 1;
                x[i] = strata[i] == 0 ? -2 * rng.NextDouble() : 2 * rng.NextDouble();
                y[i] = Truth(x[i], 0.0, rng) * (strata[i] == 0 ? 0.2 : 3.0);
            }
            int h = r % 2;
            double xt = h == 0 ? -2 * rng.NextDouble() : 2 * rng.NextDouble(), yt = Truth(xt, 0.0, rng) * (h == 0 ? 0.2 : 3.0);
            double[] beta = ExactDraws(x, y, s, rng);
            var ll = new double[s, n];
            for (int d = 0; d < s; d++)
                for (int i = 0; i < n; i++) ll[d, i] = LogNormal(y[i], beta[d] * x[i]);
            var training = new PointwiseLogDensity(ll, new UnitPartition([.. Enumerable.Range(0, n)], strata), Theta);
            double[] at = r < 200 ? [.. Enumerable.Range(0, 401).Select(g => -12.0 + g * 0.06)] : [yt];
            var lc = new double[s, at.Length];
            for (int d = 0; d < s; d++)
                for (int g = 0; g < at.Length; g++) lc[d, g] = LogNormal(at[g], beta[d] * xt);
            ConformalSet set = ConformalBayes.FullStep(training, new CandidateLogDensity(lc, [at], Theta, [h]), new ConformalOptions { Alpha = 0.2 }).Sets[0];
            Assert.Equal(13, set.Units);
            Assert.Equal(2, set.Level);
            if (r < 200) { widthSum[h] += set.Width; continue; }
            tries[h]++;
            if (set.Accepts(0)) covered[h]++;
            if (set.EffectiveSampleSizes![0] < 30) lowEss[h]++;
        }
        double nominal = 1 - 2.0 / 13;
        for (int h = 0; h < 2; h++)
        {
            double coverage = covered[h] / (double)tries[h], se = Math.Sqrt(nominal * (1 - nominal) / tries[h]);
            output.WriteLine($"stratum {h}: coverage {coverage:F4} (nominal {nominal:F4}, se {se:F4}), mean width {widthSum[h] / 100:F3}, ESS below 30 at {lowEss[h]} of {tries[h]} true values");
            Assert.True(Math.Abs(coverage - nominal) <= 3.5 * se, $"stratum {h}: coverage {coverage:F4}");
        }
        Assert.True(widthSum[1] > 2 * widthSum[0], "the noisy stratum's sets are the wider");
    }

    // ---- the scores: chosen cases and a seeded sweep against a reference in logs --------------------------------------------------

    private static double LogSumExp(IEnumerable<double> values)
    {
        double[] v = [.. values];
        double max = v.Max();
        return double.IsNegativeInfinity(max) ? max : max + Math.Log(v.Sum(x => Math.Exp(x - max)));
    }

    // Reference: w̃_s = W_s p(y | θ_s) normalized, every unit's ln Σ_s w̃_s p(z | θ_s), the candidate's own, and 1/Σ w̃², in logs.
    private static (double[] Units, double Own, double Ess) Reference(double[,] ll, int[] members, double[,] lc, int column, double[]? weights)
    {
        int s = ll.GetLength(0);
        int[] rows = [.. Enumerable.Range(0, s).Where(d => weights is null || weights[d] > 0)];
        double[] logW = [.. rows.Select(d => (weights is null ? 0.0 : Math.Log(weights[d])) + lc[d, column])];
        double norm = LogSumExp(logW);
        double[] logWt = [.. logW.Select(v => v - norm)];
        double[] units = [.. members.Select(u => LogSumExp(rows.Select((d, j) => logWt[j] + ll[d, u])))];
        double own = LogSumExp(rows.Select((d, j) => logWt[j] + lc[d, column]));
        double ess = 1.0 / logWt.Sum(v => Math.Exp(2 * v));
        return (units, own, ess);
    }

    [Fact]
    public void Scores_SeededSweep_FarRanges_Ties_StrataOfOne_AgainstReference()
    {
        var rng = new Xoshiro256PlusPlus(2026_11_06);
        int farBelow = 0, ties = 0;
        for (int trial = 0; trial < 120; trial++)
        {
            bool far = trial % 3 == 2;
            int s = 1 + rng.NextInt(far ? 40 : 300), n = 1 + rng.NextInt(20), cases = 1 + rng.NextInt(3), strataCount = 1 + rng.NextInt(3);
            var ll = new double[s, n];
            for (int d = 0; d < s; d++)
                for (int i = 0; i < n; i++)
                    ll[d, i] = far ? -1e4 + 900 * rng.NextDouble() * (i % 2 == 0 ? 1 : 0) - 50 * rng.NextDouble() : -10 * rng.NextDouble();
            int[] unitStrata = [.. Enumerable.Range(0, n).Select(_ => rng.NextInt(strataCount))];
            double[]? weights = trial % 4 == 1 ? [.. Enumerable.Range(0, s).Select(d => d % 5 == 0 && s > 1 ? 0.0 : Math.Exp(-700 * rng.NextDouble()))] : null;
            if (weights is not null && weights.All(w => w == 0)) weights[0] = 1.0;
            int[]? chains = weights is null && s >= 2 && trial % 4 == 3 ? [s / 2, s - s / 2] : null;
            double[][] grids = [.. Enumerable.Range(0, cases).Select(_ => Enumerable.Range(0, 1 + rng.NextInt(12)).Select(g => (double)g).ToArray())];
            int total = grids.Sum(g => g.Length);
            var lc = new double[s, total];
            for (int c = 0; c < total; c++)
            {
                int kind = rng.NextInt(6);
                for (int d = 0; d < s; d++)
                    lc[d, c] = kind == 0 ? ll[d, c % n]                                                // a unit's own column: ties
                             : kind == 1 && far ? -5000 - 800 * rng.NextDouble()                        // far in a tail
                             : kind == 2 && c == 0 && trial % 7 == 0 ? double.NegativeInfinity          // outside every draw's support
                             : far ? -1e4 + 900 * rng.NextDouble() : -10 * rng.NextDouble();
            }
            int[] caseStrata = [.. Enumerable.Range(0, cases).Select(_ => rng.NextInt(strataCount + 1))];   // a stratum with no unit too
            var training = new PointwiseLogDensity(ll, new UnitPartition([.. Enumerable.Range(0, n)], unitStrata), Theta) { ChainLengths = chains, Weights = weights };
            var candidates = new CandidateLogDensity(lc, grids, Theta, caseStrata, chains) { Weights = weights };
            var before = (double[,])ll.Clone();
            ConformalPrediction prediction = ConformalBayes.FullStep(training, candidates, new ConformalOptions { Alpha = 0.3 });
            Assert.Equal(before, ll);   // the layer reads; nothing upstream moves

            for (int t = 0; t < cases; t++)
            {
                ConformalSet set = prediction.Sets[t];
                int[] members = [.. Enumerable.Range(0, n).Where(i => unitStrata[i] == caseStrata[t])];
                Assert.Equal(members.Length + 1, set.Units);
                Assert.Equal((int)Math.Floor(0.3 * set.Units), set.Level);
                if (set.WholeSpace) continue;
                for (int g = 0; g < grids[t].Length; g++)
                {
                    int column = candidates.Offset(t) + g;
                    ConformalScores scores = ConformalBayes.Scores(training, candidates, column);
                    Assert.Equal(members, scores.Units);
                    Assert.Equal(BitConverter.DoubleToInt64Bits(set.CandidateScores![g]), BitConverter.DoubleToInt64Bits(scores.Candidate));
                    if (double.IsNegativeInfinity(scores.Candidate))
                    {
                        Assert.Equal(1.0 / set.Units, set.PValues![g]);
                        Assert.Equal(0.0, set.EffectiveSampleSizes![g]);
                        continue;
                    }
                    (double[] unitRef, double ownRef, double essRef) = Reference(ll, members, lc, column, weights);
                    double scale = Math.Max(1.0, Math.Abs(ownRef));
                    Assert.True(Math.Abs(scores.Candidate - ownRef) <= 1e-12 * scale, $"own {scores.Candidate:R} against {ownRef:R}");
                    Assert.True(Math.Abs(set.EffectiveSampleSizes![g] - essRef) <= 1e-9 * essRef);
                    int below = 0, tied = 0;
                    for (int i = 0; i < members.Length; i++)
                    {
                        Assert.True(Math.Abs(scores.Scores[i] - unitRef[i]) <= 1e-12 * Math.Max(1.0, Math.Abs(unitRef[i])), $"unit {members[i]}: {scores.Scores[i]:R} against {unitRef[i]:R}");
                        if (scores.Scores[i] < scores.Candidate) below++;
                        else if (scores.Scores[i] == scores.Candidate) tied++;
                        if (unitRef[i] - ownRef < -400) farBelow++;
                        // A candidate column equal to a unit's ties with it exactly.
                        bool equal = true;
                        for (int d = 0; d < s && equal; d++) equal = lc[d, column] == ll[d, members[i]];
                        if (equal)
                        {
                            Assert.Equal(scores.Candidate, scores.Scores[i]);
                            ties++;
                        }
                    }
                    Assert.Equal(below, set.Below![g]);
                    Assert.Equal(tied, set.Tied![g]);
                    Assert.Equal((below + tied + 1) / (double)set.Units, set.PValues![g]);
                }
            }
        }
        output.WriteLine($"{ties} exact ties, {farBelow} unit scores 400 nats or more below the candidate's");
        Assert.True(ties > 20);
    }

    [Fact]
    public void Scores_FarBelowTheFloor_AreSummedInLogs()
    {
        // A unit whose likelihood falls 800 nats across the draws, scored by a candidate whose weights sit where it is lowest: the
        // product's entry is below 10⁻²⁰⁰ of the unit's largest likelihood, and the score is the reference's in logs.
        const int s = 50;
        var ll = new double[s, 2];
        var lc = new double[s, 1];
        for (int d = 0; d < s; d++)
        {
            ll[d, 0] = -800.0 * d / (s - 1);
            ll[d, 1] = -1.0;
            lc[d, 0] = -20.0 * (s - 1 - d);   // weight on the last draws, where unit 0 is lowest
        }
        var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(2), Theta);
        var candidates = new CandidateLogDensity(lc, [[0.0]], Theta);
        ConformalScores scores = ConformalBayes.Scores(training, candidates, 0);
        (double[] unitRef, double ownRef, _) = Reference(ll, [0, 1], lc, 0, null);
        Assert.True(unitRef[0] < -790);
        Assert.Equal(unitRef[0], scores.Scores[0], 9);
        Assert.Equal(unitRef[1], scores.Scores[1], 12);
        Assert.Equal(ownRef, scores.Candidate, 12);
    }

    [Fact]
    public void Set_OnABimodalPredictive_IsAUnionOfIntervals()
    {
        // One draw whose density is an even mixture of N(−5, 1) and N(5, 1): each unit scores its own density, so the set is
        // {y : p(y) above the k-th lowest unit's}, two intervals around ±5 on the grid.
        static double LogMix(double y) => Math.Log(0.5 * Math.Exp(-0.5 * (y + 5) * (y + 5)) + 0.5 * Math.Exp(-0.5 * (y - 5) * (y - 5))) - 0.5 * Math.Log(2 * Math.PI);
        double[] data = [-6.2, -5.4, -5.1, -4.9, -4.0, -3.3, 3.6, 4.2, 4.8, 5.0, 5.3, 6.1, 6.8, -5.8, 5.6, -4.4, 4.5, -6.6, 6.4];
        var ll = new double[1, data.Length];
        for (int i = 0; i < data.Length; i++) ll[0, i] = LogMix(data[i]);
        double[] grid = [.. Enumerable.Range(0, 201).Select(g => -10.0 + 0.1 * g)];
        var lc = new double[1, grid.Length];
        for (int g = 0; g < grid.Length; g++) lc[0, g] = LogMix(grid[g]);
        ConformalSet set = ConformalBayes.FullStep(new PointwiseLogDensity(ll, UnitPartition.EachObservation(data.Length), Theta),
            new CandidateLogDensity(lc, [grid], Theta), new ConformalOptions { Alpha = 0.1 }).Sets[0];
        Assert.Equal(20, set.Units);
        Assert.Equal(2, set.Level);
        IReadOnlyList<GridInterval> intervals = set.Intervals;
        Assert.Equal(2, intervals.Count);
        Assert.True(intervals[0].Upper < 0 && intervals[1].Lower > 0);
        Assert.False(set.OpenBelow || set.OpenAbove);
        Assert.Equal(0.1, set.Resolution, 12);
        // The threshold is the third-lowest unit's density (|y| = 3.3 from the nearer mode, 1.7): the set is |y ∓ 5| < 1.7 within a grid step.
        Assert.InRange(intervals[1].Lower, 3.3 - 0.1, 3.3 + 0.1);
        Assert.InRange(set.Width, 2 * 2 * 1.7 - 0.4, 2 * 2 * 1.7 + 0.4);
        Assert.Equal(1.0, set.EffectiveSampleSizes![0]);
    }
}
