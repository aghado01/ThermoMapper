using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// CB-3 in the conformal layer (E16, slice 2, and the per-unit bounds' slice) against <c>lean/Enthymemata/ConformalValidity.lean</c>,
/// each test naming its declarations: the estimated band against the literal Lean count; the deterministic statements on perturbed
/// scores, with the adversary that attains them; the layer's errors, inflated ranks, bounds per position and carried over the
/// orderings, and bands against a reference read from its own scores; and the inflated set's coverage on exact draws of a model
/// the data do not share, under per-position and uniform bounds, with the margin's readings beside it.
/// </summary>
public sealed class ConformalMarginTests(ITestOutputHelper output)
{
    // ---- the Lean definitions, literally (Margin section) ---------------------------------------------------------------------

    private static int RankLe(double[] s, int u) => s.Count(v => v <= s[u]);

    private static bool Accepts(int k, double[] s, int u) => k < RankLe(s, u);

    private static int RankShift(double[] s, Func<int, double> delta, int u)
        => Enumerable.Range(0, s.Length).Count(j => j == u || s[j] + delta(j) <= s[u]);

    private static bool Fragile(int k, double[] s, double[] e, int u) => Accepts(k, s, u) && RankShift(s, j => e[j] + e[u], u) <= k;

    private static int InflatedRank(double[] s, double[] e, int u) => Enumerable.Range(0, s.Length).Count(j => s[j] - e[j] <= s[u] + e[u]);

    private static bool AcceptsInflated(int k, double[] s, double[] e, int u) => k < InflatedRank(s, e, u);

    private static bool FragileEstimated(int k, double[] s, double[] e, int u)
        => AcceptsInflated(k, s, e, u) && RankShift(s, j => 2 * e[j] + 2 * e[u], u) <= k;

    private static int Count(int m, Func<int, bool> holds) => Enumerable.Range(0, m).Count(holds);

    // Scores with ties: values drawn from a few levels half the time.
    private static double[] Scores(Xoshiro256PlusPlus rng, int m)
    {
        bool tied = rng.NextDouble() < 0.5;
        return [.. Enumerable.Range(0, m).Select(_ => tied ? rng.NextInt(4) * 0.5 : 3 * rng.NextDouble())];
    }

    [Fact]
    public void EstimatedBand_IsTheLeanCount_UnderAUniformBound()
    {
        // FragileEstimated, AcceptsInflated, rankShift: the layer's sorted reading equals the literal count, ties, a zero bound and
        // an infinite one included.
        var rng = new Xoshiro256PlusPlus(2026_1106);
        for (int trial = 0; trial < 3000; trial++)
        {
            int m = 2 + rng.NextInt(12), k = 1 + rng.NextInt(m - 1);
            double[] s = Scores(rng, m);
            double e = (trial % 10) switch { 0 => 0.0, 1 => double.PositiveInfinity, _ => rng.NextDouble() * 0.6 };
            double[] uniform = [.. Enumerable.Repeat(e, m)];
            int expected = Count(m, u => FragileEstimated(k, s, uniform, u));
            Assert.Equal(expected, ConformalBayes.EstimatedBand(s.ToArray(), e, k));
        }
    }

    [Fact]
    public void EstimatedBand_PerUnit_IsTheLeanCount()
    {
        // FragileEstimated under per-unit bounds, the carried reading: on continuous scores, where no comparison sits at rounding
        // distance, the sorted reading equals the literal count; with every bound equal it is the uniform reading, ties included.
        var rng = new Xoshiro256PlusPlus(707);
        for (int trial = 0; trial < 3000; trial++)
        {
            int m = 2 + rng.NextInt(12), k = 1 + rng.NextInt(m - 1);
            double[] s = [.. Enumerable.Range(0, m).Select(_ => 3 * rng.NextDouble())];
            double[] e = [.. Enumerable.Range(0, m).Select(_ => rng.NextDouble() < 0.2 ? 0.0 : 0.5 * rng.NextDouble())];
            Assert.Equal(Count(m, u => FragileEstimated(k, s, e, u)), ConformalBayes.EstimatedBand(s, e, k));
            double[] tied = Scores(rng, m);
            double one = rng.NextDouble() < 0.3 ? 0.0 : 0.4 * rng.NextDouble();
            Assert.Equal(ConformalBayes.EstimatedBand(tied.ToArray(), one, k), ConformalBayes.EstimatedBand(tied, [.. Enumerable.Repeat(one, m)], k));
        }
    }

    [Fact]
    public void Cb3_OnPerturbedScores_TheStatementsHold_AndTheAdversaryAttainsTheBand()
    {
        // accepts_or_fragile, acceptsInflated_of_accepts, accepts_shift_of_accepts_close, fragileEstimated_of_fragile,
        // card_fragile_le_estimated (one ordering) and fragileEstimated_mono (the uniform bound the layer reads): exact scores s,
        // per-unit bounds, estimates anywhere within them. The adversary, the unit down by its bound and every other up by its
        // own, refuses exactly the fragile units of s.
        var rng = new Xoshiro256PlusPlus(611);
        int lost = 0, fragileSeen = 0;
        double bandOverM = 0.0;
        for (int trial = 0; trial < 4000; trial++)
        {
            int m = 2 + rng.NextInt(10), k = 1 + rng.NextInt(m - 1);
            double[] s = Scores(rng, m);
            double[] e = [.. Enumerable.Range(0, m).Select(_ => rng.NextDouble() < 0.2 ? 0.0 : 0.4 * rng.NextDouble())];
            double[] est = [.. Enumerable.Range(0, m).Select(j => s[j] + e[j] * (2 * rng.NextDouble() - 1))];
            for (int u = 0; u < m; u++)
            {
                if (Accepts(k, s, u) && !Accepts(k, est, u)) { Assert.True(Fragile(k, s, e, u), "accepts_or_fragile"); lost++; }
                if (Accepts(k, s, u)) Assert.True(AcceptsInflated(k, est, e, u), "acceptsInflated_of_accepts");
                if (Accepts(k, est, u)) Assert.True(k < Count(m, j => s[j] <= s[u] + (e[j] + e[u])), "accepts_shift_of_accepts_close");
                if (Fragile(k, s, e, u)) Assert.True(FragileEstimated(k, est, e, u), "fragileEstimated_of_fragile");
            }
            int fragile = Count(m, u => Fragile(k, s, e, u)), band = Count(m, u => FragileEstimated(k, est, e, u));
            Assert.True(fragile <= band, "card_fragile_le_estimated");
            int uniformBand = ConformalBayes.EstimatedBand(est.ToArray(), e.Max(), k);
            Assert.True(band <= uniformBand, "fragileEstimated_mono");
            fragileSeen += fragile;
            bandOverM += (double)uniformBand / m;
            // The adversary: unit u down by its bound, the rest up by theirs.
            for (int u = 0; u < m; u++)
            {
                if (!Accepts(k, s, u)) continue;
                double[] adverse = [.. Enumerable.Range(0, m).Select(j => j == u ? s[j] - e[j] : s[j] + e[j])];
                Assert.Equal(Fragile(k, s, e, u), !Accepts(k, adverse, u));
            }
        }
        output.WriteLine($"{lost} decisions lost under random estimates, every one fragile; {fragileSeen} fragile units in all; the uniform band averaged {bandOverM / 4000:F3} of m");
        Assert.True(lost > 100);
    }

    // ---- the layer's reading against a reference from its own scores -------------------------------------------------------------

    [Fact]
    public void Layer_ErrorsRanksBoundsAndBands_AgainstAReference_SeededSweep()
    {
        // The layer's Monte Carlo errors against Σ_s ψ_s² over r_eff read directly in logs; its inflated rank, uniform bound and band
        // against the literal definitions on its own scores and bounds (z times the errors), far ranges, weights and chains included.
        var rng = new Xoshiro256PlusPlus(2026_1107);
        LatentReport[] theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];
        int checkedCandidates = 0, worstPath = 0;
        double worstError = 0.0;
        for (int trial = 0; trial < 80; trial++)
        {
            bool far = trial % 4 == 3;
            int s = 4 + rng.NextInt(200), n = 3 + rng.NextInt(25);
            var ll = new double[s, n];
            for (int d = 0; d < s; d++)
                for (int i = 0; i < n; i++) ll[d, i] = far ? -1e4 + 900 * rng.NextDouble() * (i % 2) - 30 * rng.NextDouble() : -6 * rng.NextDouble();
            double[]? weights = trial % 3 == 1 ? [.. Enumerable.Range(0, s).Select(d => d % 7 == 0 ? 0.0 : Math.Exp(-5 * rng.NextDouble()))] : null;
            int[]? chains = weights is null && trial % 3 == 2 ? [s / 2, s - s / 2] : null;
            int g = 1 + rng.NextInt(10);
            var lc = new double[s, g];
            for (int d = 0; d < s; d++)
                for (int c = 0; c < g; c++) lc[d, c] = c == 0 ? ll[d, 0] : far ? -1e4 + 900 * rng.NextDouble() : -6 * rng.NextDouble() * (1 + c);
            var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), theta) { ChainLengths = chains, Weights = weights };
            var candidates = new CandidateLogDensity(lc, [[.. Enumerable.Range(0, g).Select(c => (double)c)]], theta, chainLengths: chains) { Weights = weights };
            var options = new ConformalOptions { Alpha = 0.3, MarginBounds = MarginBounds.Carried };
            ConformalSet set = ConformalBayes.FullStep(training, candidates, options).Sets[0];
            if (set.WholeSpace) continue;
            ConformalSet uniformSet = ConformalBayes.FullStep(training, candidates,
                options with { MarginBounds = MarginBounds.Uniform, InflationBounds = InflationBounds.Uniform }).Sets[0];
            ConformalMargin margin = set.Margin!, uniformMargin = uniformSet.Margin!;
            Assert.Equal((MarginBounds.Carried, InflationBounds.PerPosition), (margin.Bounds, margin.Inflation));
            Assert.Equal((MarginBounds.Uniform, InflationBounds.Uniform), (uniformMargin.Bounds, uniformMargin.Inflation));
            double z = margin.Multiple;
            int[] rows = [.. Enumerable.Range(0, s).Where(d => weights is null || weights[d] > 0)];
            // The orderings in which a training unit is the candidate: its column read as a candidate's, the grid's candidates read
            // as extra units of the stratum so that their errors under its weights are read too.
            var aug = new double[s, n + g];
            for (int d = 0; d < s; d++)
            {
                for (int i = 0; i < n; i++) aug[d, i] = ll[d, i];
                for (int c = 0; c < g; c++) aug[d, n + c] = lc[d, c];
            }
            var augmented = new PointwiseLogDensity(aug, UnitPartition.EachObservation(n + g), theta) { ChainLengths = chains, Weights = weights };
            var unitWorst = new double[n];
            var candidateWorst = new double[g];
            for (int i = 0; i < n; i++)
            {
                var own = new double[s, 1];
                for (int d = 0; d < s; d++) own[d, 0] = ll[d, i];
                ConformalScores turn = ConformalBayes.Scores(augmented, new CandidateLogDensity(own, [[0.0]], theta, chainLengths: chains) { Weights = weights }, 0);
                for (int j = 0; j < n; j++) unitWorst[j] = Math.Max(unitWorst[j], z * Unbounded(turn.Errors[j]));
                for (int c = 0; c < g; c++) candidateWorst[c] = Math.Max(candidateWorst[c], z * Unbounded(turn.Errors[n + c]));
            }
            for (int j = 0; j < n; j++) AssertClose(unitWorst[j], margin.UnitCarriedErrors![j], $"unit {j}'s bound over the orderings");
            for (int c = 0; c < g; c++)
            {
                ConformalScores scores = ConformalBayes.Scores(training, candidates, c);
                // The reference error: w̃ ∝ W·p(y | θ), each estimate's Σ_s ψ_s², r_eff of the candidate's density along the chains.
                double[] logW = [.. rows.Select(d => (weights is null ? 0.0 : Math.Log(weights[d])) + lc[d, c])];
                double norm = LogSumExp(logW);
                double[] lw = [.. logW.Select(v => v - norm)];
                double top = rows.Max(d => lc[d, c]);
                double rEff = chains is null ? 1.0 : PredictiveChecks.RelativeEfficiency([.. rows.Select(d => Math.Exp(lc[d, c] - top))], chains);
                double ErrorOf(Func<int, double> logTerm)
                {
                    double est = LogSumExp(rows.Select((d, j) => lw[j] + logTerm(d)));
                    double relative = rows.Select((d, j) => PredictiveChecks.Influence(lw[j], logTerm(d), est)).Sum(t => t * t);
                    return PredictiveChecks.LogScaleError(relative, rEff);
                }
                double candidateError = ErrorOf(d => lc[d, c]);
                AssertClose(candidateError, scores.CandidateError, "the candidate's error");
                for (int i = 0; i < n; i++) AssertClose(ErrorOf(d => ll[d, i]), scores.Errors[i], $"unit {i}'s error");
                // The held ordering's per-position bounds: the inflated rule's, each its own or the largest for all.
                double[] all = [.. scores.Scores, scores.Candidate];
                double[] held = [.. scores.Errors.Select(x => z * x), z * scores.CandidateError];
                Assert.Equal(InflatedRank(all, held, n), margin.InflatedRanks[c]);
                Assert.Equal(held.Max(), margin.UniformErrors[c]);
                Assert.Equal(z * scores.CandidateError, margin.CandidateErrors[c]);
                Assert.Equal(InflatedRank(all, [.. Enumerable.Repeat(held.Max(), n + 1)], n), uniformMargin.InflatedRanks[c]);
                Assert.True(margin.InflatedRanks[c] <= uniformMargin.InflatedRanks[c], "no per-position set is wider than the uniform rule's");
                // The carried bounds: each unit's the larger of its held bound and its largest over the orderings, the candidate's
                // the larger of its own and its largest under the units' weights; the uniform bound the largest of them.
                AssertClose(Math.Max(held[n], candidateWorst[c]), margin.CandidateCarriedErrors![c], "the candidate's carried bound");
                double[] carried = [.. held[..n].Select((x, j) => Math.Max(x, margin.UnitCarriedErrors![j])), margin.CandidateCarriedErrors![c]];
                Assert.Equal(carried.Max(), margin.CarriedUniformErrors![c]);
                Assert.Equal(Count(n + 1, u => FragileEstimated(set.Level, all, carried, u)), margin.Fragile![c]);
                Assert.Equal(ConformalBayes.EstimatedBand(all.ToArray(), carried.Max(), set.Level), uniformMargin.Fragile![c]);
                Assert.Equal(Count(n + 1, u => FragileEstimated(set.Level, all, [.. Enumerable.Repeat(carried.Max(), n + 1)], u)), uniformMargin.Fragile![c]);
                Assert.True(margin.Fragile![c] <= uniformMargin.Fragile![c], "fragileEstimated_mono");
                Assert.Equal(set.Accepts(c), set.Below![c] + set.Tied![c] + 1 > set.Level);
                checkedCandidates++;
                worstError = Math.Max(worstError, Math.Abs(candidateError - scores.CandidateError));
                if (far) worstPath++;
            }
            Assert.Equal(set.Margin!.Fragile!.Max(), set.Margin.Count!.Value);
            Assert.Equal(1.0 - (double)(set.Level + set.Margin.Count!.Value) / set.Units - 0.01, set.CoverageUnderError, 12);
        }
        output.WriteLine($"{checkedCandidates} candidates checked, {worstPath} in far ranges; the largest error gap {worstError:E2}");
        Assert.True(checkedCandidates > 200);
    }

    // An error that does not exist, every weighted term at −∞, is unbounded, as the layer takes it.
    private static double Unbounded(double error) => double.IsNaN(error) ? double.PositiveInfinity : error;

    private static void AssertClose(double expected, double actual, string what)
        => Assert.True(expected == actual || Math.Abs(expected - actual) <= 1e-7 + 1e-6 * Math.Abs(expected), $"{what}: {actual:R} against {expected:R}");

    private static double LogSumExp(IEnumerable<double> values)
    {
        double[] v = [.. values];
        double max = v.Max();
        return double.IsNegativeInfinity(max) ? max : max + Math.Log(v.Sum(x => Math.Exp(x - max)));
    }

    [Fact]
    public void FallBack_IsOptIn_AndReadsTheGapBetweenTheInflatedAndPlainSets()
    {
        // The rule fires nothing unless a share is stated (the user's decision of 2026-11-07). Asked for, it reads the values whose
        // decision rests on the scores' error, those only the inflated rule accepts: it fires when they exceed the share of the plain
        // set (Gap), or when a candidate in the inflated set has weights whose k̂ is untrusted. Draws from two chains that barely move
        // fire on the gap; six draws, too few for k̂'s tail, fire on trust with no gap; slice 0's conjugate fixture on 2,000 exact
        // draws, with Gaussian noise, does not fire. Every set's signal is the rule's at each share.
        LatentReport[] theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];
        var rng = new Xoshiro256PlusPlus(5);
        const int s = 6, n = 29;
        var ll = new double[s, n];
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++) ll[d, i] = -1 - 8 * rng.NextDouble();
        double[] grid = [.. Enumerable.Range(0, 41).Select(g => 0.25 * g)];
        var lc = new double[s, grid.Length];
        for (int d = 0; d < s; d++)
            for (int c = 0; c < grid.Length; c++) lc[d, c] = -1 - 8 * rng.NextDouble();
        var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), theta);
        var candidates = new CandidateLogDensity(lc, [grid], theta);
        var asked = new ConformalOptions { FallbackShare = 0.5 };
        Assert.False(ConformalBayes.FullStep(training, candidates).Sets[0].FallBack);
        ConformalSet strained = ConformalBayes.FullStep(training, candidates, asked).Sets[0];
        output.WriteLine($"six draws: plain {strained.Width:F2}, inflated {(strained with { Inflated = true }).Width:F2}, gap {strained.Gap:F3}, "
                         + $"untrusted {strained.Margin!.Untrusted(strained.Level)}, fallback {strained.FallBack}");
        Assert.True(strained.FallBack);
        Assert.True(strained.Margin!.Untrusted(strained.Level) > 0);

        // Two chains of 500 whose β moves as an AR(1) with ρ 0.995 about the exact posterior: the relative efficiency is small, so
        // the bounds are wide, and the inflated set holds far more than the plain one.
        var cx = new double[19];
        var cy = new double[19];
        var chainRng = new Xoshiro256PlusPlus(23);
        for (int i = 0; i < 19; i++) { cx[i] = 4 * chainRng.NextDouble() - 2; cy[i] = 0.7 * cx[i] + Maths.Distributions.Normal.SampleStandard(chainRng); }
        double precision = 0.1 + cx.Sum(v => v * v), center = cx.Zip(cy, (a, b) => a * b).Sum() / precision, spread = 1.0 / Math.Sqrt(precision);
        const int perChain = 500;
        var chain = new double[2 * perChain];
        for (int c = 0; c < 2; c++)
        {
            double z = Maths.Distributions.Normal.SampleStandard(chainRng);
            for (int d = 0; d < perChain; d++)
            {
                z = 0.995 * z + Math.Sqrt(1 - 0.995 * 0.995) * Maths.Distributions.Normal.SampleStandard(chainRng);
                chain[c * perChain + d] = center + spread * z;
            }
        }
        var chainLl = new double[2 * perChain, 19];
        for (int d = 0; d < 2 * perChain; d++)
            for (int i = 0; i < 19; i++) chainLl[d, i] = ConformalBayesTests.LogNormal(cy[i], chain[d] * cx[i]);
        double[] chainGrid = [.. Enumerable.Range(0, 301).Select(g => -6.0 + g / 25.0)];
        var chainLc = new double[2 * perChain, chainGrid.Length];
        for (int d = 0; d < 2 * perChain; d++)
            for (int g = 0; g < chainGrid.Length; g++) chainLc[d, g] = ConformalBayesTests.LogNormal(chainGrid[g], chain[d] * 1.8);
        var chainTraining = new PointwiseLogDensity(chainLl, UnitPartition.EachObservation(19), theta, [perChain, perChain]);
        var chainCandidates = new CandidateLogDensity(chainLc, [chainGrid], theta, chainLengths: [perChain, perChain]);
        Assert.False(ConformalBayes.FullStep(chainTraining, chainCandidates).Sets[0].FallBack);
        ConformalSet slow = ConformalBayes.FullStep(chainTraining, chainCandidates, asked).Sets[0];
        output.WriteLine($"two slow chains: plain {slow.Width:F2}, inflated {(slow with { Inflated = true }).Width:F2}, gap {slow.Gap:F3}, "
                         + $"untrusted {slow.Margin!.Untrusted(slow.Level)}, fallback {slow.FallBack}");
        Assert.True(slow.Gap > 0.5);
        Assert.True(slow.FallBack);
        foreach (ConformalSet set in new[] { strained, slow })
            for (int g = 0; g < set.Grid.Length; g++)
                if (set.AcceptsPlain(g)) Assert.True(set.AcceptsInflated(g), "the inflated set holds the plain one");

        var x = new double[19];
        var y = new double[19];
        var draw = new Xoshiro256PlusPlus(17);
        for (int i = 0; i < 19; i++) { x[i] = 4 * draw.NextDouble() - 2; y[i] = 0.7 * x[i] + Maths.Distributions.Normal.SampleStandard(draw); }
        double[] beta = ConformalBayesTests.ExactDraws(x, y, 2000, draw);
        var exactLl = new double[2000, 19];
        for (int d = 0; d < 2000; d++)
            for (int i = 0; i < 19; i++) exactLl[d, i] = ConformalBayesTests.LogNormal(y[i], beta[d] * x[i]);
        double[] wide = [.. Enumerable.Range(0, 301).Select(g => -6.0 + g / 25.0)];
        var exactLc = new double[2000, wide.Length];
        for (int d = 0; d < 2000; d++)
            for (int g = 0; g < wide.Length; g++) exactLc[d, g] = ConformalBayesTests.LogNormal(wide[g], beta[d] * 0.5);
        var healthyTraining = new PointwiseLogDensity(exactLl, UnitPartition.EachObservation(19), theta);
        var healthyCandidates = new CandidateLogDensity(exactLc, [wide], theta);
        ConformalSet healthy = ConformalBayes.FullStep(healthyTraining, healthyCandidates, asked).Sets[0];
        output.WriteLine($"2,000 exact draws: plain {healthy.Width:F2}, inflated {(healthy with { Inflated = true }).Width:F2}, gap {healthy.Gap:F3}, "
                         + $"untrusted {healthy.Margin!.Untrusted(healthy.Level)}, fallback {healthy.FallBack}");
        Assert.False(healthy.FallBack);
        foreach (double share in new[] { 0.0, 0.01, 0.25, 1.0 })
            foreach ((PointwiseLogDensity t, CandidateLogDensity c) in new[] { (training, candidates), (chainTraining, chainCandidates), (healthyTraining, healthyCandidates) })
            {
                ConformalSet set = ConformalBayes.FullStep(t, c, new ConformalOptions { FallbackShare = share }).Sets[0];
                Assert.Equal(set.Gap > share || set.Margin!.Untrusted(set.Level) > 0, set.FallBack);
                // The rule reads the gap and k̂, not the margin's band: asking for the margin moves nothing it reads.
                Assert.Equal(set.FallBack, ConformalBayes.FullStep(t, c, new ConformalOptions { FallbackShare = share, MarginBounds = MarginBounds.Carried }).Sets[0].FallBack);
            }
        Assert.Throws<ArgumentOutOfRangeException>(() => ConformalBayes.FullStep(training, candidates, new ConformalOptions { ErrorLevel = 0.0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => ConformalBayes.FullStep(training, candidates, new ConformalOptions { FallbackShare = -1 }));
    }

    [Fact]
    public void SplitHalf_IsTheReferencesReading_FromTwoHalfSteps()
    {
        // The layer's split-half diagnostic, from two half-size score products, against the same reading built from two public steps
        // on the half densities (ConformalFallbackComparisonTests.SplitHalfMass): exact draws split by alternate rows and chains by
        // halves, one case of two cases' grids each. The halves' decisions agree bit for bit, so the masses agree to rounding.
        LatentReport[] theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];
        var rng = new Xoshiro256PlusPlus(808);
        int nonzero = 0;
        for (int trial = 0; trial < 40; trial++)
        {
            int s = 20 + rng.NextInt(300), n = 9 + rng.NextInt(20);
            bool chains = trial % 2 == 1;
            var x = new double[n];
            var y = new double[n];
            for (int i = 0; i < n; i++) { x[i] = 4 * rng.NextDouble() - 2; y[i] = ConformalBayesTests.Truth(x[i], 0.5, rng); }
            double[] beta = ConformalBayesTests.ExactDraws(x, y, s, rng);
            var ll = new double[s, n];
            for (int d = 0; d < s; d++)
                for (int i = 0; i < n; i++) ll[d, i] = ConformalBayesTests.LogNormal(y[i], beta[d] * x[i]);
            double[] grid = [.. Enumerable.Range(0, 121).Select(g => -6.0 + g / 10.0)];
            double[] xs = [2 * rng.NextDouble() - 1, 2 * rng.NextDouble() - 1];
            var lc = new double[s, 2 * grid.Length];
            for (int d = 0; d < s; d++)
                for (int c = 0; c < 2; c++)
                    for (int g = 0; g < grid.Length; g++) lc[d, c * grid.Length + g] = ConformalBayesTests.LogNormal(grid[g], beta[d] * xs[c]);
            int[]? lengths = chains ? [s / 2, s - s / 2] : null;
            var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), theta, lengths);
            var candidates = new CandidateLogDensity(lc, [grid, grid], theta, chainLengths: lengths);
            ConformalPrediction p = ConformalBayes.FullStep(training, candidates);
            double[] reference = ConformalFallbackComparisonTests.SplitHalfMass(training, candidates);
            for (int c = 0; c < 2; c++)
            {
                ConformalSet set = p.Sets[c];
                if (set.WholeSpace) continue;
                Assert.True(Math.Abs(set.SplitHalfMass - reference[c]) <= 1e-12 + 1e-9 * reference[c], $"{set.SplitHalfMass:R} against {reference[c]:R}");
                Assert.Equal(set.SplitHalfDisagreements > 0, reference[c] > 0);
                if (reference[c] > 0) nonzero++;
            }
        }
        output.WriteLine($"{nonzero} cases with halves that disagree");
        Assert.True(nonzero > 10);
        // A single value has no width: its mass is NaN, its disagreement a count.
        var one = new PointwiseLogDensity(new double[,] { { -1, -2 }, { -1.5, -1 }, { -2, -1.2 }, { -1.1, -2.2 } }, UnitPartition.EachObservation(2), theta);
        ConformalSet single = ConformalBayes.FullStep(one, new CandidateLogDensity(new double[,] { { -1 }, { -2 }, { -1 }, { -2 } }, [[0.0]], theta),
                                                      new ConformalOptions { Alpha = 0.4 }).Sets[0];
        Assert.True(double.IsNaN(single.SplitHalfMass));
        // A half in which some unit has zero density under every draw does not split: the diagnostic is NaN, and the step stands.
        var split = new PointwiseLogDensity(new double[,] { { -1, double.NegativeInfinity }, { -1.5, -1 }, { -2, double.NegativeInfinity }, { -1.1, -2.2 } },
                                            UnitPartition.EachObservation(2), theta);
        ConformalSet stands = ConformalBayes.FullStep(split, new CandidateLogDensity(new double[,] { { -1, -2 }, { -2, -1 }, { -1, -1.5 }, { -2, -1 } }, [[0.0, 1.0]], theta),
                                                      new ConformalOptions { Alpha = 0.4 }).Sets[0];
        Assert.True(double.IsNaN(stands.SplitHalfMass));
        Assert.NotNull(stands.PValues);
    }

    [Fact]
    public void OneExactDraw_HasNoError_AndNoBand()
    {
        // One draw: the estimate is exact, every error 0, and the band empty, since no score sits at the threshold.
        LatentReport[] theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];
        var exact = new PointwiseLogDensity(new double[,] { { -1, -2, -3, -4, -5, -6, -7, -8, -9, -10 } }, UnitPartition.EachObservation(10), theta);
        ConformalSet sure = ConformalBayes.FullStep(exact, new CandidateLogDensity(new double[,] { { -2.5 } }, [[0.0]], theta),
                                                    new ConformalOptions { MarginBounds = MarginBounds.Carried }).Sets[0];
        Assert.Equal(0.0, sure.Margin!.CandidateErrors[0]);
        Assert.Equal(0.0, sure.Margin.UniformErrors[0]);
        Assert.Equal(0, sure.Margin.Count);
        Assert.Equal(0.0, sure.Gap);
    }

    // ---- the margin, opt-in -----------------------------------------------------------------------------------------------------

    [Fact]
    public void TheMargin_IsOptIn_AndEveryOtherReadingStandsWithoutIt()
    {
        // The user's decision of 2026-11-08, "Make it opt-in". By default the margin is not computed, and says so: no band, no count,
        // no carried bounds, and the plain set's coverage under error NaN, not a number read as clear. Asked for, under either bound,
        // it adds those and moves nothing else: every other reading, the inflated rule's, the gap, k̂, the split-half diagnostic and
        // the fallback signal included, is the same bit for bit, over seeded fixtures with weights, chains, far ranges and a
        // candidate outside every draw's support.
        var rng = new Xoshiro256PlusPlus(2026_1111);
        LatentReport[] theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];
        int compared = 0;
        for (int trial = 0; trial < 60; trial++)
        {
            bool far = trial % 4 == 3;
            int s = 4 + rng.NextInt(200), n = 3 + rng.NextInt(25), g = 2 + rng.NextInt(12);
            var ll = new double[s, n];
            for (int d = 0; d < s; d++)
                for (int i = 0; i < n; i++) ll[d, i] = far ? -1e4 + 900 * rng.NextDouble() * (i % 2) - 30 * rng.NextDouble() : -6 * rng.NextDouble();
            double[]? weights = trial % 3 == 1 ? [.. Enumerable.Range(0, s).Select(d => d % 7 == 0 ? 0.0 : Math.Exp(-5 * rng.NextDouble()))] : null;
            int[]? chains = weights is null && trial % 3 == 2 ? [s / 2, s - s / 2] : null;
            var lc = new double[s, g];
            for (int d = 0; d < s; d++)
                for (int c = 0; c < g; c++)
                    lc[d, c] = c == g - 1 && trial % 5 == 0 ? double.NegativeInfinity : c == 0 ? ll[d, 0] : far ? -1e4 + 900 * rng.NextDouble() : -6 * rng.NextDouble() * (1 + c);
            var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), theta) { ChainLengths = chains, Weights = weights };
            var candidates = new CandidateLogDensity(lc, [[.. Enumerable.Range(0, g).Select(c => (double)c)]], theta, chainLengths: chains) { Weights = weights };
            foreach (InflationBounds inflation in new[] { InflationBounds.PerPosition, InflationBounds.Uniform })
            {
                var asked = new ConformalOptions { Alpha = 0.3, FallbackShare = 0.25, InflationBounds = inflation };
                ConformalSet off = ConformalBayes.FullStep(training, candidates, asked).Sets[0];
                if (off.WholeSpace) continue;
                ConformalMargin none = off.Margin!;
                Assert.False(none.Computed);
                Assert.Null(none.Bounds);
                Assert.Null(none.Fragile);
                Assert.Null(none.Count);
                Assert.Null(none.CandidateCarriedErrors);
                Assert.Null(none.CarriedUniformErrors);
                Assert.Null(none.UnitCarriedErrors);
                Assert.True(double.IsNaN(off.CoverageUnderError), "not read, rather than read as clear");
                Assert.Equal(off.Coverage - off.ErrorLevel, (off with { Inflated = true }).CoverageUnderError);
                foreach (MarginBounds bounds in new[] { MarginBounds.Carried, MarginBounds.Uniform })
                {
                    ConformalSet on = ConformalBayes.FullStep(training, candidates, asked with { MarginBounds = bounds }).Sets[0];
                    ConformalMargin margin = on.Margin!;
                    Assert.True(margin.Computed);
                    Assert.Equal(bounds, margin.Bounds);
                    Assert.NotNull(margin.Fragile);
                    Assert.Equal(margin.Fragile!.Max(), margin.Count);
                    Assert.False(double.IsNaN(on.CoverageUnderError));
                    Assert.Equal((off.Level, off.Units), (on.Level, on.Units));
                    Assert.Equal(off.PValues, on.PValues);
                    Assert.Equal(off.Below, on.Below);
                    Assert.Equal(off.Tied, on.Tied);
                    Assert.Equal(off.CandidateScores, on.CandidateScores);
                    Assert.Equal(off.EffectiveSampleSizes, on.EffectiveSampleSizes);
                    Assert.Equal(none.Multiple, margin.Multiple);
                    Assert.Equal(none.CandidateErrors, margin.CandidateErrors);
                    Assert.Equal(none.UniformErrors, margin.UniformErrors);
                    Assert.Equal(none.InflatedRanks, margin.InflatedRanks);
                    Assert.Equal(none.ParetoK, margin.ParetoK);
                    Assert.Equal((none.ParetoKThreshold, none.Inflation), (margin.ParetoKThreshold, margin.Inflation));
                    Assert.Equal(off.Gap, on.Gap);
                    Assert.Equal(off.SplitHalfMass, on.SplitHalfMass);
                    Assert.Equal(off.SplitHalfDisagreements, on.SplitHalfDisagreements);
                    Assert.Equal(off.FallBack, on.FallBack);
                    Assert.Equal(off.Width, on.Width);
                    Assert.Equal((off with { Inflated = true }).Width, (on with { Inflated = true }).Width);
                    Assert.Equal((off with { Inflated = true }).CoverageUnderError, (on with { Inflated = true }).CoverageUnderError);
                    compared++;
                }
            }
        }
        output.WriteLine($"{compared} sets compared with the margin and without it, every other reading the same");
        Assert.True(compared > 150);
    }

    // ---- inflated coverage on exact draws ---------------------------------------------------------------------------------------

    [Fact]
    public void InflatedCoverage_OnExactDraws_UnderNoiseTheModelDoesNotShare()
    {
        // coverage_inflated and coverage_close: exact posterior draws, few of them (S 40), so the scores carry real Monte Carlo
        // error. n 19, α 0.1, k 2: the inflated set covers at least 1 − k/m − ν = 0.89 within its binomial error under per-position
        // and uniform bounds, never wider per position; the plain set at least the margin's reading, 1 − (k + m̂)/m − ν, averaged
        // over the replications, under carried and uniform bounds, the carried band never the larger.
        const int n = 19, s = 40, replications = 4000, widths = 150;
        int[] seeds = SeedTree.Derive(2026_1108, replications);
        LatentReport[] theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];
        int plain = 0, inflated = 0, inflatedUniform = 0, fallBacks = 0;
        double underError = 0.0, underUniform = 0.0, plainWidth = 0.0, inflatedWidth = 0.0, uniformWidth = 0.0, marginSum = 0.0, uniformSum = 0.0, untrustedSets = 0;
        double[] grid = [.. Enumerable.Range(0, 601).Select(g => -10.0 + g / 20.0)];
        var uniform = new ConformalOptions { MarginBounds = MarginBounds.Uniform, InflationBounds = InflationBounds.Uniform };
        for (int r = 0; r < replications; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            var x = new double[n];
            var y = new double[n];
            for (int i = 0; i < n; i++) { x[i] = 4 * rng.NextDouble() - 2; y[i] = ConformalBayesTests.Truth(x[i], 0.5, rng); }
            double xt = 4 * rng.NextDouble() - 2, yt = ConformalBayesTests.Truth(xt, 0.5, rng);
            double[] beta = ConformalBayesTests.ExactDraws(x, y, s, rng);
            var ll = new double[s, n];
            for (int d = 0; d < s; d++)
                for (int i = 0; i < n; i++) ll[d, i] = ConformalBayesTests.LogNormal(y[i], beta[d] * x[i]);
            var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), theta);
            double[] at = r < widths ? grid : [yt];
            var lc = new double[s, at.Length];
            for (int d = 0; d < s; d++)
                for (int g = 0; g < at.Length; g++) lc[d, g] = ConformalBayesTests.LogNormal(at[g], beta[d] * xt);
            var candidates = new CandidateLogDensity(lc, [at], theta);
            ConformalSet set = ConformalBayes.FullStep(training, candidates, new ConformalOptions { MarginBounds = MarginBounds.Carried }).Sets[0];
            ConformalSet same = ConformalBayes.FullStep(training, candidates, uniform).Sets[0];
            ConformalSet off = ConformalBayes.FullStep(training, candidates).Sets[0];
            Assert.False(off.Margin!.Computed);
            ConformalSet wide = set with { Inflated = true }, wider = same with { Inflated = true };
            for (int g = 0; g < at.Length; g++)
            {
                Assert.Equal(set.AcceptsPlain(g), same.AcceptsPlain(g));
                Assert.Equal((set.AcceptsPlain(g), set.AcceptsInflated(g)), (off.AcceptsPlain(g), off.AcceptsInflated(g)));
                if (wide.Accepts(g)) Assert.True(wider.Accepts(g), "no per-position set is wider than the uniform rule's");
                if (set.Margin!.Fragile![g] >= 0) Assert.True(set.Margin.Fragile[g] <= same.Margin!.Fragile![g], "fragileEstimated_mono");
            }
            if (r < widths)
            {
                plainWidth += set.Width;
                inflatedWidth += wide.Width;
                uniformWidth += wider.Width;
                marginSum += set.Margin!.Count!.Value;
                uniformSum += same.Margin!.Count!.Value;
                if (set.FallBack) fallBacks++;
                if (set.Margin.Untrusted(set.Level) > 0) untrustedSets++;
                continue;
            }
            if (set.Accepts(0)) plain++;
            if (wide.Accepts(0)) inflated++;
            if (wider.Accepts(0)) inflatedUniform++;
            underError += set.CoverageUnderError;
            underUniform += same.CoverageUnderError;
        }
        int count = replications - widths;
        double se = Math.Sqrt(0.09 / count), plainRate = plain / (double)count, inflatedRate = inflated / (double)count, uniformRate = inflatedUniform / (double)count;
        double reading = underError / count, uniformReading = underUniform / count;
        output.WriteLine($"plain coverage {plainRate:F4}; inflated per position {inflatedRate:F4}, uniform {uniformRate:F4} (1 − k/m 0.9000, ν 0.01, se {se:F4}); "
                         + $"the margin's reading at the true value averaged {reading:F4} carried, {uniformReading:F4} uniform");
        output.WriteLine($"widths: plain {plainWidth / widths:F3}, inflated per position {inflatedWidth / widths:F3}, uniform {uniformWidth / widths:F3}; "
                         + $"over the grids m̂ averaged {marginSum / widths:F2} carried, {uniformSum / widths:F2} uniform, of 20; fallback in {fallBacks} of {widths}, "
                         + $"untrusted candidates in the inflated set in {untrustedSets}");
        Assert.True(inflatedRate >= 0.9 - 0.01 - 3.5 * se, $"coverage_inflated: {inflatedRate:F4}");
        Assert.True(uniformRate >= 0.9 - 0.01 - 3.5 * se, $"coverage_inflated, uniform: {uniformRate:F4}");
        Assert.True(plainRate >= reading - 3.5 * se, $"coverage_close: {plainRate:F4} against the margin's {reading:F4}");
        Assert.True(inflatedWidth >= plainWidth && uniformWidth >= inflatedWidth);
        Assert.True(reading >= uniformReading);
    }
}
