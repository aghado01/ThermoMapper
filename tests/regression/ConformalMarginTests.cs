using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// CB-3 in the conformal layer (E16, slice 2) against <c>lean/Enthymemata/ConformalValidity.lean</c>, each test naming its
/// declarations: the estimated band against the literal Lean count; the deterministic statements on perturbed scores, with the
/// adversary that attains them; the layer's errors, inflated ranks, uniform bounds and bands against a reference read from its
/// own scores; and the inflated set's coverage on exact draws of a model the data do not share, with the margin's reading beside
/// it.
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
            ConformalSet set = ConformalBayes.FullStep(training, candidates, new ConformalOptions { Alpha = 0.3 }).Sets[0];
            if (set.WholeSpace) continue;
            ConformalMargin margin = set.Margin!;
            double z = margin.Multiple;
            int[] rows = [.. Enumerable.Range(0, s).Where(d => weights is null || weights[d] > 0)];
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
                // The bounds, ranks and band from the layer's own scores and errors.
                double[] all = [.. scores.Scores, scores.Candidate];
                double[] bounds = [.. scores.Errors.Select(x => z * x), z * scores.CandidateError];
                Assert.Equal(InflatedRank(all, bounds, n), margin.InflatedRanks[c]);
                Assert.Equal(bounds.Max(), margin.UniformErrors[c]);
                Assert.Equal(z * scores.CandidateError, margin.CandidateErrors[c]);
                Assert.Equal(ConformalBayes.EstimatedBand(all.ToArray(), bounds.Max(), set.Level), margin.Fragile[c]);
                double[] uniform = [.. Enumerable.Repeat(bounds.Max(), n + 1)];
                Assert.Equal(Count(n + 1, u => FragileEstimated(set.Level, all, uniform, u)), margin.Fragile[c]);
                Assert.Equal(set.Accepts(c), set.Below![c] + set.Tied![c] + 1 > set.Level);
                checkedCandidates++;
                worstError = Math.Max(worstError, Math.Abs(candidateError - scores.CandidateError));
                if (far) worstPath++;
            }
            if (trial % 2 == 1)
            {
                // The carried bounds: each unit's the larger of its bound here and its bound as the candidate, read as a candidate
                // whose column is the unit's own.
                ConformalSet carried = ConformalBayes.FullStep(training, candidates, new ConformalOptions { Alpha = 0.3, MarginBounds = MarginBounds.Carried }).Sets[0];
                Assert.Equal(MarginBounds.Carried, carried.Margin!.Bounds);
                var self = new double[n];
                for (int i = 0; i < n; i++)
                {
                    var own = new double[s, 1];
                    for (int d = 0; d < s; d++) own[d, 0] = ll[d, i];
                    self[i] = ConformalBayes.Scores(training, new CandidateLogDensity(own, [[0.0]], theta, chainLengths: chains) { Weights = weights }, 0).CandidateError;
                }
                for (int c = 0; c < g; c++)
                {
                    ConformalScores scores = ConformalBayes.Scores(training, candidates, c);
                    double[] all = [.. scores.Scores, scores.Candidate];
                    double[] bounds = [.. scores.Errors.Select((x, i) => Math.Max(z * x, z * self[i])), z * scores.CandidateError];
                    Assert.Equal(Count(n + 1, u => FragileEstimated(set.Level, all, bounds, u)), carried.Margin.Fragile[c]);
                    Assert.True(carried.Margin.Fragile[c] <= set.Margin.Fragile[c] || bounds.Max() > set.Margin.UniformErrors[c]);
                }
            }
            Assert.Equal(set.Margin!.Count, set.Margin.Fragile.Max());
            Assert.Equal(1.0 - (double)(set.Level + set.Margin.Count) / set.Units - 0.01, set.CoverageUnderError, 12);
        }
        output.WriteLine($"{checkedCandidates} candidates checked, {worstPath} in far ranges; the largest error gap {worstError:E2}");
        Assert.True(checkedCandidates > 200);
    }

    private static void AssertClose(double expected, double actual, string what)
        => Assert.True(Math.Abs(expected - actual) <= 1e-7 + 1e-6 * Math.Abs(expected), $"{what}: {actual:R} against {expected:R}");

    private static double LogSumExp(IEnumerable<double> values)
    {
        double[] v = [.. values];
        double max = v.Max();
        return double.IsNegativeInfinity(max) ? max : max + Math.Log(v.Sum(x => Math.Exp(x - max)));
    }

    [Fact]
    public void FallBack_FiresOnAWideMargin_AndTheOptionsAreChecked()
    {
        // A degenerate stratum: every unit's likelihood swings across a few draws, so the bounds are wide and the band takes most
        // units; the plain set's guarantee under error falls far below 1 − α and the signal fires. Exact weights give no error.
        LatentReport[] theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];
        var rng = new Xoshiro256PlusPlus(5);
        const int s = 6, n = 29;
        var ll = new double[s, n];
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++) ll[d, i] = -1 - 8 * rng.NextDouble();
        var lc = new double[s, 3];
        for (int d = 0; d < s; d++)
            for (int c = 0; c < 3; c++) lc[d, c] = -1 - 8 * rng.NextDouble();
        var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), theta);
        var candidates = new CandidateLogDensity(lc, [[0.0, 1.0, 2.0]], theta);
        ConformalSet wide = ConformalBayes.FullStep(training, candidates).Sets[0];
        output.WriteLine($"m̂ {wide.Margin!.Count} of {wide.Units}, coverage under error {wide.CoverageUnderError:F3}, fallback {wide.FallBack}");
        Assert.True(wide.FallBack);
        Assert.True(wide.CoverageUnderError < 0.9 - 0.1);
        // One draw: the estimate is exact, every error 0, the band empty unless scores sit at the threshold.
        var exact = new PointwiseLogDensity(new double[,] { { -1, -2, -3, -4, -5, -6, -7, -8, -9, -10 } }, UnitPartition.EachObservation(10), theta);
        ConformalSet sure = ConformalBayes.FullStep(exact, new CandidateLogDensity(new double[,] { { -2.5 } }, [[0.0]], theta)).Sets[0];
        Assert.Equal(0.0, sure.Margin!.CandidateErrors[0]);
        Assert.Equal(0.0, sure.Margin.UniformErrors[0]);
        Assert.Equal(0, sure.Margin.Count);
        Assert.Throws<ArgumentOutOfRangeException>(() => ConformalBayes.FullStep(training, candidates, new ConformalOptions { ErrorLevel = 0.0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => ConformalBayes.FullStep(training, candidates, new ConformalOptions { FallbackShare = -1 }));
    }

    // ---- inflated coverage on exact draws ---------------------------------------------------------------------------------------

    [Fact]
    public void InflatedCoverage_OnExactDraws_UnderNoiseTheModelDoesNotShare()
    {
        // coverage_inflated and coverage_close: exact posterior draws, few of them (S 40), so the scores carry real Monte Carlo
        // error. n 19, α 0.1, k 2: the inflated set covers at least 1 − k/m − ν = 0.89 within its binomial error, and the plain set
        // at least the margin's reading, 1 − (k + m̂)/m − ν, averaged over the replications.
        const int n = 19, s = 40, replications = 4000, widths = 150;
        int[] seeds = SeedTree.Derive(2026_1108, replications);
        LatentReport[] theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];
        int plain = 0, inflated = 0, fallBacks = 0;
        double underError = 0.0, plainWidth = 0.0, inflatedWidth = 0.0, marginSum = 0.0, untrustedSets = 0;
        double[] grid = [.. Enumerable.Range(0, 601).Select(g => -10.0 + g / 20.0)];
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
            ConformalSet set = ConformalBayes.FullStep(training, candidates).Sets[0];
            ConformalSet wide = set with { Inflated = true };
            if (r < widths)
            {
                plainWidth += set.Width;
                inflatedWidth += wide.Width;
                marginSum += set.Margin!.Count;
                if (set.FallBack) fallBacks++;
                if (set.Margin.Untrusted(set.Level) > 0) untrustedSets++;
                continue;
            }
            if (set.Accepts(0)) plain++;
            if (wide.Accepts(0)) inflated++;
            underError += set.CoverageUnderError;
        }
        int count = replications - widths;
        double se = Math.Sqrt(0.09 / count), plainRate = plain / (double)count, inflatedRate = inflated / (double)count, reading = underError / count;
        output.WriteLine($"plain coverage {plainRate:F4}, inflated {inflatedRate:F4} (1 − k/m 0.9000, ν 0.01, se {se:F4}); the margin's reading at the true value averaged {reading:F4}");
        output.WriteLine($"widths: plain {plainWidth / widths:F3}, inflated {inflatedWidth / widths:F3}; over the grids m̂ averaged {marginSum / widths:F2} of 20, fallback in {fallBacks} of {widths}, untrusted candidates in the inflated set in {untrustedSets}");
        Assert.True(inflatedRate >= 0.9 - 0.01 - 3.5 * se, $"coverage_inflated: {inflatedRate:F4}");
        Assert.True(plainRate >= reading - 3.5 * se, $"coverage_close: {plainRate:F4} against the margin's {reading:F4}");
        Assert.True(inflatedWidth >= plainWidth);
    }
}
