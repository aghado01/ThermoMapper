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
/// The conformal layer's split form (E16, slice 3) against its Lean contract (CB-2 in <c>lean/Enthymemata/ConformalValidity.lean</c>),
/// each test naming its declarations: the accept rule against the count, ties counted for the test unit; CB-2 by exact enumeration
/// over a law on fits and units, whatever the fit; the density score against a reference in logs and against the full step's own
/// predictive, bit for bit; CQR's set an interval, and the residual's, over a seeded sweep with crossing quantiles, ties and far
/// ranges; the refusals; and coverage on the conjugate model under noise it does not share, for the exact posterior, a deliberately
/// wrong fit and a chain stopped early, under each score.
/// </summary>
public sealed class SplitConformalTests(ITestOutputHelper output)
{
    private static readonly LatentReport[] Theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];

    [Fact]
    public void Accepts_IsTheRankRule_TiesCountedForTheTestUnit()
    {
        // Accepts k C s t, k < rankLe C s t (CB-1's rule, which CB-2 reads; card_rankLe_le its count): the threshold form, a score at
        // least the k-th smallest calibration score, against #{calibration scores ≤ the score} + 1 > k, over a seeded sweep with
        // integer scores that tie, infinite ones, strata of one, strata with none and levels from 0 up.
        var rng = new Xoshiro256PlusPlus(303);
        double[] probes = [double.NegativeInfinity, -3.5, -3, -2, 0, 0.5, 3, 4, double.PositiveInfinity];
        for (int trial = 0; trial < 2000; trial++)
        {
            int n = rng.NextInt(40), strata = 1 + rng.NextInt(4);
            double alpha = 0.02 + 0.6 * rng.NextDouble();
            var scores = new double[n];
            var stratumOf = new int[n];
            for (int i = 0; i < n; i++)
            {
                double u = rng.NextDouble();
                scores[i] = u < 0.05 ? double.NegativeInfinity : u < 0.08 ? double.PositiveInfinity : rng.NextInt(7) - 3;
                stratumOf[i] = rng.NextInt(strata);
            }
            SplitCalibration fit = SplitConformal.Calibrate(scores, stratumOf, alpha);
            for (int h = 0; h <= strata; h++)   // stratum `strata` holds no unit
            {
                double[] own = [.. Enumerable.Range(0, n).Where(i => stratumOf[i] == h).Select(i => scores[i])];
                int m = own.Length + 1, k = (int)Math.Floor(alpha * m);
                Assert.Equal(m, fit.Units(h));
                Assert.Equal(k, fit.Level(h));
                foreach (double s in probes)
                {
                    int atMost = own.Count(c => c <= s), below = own.Count(c => c < s);
                    Assert.Equal((below, atMost - below), fit.Rank(h, s));
                    Assert.Equal((atMost + 1) / (double)m, fit.PValue(h, s));
                    Assert.Equal(atMost + 1 > k, fit.Accepts(h, s));
                }
            }
        }
        // Every score tied: every test unit tied with them is accepted.
        Assert.True(SplitConformal.Calibrate([.. Enumerable.Repeat(1.0, 30)], alpha: 0.5).Accepts(0, 1.0));
    }

    [Fact]
    public void SplitCoverage_ByExactEnumeration_ForAnyFit()
    {
        // split_coverage: W the training split with the fit's own randomness, under any law; given w the calibration and test units
        // are exchangeable, here independent given w (exchangeable_given_training's hypothesis); then the test unit, the last of m,
        // is accepted with chance at least 1 − k/m whatever map the fit makes of w (splitScore_equivariant). Enumerated exactly:
        // three training outcomes w, units on Z = {0, …, 4}, m of 2 to 5 units, every sequence of their values, and fits drawn at
        // random, scores tying across z or not: any table is a fit, a wrong one included. Restricted to sequences of distinct values,
        // still exchangeable given w, an injective fit covers exactly 1 − k/m (coverage_orderings_eq).
        var rng = new Xoshiro256PlusPlus(2026_1107);
        const int z = 5, ws = 3;
        for (int trial = 0; trial < 200; trial++)
        {
            int m = 2 + rng.NextInt(4), k = 1 + rng.NextInt(m - 1);
            double alpha = (k + 0.1 + 0.8 * rng.NextDouble()) / m;
            bool injective = trial % 2 == 1;
            var pw = new double[ws];
            var pz = new double[ws, z];
            var score = new double[ws, z];
            for (int w = 0; w < ws; w++)
            {
                pw[w] = 0.1 + rng.NextDouble();
                double total = 0.0;
                for (int v = 0; v < z; v++) total += pz[w, v] = 0.05 + rng.NextDouble();
                for (int v = 0; v < z; v++) pz[w, v] /= total;
                int[] order = [.. Enumerable.Range(0, z).OrderBy(_ => rng.NextDouble())];
                for (int v = 0; v < z; v++) score[w, v] = injective ? order[v] : rng.NextInt(3);
            }
            double wTotal = pw.Sum();
            double covered = 0.0, distinctMass = 0.0, distinctCovered = 0.0;
            var seq = new int[m];
            int count = (int)Math.Pow(z, m);
            for (int w = 0; w < ws; w++)
                for (int code = 0; code < count; code++)
                {
                    double p = pw[w] / wTotal;
                    for (int j = 0, c = code; j < m; j++, c /= z) { seq[j] = c % z; p *= pz[w, seq[j]]; }
                    SplitCalibration fit = SplitConformal.Calibrate([.. seq[..^1].Select(v => score[w, v])], alpha: alpha);
                    Assert.Equal(k, fit.Level(0));
                    bool accepted = fit.Accepts(0, score[w, seq[^1]]);
                    if (accepted) covered += p;
                    if (seq.Distinct().Count() == m)
                    {
                        distinctMass += p;
                        if (accepted) distinctCovered += p;
                    }
                }
            double nominal = 1.0 - (double)k / m;
            Assert.True(covered >= nominal - 1e-12, $"m {m}, k {k}: {covered:R} against {nominal:R}");
            if (injective && distinctMass > 0)
                Assert.True(Math.Abs(distinctCovered / distinctMass - nominal) <= 1e-12, $"m {m}, k {k}: {distinctCovered / distinctMass:R} against {nominal:R}");
        }
    }

    [Fact]
    public void Density_IsTheTrainingSplitsPredictive_RankedAmongTheCalibrationUnits()
    {
        // The score is ln Σ_s W_s p(y | θ_s): against a reference in logs, and against the full step's own predictive bit for bit (the
        // same arithmetic, ConformalBayes' weighing), over a seeded sweep of equal and weighted draws, chains, terms far in a tail, −∞
        // entries and whole −∞ columns, strata, and a test stratum no calibration unit holds (its whole space). The ranks are the
        // count's among the stratum's calibration units (coverage_stratum), and a test value whose column is a calibration unit's ties
        // with it exactly, counted for the test value (card_rankLe_le).
        var rng = new Xoshiro256PlusPlus(4040);
        int ties = 0;
        for (int trial = 0; trial < 300; trial++)
        {
            int s = 1 + rng.NextInt(60), n = 1 + rng.NextInt(25), cases = 1 + rng.NextInt(4), strata = 1 + rng.NextInt(3);
            bool weighted = trial % 3 == 1, chained = trial % 3 == 2 && s >= 2;
            double shift = trial % 5 == 0 ? -800.0 : trial % 5 == 1 ? 600.0 : 0.0;
            double alpha = 0.05 + 0.4 * rng.NextDouble();
            double Term() => rng.NextDouble() < 0.02 ? double.NegativeInfinity : shift - 6 * rng.NextDouble();
            var lcal = new double[s, n];
            for (int d = 0; d < s; d++)
                for (int i = 0; i < n; i++) lcal[d, i] = Term();
            if (trial % 7 == 0)
                for (int d = 0; d < s; d++) lcal[d, 0] = double.NegativeInfinity;
            var grids = new double[cases][];
            for (int t = 0; t < cases; t++) grids[t] = [.. Enumerable.Range(0, 1 + rng.NextInt(30)).Select(g => g * 0.5)];
            int columns = grids.Sum(g => g.Length);
            var lc = new double[s, columns];
            for (int d = 0; d < s; d++)
                for (int c = 0; c < columns; c++) lc[d, c] = Term();
            int copied = rng.NextInt(n);
            for (int d = 0; d < s; d++) lc[d, 0] = lcal[d, copied];
            int[] calStrata = [.. Enumerable.Range(0, n).Select(_ => rng.NextInt(strata))];
            int[] testStrata = [.. Enumerable.Range(0, cases).Select(_ => rng.NextInt(strata + 1))];
            testStrata[0] = calStrata[copied];
            double[]? weights = weighted ? [.. Enumerable.Range(0, s).Select(_ => rng.NextDouble() < 0.1 ? 0.0 : 0.1 + rng.NextDouble())] : null;
            if (weights is not null && weights.All(w => w == 0.0)) weights[0] = 1.0;
            int[]? lengths = chained ? [s / 2, s - s / 2] : null;
            var calibration = new CandidateLogDensity(lcal, [.. Enumerable.Range(0, n).Select(i => new[] { (double)i })], Theta, calStrata, lengths) { Weights = weights };
            var candidates = new CandidateLogDensity(lc, grids, Theta, testStrata, lengths) { Weights = weights };
            ConformalPrediction p = SplitConformal.Density(calibration, candidates, alpha);
            Assert.Equal(ConformalForm.Split, p.Form);
            Assert.Equal(alpha, p.Alpha);

            (int[] rows, double[]? logWeights) = PredictiveChecks.Rows(weights, s);
            var calScores = new double[n];
            SplitConformal.LogPredictive(lcal, rows, logWeights, 0, n, calScores);
            for (int i = 0; i < n; i++) Close(Reference(lcal, i, weights), calScores[i]);
            // The full step's weighing of the same columns, over a training density of finite terms.
            var ll = new double[s, 2];
            for (int d = 0; d < s; d++) { ll[d, 0] = -1 - rng.NextDouble(); ll[d, 1] = -1 - rng.NextDouble(); }
            var scorer = new ConformalBayes.StratumScorer(ll, rows, logWeights, [0, 1], lengths, 1, 1.0);
            for (int t = 0; t < cases; t++)
            {
                ConformalSet set = p.Sets[t];
                int h = testStrata[t];
                double[] own = [.. Enumerable.Range(0, n).Where(i => calStrata[i] == h).Select(i => calScores[i])];
                int m = own.Length + 1, k = (int)Math.Floor(alpha * m);
                Assert.Equal((t, h, m, k), (set.Case, set.Stratum, set.Units, set.Level));
                Assert.Null(set.Margin);
                Assert.Null(set.EffectiveSampleSizes);
                Assert.True(double.IsNaN(set.SplitHalfMass));
                if (k == 0)
                {
                    Assert.True(set.WholeSpace);
                    Assert.Equal(1.0, set.CoverageUnderError);
                    continue;
                }
                Assert.Equal(set.Coverage, set.CoverageUnderError);
                Assert.Equal(0.0, set.Gap);
                int g = grids[t].Length, offset = candidates.Offset(t);
                var reading = new ConformalBayes.Reading(g);
                scorer.Weigh(lc, offset, g, reading);
                for (int j = 0; j < g; j++)
                {
                    double score = set.CandidateScores![j];
                    Close(Reference(lc, offset + j, weights), score);
                    Assert.Equal(reading.LogPredictive[j], score);
                    int below = own.Count(c => c < score), atMost = own.Count(c => c <= score);
                    Assert.Equal(below, set.Below![j]);
                    Assert.Equal(atMost - below, set.Tied![j]);
                    Assert.Equal((atMost + 1) / (double)m, set.PValues![j]);
                    Assert.Equal(atMost + 1 > k, set.Accepts(j));
                    Assert.Equal(set.AcceptsPlain(j), set.AcceptsInflated(j));
                }
                if (t == 0)
                {
                    Assert.Equal(calScores[copied], set.CandidateScores![0]);
                    Assert.True(set.Tied![0] >= 1);
                    ties++;
                }
            }
        }
        Assert.True(ties > 60, $"{ties} copied columns read");

        // ln Σ_s W_s exp(ℓ_s)/Σ W, summed in logs from the largest term.
        static double Reference(double[,] values, int column, double[]? weights)
        {
            int s = values.GetLength(0);
            double total = weights?.Sum() ?? s, top = double.NegativeInfinity;
            for (int d = 0; d < s; d++)
                if ((weights is null || weights[d] > 0) && values[d, column] > top) top = values[d, column];
            if (double.IsNegativeInfinity(top)) return top;
            double sum = 0.0;
            for (int d = 0; d < s; d++)
                if (weights is null || weights[d] > 0) sum += (weights?[d] ?? 1.0) * Math.Exp(values[d, column] - top);
            return top + Math.Log(sum / total);
        }

        static void Close(double expected, double actual)
        {
            if (double.IsNegativeInfinity(expected)) Assert.Equal(expected, actual);
            else Assert.True(Math.Abs(expected - actual) <= 1e-12 * Math.Max(1.0, Math.Abs(expected)), $"{actual:R} against {expected:R}");
        }
    }

    [Fact]
    public void Refuses_WhatIsNotOneFunctionOfTheUnit_AndBadReadings()
    {
        var lcal = new double[,] { { -1, -2 }, { -1.5, -1 } };
        var lc = new double[,] { { -1, -2, -3 }, { -2, -1, -3 } };
        double[][] calGrids = [[0.0], [1.0]];
        var calibration = new CandidateLogDensity(lcal, calGrids, Theta);
        var candidates = new CandidateLogDensity(lc, [[0.0, 1.0, 2.0]], Theta);
        Assert.Single(SplitConformal.Density(calibration, candidates).Sets);
        // A calibration unit read on a grid, not at its response.
        Assert.Throws<ArgumentException>(() => SplitConformal.Density(new CandidateLogDensity(lc, [[0.0, 1.0], [2.0]], Theta), candidates));
        // Other draws, chains, weights or latents: not one function of the unit.
        Assert.Throws<ArgumentException>(() => SplitConformal.Density(calibration, new CandidateLogDensity(new double[,] { { -1 }, { -2 }, { -3 } }, [[0.0]], Theta)));
        Assert.Throws<ArgumentException>(() => SplitConformal.Density(calibration, new CandidateLogDensity(lc, [[0.0, 1.0, 2.0]], Theta, chainLengths: [1, 1])));
        Assert.Throws<ArgumentException>(() => SplitConformal.Density(calibration, new CandidateLogDensity(lc, [[0.0, 1.0, 2.0]], Theta) { Weights = [1.0, 2.0] }));
        Assert.Throws<ArgumentException>(() => SplitConformal.Density(calibration,
            new CandidateLogDensity(lc, [[0.0, 1.0, 2.0]], [new LatentReport("v", LatentScope.Observation, LatentTreatment.Integrated)])));
        Assert.Throws<ArgumentOutOfRangeException>(() => SplitConformal.Density(calibration, candidates, 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => SplitConformal.Calibrate([1.0], alpha: 0.0));
        Assert.Throws<ArgumentException>(() => SplitConformal.Calibrate([double.NaN]));
        Assert.Throws<ArgumentException>(() => SplitConformal.Calibrate([1.0, 2.0], [0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => SplitConformal.Calibrate([1.0], [-1]));
        // The readings: one per case, finite, a positive scale, a response per calibration case.
        var q = new QuantileReadings([0.0, 1.0], [1.0, 2.0]);
        Assert.Equal(2, SplitConformal.Quantiles(q, [0.5, 1.5], q).Intervals.Count);
        Assert.Throws<ArgumentException>(() => SplitConformal.Quantiles(new QuantileReadings([0.0], [1.0, 2.0]), [0.5], q));
        Assert.Throws<ArgumentException>(() => SplitConformal.Quantiles(new QuantileReadings([0.0, double.NaN], [1.0, 2.0]), [0.5, 1.5], q));
        Assert.Throws<ArgumentException>(() => SplitConformal.Quantiles(q, [0.5], q));
        Assert.Throws<ArgumentException>(() => SplitConformal.Quantiles(q, [0.5, double.PositiveInfinity], q));
        Assert.Throws<ArgumentException>(() => SplitConformal.Quantiles(q with { StratumOf = [0] }, [0.5, 1.5], q));
        var r = new ResidualReadings([0.0, 1.0], [1.0, 2.0]);
        Assert.Equal(2, SplitConformal.Residuals(r, [0.5, 1.5], r).Intervals.Count);
        Assert.Throws<ArgumentException>(() => SplitConformal.Residuals(r with { Scale = [1.0, 0.0] }, [0.5, 1.5], r));
        Assert.Throws<ArgumentException>(() => SplitConformal.Residuals(r with { Scale = [1.0, double.PositiveInfinity] }, [0.5, 1.5], r));
        Assert.Throws<ArgumentException>(() => SplitConformal.Residuals(r with { Center = [0.0, double.NaN] }, [0.5, 1.5], r));
        Assert.Throws<ArgumentOutOfRangeException>(() => SplitConformal.Residuals(r with { StratumOf = [0, -1] }, [0.5, 1.5], r));
    }

    [Fact]
    public void CqrSet_IsAnInterval_AndTheResidualsToo_TheirEndsTheRankRules()
    {
        // cqr_ordConnected: with the calibration held, the values y the rank rule accepts at a test case form an interval. The layer
        // gives it in closed form, [q_lo + c, q_hi − c] with c the threshold, since the score min(y − q_lo, q_hi − y) (cqrScore)
        // reaches c exactly there; the standardized residual's, μ ± Qρ, by the same argument. Over a seeded sweep (crossing
        // quantiles, responses rounded so the calibration ties, scales from 10⁻¹⁵⁰ to 10¹⁵⁰, strata, calibration bands wide enough to
        // empty a narrow test band), on a grid through and past each interval, its ends and their neighbours: the rank rule read by
        // the count agrees with the closed form except within rounding of an end, and its accepted values run unbroken.
        var rng = new Xoshiro256PlusPlus(5150);
        int empties = 0, read = 0;
        for (int trial = 0; trial < 1500; trial++)
        {
            double scale = trial % 4 == 0 ? Math.Pow(10, -150 + 300 * rng.NextDouble()) : 1.0;
            bool wide = trial % 7 == 3, rounded = trial % 2 == 0;
            int n = rng.NextInt(60), cases = 1 + rng.NextInt(5), strata = 1 + rng.NextInt(2);
            double alpha = 0.05 + 0.5 * rng.NextDouble();
            (double[] lo, double[] hi, double[] mid, double[] spread) Band(int count, bool calibrating)
            {
                var l = new double[count];
                var u = new double[count];
                var c = new double[count];
                var w = new double[count];
                for (int i = 0; i < count; i++)
                {
                    c[i] = scale * (2 * rng.NextDouble() - 1);
                    double half = calibrating && wide ? 5.0 : rng.NextDouble() - 0.25;
                    l[i] = c[i] - scale * half;
                    u[i] = c[i] + scale * half;
                    w[i] = scale * (Math.Abs(half) + 0.1);
                }
                return (l, u, c, w);
            }
            (double[] lo, double[] hi, double[] mid, double[] spread) = Band(n, true);
            var y = new double[n];
            for (int i = 0; i < n; i++)
            {
                double e = 1.2 * Normal.SampleStandard(rng);
                y[i] = mid[i] + scale * (rounded ? Math.Round(4 * e) / 4 : e);
            }
            (double[] tlo, double[] thi, double[] tmid, double[] tspread) = Band(cases, false);
            int[] calStrata = [.. Enumerable.Range(0, n).Select(_ => rng.NextInt(strata))];
            int[] testStrata = [.. Enumerable.Range(0, cases).Select(_ => rng.NextInt(strata))];
            ConformalIntervals cqr = SplitConformal.Quantiles(new QuantileReadings(lo, hi, calStrata), y, new QuantileReadings(tlo, thi, testStrata), alpha);
            ConformalIntervals res = SplitConformal.Residuals(new ResidualReadings(mid, spread, calStrata), y, new ResidualReadings(tmid, tspread, testStrata), alpha);
            Assert.Equal(alpha, cqr.Alpha);
            for (int t = 0; t < cases; t++)
            {
                int h = testStrata[t];
                int[] members = [.. Enumerable.Range(0, n).Where(i => calStrata[i] == h)];
                int m = members.Length + 1, k = (int)Math.Floor(alpha * m);
                double[] cqrScores = [.. members.Select(i => SplitConformal.CqrScore(lo[i], hi[i], y[i]))];
                double[] resScores = [.. members.Select(i => SplitConformal.ResidualScore(mid[i], spread[i], y[i]))];
                Check(cqr.Intervals[t], cqrScores, v => SplitConformal.CqrScore(tlo[t], thi[t], v), Math.Abs(tlo[t]) + Math.Abs(thi[t]));
                Check(res.Intervals[t], resScores, v => SplitConformal.ResidualScore(tmid[t], tspread[t], v), Math.Abs(tmid[t]) + tspread[t]);

                void Check(ConformalInterval interval, double[] own, Func<double, double> score, double size)
                {
                    Assert.Equal((t, h, m, k), (interval.Case, interval.Stratum, interval.Units, interval.Level));
                    if (k == 0)
                    {
                        Assert.True(interval.WholeSpace);
                        Assert.True(interval.Contains(double.MaxValue) && interval.Contains(-double.MaxValue));
                        return;
                    }
                    double threshold = own.Order().ElementAt(k - 1);
                    Assert.Equal(threshold, interval.Threshold);
                    if (interval.Empty) empties++;
                    double centre = 0.5 * (interval.Lower + interval.Upper);
                    double reach = Math.Max(Math.Abs(interval.Upper - interval.Lower), size) * 1.5 + scale;
                    var probes = new List<double>();
                    for (int g = 0; g <= 400; g++) probes.Add(centre - reach + 2 * reach * g / 400);
                    foreach (double end in (double[])[interval.Lower, interval.Upper])
                        probes.AddRange([Math.BitDecrement(end), end, Math.BitIncrement(end)]);
                    probes.Sort();
                    double tolerance = 1e-12 * (size + reach);
                    int first = -1, last = -1, accepted = 0;
                    for (int g = 0; g < probes.Count; g++)
                    {
                        double v = probes[g], sv = score(v);
                        bool rank = own.Count(c => c <= sv) + 1 > k;
                        if (rank)
                        {
                            if (first < 0) first = g;
                            last = g;
                            accepted++;
                        }
                        if (rank != interval.Contains(v))
                            Assert.True(Math.Min(Math.Abs(v - interval.Lower), Math.Abs(v - interval.Upper)) <= tolerance,
                                        $"y {v:R}: the rule {rank}, the interval [{interval.Lower:R}, {interval.Upper:R}]");
                    }
                    if (accepted > 0) Assert.Equal(last - first + 1, accepted);
                    read++;
                }
            }
        }
        output.WriteLine($"{read} sets read, {empties} empty");
        Assert.True(empties > 5);
    }

    // ---- coverage on the conjugate model -----------------------------------------------------------------------------------------

    // A fit on the training split: draws of the slope β of y = βx + N(0, sd²), with the noise sd the fit reads and its chains.
    private sealed record Fit(string Name, double[] Beta, double NoiseSd, int[]? Chains);

    // A chain stopped early: random-walk Metropolis on the conjugate model's posterior from β = 5, steps of sd 0.02, forty iterations
    // all kept, nowhere near the posterior.
    internal static double[] StoppedEarly(double[] x, double[] y, Xoshiro256PlusPlus rng)
    {
        double LogPost(double b)
        {
            double lp = -0.05 * b * b;
            for (int i = 0; i < x.Length; i++) lp -= 0.5 * (y[i] - b * x[i]) * (y[i] - b * x[i]);
            return lp;
        }
        var chain = new double[40];
        double beta = 5.0, current = LogPost(beta);
        for (int it = 0; it < chain.Length; it++)
        {
            double proposal = beta + 0.02 * Normal.SampleStandard(rng), next = LogPost(proposal);
            if (Math.Log(rng.NextDouble()) < next - current) { beta = proposal; current = next; }
            chain[it] = beta;
        }
        return chain;
    }

    private static double LogNormal(double y, double mean, double sd) => -0.5 * Math.Log(2 * Math.PI) - Math.Log(sd) - 0.5 * (y - mean) * (y - mean) / (sd * sd);

    // The coverage's standard error over r replications of t test cases sharing each calibration: given the calibration the
    // coverage is Beta(m − k, k), so its variance (m − k)k/(m²(m + 1)) over r, plus the cases' binomial remainder over r·t.
    internal static double CoverageError(int m, int k, int r, int t)
    {
        double mean = 1.0 - (double)k / m, between = (double)(m - k) * k / ((double)m * m * (m + 1));
        return Math.Sqrt(between / r + (mean * (1 - mean) - between) / ((double)r * t));
    }

    [Fact]
    public void Coverage_ForAnyFit_TheExactPosterior_AWrongFit_AndAChainStoppedEarly()
    {
        // split_coverage on a continuous space: 30 training units, 19 calibration units, so m 20 and k 2 at α 0.1, exactly 0.9 with
        // the distinct scores a continuous space gives, under noise the model does not share (ConformalBayesTests.Truth: skewed,
        // heteroscedastic, curved). Three fits on the training split: the conjugate posterior's exact draws (a wrong model, exact
        // draws); a deliberately wrong fit, the slope's sign flipped with an overconfident noise sd 0.25, CQR's quantiles crossed and
        // the residual's scale shaped against the noise; and a chain stopped early, forty Metropolis steps from β = 5. Each under the
        // density score (the seam), CQR's (the fit's normal predictive quantiles at 0.05 and 0.95) and the standardized residual's.
        // Coverage at the true y over 2,000 replications of 20 test cases, within four standard errors of the Beta law's; the width
        // beside it (the density's read on a grid, so a set past its ends is read to them, and counted), and each fit's own 90%
        // interval, uncalibrated.
        const int train = 30, calibrationUnits = 19, tests = 20, replications = 2000, widthReps = 300;
        const double z95 = 1.6448536269514722;
        double[] grid = [.. Enumerable.Range(0, 801).Select(g => -50.0 + 0.125 * g)];
        string[] scores = ["density", "CQR", "residual"];
        var covered = new int[3, 3];
        var width = new double[3, 3];
        var own = new int[3];
        int open = 0;
        int[] seeds = SeedTree.Derive(2026_1107, replications);
        for (int r = 0; r < replications; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            (double[] x, double[] y) Units(int count)
            {
                var xs = new double[count];
                var ys = new double[count];
                for (int i = 0; i < count; i++) { xs[i] = 4 * rng.NextDouble() - 2; ys[i] = ConformalBayesTests.Truth(xs[i], 0.5, rng); }
                return (xs, ys);
            }
            (double[] xTrain, double[] yTrain) = Units(train);
            (double[] xCal, double[] yCal) = Units(calibrationUnits);
            (double[] xTest, double[] yTest) = Units(tests);
            double[] exact = ConformalBayesTests.ExactDraws(xTrain, yTrain, 200, rng);
            Fit[] fits =
            [
                new("the exact posterior", exact, 1.0, null),
                new("a wrong fit", [.. exact.Select(b => -b)], 0.25, null),
                new("a chain stopped early", StoppedEarly(xTrain, yTrain, rng), 1.0, [40]),
            ];
            for (int f = 0; f < fits.Length; f++)
            {
                Fit fit = fits[f];
                bool wrong = f == 1;
                int s = fit.Beta.Length;
                double mean = fit.Beta.Average(), variance = fit.Beta.Sum(b => (b - mean) * (b - mean)) / s;
                double PredictiveSd(double x) => Math.Sqrt(fit.NoiseSd * fit.NoiseSd + x * x * variance);

                // The density score through the seam: the calibration units at their responses, the test cases at theirs.
                CandidateLogDensity Seam(double[] xs, double[][] grids)
                {
                    var values = new double[s, grids.Sum(g => g.Length)];
                    for (int d = 0; d < s; d++)
                        for (int t = 0, c = 0; t < xs.Length; t++)
                            foreach (double v in grids[t]) values[d, c++] = LogNormal(v, fit.Beta[d] * xs[t], fit.NoiseSd);
                    return new CandidateLogDensity(values, grids, Theta, chainLengths: fit.Chains);
                }
                CandidateLogDensity calibration = Seam(xCal, [.. yCal.Select(v => new[] { v })]);
                ConformalPrediction density = SplitConformal.Density(calibration, Seam(xTest, [.. yTest.Select(v => new[] { v })]));
                if (r < widthReps)
                    foreach (ConformalSet set in SplitConformal.Density(calibration, Seam(xTest[..2], [grid, grid])).Sets)
                    {
                        width[f, 0] += set.Width;
                        if (set.OpenBelow || set.OpenAbove) open++;
                    }

                // CQR: the fit's predictive quantiles, crossed for the wrong fit.
                QuantileReadings Quantiles(double[] xs) => wrong
                    ? new([.. xs.Select(x => mean * x + 0.3)], [.. xs.Select(x => mean * x - 0.3)])
                    : new([.. xs.Select(x => mean * x - z95 * PredictiveSd(x))], [.. xs.Select(x => mean * x + z95 * PredictiveSd(x))]);
                ConformalIntervals cqr = SplitConformal.Quantiles(Quantiles(xCal), yCal, Quantiles(xTest));

                // The standardized residual: the fit's centre and predictive sd, for the wrong fit a scale shaped against the noise.
                ResidualReadings Residuals(double[] xs) => new([.. xs.Select(x => mean * x)], [.. xs.Select(x => wrong ? 1.0 / (0.2 + Math.Abs(x)) : PredictiveSd(x))]);
                ConformalIntervals residual = SplitConformal.Residuals(Residuals(xCal), yCal, Residuals(xTest));

                for (int t = 0; t < tests; t++)
                {
                    Assert.Equal((20, 2), (density.Sets[t].Units, density.Sets[t].Level));
                    if (density.Sets[t].AcceptsPlain(0)) covered[f, 0]++;
                    if (cqr.Intervals[t].Contains(yTest[t])) covered[f, 1]++;
                    if (residual.Intervals[t].Contains(yTest[t])) covered[f, 2]++;
                    width[f, 1] += cqr.Intervals[t].Width;
                    width[f, 2] += residual.Intervals[t].Width;
                    if (Math.Abs(yTest[t] - mean * xTest[t]) <= z95 * PredictiveSd(xTest[t])) own[f]++;
                }
            }
        }
        int count = replications * tests;
        double se = CoverageError(20, 2, replications, tests);
        string[] names = ["the exact posterior", "a wrong fit", "a chain stopped early"];
        output.WriteLine($"density sets open at the grid's ends, their widths read to the grid: {open} of {3 * 2 * widthReps}");
        for (int f = 0; f < 3; f++)
        {
            output.WriteLine($"{names[f]}: its own 90% interval covers {own[f] / (double)count:F4}");
            for (int c = 0; c < 3; c++)
            {
                double coverage = covered[f, c] / (double)count;
                double meanWidth = c == 0 ? width[f, 0] / (2 * widthReps) : width[f, c] / count;
                output.WriteLine($"  {scores[c]}: coverage {coverage:F4} (nominal 0.9000, se {se:F4}), mean width {meanWidth:F2}");
                Assert.True(Math.Abs(coverage - 0.9) <= 4 * se, $"{names[f]}, {scores[c]}: {coverage:F4}");
            }
        }
    }
}
