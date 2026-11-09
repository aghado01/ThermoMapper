using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The smoothed conformal p-values of a stream (E35's slice 6a, M-1). The order-statistic tree ranks as a brute-force count does. M-1's
/// finite core holds by enumeration, ties included (<c>seqRank_uniform</c>): over every ordering of a data set, each rank sequence has
/// probability 1. The p-values are i.i.d. uniform on exchangeable streams, scores tied or not. Each premise the type enforces is broken in
/// a reference implementation beside it, with the failure each causes: a fixed calibration set, CB-1's tie rule, and scores kept from
/// earlier fits. A fixed fit's three scores are the split form's. With the qualification switch, the cost per arrival.
/// </summary>
public sealed class SmoothedPValueTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    // ---- The tree ----

    [Fact]
    public void TheTree_RanksAsABruteForceCount_TiesAndFarValuesIncluded()
    {
        double[] far = [double.NegativeInfinity, -1e300, -double.Epsilon, -0.0, 0.0, double.Epsilon, 1e300, double.PositiveInfinity];
        var rng = new Xoshiro256PlusPlus(seed: 20261108);
        int worstHeight = 0;
        for (int trial = 0; trial < 12; trial++)
        {
            var tree = new OrderStatisticTree<double>();
            var held = new List<double>();
            int n = trial < 4 ? 300 : 3000;
            for (int i = 0; i < n; i++)
            {
                double u = rng.NextDouble();
                double key = (trial % 4) switch
                {
                    0 => Math.Floor(5.0 * u),                       // five values, ties everywhere
                    1 => Gaussian(rng),                             // continuous
                    2 => far[(int)(far.Length * u)],                // the far values and both zeros
                    _ => u < 0.5 ? Math.Floor(3.0 * Gaussian(rng)) : 1e6 * Gaussian(rng),
                };
                if (trial >= 8) key = i;                            // ascending: the rotations' worst case
                (int Below, int Tied) before = tree.Add(key);
                if (i % 37 == 0 || i == n - 1) Assert.Equal((held.Count(v => v < key), held.Count(v => v == key)), before);
                held.Add(key);
                if (i % 37 == 0 || i == n - 1)
                    foreach (double query in new[] { key, Gaussian(rng), far[i % far.Length], held[(int)(held.Count * rng.NextDouble())] })
                    {
                        int below = held.Count(v => v < query), tied = held.Count(v => v == query);
                        Assert.Equal((below, tied), tree.Rank(query));
                    }
            }
            Assert.Equal(n, tree.Count);
            Assert.Equal(held.Distinct().Count(), tree.Distinct);   // 0 and −0 are one key, as == holds
            Assert.True(tree.Height <= 1.45 * Math.Log2(tree.Distinct + 2), $"height {tree.Height} for {tree.Distinct} keys");
            worstHeight = Math.Max(worstHeight, tree.Height);
        }
        output.WriteLine($"the tallest tree: {worstHeight}");
    }

    // ---- M-1's finite core, by enumeration ----

    private static IEnumerable<int[]> Permutations(int n)
    {
        int[] p = [.. Enumerable.Range(0, n)];
        while (true)
        {
            yield return (int[])p.Clone();
            int i = n - 2;
            while (i >= 0 && p[i] >= p[i + 1]) i--;
            if (i < 0) yield break;
            int j = n - 1;
            while (p[j] <= p[i]) j--;
            (p[i], p[j]) = (p[j], p[i]);
            Array.Reverse(p, i + 1, n - i - 1);
        }
    }

    public static TheoryData<double[]> DataSets =>
    [
        [1, 2, 3, 4], [1, 1, 2, 2], [3, 3, 3, 3], [1, 2, 2, 3, 3], [5, 5, 5, 7, 9], [0, -0.0, 2, 2, 2], [4, 1, 4, 1, 4],
    ];

    [Theory]
    [MemberData(nameof(DataSets))]
    public void M1_OverEveryOrdering_EachRankSequenceHasProbabilityOne(double[] data)
    {
        // At step n the tie-broken rank Rₙ = Lₙ + ⌈θₙTₙ⌉ is uniform on Lₙ + 1, …, Lₙ + Tₙ. Over the n! orderings, each sequence of
        // ranks rₙ ∈ 1, …, n has total probability Σ_orderings Π_n 1[rₙ ∈ (Lₙ, Lₙ + Tₙ]]/Tₙ = 1, so 1/n! under a uniform ordering.
        int n = data.Length;
        var total = new Dictionary<string, double>();
        foreach (int[] order in Permutations(n))
        {
            SmoothedPValue[] steps = new SmoothedPValues(seed: 1).AddRange(order.Select(i => data[i]));
            // Every rank sequence this ordering reaches, with its probability.
            var paths = new List<(string Key, double Probability)> { ("", 1.0) };
            foreach (SmoothedPValue step in steps)
            {
                var next = new List<(string, double)>();
                foreach ((string key, double probability) in paths)
                    for (int r = step.Below + 1; r <= step.Below + step.Tied; r++) next.Add((key + r + ",", probability / step.Tied));
                paths = next;
            }
            foreach ((string key, double probability) in paths) total[key] = total.GetValueOrDefault(key) + probability;
        }
        long sequences = 1;
        for (int k = 2; k <= n; k++) sequences *= k;
        Assert.Equal(sequences, total.Count);
        Assert.All(total.Values, v => Assert.Equal(1.0, v, 12));
    }

    // ---- i.i.d. uniform on exchangeable streams ----

    // The Kolmogorov–Smirnov distance of the values from U(0, 1).
    private static double KolmogorovSmirnov(IEnumerable<double> values)
    {
        double[] sorted = [.. values.Order()];
        double d = 0.0;
        for (int i = 0; i < sorted.Length; i++) d = Math.Max(d, Math.Max((i + 1.0) / sorted.Length - sorted[i], sorted[i] - (double)i / sorted.Length));
        return d;
    }

    private static double Correlation(IReadOnlyList<double> a, IReadOnlyList<double> b)
    {
        double ma = a.Average(), mb = b.Average(), sab = 0.0, saa = 0.0, sbb = 0.0;
        for (int i = 0; i < a.Count; i++)
        {
            sab += (a[i] - ma) * (b[i] - mb);
            saa += (a[i] - ma) * (a[i] - ma);
            sbb += (b[i] - mb) * (b[i] - mb);
        }
        return sab / Math.Sqrt(saa * sbb);
    }

    [Theory]
    [InlineData("continuous")]
    [InlineData("three values")]
    [InlineData("all tied")]
    public void ThePValues_AreIidUniform_OnExchangeableStreams(string scores)
    {
        const int streams = 400, length = 50;
        var rng = new Xoshiro256PlusPlus(seed: 20261108);
        var at = new List<double>[length];
        for (int n = 0; n < length; n++) at[n] = [];
        var first = new List<double>();
        var second = new List<double>();
        for (int r = 0; r < streams; r++)
        {
            var stream = new SmoothedPValues(seed: 1000 + r);
            for (int n = 0; n < length; n++)
            {
                double u = rng.NextDouble();
                double score = scores switch { "continuous" => Gaussian(rng), "three values" => u < 0.2 ? 0 : u < 0.7 ? 1 : 2, _ => 7.0 };
                double p = stream.Add(score).P;
                Assert.InRange(p, double.Epsilon, 1.0);
                at[n].Add(p);
                if (n >= 1) { first.Add(at[n - 1][^1]); second.Add(p); }
            }
        }
        double pooled = KolmogorovSmirnov(at.SelectMany(v => v)), worstPosition = at.Max(KolmogorovSmirnov);
        double lag = Correlation(first, second);
        output.WriteLine($"{scores}: KS pooled {pooled:F4} (1% critical {1.63 / Math.Sqrt(streams * length):F4}), worst position {worstPosition:F4} "
                         + $"(1% critical over {length} positions {1.95 / Math.Sqrt(streams):F4}); lag-one correlation {lag:F4}");
        Assert.True(pooled < 1.63 / Math.Sqrt(streams * length));
        Assert.True(worstPosition < 1.95 / Math.Sqrt(streams));   // 1% family-wise over 50 positions, Bonferroni
        Assert.True(Math.Abs(lag) < 4.0 / Math.Sqrt(first.Count));
    }

    // ---- Each premise broken, beside the type that enforces it ----

    [Fact]
    public void AFixedCalibrationSet_MakesThePValuesDependent_TheStreamsAreNot()
    {
        // Three calibration units, then two tests, of five distinct scores in every order. Against the fixed calibration set each test's
        // rank among four is uniform, yet P(both lowest) = 1/10, not 1/16. Ranked against the stream, the fourth and fifth ranks are
        // independent: P(both lowest) = 1/4 · 1/5.
        double[] data = [1, 2, 3, 4, 5];
        int splitBoth = 0, streamBoth = 0, orders = 0;
        foreach (int[] order in Permutations(5))
        {
            double[] v = [.. order.Select(i => data[i])];
            double[] calibration = v[..3];
            if (calibration.All(c => v[3] < c) && calibration.All(c => v[4] < c)) splitBoth++;
            SmoothedPValue[] steps = new SmoothedPValues(seed: 1).AddRange(v);
            if (steps[3].Below == 0 && steps[4].Below == 0) streamBoth++;
            orders++;
        }
        Assert.Equal(1.0 / 10.0, splitBoth / (double)orders, 12);
        Assert.Equal(1.0 / 20.0, streamBoth / (double)orders, 12);
    }

    [Fact]
    public void CB1sTieRule_IsNotUniform_TheRandomizedTieIs()
    {
        // An all-tied stream: CB-1's rule, ties counted for the unit, reads (L + T)/n = 1 at every step; the smoothed p is θ, uniform.
        SmoothedPValue[] steps = new SmoothedPValues(seed: 5).AddRange(Enumerable.Repeat(3.0, 4000));
        Assert.All(steps, s => Assert.Equal(1.0, (s.Below + s.Tied) / (double)s.Count));
        Assert.True(KolmogorovSmirnov(steps.Select(s => s.P)) < 1.63 / Math.Sqrt(steps.Length));
    }

    [Fact]
    public void ScoresKeptFromEarlierFits_Drift_AFixedScoreDoesNot()
    {
        // Each unit scored by a fit to the stream up to it, −|yᵢ − ȳᵢ| with ȳᵢ the mean of y₁, …, yᵢ, and the score kept: the early fits
        // sit nearer their own units, so later units rank low and p falls. One fixed score, −|yᵢ − μ|, keeps p uniform.
        const int streams = 500, length = 60;
        var rng = new Xoshiro256PlusPlus(seed: 9);
        var kept = new List<double>();
        var fixedScore = new List<double>();
        for (int r = 0; r < streams; r++)
        {
            var drifting = new SmoothedPValues(seed: r);
            var steady = new SmoothedPValues(seed: r + 100000);
            double sum = 0.0;
            for (int n = 1; n <= length; n++)
            {
                double y = Gaussian(rng);
                sum += y;
                double p = drifting.Add(-Math.Abs(y - sum / n)).P;
                double q = steady.Add(-Math.Abs(y)).P;
                kept.Add(p);
                fixedScore.Add(q);
            }
        }
        double critical = 1.63 / Math.Sqrt(kept.Count);
        output.WriteLine($"KS distance over every position: scores kept {KolmogorovSmirnov(kept):F4}, mean p {kept.Average():F4}; "
                         + $"a fixed score {KolmogorovSmirnov(fixedScore):F4}, mean p {fixedScore.Average():F4} (1% critical {critical:F4})");
        Assert.True(KolmogorovSmirnov(kept) > 2.0 * critical);
        Assert.True(KolmogorovSmirnov(fixedScore) < critical);
    }

    [Fact]
    public void ANaNScore_IsRefused()
        => Assert.Throws<ArgumentException>(() => new SmoothedPValues(seed: 1).Add(double.NaN));

    // ---- A fixed fit's scores ----

    // Exact draws from a conjugate normal posterior for a mean (σ² 1, prior N(0, 10²)) fitted to a training set proper, read at the
    // stream's units: the candidate seam's matrix and log CDF, by hand.
    private static CandidateLogDensity FixedFit(double[] training, double[] stream, int draws, Xoshiro256PlusPlus rng)
    {
        double precision = 1.0 / 100.0 + training.Length, mean = training.Sum() / precision;
        var values = new double[draws, stream.Length];
        var logCdf = new double[draws, stream.Length];
        for (int s = 0; s < draws; s++)
        {
            double mu = mean + Gaussian(rng) / Math.Sqrt(precision);
            for (int t = 0; t < stream.Length; t++)
            {
                double e = stream[t] - mu;
                values[s, t] = -0.5 * Math.Log(2.0 * Math.PI) - 0.5 * e * e;
                logCdf[s, t] = Maths.Distributions.Normal.StandardLogCdf(e);
            }
        }
        return new CandidateLogDensity(values, [.. stream.Select(y => new[] { y })], [new LatentReport("mean", LatentScope.Global, LatentTreatment.Conditioned)])
        {
            LogCdf = logCdf,
        };
    }

    [Fact]
    public void AFixedFitsScores_AreTheSplitFormsArithmetic()
    {
        var rng = new Xoshiro256PlusPlus(seed: 21);
        double[] training = [.. Enumerable.Range(0, 30).Select(_ => 0.4 + Gaussian(rng))];
        double[] stream = [.. Enumerable.Range(0, 40).Select(_ => 0.4 + 1.3 * Gaussian(rng)), -8.0, 9.0];
        CandidateLogDensity units = FixedFit(training, stream, 200, rng);
        double[] density = StreamScores.Density(units), centrality = StreamScores.Centrality(units);
        for (int t = 0; t < stream.Length; t++)
        {
            double top = Enumerable.Range(0, 200).Max(s => units.Values[s, t]);
            double expected = top + Math.Log(Enumerable.Range(0, 200).Sum(s => Math.Exp(units.Values[s, t] - top)) / 200.0);
            Assert.Equal(expected, density[t], 12);
            double f = Enumerable.Range(0, 200).Average(s => Math.Exp(units.LogCdf![s, t]));
            Assert.Equal(Math.Min(f, 1.0 - f), centrality[t], 12);
        }
        // Far out, both tails keep their digits: at 8.6 sds above, 1 − F is the tail's own, not 1 − F rounded to 0.
        Assert.True(centrality[^1] > 0.0 && centrality[^2] > 0.0);
        double[] residual = StreamScores.Residuals(new ResidualReadings([.. stream.Select(_ => 0.4)], [.. stream.Select(_ => 2.0)]), stream);
        for (int t = 0; t < stream.Length; t++) Assert.Equal(-Math.Abs(stream[t] - 0.4) / 2.0, residual[t]);
        Assert.Throws<ArgumentException>(() => StreamScores.Density(new CandidateLogDensity(new double[1, 2], [[1.0, 2.0]], [])));
    }

    [Theory]
    [InlineData("density")]
    [InlineData("centrality")]
    public void AFixedFitsScores_GiveUniformPValues_OnAnExchangeableStream_WhateverTheFit(string score)
    {
        // The fit assumes N(μ, 1); the stream is t₃, scaled and shifted, which the fit does not share. The fit is fixed, so the scores are
        // i.i.d. and the p-values uniform all the same.
        const int streams = 200, length = 40;
        var rng = new Xoshiro256PlusPlus(seed: 33);
        var all = new List<double>();
        for (int r = 0; r < streams; r++)
        {
            double[] training = [.. Enumerable.Range(0, 20).Select(_ => 1.0 + Gaussian(rng))];
            double[] stream = [.. Enumerable.Range(0, length).Select(_ => 1.0 + 0.7 * Gaussian(rng) / Math.Sqrt((Math.Pow(Gaussian(rng), 2) + Math.Pow(Gaussian(rng), 2) + Math.Pow(Gaussian(rng), 2)) / 3.0))];
            CandidateLogDensity units = FixedFit(training, stream, 100, rng);
            double[] scores = score == "density" ? StreamScores.Density(units) : StreamScores.Centrality(units);
            all.AddRange(new SmoothedPValues(seed: r).AddRange(scores).Select(s => s.P));
        }
        double d = KolmogorovSmirnov(all), critical = 1.63 / Math.Sqrt(all.Count);
        output.WriteLine($"{score}: KS {d:F4} (1% critical {critical:F4})");
        Assert.True(d < critical);
    }

    // ---- The cost per arrival ----

    [Fact]
    public void CostPerArrival()
    {
        // One arrival's rank and insert, at stream lengths 10³ to 10⁶, continuous scores and three-valued ones, after a warm-up stream
        // that leaves the JIT's optimized code in place. Cycle clock, indicative.
        if (!BassQualification.Enabled) return;
        var warm = new SmoothedPValues(seed: 2);
        var warmRng = new Xoshiro256PlusPlus(seed: 3);
        for (int i = 0; i < 400_000; i++) warm.Add(Gaussian(warmRng));
        foreach (string kind in new[] { "continuous", "three values" })
        {
            var rng = new Xoshiro256PlusPlus(seed: 4);
            double[] scores = [.. Enumerable.Range(0, 1_000_000).Select(_ => kind == "continuous" ? Gaussian(rng) : Math.Floor(3.0 * rng.NextDouble()))];
            var stream = new SmoothedPValues(seed: 1);
            var line = new List<string>();
            int done = 0;
            foreach (int upTo in new[] { 1_000, 10_000, 100_000, 1_000_000 })
            {
                int from = done;
                long t0 = BassCost.Tick();
                for (; done < upTo; done++) stream.Add(scores[done]);
                double perArrival = (BassCost.Tick() - t0) / BassCost.TicksPerSecond * 1e9 / (upTo - from);
                line.Add($"{perArrival:F0} ns per arrival from n {from:N0} to {upTo:N0}");
            }
            output.WriteLine($"{kind} ({(BassCost.Cycles ? "cycles" : "wall")} clock, indicative): {string.Join("; ", line)}");
        }
    }
}
