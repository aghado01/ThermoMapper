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
/// Batch predictive inference (E16, the conformal layer's slice 16; Lee, Tchetgen Tchetgen and Dobriban, arXiv 2409.13990): the rank law
/// and the box counts are enumerations of H (<c>card_testRank_mem</c>, <c>card_rankSeqs</c>); at one test unit the quantile set is split
/// conformal's; on exchangeable batches the quantile, simultaneous and function sets cover (<c>batch_interval</c>, <c>batch_simultaneous</c>,
/// <c>batch_coverage</c>), narrower than Bonferroni over single sets; and selection keeps its false claims within η (<c>batch_exceed</c>).
/// </summary>
public sealed class BatchConformalTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    // Every sorted rank sequence of H.
    private static List<int[]> H(int n, int m)
    {
        var all = new List<int[]>();
        var r = new int[m];
        void Fill(int j, int start)
        {
            if (j == m) { all.Add((int[])r.Clone()); return; }
            for (int v = start; v <= n + 1; v++) { r[j] = v; Fill(j + 1, v); }
        }
        Fill(0, 1);
        return all;
    }

    [Fact]
    public void TheCounts_AreEnumerationsOfH()
    {
        // For n ≤ 6 and m ≤ 4: the law of r_(ζ) against the share of H with r_(ζ) = k, and 300 random boxes against the share of H inside.
        var rng = new Xoshiro256PlusPlus(1801);
        for (int n = 0; n <= 6; n++)
            for (int m = 1; m <= 4; m++)
            {
                List<int[]> h = H(n, m);
                for (int zeta = 1; zeta <= m; zeta++)
                {
                    double[] law = BatchConformal.RankLaw(n, m, zeta);
                    for (int k = 1; k <= n + 1; k++)
                        Assert.Equal(h.Count(r => r[zeta - 1] == k) / (double)h.Count, law[k - 1], 12);
                }
            }
        for (int c = 0; c < 300; c++)
        {
            int n = rng.NextInt(7), m = 1 + rng.NextInt(4);
            int[] orders = [.. Enumerable.Range(1, m).Where(_ => rng.NextDouble() < 0.6)];
            if (orders.Length == 0) orders = [1 + rng.NextInt(m)];
            int[] w = [.. orders.Select(_ => 1 + rng.NextInt(n + 1))];
            int[] q = [.. w.Select(v => v + rng.NextInt(n + 2 - v))];
            double share = H(n, m).Count(r => orders.Select((t, j) => r[t - 1] >= w[j] && r[t - 1] <= q[j]).All(x => x)) / (double)H(n, m).Count;
            Assert.Equal(share, BatchConformal.BoxCoverage(n, m, orders, w, q), 12);
        }
    }

    [Fact]
    public void AtOneTestUnit_TheQuantileSet_IsSplitConformals()
    {
        // m = 1: q_U = ⌈(1 − γ)(n + 1)⌉, exactly, at levels where (1 − γ)(n + 1) is whole; q_L = ⌊β(n + 1)⌋ + 1 with β > 0.
        var rng = new Xoshiro256PlusPlus(1802);
        foreach ((int n, double gamma) in new[] { (19, 0.1), (99, 0.1), (37, 0.05), (9, 0.2), (3, 0.5) })
        {
            double[] s = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))];
            double[] sorted = [.. s.Order()];
            BatchInterval set = BatchConformal.Quantile(s, 1, 1, 0.0, gamma);
            int k = (int)Math.Ceiling((1 - gamma) * (n + 1) - 1e-9);
            Assert.Equal(k, set.UpperRank);
            Assert.Equal(k > n ? double.PositiveInfinity : sorted[k - 1], set.Upper);
            Assert.Equal(double.NegativeInfinity, set.Lower);
            Assert.Equal((double)k / (n + 1), set.Coverage, 12);
        }
    }

    [Fact]
    public void TheQuantileSet_Covers_NarrowerThanBonferroni()
    {
        // n 499 calibration and m 20 test scores, exchangeable N(0, 1); the median of the batch (ζ 10), β = γ = 0.05; 4,000 batches. Beside
        // it Bonferroni over the single sets: every test score in [S_(⌊(β/m)(n + 1)⌋), S_(⌈(1 − γ/m)(n + 1)⌉)].
        const int n = 499, m = 20, zeta = 10, replications = 4000;
        const double beta = 0.05, gamma = 0.05;
        var rng = new Xoshiro256PlusPlus(2026_1126);
        int covered = 0, bonferroni = 0;
        double width = 0.0, bonferroniWidth = 0.0, top = BatchConformal.RankLaw(n, m, zeta).Max();
        for (int r = 0; r < replications; r++)
        {
            double[] s = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))];
            double[] test = [.. Enumerable.Range(0, m).Select(_ => Gaussian(rng)).Order()];
            BatchInterval set = BatchConformal.Quantile(s, m, zeta, beta, gamma);
            if (set.Contains(test[zeta - 1])) covered++;
            width += set.Width;
            double[] sorted = [.. s.Order()];
            int lo = (int)Math.Floor(beta / m * (n + 1)), hi = (int)Math.Ceiling((1 - gamma / m) * (n + 1));
            double low = lo == 0 ? double.NegativeInfinity : sorted[lo - 1], high = hi > n ? double.PositiveInfinity : sorted[hi - 1];
            if (test[zeta - 1] >= low && test[zeta - 1] <= high) bonferroni++;
            bonferroniWidth += high - low;
        }
        double rate = covered / (double)replications, se = Math.Sqrt(0.09 / replications);
        output.WriteLine($"the batch median's set covers {rate:F4} (nominal 0.9, at most 0.9 + {top:F4} with distinct scores), mean width {width / replications:F3}; "
                         + $"Bonferroni covers {bonferroni / (double)replications:F4}, mean width {bonferroniWidth / replications:F3}; se {se:F4}");
        Assert.True(rate >= 0.9 - 3.5 * se);
        Assert.True(rate <= 0.9 + top + 3.5 * se);
        Assert.True(width < 0.5 * bonferroniWidth);
    }

    [Fact]
    public void TheSimultaneousBounds_Cover()
    {
        // The batch's quartiles at once: m 40 (orders 10 and 30), n 200, α 0.1, 3,000 exchangeable batches.
        const int n = 200, m = 40, replications = 3000;
        var rng = new Xoshiro256PlusPlus(2026_1127);
        int covered = 0;
        // The centred search's ranks depend on n, m, the orders and α alone: found once, then each batch read through its box.
        BatchBounds first = BatchConformal.Simultaneous([.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))], m, [10, 30], 0.1);
        double coverage = first.Coverage;
        for (int r = 0; r < replications; r++)
        {
            double[] s = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))];
            double[] test = [.. Enumerable.Range(0, m).Select(_ => Gaussian(rng))];
            BatchBounds bounds = BatchConformal.Box(s, m, first.Orders, first.LowerRanks, first.UpperRanks);
            if (bounds.Contains(test)) covered++;
        }
        double rate = covered / (double)replications, se = Math.Sqrt(0.09 / replications);
        output.WriteLine($"the quartiles at once: covered {rate:F4}, the box's share of H {coverage:F4}, ranks [{first.LowerRanks[0]}, {first.UpperRanks[0]}] "
                         + $"and [{first.LowerRanks[1]}, {first.UpperRanks[1]}] of {n + 1}; se {se:F4}");
        Assert.True(coverage >= 0.9);
        Assert.True(rate >= 0.9 - 3.5 * se);
    }

    [Fact]
    public void TheFunctionSet_ForTheMean_Covers_AndForAQuantile_IsTheQuantileSet()
    {
        // n 12, m 3 (|H| = 455): the batch mean of scores in [0, 1] (S_(0) = 0, S_(n+1) = 1; with unbounded scores the mean's set is the
        // whole line, a sequence holding S_(0) = −∞ always in reach), h̃ the mean of the ranks, β = γ = 0.05, 4,000 exchangeable batches.
        const int n = 12, m = 3, replications = 4000;
        var rng = new Xoshiro256PlusPlus(2026_1128);
        int covered = 0;
        double width = 0.0;
        for (int r = 0; r < replications; r++)
        {
            double[] s = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            double[] test = [.. Enumerable.Range(0, m).Select(_ => rng.NextDouble())];
            BatchInterval set = BatchConformal.Function(s, m, v => v.Average(), 0.05, 0.05, floor: 0.0, ceiling: 1.0);
            if (set.Contains(test.Average())) covered++;
            width += set.Width;
        }
        double rate = covered / (double)replications, se = Math.Sqrt(0.09 / replications);
        output.WriteLine($"the batch mean's set covers {rate:F4}, mean width {width / replications:F3} of the scores' range 1; se {se:F4}");
        Assert.True(rate >= 0.9 - 3.5 * se);
        Assert.True(width / replications < 1.0);

        // h the second smallest, h̃ = r_(2): Algorithm 1 is the quantile set.
        double[] scores = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))];
        BatchInterval byFunction = BatchConformal.Function(scores, m, v => v[1], 0.05, 0.05, r => r[1]);
        BatchInterval byQuantile = BatchConformal.Quantile(scores, m, 2, 0.05, 0.05);
        Assert.Equal(byQuantile.Lower, byFunction.Lower);
        Assert.Equal(byQuantile.Upper, byFunction.Upper);
        Assert.Equal(byQuantile.Coverage, byFunction.Coverage, 12);
    }

    [Fact]
    public void TheSelection_KeepsItsFalseClaimsWithinEta()
    {
        // x uniform, y = x + 0.3ε, μ̂(x) = x, c = 0.7: units declared above c, at most η = 2 false claims with chance 0.9; n 200, m 50,
        // 3,000 replications.
        const int n = 200, m = 50, eta = 2, replications = 3000;
        const double cutoff = 0.7, alpha = 0.1;
        var rng = new Xoshiro256PlusPlus(2026_1129);
        int kept = 0;
        double selected = 0.0;
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            double[] y = [.. x.Select(v => v + 0.3 * Gaussian(rng))];
            double[] xt = [.. Enumerable.Range(0, m).Select(_ => rng.NextDouble())];
            double[] yt = [.. xt.Select(v => v + 0.3 * Gaussian(rng))];
            BatchSelection pick = BatchConformal.Select(x, y, cutoff, xt, eta, alpha);
            int falseClaims = pick.Selected.Count(j => yt[j] <= cutoff);
            if (falseClaims <= eta) kept++;
            selected += pick.Selected.Count;
        }
        double rate = kept / (double)replications, se = Math.Sqrt(alpha * (1 - alpha) / replications);
        output.WriteLine($"false claims within {eta} in {rate:F4} of the batches, {selected / replications:F1} selected on average of {m}; se {se:F4}");
        Assert.True(rate >= 1 - alpha - 3.5 * se);
    }

    [Fact]
    public void WhatTheBatchDoesNotServe_IsRefused()
    {
        double[] s = [0.1, 0.5, 0.3];
        Assert.Throws<ArgumentOutOfRangeException>(() => BatchConformal.Quantile(s, 0, 1, 0.05, 0.05));
        Assert.Throws<ArgumentOutOfRangeException>(() => BatchConformal.Quantile(s, 3, 4, 0.05, 0.05));
        Assert.Throws<ArgumentOutOfRangeException>(() => BatchConformal.Quantile(s, 3, 2, 0.6, 0.5));
        Assert.Throws<ArgumentException>(() => BatchConformal.Quantile([0.1, double.NaN], 2, 1, 0.05, 0.05));
        Assert.Throws<ArgumentException>(() => BatchConformal.Quantile([-0.2, 0.3], 2, 1, 0.05, 0.05, floor: 0.0));
        Assert.Throws<ArgumentException>(() => BatchConformal.Simultaneous(s, 4, [2, 2], 0.1));
        Assert.Throws<ArgumentException>(() => BatchConformal.BoxCoverage(3, 4, [1, 3], [1, 3], [2, 2]));
        Assert.Throws<ArgumentException>(() => BatchConformal.Function([.. Enumerable.Range(0, 60).Select(i => (double)i)], 8, v => v.Sum(), 0.05, 0.05));
        Assert.Throws<ArgumentException>(() => BatchConformal.Select([0.4, -0.1], [0.2, 0.9], 0.5, [0.3], 0, 0.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => BatchConformal.Select([0.4, 0.1], [0.2, 0.9], 0.5, [0.3, 0.7], 2, 0.1));
        // A floor of 0 reads S_(0) = 0, where the lower end has no level.
        Assert.Equal(0.0, BatchConformal.Quantile([0.2, 0.4, 0.9], 2, 1, 0.0, 0.2, floor: 0.0).Lower);
    }
}
