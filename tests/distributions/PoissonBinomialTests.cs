using System;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// The Poisson-binomial pmf against brute force over the 2ᵖ inclusion vectors, at chosen coins and in a seeded sweep whose
/// probabilities span twelve decades and reach 0 and 1; against the binomial when the coins agree, at two thousand coins; the
/// truncated recursion against the full one; the interval probabilities against the pmf; and the refusals.
/// </summary>
public sealed class PoissonBinomialTests
{
    // Σ over χ ∈ {0, 1}ᵖ of Π pⱼ^χⱼ (1 − pⱼ)^(1 − χⱼ), by its count.
    private static double[] BruteForce(double[] p)
    {
        var pmf = new double[p.Length + 1];
        for (long mask = 0; mask < 1L << p.Length; mask++)
        {
            double product = 1.0;
            int count = 0;
            for (int j = 0; j < p.Length; j++)
            {
                bool on = ((mask >> j) & 1) != 0;
                product *= on ? p[j] : 1.0 - p[j];
                count += on ? 1 : 0;
            }
            pmf[count] += product;
        }
        return pmf;
    }

    private static double[] Pmf(double[] p, int length)
    {
        var pmf = new double[length];
        PoissonBinomial.Pmf(p, pmf);
        return pmf;
    }

    [Theory]
    [InlineData(new[] { 0.5 })]
    [InlineData(new[] { 0.1, 0.9 })]
    [InlineData(new[] { 0.0, 1.0, 0.3 })]
    [InlineData(new[] { 0.2, 0.2, 0.2, 0.2, 0.2 })]
    [InlineData(new[] { 1e-9, 0.5, 1.0 - 1e-9, 0.25, 0.75, 0.01, 0.99, 0.6 })]
    public void Pmf_MatchesBruteForce(double[] p)
    {
        double[] exact = BruteForce(p), pmf = Pmf(p, p.Length + 1);
        for (int k = 0; k <= p.Length; k++) Assert.Equal(exact[k], pmf[k], 1e-15);
        Assert.Equal(1.0, pmf.Sum(), 1e-14);
    }

    [Fact]
    public void SeededSweep_MatchesBruteForce_AcrossTwelveDecades()
    {
        var rng = new Xoshiro256PlusPlus(20261107);
        double worst = 0.0;
        for (int c = 0; c < 400; c++)
        {
            int n = 1 + (int)(rng.NextDouble() * 14);
            var p = new double[n];
            for (int j = 0; j < n; j++)
            {
                double u = rng.NextDouble();
                p[j] = (c % 4) switch
                {
                    0 => u,
                    1 => Math.Pow(10.0, -12.0 * u),             // small coins, log-uniform
                    2 => 1.0 - Math.Pow(10.0, -12.0 * u),       // near-sure coins
                    _ => u < 0.15 ? 0.0 : u > 0.85 ? 1.0 : u,   // certain coins among the rest
                };
            }
            double[] exact = BruteForce(p), pmf = Pmf(p, n + 1);
            for (int k = 0; k <= n; k++)
            {
                // Each entry is a sum of nonnegative products, so it holds its relative precision.
                double error = Math.Abs(pmf[k] - exact[k]) / Math.Max(exact[k], 1e-300);
                if (exact[k] > 1e-290) worst = Math.Max(worst, error);
            }
        }
        Assert.True(worst < 1e-12, $"worst relative error {worst:E2}");
    }

    [Fact]
    public void EqualCoins_AreTheBinomial_AtTwoThousandCoins()
    {
        const int n = 2000;
        foreach (double q in new[] { 1e-4, 0.03, 0.5, 0.97 })
        {
            double[] p = Enumerable.Repeat(q, n).ToArray();
            double[] pmf = Pmf(p, n + 1);
            Assert.Equal(1.0, pmf.Sum(), 1e-12);
            double worst = 0.0;
            for (int k = 0; k <= n; k++)
            {
                double log = SpecialFunctions.LogGamma(n + 1.0) - SpecialFunctions.LogGamma(k + 1.0) - SpecialFunctions.LogGamma(n - k + 1.0)
                             + k * Math.Log(q) + (n - k) * SpecialFunctions.Log1p(-q);
                if (log < -600.0) continue;
                worst = Math.Max(worst, Math.Abs(pmf[k] / Math.Exp(log) - 1.0));
            }
            Assert.True(worst < 1e-9, $"q {q}: worst relative error {worst:E2}");
        }
    }

    [Fact]
    public void TruncatedRecursion_KeepsTheLowCountsBitForBit()
    {
        var rng = new Xoshiro256PlusPlus(7);
        double[] p = [.. Enumerable.Range(0, 60).Select(_ => rng.NextDouble() * 0.2)];
        double[] full = Pmf(p, p.Length + 1);
        for (int top = 0; top <= 8; top++)
        {
            double[] truncated = Pmf(p, top + 1);
            for (int k = 0; k <= top; k++) Assert.Equal(full[k], truncated[k]);
        }
    }

    [Fact]
    public void IntervalProbabilities_SumThePmf()
    {
        double[] p = [0.3, 0.05, 0.6, 0.2, 0.9, 0.1, 0.45];
        double[] pmf = BruteForce(p);
        for (int lower = 0; lower <= p.Length + 1; lower++)
            for (int upper = lower - 1; upper <= p.Length + 2; upper++)
            {
                double expected = 0.0;
                for (int k = Math.Max(lower, 0); k <= Math.Min(upper, p.Length); k++) expected += pmf[k];
                Assert.Equal(expected, PoissonBinomial.Probability(p, lower, upper), 1e-15);
            }
        Assert.Equal(pmf[0] + pmf[1] + pmf[2], PoissonBinomial.Cdf(p, 2), 1e-15);
        Assert.Equal(0.0, PoissonBinomial.Cdf(p, -1));
        Assert.Equal(1.0, PoissonBinomial.Probability([], 0, 3));
    }

    [Fact]
    public void Refuses_ProbabilitiesOutsideTheUnitInterval()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => PoissonBinomial.Probability([0.2, 1.5], 0, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => PoissonBinomial.Probability([-0.1], 0, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => PoissonBinomial.Probability([double.NaN], 0, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => PoissonBinomial.Probability([0.5, 2.0], 5, 7));
        Assert.Throws<ArgumentOutOfRangeException>(() => PoissonBinomial.Probability([0.5], -1, 0));
    }
}
