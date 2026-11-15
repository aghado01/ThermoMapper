using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// <see cref="SpecialFunctions.LogGammaRatio"/> against references that subtract no large numbers: at an integer step s = m, the sum
/// Σ ln(x + i); at s = ½ and an integer x = n, ½ ln π + ln n + Σₖ ln(1 − 1/(2k)), from Γ(n + ½)/Γ(n) = √π·n·Πₖ (2k − 1)/(2k); and against
/// the difference of two ln Γ where both are small enough to keep their digits. A seeded sweep from x = 10⁻³⁰⁰ to 10¹⁵ checks
/// additivity, ln Γ(x + s + r)/Γ(x) = ln Γ(x + s)/Γ(x) + ln Γ(x + s + r)/Γ(x + s), the property no single reference spans.
/// </summary>
public sealed class LogGammaRatioTests
{
    private const double Eps = 2.220446049250313e-16;

    [Theory]
    [InlineData(1e-6)]
    [InlineData(0.37)]
    [InlineData(1.0)]
    [InlineData(15.5)]
    [InlineData(16.0)]
    [InlineData(123.25)]
    [InlineData(1e5)]
    [InlineData(1e9)]
    public void IntegerSteps_MatchTheProduct(double x)
    {
        for (int m = 1; m <= 40; m++)
        {
            var (sum, scale) = (new Neumaier(), 0.0);
            for (int i = 0; i < m; i++)
            {
                double term = Math.Log(x + i);
                sum.Add(term);
                scale += Math.Abs(term);
            }
            double got = SpecialFunctions.LogGammaRatio(x, m), error = Math.Abs(got - sum.Value);
            Assert.True(error <= 8 * Eps * (1.0 + scale), $"x = {x:R}, s = {m}: {got:R} against {sum.Value:R}");
        }
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(15)]
    [InlineData(16)]
    [InlineData(1000)]
    [InlineData(100_000)]
    [InlineData(1_000_000)]
    public void HalfStep_MatchesTheCentralBinomial(int n)
    {
        var sum = new Neumaier();
        sum.Add(0.5 * Math.Log(Math.PI));
        sum.Add(Math.Log(n));
        for (int k = 1; k <= n; k++) sum.Add(SpecialFunctions.Log1p(-0.5 / k));
        double got = SpecialFunctions.LogGammaRatio(n, 0.5);
        Assert.True(Math.Abs(got - sum.Value) <= 1e-14 * (1.0 + Math.Abs(sum.Value)), $"n = {n}: {got:R} against {sum.Value:R}");
    }

    [Fact]
    public void ModerateArguments_MatchTheLogGammaDifference()
    {
        var rng = new Xoshiro256PlusPlus(20261114);
        for (int c = 0; c < 2000; c++)
        {
            double x = Math.Pow(10.0, -3.0 + 4.6 * rng.NextDouble()), s = 30.0 * rng.NextDouble();   // x below 40, x + s below 70
            double reference = SpecialFunctions.LogGamma(x + s) - SpecialFunctions.LogGamma(x);
            double scale = Math.Abs(SpecialFunctions.LogGamma(x + s)) + Math.Abs(SpecialFunctions.LogGamma(x));
            double got = SpecialFunctions.LogGammaRatio(x, s);
            Assert.True(Math.Abs(got - reference) <= 64 * Eps * (1.0 + scale), $"x = {x:R}, s = {s:R}: {got:R} against {reference:R}");
        }
    }

    [Fact]
    public void Additivity_SeededSweep()
    {
        var rng = new Xoshiro256PlusPlus(20261115);
        var failures = new List<string>();
        double worst = 0.0;
        for (int c = 0; c < 20_000; c++)
        {
            double x = Math.Pow(10.0, -300.0 + 315.0 * rng.NextDouble());
            double s = rng.NextDouble() < 0.1 ? rng.NextInt(8) * 0.5 : Math.Pow(10.0, -6.0 + 12.0 * rng.NextDouble());
            double r = rng.NextDouble() < 0.1 ? rng.NextInt(8) * 0.5 : Math.Pow(10.0, -6.0 + 12.0 * rng.NextDouble());
            if (x + s == x) continue;   // the second leg's x would equal the first's: nothing to add
            double whole = SpecialFunctions.LogGammaRatio(x, s + r);
            double first = SpecialFunctions.LogGammaRatio(x, s), second = SpecialFunctions.LogGammaRatio(x + s, r);
            Assert.True(double.IsFinite(whole) && double.IsFinite(first) && double.IsFinite(second), $"x = {x:R}, s = {s:R}, r = {r:R}");
            // Rounding s + r and x + s moves each leg by its slope, ψ(·) ≈ ln(·), times the rounding.
            double slope = Math.Abs(Math.Log(x + s + r)) + Math.Abs(Math.Log(x + s)) + 1.0;
            double tolerance = 64 * Eps * (Math.Abs(first) + Math.Abs(second) + 1.0) + 4 * Eps * (x + s + r) * slope;
            double error = Math.Abs(whole - (first + second));
            worst = Math.Max(worst, error / tolerance);
            if (!(error <= tolerance)) failures.Add($"x = {x:R}, s = {s:R}, r = {r:R}: {whole:R} against {first + second:R}");
        }
        Assert.True(failures.Count == 0, string.Join("\n", failures.GetRange(0, Math.Min(20, failures.Count))));
    }

    [Fact]
    public void Boundaries()
    {
        Assert.Equal(0.0, SpecialFunctions.LogGammaRatio(3.5, 0.0));
        Assert.Equal(Math.Log(4.0), SpecialFunctions.LogGammaRatio(4.0, 1.0), 14);
        // At a subnormal x, Γ(x + 1)/Γ(x) = x: ln x itself.
        Assert.Equal(Math.Log(1e-310), SpecialFunctions.LogGammaRatio(1e-310, 1.0), 12);
        Assert.Throws<ArgumentOutOfRangeException>(() => SpecialFunctions.LogGammaRatio(0.0, 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => SpecialFunctions.LogGammaRatio(-1.0, 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => SpecialFunctions.LogGammaRatio(double.NaN, 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => SpecialFunctions.LogGammaRatio(1.0, -0.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => SpecialFunctions.LogGammaRatio(1.0, double.PositiveInfinity));
    }

    private sealed class Neumaier
    {
        private double _sum, _carry;

        public double Value => _sum + _carry;

        public void Add(double v)
        {
            double t = _sum + v;
            _carry += Math.Abs(_sum) >= Math.Abs(v) ? _sum - t + v : v - t + _sum;
            _sum = t;
        }
    }
}
