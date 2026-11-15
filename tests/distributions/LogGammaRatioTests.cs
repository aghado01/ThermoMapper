using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Distributions.Tests;

/// <summary>
/// <see cref="SpecialFunctions.LogGammaRatio"/> against references that subtract no large numbers: at an integer step s = m, the sum
/// Σ ln(x + i); at s = ½ and an integer x = n, ½ ln π + ln n + Σₖ ln(1 − 1/(2k)), from Γ(n + ½)/Γ(n) = √π·n·Πₖ (2k − 1)/(2k); at a small
/// step, Taylor's series in s through the polygammas; and against the difference of two ln Γ where both are small enough to keep their
/// digits. A seeded sweep from x = 10⁻³⁰⁰ to 10¹⁵ checks additivity, ln Γ(x + s + r)/Γ(x) = ln Γ(x + s)/Γ(x) + ln Γ(x + s + r)/Γ(x + s),
/// the property no single reference spans.
/// </summary>
public sealed class LogGammaRatioTests(ITestOutputHelper output)
{
    private const double Eps = 2.220446049250313e-16;

    /// <summary>
    /// A small step keeps its relative digits: against s·ψ(x) + s²/2·ψ′(x) + s³/6·ψ″(x), whose next term is below 2·10⁻¹⁷ of the value
    /// at s ≤ 10⁻⁶ and x ≥ ¼ (1.5·10⁻¹⁷ at x = ¼). The polygammas come from outside the function: from their values at ¼, ½ or 1 by their
    /// recurrences, summed in decimal, or at 10⁶ and 10⁹ by their asymptotic series with ln x = m·ln 10. Until 2026-11-14 the function
    /// subtracted Stirling's remainder at x + s from the one at x, two rounded values near 1/(12x), for an absolute error up to about
    /// 2·10⁻¹⁸ however small s was: up to 4·10⁻⁶ of the value at s = 10⁻¹² in the slice 0 review's sweep, and 2.5·10⁻¹⁴ to 3.3·10⁻¹² here.
    /// </summary>
    [Theory]
    [InlineData(0.25)]
    [InlineData(0.5)]
    [InlineData(1.0)]
    [InlineData(2.0)]
    [InlineData(2.5)]
    [InlineData(10.0)]
    [InlineData(15.5)]
    [InlineData(16.0)]
    [InlineData(17.0)]
    [InlineData(99.5)]
    [InlineData(12345.0)]
    [InlineData(1e5)]
    [InlineData(1e6)]
    [InlineData(1e9)]
    public void SmallSteps_KeepTheirRelativeDigits(double x)
    {
        var (psi, trigamma, tetragamma) = Polygammas(x);
        double worst = 0.0;
        foreach (double s in new[] { 1e-6, 1e-8, 1e-10, 1e-12, 1e-15 })
        {
            double truth = s * (psi + 0.5 * s * (trigamma + s / 3.0 * tetragamma));
            double got = SpecialFunctions.LogGammaRatio(x, s), relative = Math.Abs(got - truth) / Math.Abs(truth);
            worst = Math.Max(worst, relative);
            Assert.True(relative <= 1e-14, $"x = {x:R}, s = {s:R}: {got:R} against {truth:R}, relative {relative:E2}");
        }
        output.WriteLine($"x = {x:R}: ψ = {psi:R}, worst relative error {worst:E2}");
    }

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

    // ψ, ψ′ and ψ″ at x = ¼, at an integer or half-integer x ≥ ½, or at x = 10ᵐ ≥ 10⁶, in decimal, then rounded once.
    private static (double Psi, double Trigamma, double Tetragamma) Polygammas(double x)
    {
        const decimal EulerGamma = 0.5772156649015328606065120901m, Ln2 = 0.6931471805599453094172321215m;
        const decimal Ln10 = 2.3025850929940456840179914547m, Pi = 3.1415926535897932384626433833m;
        const decimal Zeta3 = 1.2020569031595942853997381615m, Catalan = 0.9159655941772190150546035149m;
        decimal psi, trigamma, tetragamma;
        if (x >= 1e6)
        {
            // ψ ~ ln x − 1/(2x) − 1/(12x²) + 1/(120x⁴), ψ′ ~ 1/x + 1/(2x²) + 1/(6x³), ψ″ ~ −1/x² − 1/x³ − 1/(2x⁴): the next terms are below 10⁻³⁰.
            int m = (int)Math.Round(Math.Log10(x));
            Assert.Equal(Math.Pow(10.0, m), x);
            decimal u = 1m / (decimal)x, u2 = u * u;
            psi = m * Ln10 - u / 2m - u2 / 12m + u2 * u2 / 120m;
            trigamma = u + u2 / 2m + u2 * u / 6m;
            tetragamma = -u2 - u2 * u - u2 * u2 / 2m;
        }
        else
        {
            decimal t = (decimal)x - decimal.Floor((decimal)x);
            if (x == 0.25)
            {
                psi = -EulerGamma - Pi / 2m - 3m * Ln2;
                trigamma = Pi * Pi + 8m * Catalan;
                tetragamma = -2m * Pi * Pi * Pi - 56m * Zeta3;
            }
            else if (t == 0.5m)
            {
                psi = -EulerGamma - 2m * Ln2;
                trigamma = Pi * Pi / 2m;
                tetragamma = -14m * Zeta3;
            }
            else
            {
                Assert.Equal(0m, t);
                t = 1m;
                psi = -EulerGamma;
                trigamma = Pi * Pi / 6m;
                tetragamma = -2m * Zeta3;
            }
            // ψ(t + 1) = ψ(t) + 1/t, ψ′(t + 1) = ψ′(t) − 1/t², ψ″(t + 1) = ψ″(t) + 2/t³.
            for (; t < (decimal)x; t += 1m)
            {
                decimal inverse = 1m / t;
                psi += inverse;
                trigamma -= inverse * inverse;
                tetragamma += 2m * inverse * inverse * inverse;
            }
        }
        return ((double)psi, (double)trigamma, (double)tetragamma);
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
