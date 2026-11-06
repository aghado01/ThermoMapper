using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Distributions.Tests;

/// <summary>
/// Ratio-of-uniforms (E1, after slice 3): the primitive on laws with known rectangles, and Hörmann and Leydold's GIG draw forced, on chosen
/// cases in each of its three regions and a seeded sweep over far ranges, against the GIG's law by quadrature in ln x
/// (<see cref="LogConcaveEnvelopeTests.Law{T}"/>), Kolmogorov–Smirnov at 2,000 draws a case. The default draw meets
/// <c>gbass_primitives_oracle.R</c>'s law in the oracle tests.
/// </summary>
public sealed class GigRatioOfUniformsTests(ITestOutputHelper output)
{
    private const int Draws = 2_000;
    private const double KolmogorovCut = 2.6;

    private readonly struct StandardNormal : IRelativeLogDensity
    {
        public double Value(double x) => -0.5 * x * x;
    }

    private readonly struct UnitGamma(double shape) : IRelativeLogDensity
    {
        private readonly double _mode = shape - 1.0;
        public double Value(double x) => x > 0.0 ? (shape - 1.0) * Math.Log(x / _mode) - (x - _mode) : double.NegativeInfinity;
    }

    [Fact]
    public void ThePrimitive_FindsItsRectangle_AndDrawsTheLaw()
    {
        // N(0, 1): x·e^(−x²/4) is greatest at x = √2, so u± = ±√(2/e).
        (double minus, double plus) = RatioOfUniforms.Bounds(new StandardNormal(), 0.0, double.NegativeInfinity, double.PositiveInfinity);
        Assert.Equal(Math.Sqrt(2.0 / Math.E), plus, 7);
        Assert.Equal(-Math.Sqrt(2.0 / Math.E), minus, 7);
        var rng = new Xoshiro256PlusPlus(20261115);
        var normal = Normal.MeanSd(0.0, 1.0);
        var xs = new double[20_000];
        for (int i = 0; i < xs.Length; i++) xs[i] = RatioOfUniforms.Draw(rng, new StandardNormal(), 0.0, minus, plus);
        Assert.True(Kolmogorov(xs, normal.Cdf) <= KolmogorovCut);
        // Gamma(3.5): shifted by its mode, its rectangle by search.
        var gamma = new UnitGamma(3.5);
        (double gm, double gp) = RatioOfUniforms.Bounds(gamma, 2.5, 0.0, double.PositiveInfinity);
        var law = Gamma.ShapeRate(3.5, 1.0);
        for (int i = 0; i < xs.Length; i++) xs[i] = RatioOfUniforms.Draw(rng, gamma, 2.5, gm, gp);
        Assert.True(Kolmogorov(xs, law.Cdf) <= KolmogorovCut);
    }

    public static IEnumerable<object[]> Chosen() =>
    [
        // Without the shift: T₋½-concave, λ ≥ 1 − 9ω²/4, or ω > 0.2.
        [0.5, 2.0, 0.5], [0.5, 2.0, 0.03], [1.0, 3.0, 1.2], [-1.0, 4.0, 1.5], [0.0, 1.0, 1.0], [1.5, 0.5, 0.5],
        // With the shift: λ > 2 or ω > 3.
        [2.5, 1.0, 1.0], [-5.0, 2.0, 30.0], [0.5, 40.0, 40.0], [40.0, 1e-4, 1e4], [-3.0, 1e3, 1e-6], [200.0, 3.0, 0.1],
        // The dominated hat: λ < 1 and ω ≤ 0.2 below the T-concave line.
        [0.5, 2.0, 1e-4], [0.0107, 2e-6, 1.5e-7], [-0.3, 1e-3, 1e-3], [0.9, 0.01, 0.01], [0.0, 1e-6, 1e-6], [-0.99, 0.02, 0.5],
    ];

    [Theory]
    [MemberData(nameof(Chosen))]
    public void TheRatioDraw_CarriesTheLaw_OnChosenCases(double p, double a, double b)
    {
        double statistic = Statistic(p, a, b, new Xoshiro256PlusPlus(20261116));
        output.WriteLine($"GIG({p}, {a}, {b}): ω {Math.Sqrt(a * b):G3}, √n·D {statistic:F2}");
        Assert.True(statistic <= KolmogorovCut, $"√n·D = {statistic:F2}");
    }

    [Fact]
    public void TheRatioDraw_CarriesTheLaw_OverFarRanges()
    {
        var rng = new Xoshiro256PlusPlus(20261117);
        var failures = new List<string>();
        double worst = 0.0;
        for (int c = 0; c < 300; c++)
        {
            double p = (rng.NextDouble() < 0.5 ? -1.0 : 1.0) * (rng.NextDouble() < 0.1 ? 0.0 : Math.Pow(10.0, -3.0 + 6.0 * rng.NextDouble()));
            double a = Math.Pow(10.0, -8.0 + 16.0 * rng.NextDouble()), b = Math.Pow(10.0, -8.0 + 16.0 * rng.NextDouble());
            double statistic;
            try { statistic = Statistic(p, a, b, rng); }
            catch (ArithmeticException e) { failures.Add($"GIG({p:R}, {a:R}, {b:R}): {e.Message}"); continue; }
            worst = Math.Max(worst, statistic);
            if (!(statistic <= KolmogorovCut)) failures.Add($"GIG({p:R}, {a:R}, {b:R}), ω {Math.Sqrt(a * b):G3}: √n·D = {statistic:F2}");
        }
        output.WriteLine($"300 cases, worst √n·D {worst:F2}");
        Assert.True(failures.Count == 0, string.Join("\n", failures.Take(20)));
    }

    // √n·D of the ratio draw's ln x against the GIG's law in u = ln x about its peak, the hat's own log ratio integrated by quadrature.
    private static double Statistic(double p, double a, double b, Xoshiro256PlusPlus rng)
    {
        var gig = GeneralizedInverseGaussian.Of(p, a, b);
        double omega = Math.Sqrt(a) * Math.Sqrt(b), radius = Math.Sqrt(p * p + omega * omega), small = omega * omega / (2.0 * (radius + Math.Abs(p)));
        double up = p >= 0.0 ? 0.5 * (radius + p) : small, down = p >= 0.0 ? small : 0.5 * (radius - p);
        double peak = Math.Log(Math.Sqrt(b)) - Math.Log(Math.Sqrt(a)) + Math.Asinh(p / omega);
        var law = new LogConcaveEnvelopeTests.Law<GeneralizedInverseGaussian.LogRatio>(new GeneralizedInverseGaussian.LogRatio(up, down, 0.0),
                                                                                         double.PositiveInfinity, double.PositiveInfinity);
        var deltas = new double[Draws];
        for (int i = 0; i < Draws; i++) deltas[i] = Math.Log(gig.SampleByRatio(rng)) - peak;
        return Kolmogorov(deltas, law.Cdf);
    }

    private static double Kolmogorov(double[] xs, Func<double, double> cdf)
    {
        double[] sorted = [.. xs.Order()];
        double gap = 0.0;
        for (int i = 0; i < sorted.Length; i++)
        {
            double f = cdf(sorted[i]);
            gap = Math.Max(gap, Math.Max(Math.Abs(f - (double)i / sorted.Length), Math.Abs(f - (i + 1.0) / sorted.Length)));
        }
        return Math.Sqrt(sorted.Length) * gap;
    }
}
