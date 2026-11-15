using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Distributions.Tests;

/// <summary>
/// The interval Poincaré constant by finite elements (E8, slice 8), against what shares nothing with it: the uniform law's discrete
/// spectrum in closed form and its O(h²) approach to π²/(b − a)²; the Gaussian constant σ² on a normal truncated far out; a bimodal
/// mixture above the uniform's constant; and a seeded sweep of truncated mixtures of normal and Student t components against the
/// two-sided comparison with the uniform law, Rayleigh quotients of random smooth functions by quadrature, the linear function's lower
/// bound, the closed bound min(σ², (b − a)²/π²) where the law is one normal, and the grids' nesting.
/// </summary>
public sealed class IntervalPoincareTests(ITestOutputHelper output)
{
    /// <summary>
    /// On a uniform law the consistent-mass P1 pencil has λ₁ʰ = (6/h²)·2sin²(θ/2)/(2 + cos θ), θ = π/n, from n = 2 to 4,096 elements;
    /// its error to π²/L² falls by four per halving; and the constant read lies under L²/π², by at most 1.5 of its stated error.
    /// </summary>
    [Fact]
    public void Uniform_DiscreteSpectrum_AndTheRate()
    {
        const double a = -0.3, b = 1.7, length = b - a;
        double Flat(double _) => -Math.Log(length);
        double worst = 0.0;
        foreach (int n in new[] { 2, 3, 8, 64, 1000, 4096 })
        {
            double h = length / n, theta = Math.PI / n, s = Math.Sin(0.5 * theta);
            double expected = 6.0 / (h * h) * 2.0 * s * s / (2.0 + Math.Cos(theta));
            worst = Math.Max(worst, Math.Abs(IntervalPoincare.FirstGap(Flat, a, b, n) - expected) / expected);
        }
        double exact = Math.PI * Math.PI / (length * length);
        double[] errors = [.. new[] { 64, 128, 256 }.Select(n => IntervalPoincare.FirstGap(Flat, a, b, n) - exact)];
        double rate1 = errors[0] / errors[1], rate2 = errors[1] / errors[2];
        IntervalPoincare of = IntervalPoincare.Of(Flat, a, b);
        double under = 1.0 - of.Constant * exact;
        output.WriteLine($"worst gap to the discrete spectrum {worst:E2}; error ratios per halving {rate1:F4}, {rate2:F4}; " +
                         $"constant under L²/π² by {under:E2} at {of.Elements} elements, stated {of.RelativeError:E2}");
        Assert.True(worst <= 1e-10, $"{worst:E2}");
        Assert.InRange(rate1, 3.99, 4.01);
        Assert.InRange(rate2, 3.99, 4.01);
        Assert.True(under >= -1e-12, $"{under:E2}");
        Assert.True(under <= 1.5 * of.RelativeError, $"{under:E2} against {of.RelativeError:E2}");
        Assert.True(of.RelativeError <= 1e-7, $"{of.RelativeError:E2}");
    }

    /// <summary>
    /// The Gaussian constant: N(0.5, 0.04²) truncated to [0, 1], 12.5σ each side, has constant σ² to the truncation's e^(−78); its
    /// first eigenfunction is linear, which the elements hold exactly, so the constant read is σ² within quadrature and bisection.
    /// </summary>
    [Fact]
    public void NormalTruncatedFarOut_IsSigmaSquared()
    {
        const double sigma = 0.04;
        var law = Truncated.Of(Normal.MeanSd(0.5, sigma), 0.0, 1.0);
        IntervalPoincare of = IntervalPoincare.Of(law.LogPdf, 0.0, 1.0);
        double gap = Math.Abs(of.Constant / (sigma * sigma) - 1.0);
        output.WriteLine($"constant/σ² − 1 = {of.Constant / (sigma * sigma) - 1.0:E2} at {of.Elements} elements, stated error {of.RelativeError:E2}");
        Assert.True(gap <= 1e-9, $"{gap:E2}");
    }

    /// <summary>
    /// Two normals N(0.25, 0.06²) and N(0.75, 0.06²), equally weighted on [0, 1], whose density between them falls to about 10⁻⁴ of its
    /// peaks: not log-concave, and its constant is above the uniform law's 1/π² on the same interval, where neither closed bound holds.
    /// </summary>
    [Fact]
    public void BimodalMixture_ExceedsTheUniformConstant()
    {
        var law = Truncated.Of(Mixture.Of([Normal.MeanSd(0.25, 0.06), Normal.MeanSd(0.75, 0.06)], [0.5, 0.5]), 0.0, 1.0);
        IntervalPoincare of = IntervalPoincare.Of(law.LogPdf, 0.0, 1.0);
        double uniform = 1.0 / (Math.PI * Math.PI);
        output.WriteLine($"constant {of.Constant:E4} against 1/π² = {uniform:E4} and σ² = {0.06 * 0.06:E4}; {of.Elements} elements, stated error {of.RelativeError:E2}");
        Assert.True(of.Constant > uniform, $"{of.Constant:E4}");
        Assert.True(of.RelativeError <= 1e-7, $"{of.RelativeError:E2}");
    }

    /// <summary>
    /// 120 seeded truncated mixtures of one to three components, each normal or Student t (ν from 1 to 30), scales from 0.03 to 1.5,
    /// on intervals inside [−0.5, 1.5]. For each, the constant C read:
    /// <list type="bullet">
    /// <item>lies within the uniform law's L²/π² times and over sup ρ/inf ρ, both ways (each Rayleigh quotient's numerator and denominator
    /// compare with the uniform's within those factors), the extremes scanned at 20,001 points;</item>
    /// <item>is at least Var(f)/E[f′²] for three random smooth f, by Gauss–Legendre quadrature on 64 panels, within its stated error;</item>
    /// <item>is at least Var_ρ(x), the linear function's quotient;</item>
    /// <item>is at most min(σ², L²/π²) where the law is one normal;</item>
    /// <item>and its grid's λ₁ʰ falls as the grid is halved, from 1,024 elements to 2,048.</item>
    /// </list>
    /// </summary>
    [Fact]
    public void SeededSweep_OfTruncatedMixtures()
    {
        var rng = new Xoshiro256PlusPlus(20261115);
        double worstComparison = double.NegativeInfinity, worstQuotient = double.NegativeInfinity, worstLinear = double.NegativeInfinity,
               worstClosed = double.NegativeInfinity, worstNesting = double.NegativeInfinity, slowest = 0.0;
        var elements = new SortedDictionary<int, int>();
        var seconds = new List<double>();
        GaussLegendre rule = GaussLegendre.Of(20);
        for (int i = 0; i < 120; i++)
        {
            int count = 1 + (int)(3 * rng.NextDouble());
            var parts = new IContinuousDistribution[count];
            var weights = new double[count];
            bool oneNormal = count == 1;
            double sigma = 0.0;
            for (int c = 0; c < count; c++)
            {
                double location = -0.5 + 2.0 * rng.NextDouble(), scale = 0.03 * Math.Pow(50.0, rng.NextDouble());
                if (rng.NextDouble() < 0.5)
                {
                    parts[c] = Normal.MeanSd(location, scale);
                    sigma = scale;
                }
                else
                {
                    parts[c] = StudentT.LocationScale(1.0 + 29.0 * rng.NextDouble(), location, scale);
                    oneNormal = false;
                }
                weights[c] = 0.1 + 0.9 * rng.NextDouble();
            }
            double a = -0.5 + 1.5 * rng.NextDouble(), b = Math.Min(1.5, a + 0.1 + 1.4 * rng.NextDouble());
            var law = Truncated.Of(Mixture.Of(parts, weights), a, b);
            Func<double, double> logDensity = law.LogPdf;
            long start = Stopwatch.GetTimestamp();
            IntervalPoincare of = IntervalPoincare.Of(logDensity, a, b);
            seconds.Add(Stopwatch.GetElapsedTime(start).TotalSeconds);
            slowest = Math.Max(slowest, seconds[^1]);
            elements[of.Elements] = elements.GetValueOrDefault(of.Elements) + 1;
            Assert.True(of.RelativeError <= 1e-7 || of.Elements == IntervalPoincare.MaxElements, $"law {i}: {of.RelativeError:E2} at {of.Elements}");
            double length = b - a, uniform = length * length / (Math.PI * Math.PI);

            double top = double.NegativeInfinity, bottom = double.PositiveInfinity;
            for (int r = 0; r <= 20000; r++)
            {
                double l = logDensity(a + length * r / 20000.0);
                top = Math.Max(top, l);
                bottom = Math.Min(bottom, l);
            }
            double spread = Math.Exp(top - bottom);
            if (double.IsFinite(spread))
            {
                double low = uniform / spread, high = uniform * spread;
                worstComparison = Math.Max(worstComparison, Math.Max(low / of.Constant - 1.0, of.Constant / high - 1.0));
            }

            double[] amp = [.. Enumerable.Range(0, 3).Select(_ => Normal.SampleStandard(rng))];
            double[] freq = [.. Enumerable.Range(0, 3).Select(_ => 20.0 * rng.NextDouble())];
            double[] phase = [.. Enumerable.Range(0, 3).Select(_ => 6.3 * rng.NextDouble())];
            (double variance, double energy, double linear) = Moments(rule, logDensity, top, a, b,
                x => amp.Select((c, j) => c * Math.Sin(freq[j] * x + phase[j])).Sum(),
                x => amp.Select((c, j) => c * freq[j] * Math.Cos(freq[j] * x + phase[j])).Sum());
            if (energy > 0.0) worstQuotient = Math.Max(worstQuotient, variance / energy / of.Constant - 1.0 - of.RelativeError);
            worstLinear = Math.Max(worstLinear, linear / of.Constant - 1.0);
            if (oneNormal) worstClosed = Math.Max(worstClosed, of.Constant / Math.Min(sigma * sigma, uniform) - 1.0);
            double coarse = IntervalPoincare.FirstGap(logDensity, a, b, 1024), fine = IntervalPoincare.FirstGap(logDensity, a, b, 2048);
            worstNesting = Math.Max(worstNesting, fine / coarse - 1.0);
        }
        output.WriteLine($"beyond the uniform comparison {worstComparison:E2}; quotient over C {worstQuotient:E2}; Var(x) over C {worstLinear:E2}; " +
                         $"over the closed bound {worstClosed:E2}; λ at 2,048 over 1,024 {worstNesting:E2}; median {seconds.Order().ElementAt(seconds.Count / 2):F3} s, slowest {slowest:F3} s; " +
                         $"grids {string.Join(", ", elements.Select(e => $"{e.Key}: {e.Value}"))}");
        Assert.True(worstComparison <= 1e-4, $"{worstComparison:E2}");
        Assert.True(worstQuotient <= 1e-9, $"{worstQuotient:E2}");
        Assert.True(worstLinear <= 1e-9, $"{worstLinear:E2}");
        Assert.True(worstClosed <= 1e-9, $"{worstClosed:E2}");
        Assert.True(worstNesting <= 1e-8, $"{worstNesting:E2}");
    }

    /// <summary>Refused: a reversed or infinite interval, a grid of one element, a density zero across the interval, and a NaN in it.</summary>
    [Fact]
    public void Refusals()
    {
        Assert.Throws<ArgumentException>(() => IntervalPoincare.Of(_ => 0.0, 1.0, 0.0));
        Assert.Throws<ArgumentException>(() => IntervalPoincare.Of(_ => 0.0, 0.0, double.PositiveInfinity));
        Assert.Throws<ArgumentOutOfRangeException>(() => IntervalPoincare.FirstGap(_ => 0.0, 0.0, 1.0, 1));
        Assert.Throws<ArgumentException>(() => IntervalPoincare.Of(_ => double.NegativeInfinity, 0.0, 1.0));
        Assert.Throws<ArgumentException>(() => IntervalPoincare.Of(x => x > 0.5 ? double.NaN : 0.0, 0.0, 1.0));
    }

    // Var(f), E[f′²] and Var(x) under the law, its density scaled by e^(−top), by Gauss–Legendre on 64 panels.
    private static (double Variance, double Energy, double Linear) Moments(GaussLegendre rule, Func<double, double> logDensity, double top, double a, double b,
                                                                           Func<double, double> f, Func<double, double> slope)
    {
        const int panels = 64;
        double step = (b - a) / panels, mass = 0.0, m1 = 0.0, m2 = 0.0, energy = 0.0, x1 = 0.0, x2 = 0.0, middle = 0.5 * (a + b);
        for (int q = 0; q < panels; q++)
        {
            double centre = a + (q + 0.5) * step;
            for (int r = 0; r < rule.Count; r++)
            {
                double x = centre + 0.5 * step * rule.Nodes[r], w = 0.5 * step * rule.Weights[r] * Math.Exp(logDensity(x) - top);
                double v = f(x), s = slope(x), d = x - middle;
                mass += w;
                m1 += w * v;
                m2 += w * v * v;
                energy += w * s * s;
                x1 += w * d;
                x2 += w * d * d;
            }
        }
        double mean = m1 / mass, xm = x1 / mass;
        return (m2 / mass - mean * mean, energy / mass, x2 / mass - xm * xm);
    }
}
