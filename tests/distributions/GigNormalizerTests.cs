using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// The GIG's normalizer, density and moments against quadrature of its kernel in u = ln x, where p·u − (a·eᵘ + b·e⁻ᵘ)/2 is concave, over a
/// seeded sweep of p, a and b across far ranges and at both limits.
/// </summary>
public sealed class GigNormalizerTests
{
    [Fact]
    public void Sweep_MeetsTheKernelsIntegral()
    {
        var rng = new Xoshiro256PlusPlus(20261111);
        var failures = new List<string>();
        for (int c = 0; c < 1_500; c++)
        {
            double p = (rng.NextDouble() < 0.5 ? -1.0 : 1.0) * (rng.NextDouble() < 0.1 ? 0.5 : LogUniform(rng, 1e-3, 100.0));
            double a = LogUniform(rng, 1e-12, 1e6), b = LogUniform(rng, 1e-12, 1e6);
            switch (rng.NextDouble())
            {
                case < 0.15 when p > 0.0: b = 0.0; break;
                case < 0.3 when p < 0.0: a = 0.0; break;
            }
            if (rng.NextDouble() < 0.05) p = 0.0;
            if ((p >= 0.0 && a == 0.0) || (p <= 0.0 && b == 0.0)) continue;
            var gig = GeneralizedInverseGaussian.Of(p, a, b);
            string at = $"GIG({p:R}, {a:R}, {b:R})";
            double l0 = LogMoment(p, a, b, 0), l1 = LogMoment(p, a, b, 1), l2 = LogMoment(p, a, b, 2);
            Check(failures, $"{at} ln Z", GeneralizedInverseGaussian.LogNormalizer(p, a, b), l0, 1e-11, log: true);
            double x = gig.Mode > 0.0 ? gig.Mode : Math.Exp(l1 - l0);
            // ln f is a difference of the kernel and the normalizer, each as large as (ax + b/x)/2: its rounding is theirs.
            Check(failures, $"{at} ln f({x:R})", gig.LogPdf(x), gig.LogKernel(x) - l0, 1e-11 * Math.Max(1.0, Math.Abs(l0)), log: true);
            double mean = double.IsFinite(l1) ? Math.Exp(l1 - l0) : double.PositiveInfinity;
            Check(failures, $"{at} mean", gig.Mean, mean, 1e-10, log: false);
            double variance = double.IsFinite(l2) ? Math.Exp(l2 - l0) - mean * mean : double.PositiveInfinity;
            // The quadrature's own variance cancels where it is far below the mean squared; read it where it holds.
            if (double.IsPositiveInfinity(variance) || variance > 1e-4 * mean * mean)
                Check(failures, $"{at} variance", gig.Variance, variance, 1e-6, log: false);
        }
        Assert.True(failures.Count == 0, string.Join("\n", failures.Take(20)));
    }

    [Fact]
    public void Limits_AreTheGammaAndInverseGamma()
    {
        Assert.Equal(SpecialFunctions.LogGamma(2.5) + 2.5 * Math.Log(2.0 / 3.0), GeneralizedInverseGaussian.LogNormalizer(2.5, 3.0, 0.0), 12);
        Assert.Equal(SpecialFunctions.LogGamma(1.5) - 1.5 * Math.Log(4.0 / 2.0), GeneralizedInverseGaussian.LogNormalizer(-1.5, 0.0, 4.0), 12);
        Assert.Equal(double.PositiveInfinity, GeneralizedInverseGaussian.LogNormalizer(-0.5, 3.0, 0.0));
        Assert.Equal(double.PositiveInfinity, GeneralizedInverseGaussian.LogNormalizer(0.0, 0.0, 4.0));
        // Inverse Gaussian, GIG(−½, γ², 1): mean 1/γ, variance 1/γ³, normalizer √(2π)e^(−γ).
        var ig = GeneralizedInverseGaussian.Of(-0.5, 9.0, 1.0);
        Assert.Equal(1.0 / 3.0, ig.Mean, 14);
        Assert.Equal(1.0 / 27.0, ig.Variance, 13);
        Assert.Equal(0.5 * Math.Log(2.0 * Math.PI) - 3.0, GeneralizedInverseGaussian.LogNormalizer(-0.5, 9.0, 1.0), 13);
    }

    private static void Check(List<string> failures, string what, double got, double expected, double tolerance, bool log)
    {
        if (double.IsPositiveInfinity(expected))
        {
            if (!double.IsPositiveInfinity(got)) failures.Add($"{what}: {got:R} against +∞");
            return;
        }
        double gap = log ? Math.Abs(got - expected) / Math.Max(1.0, Math.Abs(expected)) : Math.Abs(got - expected) / Math.Abs(expected);
        if (!(gap <= tolerance)) failures.Add($"{what}: {got:R} against {expected:R}, gap {gap:G3}");
    }

    // ln ∫ x^k·x^(p−1)e^(−(ax + b/x)/2) dx = ln ∫ exp((p + k)u − (a·eᵘ + b·e⁻ᵘ)/2) du: the peak by bisection on the slope, the ends where it
    // has fallen by 60 (+∞ where it never falls on a side, an integral that diverges), Gauss–Legendre on even panels between.
    internal static double LogMoment(double p, double a, double b, int k)
    {
        double q = p + k;
        if ((a == 0.0 && q >= 0.0) || (b == 0.0 && q <= 0.0)) return double.PositiveInfinity;
        double G(double u) => q * u - 0.5 * (a > 0.0 ? a * Math.Exp(u) : 0.0) - 0.5 * (b > 0.0 ? b * Math.Exp(-u) : 0.0);
        double Slope(double u) => q - 0.5 * (a > 0.0 ? a * Math.Exp(u) : 0.0) + 0.5 * (b > 0.0 ? b * Math.Exp(-u) : 0.0);
        double lo = -1.0, hi = 1.0;
        while (Slope(lo) < 0.0) lo *= 2.0;
        while (Slope(hi) > 0.0) hi *= 2.0;
        for (int i = 0; i < 300; i++)
        {
            double mid = 0.5 * (lo + hi);
            if (Slope(mid) > 0.0) lo = mid;
            else hi = mid;
        }
        double peak = 0.5 * (lo + hi), top = G(peak), floor = top - 60.0;
        double Edge(double direction)
        {
            double step = 1e-3, at = peak;
            while (G(at + direction * step) >= floor) step *= 2.0;
            double inner = peak, outer = peak + direction * step;
            for (int i = 0; i < 300; i++)
            {
                double mid = 0.5 * (inner + outer);
                if (G(mid) < floor) outer = mid;
                else inner = mid;
            }
            return outer;
        }
        double left = Edge(-1.0), right = Edge(1.0);
        int panels = (int)Math.Clamp((right - left) / 0.02, 2_000.0, 5_000_000.0);
        double[] nodes = [0.0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640];
        double[] weights = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891];
        double width = (right - left) / panels, sum = 0.0;
        for (int j = 0; j < panels; j++)
        {
            double mid = left + (j + 0.5) * width;
            for (int i = 0; i < nodes.Length; i++) sum += weights[i] * Math.Exp(G(mid + 0.5 * width * nodes[i]) - top);
        }
        return top + Math.Log(sum * 0.5 * width);
    }

    private static double LogUniform(Xoshiro256PlusPlus rng, double lo, double hi) => DistributionSweepTests.LogUniform(rng, lo, hi);
}
