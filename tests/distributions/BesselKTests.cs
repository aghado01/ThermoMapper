using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// The log Bessel K against what shares nothing with it: its integral eˣK_ν(x) = ∫₀^∞ e^(−x(cosh t − 1))·cosh(νt) dt by quadrature in
/// logs, over a seeded sweep of orders and arguments across their far ranges, and the half-integer orders in closed form.
/// </summary>
public sealed class BesselKTests
{
    [Fact]
    public void HalfIntegerOrders_AreClosed()
    {
        // K_{n+½}(x) = √(π/(2x))e^(−x) Σₖ (n + k)!/(k!(n − k)!)(2x)^(−k).
        foreach (double x in new[] { 1e-300, 1e-12, 1e-3, 0.5, 1.999, 2.0, 2.001, 7.0, 150.0, 1e6, 1e300 })
            for (int n = 0; n <= 6; n++)
            {
                double sum = 0.0;
                var logs = new List<double>();
                for (int k = 0; k <= n; k++)
                {
                    logs.Add(SpecialFunctions.LogGamma(n + k + 1.0) - SpecialFunctions.LogGamma(k + 1.0) - SpecialFunctions.LogGamma(n - k + 1.0) - k * Math.Log(2.0 * x));
                }
                double top = logs.Max();
                foreach (double l in logs) sum += Math.Exp(l - top);
                double expected = 0.5 * Math.Log(Math.PI / (2.0 * x)) + top + Math.Log(sum);
                double got = SpecialFunctions.LogBesselKScaled(n + 0.5, x);
                Assert.True(Math.Abs(got - expected) <= 1e-13 * Math.Max(1.0, Math.Abs(expected)), $"K_{n + 0.5}({x:R}): {got:R} against {expected:R}");
                Assert.Equal(got, SpecialFunctions.LogBesselKScaled(-(n + 0.5), x));
            }
    }

    [Fact]
    public void Sweep_MeetsTheIntegral()
    {
        var rng = new Xoshiro256PlusPlus(20261110);
        var failures = new List<string>();
        double worst = 0.0;
        for (int c = 0; c < 2_000; c++)
        {
            double nu = (rng.NextDouble()) switch
            {
                < 0.15 => 0.0,
                < 0.3 => rng.NextDouble(),
                < 0.45 => Math.Round(LogUniform(rng, 1.0, 50.0)) + (rng.NextDouble() < 0.5 ? 0.5 : 0.0),
                _ => LogUniform(rng, 1e-3, 1e4),
            };
            double x = rng.NextDouble() < 0.2 ? LogUniform(rng, 1.0, 3.0) : LogUniform(rng, 1e-300, 1e300);
            double expected = Integral(nu, x), got = SpecialFunctions.LogBesselKScaled(nu, x);
            double gap = Math.Abs(got - expected) / Math.Max(1.0, Math.Abs(expected));
            worst = Math.Max(worst, gap);
            if (!(gap <= 1e-12)) failures.Add($"K_{nu:R}({x:R}): {got:R} against {expected:R}, gap {gap:G3}");
        }
        Assert.True(failures.Count == 0, $"worst {worst:G3}\n" + string.Join("\n", failures.Take(20)));
    }

    [Fact]
    public void Edges()
    {
        Assert.Equal(double.PositiveInfinity, SpecialFunctions.LogBesselK(1.0, 0.0));
        Assert.Equal(double.NegativeInfinity, SpecialFunctions.LogBesselK(1.0, double.PositiveInfinity));
        Assert.True(double.IsNaN(SpecialFunctions.LogBesselK(double.NaN, 1.0)));
        Assert.Throws<ArgumentOutOfRangeException>(() => SpecialFunctions.LogBesselK(1.0, -1.0));
        // Continuous across the methods' seams: x = 2 (Temme to Steed) and ν = 300 (recurrence to Debye).
        Assert.True(Math.Abs(SpecialFunctions.LogBesselKScaled(3.3, Math.BitDecrement(2.0)) - SpecialFunctions.LogBesselKScaled(3.3, Math.BitIncrement(2.0))) < 1e-14);
        Assert.True(Math.Abs(SpecialFunctions.LogBesselKScaled(Math.BitDecrement(300.0), 40.0) - SpecialFunctions.LogBesselKScaled(300.0, 40.0)) < 1e-11);
    }

    // ln ∫₀^∞ exp(g(t)) dt, g(t) = −x(cosh t − 1) + ln cosh(νt): unimodal (its slope ν·tanh(νt) − x·sinh t changes sign at most once), so
    // the peak by bisection on the slope, the ends where g has fallen by 60, and Gauss–Legendre on even panels between them.
    private static double Integral(double nu, double x)
    {
        double G(double t) => -2.0 * x * Math.Pow(Math.Sinh(0.5 * t), 2.0) + LogCosh(nu * t);
        double Slope(double t) => nu * Math.Tanh(nu * t) - x * Math.Sinh(t);
        double peak = 0.0;
        if (nu * nu > x)
        {
            double lo = 0.0, hi = 1.0;
            while (Slope(hi) > 0.0) { lo = hi; hi *= 2.0; }
            for (int i = 0; i < 300; i++)
            {
                double mid = 0.5 * (lo + hi);
                if (Slope(mid) > 0.0) lo = mid;
                else hi = mid;
            }
            peak = 0.5 * (lo + hi);
        }
        double top = G(peak), floor = top - 60.0;
        double Edge(double from, double direction)
        {
            double step = Math.Max(Math.Abs(from) * 1e-12, 1e-300), at = from;
            while (true)
            {
                double next = at + direction * step;
                if (direction < 0.0 && next <= 0.0) return 0.0;
                if (G(next) < floor) { at = next; break; }
                step *= 2.0;
            }
            double inner = from, outer = at;
            for (int i = 0; i < 300; i++)
            {
                double mid = 0.5 * (inner + outer);
                if (G(mid) < floor) outer = mid;
                else inner = mid;
            }
            return outer;
        }
        double left = peak > 0.0 ? Edge(peak, -1.0) : 0.0, right = Edge(peak, 1.0);
        const int panels = 3_000;
        double[] nodes = [0.0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640];
        double[] weights = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891];
        double width = (right - left) / panels, sum = 0.0;
        for (int k = 0; k < panels; k++)
        {
            double mid = left + (k + 0.5) * width;
            for (int i = 0; i < nodes.Length; i++) sum += weights[i] * Math.Exp(G(mid + 0.5 * width * nodes[i]) - top);
        }
        return top + Math.Log(sum * 0.5 * width);
    }

    private static double LogCosh(double u)
    {
        u = Math.Abs(u);
        return u + SpecialFunctions.Log1p(Math.Exp(-2.0 * u)) - Math.Log(2.0);
    }

    private static double LogUniform(Xoshiro256PlusPlus rng, double lo, double hi) => DistributionSweepTests.LogUniform(rng, lo, hi);
}
