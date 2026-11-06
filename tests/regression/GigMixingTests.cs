using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// GIG mixing (E1's slice 3), the generalized hyperbolic: its marginal against quadrature of the mixture itself, ∫ N(r | m·v, s²·v)
/// p(v) dv over ln v with the normal written out and the GIG normalized by its own quadrature, so nothing is shared with the Bessel K
/// closed form; across the Normal-Wald, the skewed t (a = 0), the variance gamma (b = 0) and general GIGs, drifts of both signs and none.
/// The skewed t falls out of t mixing at a drift and meets the closed t at none; the conditional's draws carry its GIG's mean.
/// </summary>
public sealed class GigMixingTests
{
    [Fact]
    public void Marginal_IsTheMixturesIntegral()
    {
        var rng = new Xoshiro256PlusPlus(20261112);
        var failures = new List<string>();
        for (int c = 0; c < 600; c++)
        {
            GigParameters prior = (c % 4) switch
            {
                0 => new GigParameters(-0.5, Math.Pow(LogUniform(rng, 1e-3, 1e3), 2.0), 1.0),
                // Shapes from 0.2: below it the reference quadrature's normalizer decays too slowly to hold (the GIG's own tests read them).
                1 => GigParameters.OfInverseGamma(LogUniform(rng, 0.2, 50.0), LogUniform(rng, 0.05, 50.0)),
                2 => GigParameters.OfGamma(LogUniform(rng, 0.2, 50.0), LogUniform(rng, 0.05, 50.0)),
                _ => new GigParameters(Signed(rng) * LogUniform(rng, 1e-3, 30.0), LogUniform(rng, 1e-4, 1e3), LogUniform(rng, 1e-4, 1e3)),
            };
            var mixing = new GigMixing(prior);
            double s2 = LogUniform(rng, 1e-4, 1e4), m = rng.NextDouble() < 0.2 ? 0.0 : Signed(rng) * LogUniform(rng, 1e-3, 10.0) * Math.Sqrt(s2);
            double r = Signed(rng) * LogUniform(rng, 1e-4, 30.0) * Math.Sqrt(s2);
            if (prior.B == 0.0 && prior.P <= 0.5) r = Signed(rng) * LogUniform(rng, 1e-2, 30.0) * Math.Sqrt(s2);   // away from the pole at 0
            double logZ = Quadrature.LogIntegral(u => prior.P * u - 0.5 * (prior.A * Math.Exp(u) + prior.B * Math.Exp(-u)), 0.0);
            double Integrand(double u)
            {
                double v = Math.Exp(u), z = r - m * v;
                return -0.5 * Math.Log(2.0 * Math.PI * s2 * v) - z * z / (2.0 * s2 * v)
                       + prior.P * u - 0.5 * (prior.A * v + prior.B / v) - logZ;
            }
            double expected = Quadrature.LogIntegral(Integrand, 0.0), got = mixing.MarginalLogDensity(r, s2, m);
            double gap = Math.Abs(got - expected) / Math.Max(1.0, Math.Abs(expected));
            if (!(gap <= 1e-8)) failures.Add($"GIG({prior.P:R}, {prior.A:R}, {prior.B:R}), s² {s2:R}, m {m:R}, r {r:R}: {got:R} against {expected:R}");
        }
        Assert.True(failures.Count == 0, string.Join("\n", failures.Take(20)));
    }

    [Theory]
    [InlineData(3.0, 0.7)]
    [InlineData(0.6, 2.5)]
    [InlineData(40.0, 0.01)]
    public void SkewedT_FallsOutOfTMixing(double nu, double s2)
    {
        var t = new StudentTMixing(nu);
        var gig = new GigMixing(t.Prior);
        foreach (double r in new[] { -40.0, -1.3, 0.0, 0.2, 7.0 })
        {
            Assert.Equal(gig.MarginalLogDensity(r, s2, 0.4), t.MarginalLogDensity(r, s2, 0.4));
            // No drift: the closed t, and GIG mixing's a = 0 limit, agree.
            Assert.Equal(t.MarginalLogDensity(r, s2, 0.0), gig.MarginalLogDensity(r, s2, 0.0), 10);
            // A tiny drift moves the closed t's value by about m·r/s², its first-order term, and no more.
            double closed = t.MarginalLogDensity(r, s2, 0.0), drifted = t.MarginalLogDensity(r, s2, 1e-12);
            Assert.True(Math.Abs(drifted - closed) <= 2.0 * Math.Abs(1e-12 * r / s2) + 1e-10 * Math.Max(1.0, Math.Abs(closed)), $"{drifted:R} against {closed:R}");
        }
    }

    [Theory]
    [InlineData(-0.5, 4.0, 1.0, 0.3, 1.7, 0.6)]
    [InlineData(-3.0, 0.0, 6.0, -0.8, 0.5, 2.0)]
    [InlineData(1.5, 2.0, 0.0, 1.1, -0.4, 0.3)]
    public void Conditional_IsItsGig(double p, double a, double b, double r, double m, double s2)
    {
        var mixing = new GigMixing(new GigParameters(p, a, b));
        var law = GeneralizedInverseGaussian.Of(p - 0.5, a + m * m / s2, b + r * r / s2);
        var rng = new Xoshiro256PlusPlus(20261113);
        const int draws = 200_000;
        double sum = 0.0;
        for (int i = 0; i < draws; i++) sum += mixing.SampleConditional(r, s2, m, 1.0, rng);
        double z = (sum / draws - law.Mean) / Math.Sqrt(law.Variance / draws);
        Assert.True(Math.Abs(z) < 4.0, $"z = {z:F2}");
        Assert.Equal(GeneralizedInverseGaussian.Of(p, a, b).Mean, mixing.Mean, 12);
    }

    [Fact]
    public void NormalWald_IsTheInverseGaussian()
    {
        var mixing = GigMixing.NormalWald(2.5);
        Assert.Equal(1.0 / 2.5, mixing.Mean, 14);
        Assert.Equal(1.0 / (2.5 * 2.5 * 2.5), mixing.Variance, 14);
        Assert.Throws<ArgumentOutOfRangeException>(() => GigMixing.NormalWald(0.0));
    }

    private static double LogUniform(Xoshiro256PlusPlus rng, double lo, double hi) => Math.Exp(Math.Log(lo) + rng.NextDouble() * (Math.Log(hi) - Math.Log(lo)));

    private static double Signed(Xoshiro256PlusPlus rng) => rng.NextDouble() < 0.5 ? -1.0 : 1.0;
}
