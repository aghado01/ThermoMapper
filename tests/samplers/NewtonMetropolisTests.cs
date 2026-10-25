using System;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Xunit;

namespace Maths.Samplers.Tests;

/// <summary>
/// The Newton-proposal Metropolis step. On a correlated Gaussian target the proposal is the target itself, so every step
/// is accepted and the draws are independent with the target's moments. On a product of a Student t (whose log density
/// is convex in the tails, so the eigenvalue cap acts there) and a log-gamma variable, the chain has the target's moments.
/// </summary>
public sealed class NewtonMetropolisTests
{
    [Fact]
    public void GaussianTarget_IsProposedExactly_AndAlwaysAccepted()
    {
        // N(μ, Σ) with μ = (1, −2), Σ = [[2, 0.9], [0.9, 1]]: precision Q = Σ⁻¹.
        double det = 2.0 * 1.0 - 0.9 * 0.9;
        double q00 = 1.0 / det, q11 = 2.0 / det, q01 = -0.9 / det;
        double LogTarget(ReadOnlySpan<double> x, Span<double> g, Span<double> h)
        {
            double a = x[0] - 1.0, b = x[1] + 2.0;
            g[0] = -(q00 * a + q01 * b);
            g[1] = -(q01 * a + q11 * b);
            h[0] = -q00; h[1] = -q01; h[2] = -q01; h[3] = -q11;
            return -0.5 * (q00 * a * a + 2.0 * q01 * a * b + q11 * b * b);
        }

        var step = new NewtonMetropolis(2);
        var rng = new Xoshiro256PlusPlus(seed: 4);
        double[] x = [5.0, 5.0];
        const int n = 40_000;
        double s0 = 0, s1 = 0, s00 = 0, s11 = 0, s01 = 0;
        int accepted = 0;
        for (int i = 0; i < n; i++)
        {
            if (step.Step(x, LogTarget, rng)) accepted++;
            s0 += x[0]; s1 += x[1]; s00 += x[0] * x[0]; s11 += x[1] * x[1]; s01 += x[0] * x[1];
        }
        Assert.Equal(n, accepted);
        double m0 = s0 / n, m1 = s1 / n;
        Assert.InRange(m0, 1.0 - 0.03, 1.0 + 0.03);
        Assert.InRange(m1, -2.0 - 0.03, -2.0 + 0.03);
        Assert.InRange(s00 / n - m0 * m0, 2.0 * 0.96, 2.0 * 1.04);
        Assert.InRange(s11 / n - m1 * m1, 0.96, 1.04);
        Assert.InRange(s01 / n - m0 * m1, 0.9 - 0.05, 0.9 + 0.05);
    }

    [Fact]
    public void HeavyTailedTarget_HasItsMoments_WithTheCapActingInTheTails()
    {
        // x₀ ~ t₅ (variance 5/3; log density convex for |x₀| > √5), x₁ = ln of Gamma(3, 1) (mean ψ(3), variance ψ′(3)).
        const double nu = 5.0;
        static double LogTarget(ReadOnlySpan<double> x, Span<double> g, Span<double> h)
        {
            double t = 1.0 + x[0] * x[0] / nu;
            g[0] = -(nu + 1.0) * x[0] / (nu * t);
            h[0] = -(nu + 1.0) * (nu - x[0] * x[0]) / (nu * nu * t * t);
            g[1] = 3.0 - Math.Exp(x[1]);
            h[3] = -Math.Exp(x[1]);
            h[1] = 0.0; h[2] = 0.0;
            return -0.5 * (nu + 1.0) * Math.Log(t) + 3.0 * x[1] - Math.Exp(x[1]);
        }

        var step = new NewtonMetropolis(2);
        var rng = new Xoshiro256PlusPlus(seed: 6);
        double[] x = [0.0, 0.0];
        const int n = 400_000;
        double s0 = 0, s1 = 0, s00 = 0, s11 = 0;
        int accepted = 0, tail = 0;
        for (int i = 0; i < n; i++)
        {
            if (step.Step(x, LogTarget, rng)) accepted++;
            if (Math.Abs(x[0]) > Math.Sqrt(nu)) tail++;
            s0 += x[0]; s1 += x[1]; s00 += x[0] * x[0]; s11 += x[1] * x[1];
        }
        double m0 = s0 / n, m1 = s1 / n;
        Assert.True(tail > n / 50, "the chain visits the convex tails");
        Assert.InRange(m0, -0.03, 0.03);
        Assert.InRange(s00 / n - m0 * m0, 5.0 / 3.0 * 0.92, 5.0 / 3.0 * 1.08);
        const double digamma3 = 0.9227843350984671, trigamma3 = 0.3949340668482264;
        Assert.InRange(m1, digamma3 - 0.01, digamma3 + 0.01);
        Assert.InRange(s11 / n - m1 * m1, trigamma3 * 0.97, trigamma3 * 1.03);
        Assert.InRange(accepted / (double)n, 0.5, 1.0);
    }

    [Fact]
    public void InvalidSettings_AreRefused()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new NewtonMetropolis(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new NewtonMetropolis(2, delta: 0.0));
        var step = new NewtonMetropolis(2);
        Assert.Throws<ArgumentException>(() => step.Step(new double[3], (x, g, h) => 0.0, new Xoshiro256PlusPlus(1)));
    }
}
