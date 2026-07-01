using System;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Spline.Bars;

/// <summary>Uniform proposal kernel: ignores the center and draws on (0,1) — the Denison baseline, no locality.</summary>
public sealed class UniformKernel : IKnotKernel
{
    public double Sample(double center, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        return rng.NextDouble();
    }

    // Density 1 on the proposal support; candidates are constructed in [0,1), so the log-density is 0.
    public double LogDensity(double x, double center) => 0.0;
}

/// <summary>
/// DMGK locality kernel: a knot near center <c>c</c> is proposed from <c>Beta(cτ, (1−c)τ)</c>, concentrating
/// new knots where the curve already has them (roughness clusters). Larger <c>τ</c> is tighter; the asymmetric
/// density supplies the Hastings correction. The center is clamped off the open-interval endpoints so the Beta
/// shapes stay positive.
/// </summary>
public sealed class LocalBetaKernel : IKnotKernel
{
    private const double Eps = 1e-6;
    private readonly double _tau;

    public LocalBetaKernel(double tau = 50.0)
    {
        if (!(tau > 0.0)) throw new ArgumentOutOfRangeException(nameof(tau), "Spread τ must be positive.");
        _tau = tau;
    }

    private (double A, double B) Shapes(double center)
    {
        double c = Math.Clamp(center, Eps, 1.0 - Eps);
        return (c * _tau, (1.0 - c) * _tau);
    }

    public double Sample(double center, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        var (a, b) = Shapes(center);
        return Beta.Sample(rng, a, b);
    }

    public double LogDensity(double x, double center)
    {
        var (a, b) = Shapes(center);
        return Beta.LogPdf(x, a, b);
    }
}

/// <summary>Shared proposal-density math for the knot moves.</summary>
public static class ProposalMath
{
    /// <summary>
    /// Log of the mixture proposal density <c>(1/n) Σᵢ h(x | centerᵢ)</c> — the density of a knot proposed by
    /// first picking an existing knot uniformly, then perturbing it with <paramref name="kernel"/> (DMGK birth).
    /// Empty centers give the base uniform density (log 0). Computed by log-sum-exp for stability.
    /// </summary>
    public static double LogMixtureDensity(double x, double[] centers, IKnotKernel kernel)
    {
        ArgumentNullException.ThrowIfNull(centers);
        ArgumentNullException.ThrowIfNull(kernel);

        int n = centers.Length;
        if (n == 0) return 0.0;

        Span<double> ld = n <= 64 ? stackalloc double[n] : new double[n];
        double max = double.NegativeInfinity;
        for (int i = 0; i < n; i++)
        {
            ld[i] = kernel.LogDensity(x, centers[i]);
            if (ld[i] > max) max = ld[i];
        }
        if (double.IsNegativeInfinity(max)) return double.NegativeInfinity;

        double sum = 0.0;
        for (int i = 0; i < n; i++) sum += Math.Exp(ld[i] - max);
        return max + Math.Log(sum) - Math.Log(n);
    }
}
