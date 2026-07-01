using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The beta distribution as a sampler+density primitive, sampled as <c>G_a / (G_a + G_b)</c> from two
/// <see cref="Gamma"/> draws. The home of BARS's locality proposal: a knot at center <c>ξ</c> with spread
/// <c>τ</c> proposes from <c>Beta(ξτ, (1−ξ)τ)</c>, and the (asymmetric) <see cref="LogPdf"/> supplies the
/// Hastings correction. Static and zero-alloc for the proposal inner loop.
/// </summary>
public static class Beta
{
    /// <summary>One draw from <c>Beta(a, b)</c> on (0,1) (mean <c>a/(a+b)</c>).</summary>
    public static double Sample(Xoshiro256PlusPlus rng, double a, double b)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (!(a > 0.0)) throw new ArgumentOutOfRangeException(nameof(a), "Shape a must be positive.");
        if (!(b > 0.0)) throw new ArgumentOutOfRangeException(nameof(b), "Shape b must be positive.");

        double ga = Gamma.Sample(rng, a);
        double gb = Gamma.Sample(rng, b);
        double sum = ga + gb;
        return sum > 0.0 ? ga / sum : 0.5;   // degenerate underflow guard
    }

    /// <summary>Log density of <c>Beta(a, b)</c> at <paramref name="x"/>; −∞ outside (0,1).</summary>
    public static double LogPdf(double x, double a, double b)
    {
        if (!(a > 0.0)) throw new ArgumentOutOfRangeException(nameof(a));
        if (!(b > 0.0)) throw new ArgumentOutOfRangeException(nameof(b));
        if (x <= 0.0 || x >= 1.0) return double.NegativeInfinity;

        double logBeta = SpecialFunctions.LogGamma(a) + SpecialFunctions.LogGamma(b)
                       - SpecialFunctions.LogGamma(a + b);
        return (a - 1.0) * Math.Log(x) + (b - 1.0) * Math.Log(1.0 - x) - logBeta;
    }
}
