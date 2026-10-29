using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The half-Cauchy distribution C⁺(0, s): density <c>2/(πs(1 + (x/s)²))</c> on x ≥ 0, CDF <c>(2/π)·atan(x/s)</c>. As
/// the prior on a standard deviation (Gelman 2006) it is proper and scale-equivariant: σ ~ C⁺(0, s) gives cσ ~ C⁺(0, cs),
/// so a posterior under it moves with the data's units. Its median is s; it has no mean. A readonly value type:
/// zero-alloc in inner loops.
/// </summary>
public readonly struct HalfCauchy : IContinuousDistribution
{
    /// <summary>Scale s, the median.</summary>
    public double Scale { get; }

    private HalfCauchy(double scale) => Scale = scale;

    /// <summary>C⁺(0, <paramref name="scale"/>).</summary>
    public static HalfCauchy WithScale(double scale)
    {
        if (!(scale > 0.0) || double.IsPositiveInfinity(scale))
            throw new ArgumentOutOfRangeException(nameof(scale), "Scale must be positive and finite.");
        return new HalfCauchy(scale);
    }

    /// <summary>The median, s.</summary>
    public double Median => Scale;

    /// <summary>One draw by inversion, <c>s·tan(πU/2)</c>.</summary>
    public double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        return Quantile(rng.NextDouble());
    }

    /// <summary>The p-quantile <c>s·tan(πp/2)</c>; +∞ at p = 1.</summary>
    public double Quantile(double p)
    {
        if (!(p >= 0.0 && p <= 1.0)) throw new ArgumentOutOfRangeException(nameof(p), "A probability must lie in [0, 1].");
        if (p == 1.0) return double.PositiveInfinity;
        // tan(πp/2) = 1/tan(π(1 − p)/2), the form that keeps its digits as p approaches 1.
        return p <= 0.5 ? Scale * Math.Tan(0.5 * Math.PI * p) : Scale / Math.Tan(0.5 * Math.PI * (1.0 - p));
    }

    /// <summary>Log density; −∞ for <c>x &lt; 0</c>.</summary>
    public double LogPdf(double x)
    {
        if (x < 0.0 || double.IsNaN(x)) return double.NegativeInfinity;
        if (double.IsPositiveInfinity(x)) return double.NegativeInfinity;
        double z = x / Scale;
        return Math.Log(2.0 / (Math.PI * Scale)) - SpecialFunctions.Log1p(z * z);
    }

    public double Cdf(double x) => x <= 0.0 ? 0.0 : 2.0 / Math.PI * Math.Atan(x / Scale);

    // S(x) = (2/π)·atan(s/x), exact in the upper tail where 1 − F cancels.
    public double Survival(double x) => x <= 0.0 ? 1.0 : 2.0 / Math.PI * Math.Atan(Scale / x);

    public double LogCdf(double x) => x <= 0.0 ? double.NegativeInfinity : Math.Log(Cdf(x));
    public double LogSurvival(double x) => x <= 0.0 ? 0.0 : Math.Log(Survival(x));
    public double IntervalProbability(double lower, double upper) => Math.Exp(LogIntervalProbability(lower, upper));

    public double LogIntervalProbability(double lower, double upper)
    {
        if (Tails.IsEmpty(lower, upper) || upper <= 0.0) return double.NegativeInfinity;
        return Tails.LogInterval(LogCdf(lower), LogSurvival(lower), LogCdf(upper), LogSurvival(upper));
    }
}
