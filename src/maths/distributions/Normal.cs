using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The normal distribution N(mean, sd²), built with <see cref="MeanSd"/>, <see cref="MeanVariance"/> or
/// <see cref="Standard"/>. <see cref="SampleStandard"/> is the standard-normal draw that also backs the
/// Marsaglia–Tsang <see cref="Gamma"/> sampler. The tails come from <c>Φ(−z) = ½·Q(½, z²/2)</c> in log space, so
/// <c>ln Φ(−40)</c> is exact where <c>Φ(−40)</c> underflows, and an interval straddling the mean is the sum of
/// two error functions rather than a difference of CDFs near ½. A readonly value type: zero-alloc in inner loops.
/// </summary>
public readonly struct Normal : IContinuousDistribution
{
    private static readonly double LogSqrt2Pi = 0.5 * Math.Log(2.0 * Math.PI);

    /// <summary>The mean.</summary>
    public double Mean { get; }

    /// <summary>The standard deviation.</summary>
    public double Sd { get; }

    private Normal(double mean, double sd)
    {
        Mean = mean;
        Sd = sd;
    }

    /// <summary>N(0, 1).</summary>
    public static Normal Standard => new(0.0, 1.0);

    /// <summary>N(mean, sd²) from the standard deviation.</summary>
    public static Normal MeanSd(double mean, double sd)
    {
        if (!double.IsFinite(mean)) throw new ArgumentOutOfRangeException(nameof(mean), "Mean must be finite.");
        if (!(sd > 0.0) || double.IsPositiveInfinity(sd))
            throw new ArgumentOutOfRangeException(nameof(sd), "Standard deviation must be positive and finite.");
        return new Normal(mean, sd);
    }

    /// <summary>N(mean, variance) from the variance.</summary>
    public static Normal MeanVariance(double mean, double variance)
    {
        if (!(variance > 0.0) || double.IsPositiveInfinity(variance))
            throw new ArgumentOutOfRangeException(nameof(variance), "Variance must be positive and finite.");
        return MeanSd(mean, Math.Sqrt(variance));
    }

    /// <summary>The variance <c>sd²</c>.</summary>
    public double Variance => Sd * Sd;

    /// <summary>One standard-normal draw via the Box–Muller transform (cosine branch; two uniforms).</summary>
    public static double SampleStandard(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        double u1 = 1.0 - rng.NextDouble();   // (0,1] avoids log(0)
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    /// <summary>One draw from this distribution.</summary>
    public double Sample(Xoshiro256PlusPlus rng) => Mean + Sd * SampleStandard(rng);

    public double LogPdf(double x)
    {
        double z = (x - Mean) / Sd;
        return -0.5 * z * z - Math.Log(Sd) - LogSqrt2Pi;
    }

    public double Cdf(double x) => Math.Exp(LogCdf(x));
    public double Survival(double x) => Math.Exp(LogSurvival(x));
    public double LogCdf(double x) => StandardLogCdf((x - Mean) / Sd);
    public double LogSurvival(double x) => StandardLogCdf((Mean - x) / Sd);
    public double IntervalProbability(double lower, double upper) => Math.Exp(LogIntervalProbability(lower, upper));

    public double LogIntervalProbability(double lower, double upper)
    {
        if (Tails.IsEmpty(lower, upper)) return double.NegativeInfinity;
        double a = (lower - Mean) / Sd, b = (upper - Mean) / Sd;
        if (a >= 0.0) return SameTail(LogLowerTail(a), LogLowerTail(b));      // Φ(−a) − Φ(−b)
        if (b <= 0.0) return SameTail(LogLowerTail(-b), LogLowerTail(-a));    // Φ(b) − Φ(a)
        return Math.Log(0.5 * (Erf(-a) + Erf(b)));                             // straddles the mean
    }

    /// <summary><c>ln Φ(z)</c> of the standard normal.</summary>
    public static double StandardLogCdf(double z)
        => z < 0.0 ? LogLowerTail(-z) : SpecialFunctions.Log1mExp(LogLowerTail(z));

    // ln Φ(−t) for t ≥ 0: Φ(−t) = ½·erfc(t/√2) = ½·Q(½, t²/2).
    private static double LogLowerTail(double t)
        => double.IsPositiveInfinity(t) ? double.NegativeInfinity
         : SpecialFunctions.LogRegularizedGamma(0.5, 0.5 * t * t).LogQ - SpecialFunctions.Ln2;

    // erf(t/√2) = P(½, t²/2) for t ≥ 0, accurate near 0 where 1 − 2Φ(−t) is not.
    private static double Erf(double t)
        => double.IsPositiveInfinity(t) ? 1.0 : Math.Exp(SpecialFunctions.LogRegularizedGamma(0.5, 0.5 * t * t).LogP);

    // ln(e^hi − e^lo) for the two log tail masses hi ≥ lo.
    private static double SameTail(double logHi, double logLo)
        => double.IsNegativeInfinity(logHi) ? double.NegativeInfinity
         : logHi + SpecialFunctions.Log1mExp(logLo - logHi);
}
