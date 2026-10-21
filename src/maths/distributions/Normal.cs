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

    /// <summary>The quantile function: the x with <c>P(X ≤ x) = p</c>; −∞ at 0, +∞ at 1, NaN outside [0, 1].</summary>
    public double Quantile(double p) => Mean + Sd * StandardQuantile(p);

    /// <summary>
    /// The standard normal quantile Φ⁻¹(p), by Wichura's algorithm AS 241 (PPND16, 1988): rational approximations
    /// in three regions, accurate to about 1 part in 10¹⁶. −∞ at 0, +∞ at 1, NaN outside [0, 1].
    /// </summary>
    public static double StandardQuantile(double p)
    {
        if (!(p >= 0.0 && p <= 1.0)) return double.NaN;
        if (p == 0.0) return double.NegativeInfinity;
        if (p == 1.0) return double.PositiveInfinity;

        double q = p - 0.5;
        if (Math.Abs(q) <= 0.425)
        {
            double r = 0.180625 - q * q;
            return q * (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r + 67265.770927008700853) * r
                              + 45921.953931549871457) * r + 13731.693765509461125) * r + 1971.5909503065514427) * r
                              + 133.14166789178437745) * r + 3.387132872796366608)
                     / (((((((r * 5226.495278852545925 + 28729.085735721942674) * r + 39307.89580009271061) * r
                              + 21213.794301586595867) * r + 5394.1960214247511077) * r + 687.1870074920579083) * r
                              + 42.313330701600911252) * r + 1.0);
        }

        double t = Math.Sqrt(-Math.Log(q < 0.0 ? p : 1.0 - p));
        double value;
        if (t <= 5.0)
        {
            t -= 1.6;
            value = (((((((t * 7.7454501427834140764e-4 + 0.0227238449892691845833) * t + 0.24178072517745061177) * t
                          + 1.27045825245236838258) * t + 3.64784832476320460504) * t + 5.7694972214606914055) * t
                          + 4.6303378461565452959) * t + 1.42343711074968357734)
                  / (((((((t * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * t + 0.0151986665636164571966) * t
                          + 0.14810397642748007459) * t + 0.68976733498510000455) * t + 1.6763848301838038494) * t
                          + 2.05319162663775882187) * t + 1.0);
        }
        else
        {
            t -= 5.0;
            value = (((((((t * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * t + 0.0012426609473880784386) * t
                          + 0.026532189526576123093) * t + 0.29656057182850489123) * t + 1.7848265399172913358) * t
                          + 5.4637849111641143699) * t + 6.6579046435011037772)
                  / (((((((t * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * t + 1.8463183175100546818e-5) * t
                          + 7.868691311456132591e-4) * t + 0.0148753612908506148525) * t + 0.13692988092273580531) * t
                          + 0.59983220655588793769) * t + 1.0);
        }
        return q < 0.0 ? -value : value;
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
