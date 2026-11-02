using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// A univariate continuous distribution with log-space tails, quantiles that read them, and draws. Interval probabilities
/// are computed from whichever tail represents them without cancellation, so a far-tail interval keeps its digits as a log.
/// That is what <see cref="Truncated{T}"/> builds on: its normalizer is the parent's interval probability, and its draws are
/// the parent's <see cref="SampleOn"/>.
/// </summary>
public interface IContinuousDistribution
{
    double LogPdf(double x);
    double Cdf(double x);
    double Survival(double x);
    double LogCdf(double x);
    double LogSurvival(double x);

    /// <summary><c>P(lower &lt; X ≤ upper)</c>; 0 for an empty interval.</summary>
    double IntervalProbability(double lower, double upper);

    /// <summary><c>ln P(lower &lt; X ≤ upper)</c>; −∞ for an empty interval.</summary>
    double LogIntervalProbability(double lower, double upper);

    /// <summary>One draw.</summary>
    double Sample(Xoshiro256PlusPlus rng);

    /// <summary>
    /// One draw conditioned on [<paramref name="lower"/>, <paramref name="upper"/>], either bound possibly infinite. It
    /// needs no normalizer, so it costs little in a loop whose bounds or parameters move; an interval that is empty or
    /// misses the support is refused.
    /// </summary>
    double SampleOn(Xoshiro256PlusPlus rng, double lower, double upper);

    /// <summary>The x with <c>P(X ≤ x) = p</c>, the support's ends at 0 and 1; a p outside [0, 1] is refused.</summary>
    double Quantile(double p);

    /// <summary>The x with <c>ln P(X ≤ x) = logP</c>, for a lower tail below the double range; a positive log is refused.</summary>
    double InverseLogCdf(double logP);

    /// <summary>The x with <c>ln P(X &gt; x) = logS</c>, for an upper tail that 1 − p cannot hold; a positive log is refused.</summary>
    double InverseLogSurvival(double logS);
}

/// <summary>Interval arithmetic on log tails, shared by the distributions.</summary>
internal static class Tails
{
    /// <summary>Refuses a probability outside [0, 1].</summary>
    public static void CheckProbability(double p)
    {
        if (!(p >= 0.0 && p <= 1.0)) throw new ArgumentOutOfRangeException(nameof(p), "A probability must lie in [0, 1].");
    }

    /// <summary>Refuses a log probability above 0.</summary>
    public static void CheckLogProbability(double logP)
    {
        if (!(logP <= 0.0)) throw new ArgumentOutOfRangeException(nameof(logP), "A log probability must be at most 0.");
    }

    /// <summary>[lower, upper] clipped to a support; refuses an interval that is empty or misses it.</summary>
    public static (double Lower, double Upper) Clip(double lower, double upper, double supportLower, double supportUpper)
    {
        if (double.IsNaN(lower) || double.IsNaN(upper) || !(upper > lower))
            throw new ArgumentException($"A truncated draw needs lower < upper, got [{lower}, {upper}].");
        double lo = Math.Max(lower, supportLower), hi = Math.Min(upper, supportUpper);
        if (!(hi > lo)) throw new ArgumentException($"The truncation interval [{lower}, {upper}] has no mass.");
        return (lo, hi);
    }

    /// <summary><c>ln(eᵃ + eᵇ)</c>.</summary>
    public static double LogAdd(double a, double b)
    {
        double high = Math.Max(a, b), low = Math.Min(a, b);
        return double.IsNegativeInfinity(low) ? high : high + SpecialFunctions.Log1p(Math.Exp(low - high));
    }

    /// <summary>Rejects NaN bounds and <c>lower &gt; upper</c>; true when the interval is empty.</summary>
    public static bool IsEmpty(double lower, double upper)
    {
        if (double.IsNaN(lower) || double.IsNaN(upper) || lower > upper)
            throw new ArgumentException($"An interval needs lower ≤ upper, got [{lower}, {upper}].");
        return lower == upper;
    }

    /// <summary>
    /// <c>ln P(lo &lt; X ≤ hi)</c> from the log CDF and log survival at both ends. Above the median it is
    /// S(lo) − S(hi), below it F(hi) − F(lo), so neither subtracts from 1. An interval that straddles the median
    /// is 1 − F(lo) − S(hi); a distribution with an accurate central form (normal, Student-t) uses that instead,
    /// since this one loses digits for an interval narrow around the median.
    /// </summary>
    public static double LogInterval(double logFLo, double logSLo, double logFHi, double logSHi)
    {
        if (logFLo >= -SpecialFunctions.Ln2)
            return double.IsNegativeInfinity(logSLo) ? double.NegativeInfinity
                 : logSLo + SpecialFunctions.Log1mExp(logSHi - logSLo);
        if (logFHi <= -SpecialFunctions.Ln2)
            return double.IsNegativeInfinity(logFHi) ? double.NegativeInfinity
                 : logFHi + SpecialFunctions.Log1mExp(logFLo - logFHi);
        return SpecialFunctions.Log1p(-(Math.Exp(logFLo) + Math.Exp(logSHi)));
    }
}
