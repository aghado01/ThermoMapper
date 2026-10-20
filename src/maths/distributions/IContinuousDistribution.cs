using System;

namespace Maths.Distributions;

/// <summary>
/// A univariate continuous distribution with log-space tails. Interval probabilities are computed from whichever
/// tail represents them without cancellation, so a far-tail interval keeps its digits as a log. That is what
/// <see cref="Truncated{T}"/> builds on: its normalizer is the parent's interval probability.
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
}

/// <summary>Interval arithmetic on log tails, shared by the distributions.</summary>
internal static class Tails
{
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
