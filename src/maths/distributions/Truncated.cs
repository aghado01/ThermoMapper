using System;

namespace Maths.Distributions;

/// <summary>
/// A distribution conditioned on an interval: <c>X | lower ≤ X ≤ upper</c>, with the parent's interval probability
/// as the normalizer. Because that probability comes from the parent's log tails, a far-tail truncation such as
/// N(0, 1) on [10, ∞) has the right density, φ(x)/Q(10), where 1 − Φ(10) would be 0. Generic over any
/// <see cref="IContinuousDistribution"/> value type, so it adds no allocation; build it with
/// <see cref="Truncated.Of{T}"/>. It has no sampler yet: no caller draws from a truncated law. Its interval moments are the
/// parent's on the interval clipped to the bounds, scaled by the log mass, for a parent that has them.
/// </summary>
public readonly struct Truncated<T> : IContinuousDistribution, IIntervalMoments where T : struct, IContinuousDistribution
{
    /// <summary>The untruncated distribution.</summary>
    public T Parent { get; }

    /// <summary>Lower bound (may be −∞).</summary>
    public double Lower { get; }

    /// <summary>Upper bound (may be +∞).</summary>
    public double Upper { get; }

    /// <summary><c>ln P(lower ≤ X ≤ upper)</c> under the parent.</summary>
    public double LogMass { get; }

    internal Truncated(T parent, double lower, double upper, double logMass)
    {
        Parent = parent;
        Lower = lower;
        Upper = upper;
        LogMass = logMass;
    }

    /// <summary><c>P(lower ≤ X ≤ upper)</c> under the parent (may underflow; <see cref="LogMass"/> does not).</summary>
    public double Mass => Math.Exp(LogMass);

    public double LogPdf(double x) => x >= Lower && x <= Upper ? Parent.LogPdf(x) - LogMass : double.NegativeInfinity;
    public double Cdf(double x) => Math.Exp(LogCdf(x));
    public double Survival(double x) => Math.Exp(LogSurvival(x));

    public double LogCdf(double x)
        => x <= Lower ? double.NegativeInfinity : x >= Upper ? 0.0 : Parent.LogIntervalProbability(Lower, x) - LogMass;

    public double LogSurvival(double x)
        => x <= Lower ? 0.0 : x >= Upper ? double.NegativeInfinity : Parent.LogIntervalProbability(x, Upper) - LogMass;

    public double IntervalProbability(double lower, double upper) => Math.Exp(LogIntervalProbability(lower, upper));

    public double LogIntervalProbability(double lower, double upper)
    {
        if (Tails.IsEmpty(lower, upper)) return double.NegativeInfinity;
        double lo = Math.Max(lower, Lower), hi = Math.Min(upper, Upper);
        if (!(hi > lo)) return double.NegativeInfinity;
        return Parent.LogIntervalProbability(lo, hi) - LogMass;
    }

    /// <inheritdoc/>
    public double IntervalMoments(double lower, double upper, double centre, double logScale, Span<double> moments)
    {
        if (Parent is not IIntervalMoments parent)
            throw new NotSupportedException($"{typeof(T).Name} has no interval moments.");
        moments.Clear();
        if (Tails.IsEmpty(lower, upper)) return 0.0;
        double lo = Math.Max(lower, Lower), hi = Math.Min(upper, Upper);
        if (!(hi > lo)) return 0.0;
        return parent.IntervalMoments(lo, hi, centre, logScale + LogMass, moments);
    }
}

/// <summary>Builds <see cref="Truncated{T}"/> with the parent's type inferred.</summary>
public static class Truncated
{
    /// <summary>
    /// <paramref name="parent"/> conditioned on <c>[lower, upper]</c>. Either bound may be infinite; the interval
    /// must be non-empty and carry mass under the parent.
    /// </summary>
    public static Truncated<T> Of<T>(T parent, double lower, double upper) where T : struct, IContinuousDistribution
    {
        if (double.IsNaN(lower) || double.IsNaN(upper) || !(upper > lower))
            throw new ArgumentException($"Truncation needs lower < upper, got [{lower}, {upper}].");
        double logMass = parent.LogIntervalProbability(lower, upper);
        if (double.IsNegativeInfinity(logMass) || double.IsNaN(logMass))
            throw new ArgumentException($"The truncation interval [{lower}, {upper}] has no mass under the parent.");
        return new Truncated<T>(parent, lower, upper, logMass);
    }
}
