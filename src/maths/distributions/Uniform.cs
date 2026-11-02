using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The continuous uniform distribution on [lower, upper], built with <see cref="On"/>. Its interval moments are closed
/// forms at every power: ((b − c)ᵏ⁺¹ − (a − c)ᵏ⁺¹)/((k + 1)(upper − lower)) on the interval clipped to the support.
/// </summary>
public readonly struct Uniform : IContinuousDistribution, IIntervalMoments
{
    /// <summary>The lower end of the support.</summary>
    public double Lower { get; }

    /// <summary>The upper end of the support.</summary>
    public double Upper { get; }

    private Uniform(double lower, double upper)
    {
        Lower = lower;
        Upper = upper;
    }

    /// <summary>U(<paramref name="lower"/>, <paramref name="upper"/>), finite with lower below upper.</summary>
    public static Uniform On(double lower, double upper)
    {
        if (!double.IsFinite(lower) || !double.IsFinite(upper) || !(upper > lower))
            throw new ArgumentException($"A uniform needs finite lower < upper, got [{lower}, {upper}].");
        return new Uniform(lower, upper);
    }

    /// <summary>The width upper − lower.</summary>
    public double Width => Upper - Lower;

    /// <summary>One draw.</summary>
    public double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        return Lower + rng.NextDouble() * Width;
    }

    /// <inheritdoc/>
    public double SampleOn(Xoshiro256PlusPlus rng, double lower, double upper)
    {
        ArgumentNullException.ThrowIfNull(rng);
        (lower, upper) = Tails.Clip(lower, upper, Lower, Upper);
        return Math.Clamp(lower + rng.NextDouble() * (upper - lower), lower, upper);
    }

    /// <inheritdoc/>
    public double Quantile(double p)
    {
        Tails.CheckProbability(p);
        return Math.Min(Lower + p * Width, Upper);
    }

    /// <inheritdoc/>
    public double InverseLogCdf(double logP)
    {
        Tails.CheckLogProbability(logP);
        return Math.Min(Lower + Math.Exp(logP) * Width, Upper);
    }

    /// <inheritdoc/>
    public double InverseLogSurvival(double logS)
    {
        Tails.CheckLogProbability(logS);
        return Math.Max(Upper - Math.Exp(logS) * Width, Lower);
    }

    public double LogPdf(double x) => x >= Lower && x <= Upper ? -Math.Log(Width) : double.NegativeInfinity;
    public double Cdf(double x) => Math.Clamp((x - Lower) / Width, 0.0, 1.0);
    public double Survival(double x) => Math.Clamp((Upper - x) / Width, 0.0, 1.0);
    public double LogCdf(double x) => Math.Log(Cdf(x));
    public double LogSurvival(double x) => Math.Log(Survival(x));

    public double IntervalProbability(double lower, double upper)
    {
        if (Tails.IsEmpty(lower, upper)) return 0.0;
        return Math.Max(0.0, Math.Min(upper, Upper) - Math.Max(lower, Lower)) / Width;
    }

    public double LogIntervalProbability(double lower, double upper) => Math.Log(IntervalProbability(lower, upper));

    /// <inheritdoc/>
    public double IntervalMoments(double lower, double upper, double centre, double logScale, Span<double> moments)
    {
        moments.Clear();
        if (Tails.IsEmpty(lower, upper)) return 0.0;
        double a = Math.Max(lower, Lower) - centre, b = Math.Min(upper, Upper) - centre;
        if (!(b > a)) return 0.0;
        double density = Math.Exp(-logScale - Math.Log(Width)), pa = a, pb = b;
        for (int k = 0; k < moments.Length; k++)
        {
            moments[k] = (pb - pa) / (k + 1) * density;
            pa *= a;
            pb *= b;
        }
        return 0.0;
    }
}
