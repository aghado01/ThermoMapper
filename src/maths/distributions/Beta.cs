using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The beta distribution Beta(a, b) on (0, 1), built with <see cref="Shapes"/> and sampled as
/// <c>G_a / (G_a + G_b)</c> from two unit-scale <see cref="Gamma"/> draws. The home of BARS's locality proposal: a
/// knot at center <c>ξ</c> with spread <c>τ</c> proposes from <c>Beta(ξτ, (1−ξ)τ)</c>, and the (asymmetric)
/// <see cref="LogPdf"/> supplies the Hastings correction. The CDF is <c>I_x(a, b)</c>. A readonly value type:
/// zero-alloc in the proposal inner loop.
/// </summary>
public readonly struct Beta : IContinuousDistribution
{
    /// <summary>First shape a.</summary>
    public double A { get; }

    /// <summary>Second shape b.</summary>
    public double B { get; }

    private Beta(double a, double b)
    {
        A = a;
        B = b;
    }

    /// <summary>Beta(a, b): mean a/(a + b).</summary>
    public static Beta Shapes(double a, double b)
    {
        if (!(a > 0.0) || double.IsPositiveInfinity(a)) throw new ArgumentOutOfRangeException(nameof(a), "Shape a must be positive and finite.");
        if (!(b > 0.0) || double.IsPositiveInfinity(b)) throw new ArgumentOutOfRangeException(nameof(b), "Shape b must be positive and finite.");
        return new Beta(a, b);
    }

    public double Mean => A / (A + B);
    public double Variance => A * B / ((A + B) * (A + B) * (A + B + 1.0));

    /// <summary>One draw on (0,1). If both gamma draws underflow it returns 0.5.</summary>
    public double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        double ga = Gamma.SampleCore(rng, A, 1.0);
        double gb = Gamma.SampleCore(rng, B, 1.0);
        double sum = ga + gb;
        return sum > 0.0 ? ga / sum : 0.5;   // degenerate underflow guard
    }

    /// <summary>Log density; −∞ outside (0,1).</summary>
    public double LogPdf(double x)
    {
        if (x <= 0.0 || x >= 1.0) return double.NegativeInfinity;
        return (A - 1.0) * Math.Log(x) + (B - 1.0) * Math.Log(1.0 - x) - SpecialFunctions.LogBeta(A, B);
    }

    public double Cdf(double x) => Math.Exp(LogCdf(x));
    public double Survival(double x) => Math.Exp(LogSurvival(x));

    public double LogCdf(double x)
        => x <= 0.0 ? double.NegativeInfinity : x >= 1.0 ? 0.0 : SpecialFunctions.LogRegularizedBeta(x, A, B);

    public double LogSurvival(double x)
        => x <= 0.0 ? 0.0 : x >= 1.0 ? double.NegativeInfinity : SpecialFunctions.LogRegularizedBetaComplement(x, A, B);

    public double IntervalProbability(double lower, double upper) => Math.Exp(LogIntervalProbability(lower, upper));

    public double LogIntervalProbability(double lower, double upper)
    {
        if (Tails.IsEmpty(lower, upper) || upper <= 0.0 || lower >= 1.0) return double.NegativeInfinity;
        return Tails.LogInterval(LogCdf(lower), LogSurvival(lower), LogCdf(upper), LogSurvival(upper));
    }
}
