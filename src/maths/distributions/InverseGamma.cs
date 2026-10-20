using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The inverse-gamma distribution IG(shape α, scale β): <c>X = 1/G</c> with <c>G ~ Gamma(α, rate β)</c>, density
/// <c>βᵅ x^(−α−1) e^(−β/x) / Γ(α)</c>. β is a scale of X and the rate of 1/X, so the factory names it once, as
/// <see cref="ShapeScale"/>. The conjugate variance draw <c>σ² ~ IG(m/2, S/2)</c> is
/// <c>InverseGamma.ShapeScale(m/2, S/2)</c>; Zellner–Siow's <c>g ~ IG(½, n/2)</c> is
/// <c>InverseGamma.ShapeScale(0.5, n/2)</c>. A readonly value type: zero-alloc in inner loops.
/// </summary>
public readonly struct InverseGamma : IContinuousDistribution
{
    /// <summary>Shape α.</summary>
    public double Shape { get; }

    /// <summary>Scale β of X (the rate of 1/X).</summary>
    public double Scale { get; }

    private InverseGamma(double shape, double scale)
    {
        Shape = shape;
        Scale = scale;
    }

    /// <summary>IG(shape α, scale β): mean β/(α − 1) for α &gt; 1.</summary>
    public static InverseGamma ShapeScale(double shape, double scale)
    {
        if (!(shape > 0.0) || double.IsPositiveInfinity(shape))
            throw new ArgumentOutOfRangeException(nameof(shape), "Shape must be positive and finite.");
        if (!(scale > 0.0) || double.IsPositiveInfinity(scale))
            throw new ArgumentOutOfRangeException(nameof(scale), "Scale must be positive and finite.");
        return new InverseGamma(shape, scale);
    }

    /// <summary><c>β/(α − 1)</c> for α &gt; 1; +∞ where the mean diverges.</summary>
    public double Mean => Shape > 1.0 ? Scale / (Shape - 1.0) : double.PositiveInfinity;

    /// <summary><c>β²/((α − 1)²(α − 2))</c> for α &gt; 2; +∞ where the variance diverges.</summary>
    public double Variance => Shape > 2.0
        ? Scale * Scale / ((Shape - 1.0) * (Shape - 1.0) * (Shape - 2.0))
        : double.PositiveInfinity;

    /// <summary>The mode <c>β/(α + 1)</c>.</summary>
    public double Mode => Scale / (Shape + 1.0);

    /// <summary>One draw, as the reciprocal of a Gamma(α, rate β) draw.</summary>
    public double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        return 1.0 / Gamma.SampleCore(rng, Shape, 1.0 / Scale);
    }

    /// <summary>Log density; −∞ for <c>x ≤ 0</c>.</summary>
    public double LogPdf(double x)
    {
        if (x <= 0.0) return double.NegativeInfinity;
        return Shape * Math.Log(Scale) - SpecialFunctions.LogGamma(Shape) - (Shape + 1.0) * Math.Log(x) - Scale / x;
    }

    public double Cdf(double x) => Math.Exp(LogCdf(x));
    public double Survival(double x) => Math.Exp(LogSurvival(x));

    // P(X ≤ x) = P(1/X ≥ 1/x) = Q(α, β/x).
    public double LogCdf(double x) => x <= 0.0 ? double.NegativeInfinity : SpecialFunctions.LogRegularizedGamma(Shape, Scale / x).LogQ;
    public double LogSurvival(double x) => x <= 0.0 ? 0.0 : SpecialFunctions.LogRegularizedGamma(Shape, Scale / x).LogP;
    public double IntervalProbability(double lower, double upper) => Math.Exp(LogIntervalProbability(lower, upper));

    public double LogIntervalProbability(double lower, double upper)
    {
        if (Tails.IsEmpty(lower, upper) || upper <= 0.0) return double.NegativeInfinity;
        return Tails.LogInterval(LogCdf(lower), LogSurvival(lower), LogCdf(upper), LogSurvival(upper));
    }
}
