using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The Weibull distribution W(shape k, scale λ): density <c>(k/λ)(x/λ)^(k−1) e^(−(x/λ)^k)</c> on x &gt; 0 and survival
/// <c>e^(−(x/λ)^k)</c>. k = 1 is the exponential with mean λ. With k = ½ it is the law of τ² when τ is exponential with
/// rate 1/√λ: the scale-dependent prior on a variance (Klein &amp; Kneib 2016), which Bach &amp; Klein (2025) put on
/// smoothing variances. A readonly value type: zero-alloc in inner loops.
/// </summary>
public readonly struct Weibull : IContinuousDistribution
{
    /// <summary>Shape k.</summary>
    public double Shape { get; }

    /// <summary>Scale λ.</summary>
    public double Scale { get; }

    private Weibull(double shape, double scale)
    {
        Shape = shape;
        Scale = scale;
    }

    /// <summary>W(shape k, scale λ): mean λΓ(1 + 1/k).</summary>
    public static Weibull ShapeScale(double shape, double scale)
    {
        if (!(shape > 0.0) || double.IsPositiveInfinity(shape))
            throw new ArgumentOutOfRangeException(nameof(shape), "Shape must be positive and finite.");
        if (!(scale > 0.0) || double.IsPositiveInfinity(scale))
            throw new ArgumentOutOfRangeException(nameof(scale), "Scale must be positive and finite.");
        return new Weibull(shape, scale);
    }

    /// <summary><c>λΓ(1 + 1/k)</c>.</summary>
    public double Mean => Scale * Math.Exp(SpecialFunctions.LogGamma(1.0 + 1.0 / Shape));

    /// <summary><c>λ²[Γ(1 + 2/k) − Γ(1 + 1/k)²]</c>.</summary>
    public double Variance
    {
        get
        {
            double g1 = Math.Exp(SpecialFunctions.LogGamma(1.0 + 1.0 / Shape));
            return Scale * Scale * (Math.Exp(SpecialFunctions.LogGamma(1.0 + 2.0 / Shape)) - g1 * g1);
        }
    }

    /// <summary><c>λ(ln 2)^(1/k)</c>.</summary>
    public double Median => Scale * Math.Pow(SpecialFunctions.Ln2, 1.0 / Shape);

    /// <summary><c>λ((k − 1)/k)^(1/k)</c> for k &gt; 1; for k ≤ 1 the density is largest at 0, where it is unbounded below
    /// k = 1, and the mode is 0.</summary>
    public double Mode => Shape > 1.0 ? Scale * Math.Pow((Shape - 1.0) / Shape, 1.0 / Shape) : 0.0;

    /// <summary>One draw by inversion, <c>λ(−ln(1 − U))^(1/k)</c>.</summary>
    public double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        return Scale * Math.Pow(-Math.Log(1.0 - rng.NextDouble()), 1.0 / Shape);
    }

    /// <summary>The p-quantile <c>λ(−ln(1 − p))^(1/k)</c>, with ln(1 − p) taken accurately for small p.</summary>
    public double Quantile(double p)
    {
        if (!(p >= 0.0 && p <= 1.0)) throw new ArgumentOutOfRangeException(nameof(p), "A probability must lie in [0, 1].");
        return Scale * Math.Pow(-SpecialFunctions.Log1p(-p), 1.0 / Shape);
    }

    /// <summary>Log density; −∞ for <c>x ≤ 0</c>.</summary>
    public double LogPdf(double x)
    {
        if (x <= 0.0) return double.NegativeInfinity;
        double z = x / Scale;
        return Math.Log(Shape / Scale) + (Shape - 1.0) * Math.Log(z) - Math.Pow(z, Shape);
    }

    public double Cdf(double x) => Math.Exp(LogCdf(x));
    public double Survival(double x) => Math.Exp(LogSurvival(x));

    // ln S(x) = −(x/λ)^k exactly; ln F(x) = ln(1 − e^(−(x/λ)^k)), accurate in both tails.
    public double LogCdf(double x) => x <= 0.0 ? double.NegativeInfinity : SpecialFunctions.Log1mExp(-Math.Pow(x / Scale, Shape));
    public double LogSurvival(double x) => x <= 0.0 ? 0.0 : -Math.Pow(x / Scale, Shape);
    public double IntervalProbability(double lower, double upper) => Math.Exp(LogIntervalProbability(lower, upper));

    public double LogIntervalProbability(double lower, double upper)
    {
        if (Tails.IsEmpty(lower, upper) || upper <= 0.0) return double.NegativeInfinity;
        return Tails.LogInterval(LogCdf(lower), LogSurvival(lower), LogCdf(upper), LogSurvival(upper));
    }
}
