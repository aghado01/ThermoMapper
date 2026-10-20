using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The gamma distribution, built with an explicit parameterization: <see cref="ShapeRate"/> (density
/// <c>βᵅ xᵅ⁻¹ e^(−βx) / Γ(α)</c>, mean α/β) or <see cref="ShapeScale"/> (scale θ = 1/β, mean αθ). Conjugate
/// updates produce rates, so a caller that holds a rate never converts it at the call site; R's BASS τ prior
/// (derivations review §6) is the cautionary case of reading one as the other. Sampling is Marsaglia–Tsang
/// (2000) for shape ≥ 1, with the boost <c>G(a) = G(a+1)·U^(1/a)</c> below 1; <see cref="SampleStandard"/> is the
/// unit-scale draw behind <see cref="Beta"/>. A readonly value type: zero-alloc in inner loops.
/// </summary>
public readonly struct Gamma : IContinuousDistribution
{
    /// <summary>Shape α.</summary>
    public double Shape { get; }

    /// <summary>Rate β = 1/θ.</summary>
    public double Rate { get; }

    /// <summary>Scale θ = 1/β.</summary>
    public double Scale { get; }

    private Gamma(double shape, double rate, double scale)
    {
        Shape = shape;
        Rate = rate;
        Scale = scale;
    }

    /// <summary>Gamma(shape α, rate β): mean α/β.</summary>
    public static Gamma ShapeRate(double shape, double rate)
    {
        CheckShape(shape);
        if (!(rate > 0.0) || double.IsPositiveInfinity(rate))
            throw new ArgumentOutOfRangeException(nameof(rate), "Rate must be positive and finite.");
        return new Gamma(shape, rate, 1.0 / rate);
    }

    /// <summary>Gamma(shape α, scale θ): mean αθ.</summary>
    public static Gamma ShapeScale(double shape, double scale)
    {
        CheckShape(shape);
        if (!(scale > 0.0) || double.IsPositiveInfinity(scale))
            throw new ArgumentOutOfRangeException(nameof(scale), "Scale must be positive and finite.");
        return new Gamma(shape, 1.0 / scale, scale);
    }

    public double Mean => Shape * Scale;
    public double Variance => Shape * Scale * Scale;

    /// <summary>One draw.</summary>
    public double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        return SampleCore(rng, Shape, Scale);
    }

    /// <summary>One draw from Gamma(shape, scale 1).</summary>
    public static double SampleStandard(Xoshiro256PlusPlus rng, double shape)
    {
        ArgumentNullException.ThrowIfNull(rng);
        CheckShape(shape);
        return SampleCore(rng, shape, 1.0);
    }

    // The draw scales inside the loop (d·v·scale), as it always has, so a seeded stream is unchanged.
    internal static double SampleCore(Xoshiro256PlusPlus rng, double shape, double scale)
    {
        if (shape < 1.0)
        {
            double u = 1.0 - rng.NextDouble();   // (0,1]
            return SampleCore(rng, shape + 1.0, scale) * Math.Pow(u, 1.0 / shape);
        }

        double d = shape - 1.0 / 3.0;
        double c = 1.0 / Math.Sqrt(9.0 * d);
        while (true)
        {
            double x, v;
            do
            {
                x = Normal.SampleStandard(rng);
                v = 1.0 + c * x;
            }
            while (v <= 0.0);

            v = v * v * v;
            double u = 1.0 - rng.NextDouble();
            double x2 = x * x;
            if (u < 1.0 - 0.0331 * x2 * x2)
                return d * v * scale;
            if (Math.Log(u) < 0.5 * x2 + d * (1.0 - v + Math.Log(v)))
                return d * v * scale;
        }
    }

    /// <summary>Log density; −∞ for <c>x ≤ 0</c>.</summary>
    public double LogPdf(double x)
    {
        if (x <= 0.0) return double.NegativeInfinity;
        return (Shape - 1.0) * Math.Log(x) - Rate * x + Shape * Math.Log(Rate) - SpecialFunctions.LogGamma(Shape);
    }

    public double Cdf(double x) => Math.Exp(LogCdf(x));
    public double Survival(double x) => Math.Exp(LogSurvival(x));
    public double LogCdf(double x) => x <= 0.0 ? double.NegativeInfinity : SpecialFunctions.LogRegularizedGamma(Shape, Rate * x).LogP;
    public double LogSurvival(double x) => x <= 0.0 ? 0.0 : SpecialFunctions.LogRegularizedGamma(Shape, Rate * x).LogQ;
    public double IntervalProbability(double lower, double upper) => Math.Exp(LogIntervalProbability(lower, upper));

    public double LogIntervalProbability(double lower, double upper)
    {
        if (Tails.IsEmpty(lower, upper) || upper <= 0.0) return double.NegativeInfinity;
        return Tails.LogInterval(LogCdf(lower), LogSurvival(lower), LogCdf(upper), LogSurvival(upper));
    }

    private static void CheckShape(double shape)
    {
        if (!(shape > 0.0) || double.IsPositiveInfinity(shape))
            throw new ArgumentOutOfRangeException(nameof(shape), "Shape must be positive and finite.");
    }
}
