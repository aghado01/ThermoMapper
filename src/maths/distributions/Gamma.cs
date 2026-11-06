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
/// <para>In ln x the density is log-concave at every shape, ∝ e^(αy − βe^y), so <see cref="SampleOn"/> is rejection under
/// <see cref="LogConcaveEnvelope"/> and <see cref="InverseLogCdf"/> is Newton on the log tail there, for shape below 1 as
/// for shape above it. <see cref="InverseGamma"/> is the same coordinate reflected.</para>
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

    /// <summary>ln X of one draw, on <see cref="Sample"/>'s stream: a tiny shape's X falls below the doubles, where ln X still stands.</summary>
    public double SampleLog(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (Shape >= 1.0) return Math.Log(Sample(rng));
        double u = 1.0 - rng.NextDouble();   // (0,1]
        return Math.Log(SampleCore(rng, Shape + 1.0, Scale)) + Math.Log(u) / Shape;
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

    /// <inheritdoc/>
    public double SampleOn(Xoshiro256PlusPlus rng, double lower, double upper)
    {
        ArgumentNullException.ThrowIfNull(rng);
        (lower, upper) = Tails.Clip(lower, upper, 0.0, double.PositiveInfinity);
        double anchor = Math.Clamp(Mean, lower, upper);
        double delta = DrawLogOffset(rng, Shape, Rate * anchor, LogConcaveEnvelope.LogRatio(lower, anchor), LogConcaveEnvelope.LogRatio(anchor, upper));
        return Math.Clamp(LogConcaveEnvelope.Scaled(anchor, delta), lower, upper);
    }

    // δ = ln(g/g₀) for a unit-rate gamma g on [g₀e^−left, g₀e^right], with g₀ the point closest to the shape.
    internal static double DrawLogOffset(Xoshiro256PlusPlus rng, double shape, double g0, double left, double right)
    {
        double rightReach = LogConcaveEnvelope.Reach(Math.Max(0.0, g0 - shape), g0);
        // The quadratic overshoots an exponential: g₀(eᵗ − 1) = 1 + αt, solved once, is the closer start from the same side.
        rightReach = Math.Min(rightReach, SpecialFunctions.Log1p((1.0 + shape * rightReach) / g0));
        return LogConcaveEnvelope.Draw(rng, new LogGammaRatio(shape, g0), left, right,
                                       LogConcaveEnvelope.Reach(Math.Max(0.0, shape - g0), g0), rightReach);
    }

    // h(δ) = αδ − g₀(e^δ − 1): the log-gamma density about ln g₀.
    internal readonly struct LogGammaRatio(double shape, double g0) : IConcaveLogRatio
    {
        public double Value(double delta) => shape * delta - g0 * SpecialFunctions.Expm1(delta);
        public double Slope(double delta) => shape - g0 * Math.Exp(delta);
    }

    /// <inheritdoc/>
    public double Quantile(double p)
    {
        Tails.CheckProbability(p);
        return p <= 0.5 ? InverseLogCdf(Math.Log(p)) : InverseLogSurvival(SpecialFunctions.Log1p(-p));
    }

    /// <inheritdoc/>
    public double InverseLogCdf(double logP)
    {
        Tails.CheckLogProbability(logP);
        return Scale * UnitQuantile(Shape, logP, SpecialFunctions.Log1mExp(logP));
    }

    /// <inheritdoc/>
    public double InverseLogSurvival(double logS)
    {
        Tails.CheckLogProbability(logS);
        return Scale * UnitQuantile(Shape, SpecialFunctions.Log1mExp(logS), logS);
    }

    // The g with ln P(α, g) = logLower and ln Q(α, g) = logUpper, solved on the smaller of the two. ln P ≤ α ln g − ln Γ(α + 1),
    // so that line's root is a bound from below, and its value where g underflows, since the two then agree.
    internal static double UnitQuantile(double shape, double logLower, double logUpper)
    {
        if (double.IsNegativeInfinity(logLower)) return 0.0;
        if (double.IsNegativeInfinity(logUpper)) return double.PositiveInfinity;
        double logGamma = SpecialFunctions.LogGamma(shape);
        double power = (logLower + logGamma + Math.Log(shape)) / shape;
        if (power < -700.0) return Math.Exp(power);

        bool lower = logLower <= logUpper;
        double z = lower ? Normal.StandardInverseLogCdf(logLower) : -Normal.StandardInverseLogCdf(logUpper);
        double cube = 1.0 - 1.0 / (9.0 * shape) + z / (3.0 * Math.Sqrt(shape));   // Wilson and Hilferty 1931
        double centre = cube > 0.0 ? Math.Log(shape) + 3.0 * Math.Log(cube) : double.NegativeInfinity;
        double start = lower ? Math.Max(power, centre) : shape >= 1.0 && cube > 0.0 ? centre : power;
        if (!lower && -logUpper > Math.Max(1.0, 4.0 * shape))
        {
            // Far above the bulk, ln Q ≈ (α − 1) ln g − g − ln Γ(α). Past ε⁻¹ it is the value itself, where the slope's
            // α ln g − g − ln Q is all rounding.
            double g = Math.Max(shape, -logUpper);
            for (int i = 0; i < 2; i++) g = -logUpper - logGamma + (shape - 1.0) * Math.Log(g);
            if (-logUpper > 1e15) return g;
            if (g > Math.Max(shape, 1.0)) start = Math.Log(g);
        }

        double y = TailInversion.Solve(new LogTail(shape, logGamma, lower), lower ? logLower : logUpper, start, -745.0, 710.0, lower, scale: 1.0, out double correction);
        double value = Math.Exp(y);
        return value + value * SpecialFunctions.Expm1(correction);
    }

    // ln P or ln Q in y = ln g, whose slope is ±gf(g)/tail.
    private readonly struct LogTail(double shape, double logGamma, bool lower) : ILogTail
    {
        public double At(double y, out double slope)
        {
            double g = Math.Exp(y);
            (double logP, double logQ) = SpecialFunctions.LogRegularizedGamma(shape, g);
            double tail = lower ? logP : logQ, rate = Math.Exp(shape * y - g - logGamma - tail);
            slope = lower ? rate : -rate;
            return tail;
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

    // Where Rate·x nears the doubles' end, P(a, z) is zᵃ/Γ(a + 1) to double precision, and in logs it stays in range. A tiny shape keeps mass there.
    public double LogCdf(double x)
    {
        if (x <= 0.0) return double.NegativeInfinity;
        double z = Rate * x;
        return z < 1e-300 ? Shape * (Math.Log(Rate) + Math.Log(x)) - SpecialFunctions.LogGamma(Shape + 1.0) : SpecialFunctions.LogRegularizedGamma(Shape, z).LogP;
    }

    public double LogSurvival(double x)
    {
        if (x <= 0.0) return 0.0;
        double z = Rate * x;
        return z < 1e-300 ? SpecialFunctions.Log1mExp(LogCdf(x)) : SpecialFunctions.LogRegularizedGamma(Shape, z).LogQ;
    }

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
