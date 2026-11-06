using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The beta distribution Beta(a, b) on (0, 1), built with <see cref="Shapes"/> and sampled as
/// <c>G_a / (G_a + G_b)</c> from two unit-scale <see cref="Gamma"/> draws. The home of BARS's locality proposal: a
/// knot at center <c>ξ</c> with spread <c>τ</c> proposes from <c>Beta(ξτ, (1−ξ)τ)</c>, and the (asymmetric)
/// <see cref="LogPdf"/> supplies the Hastings correction. The CDF is <c>I_x(a, b)</c>. A readonly value type:
/// zero-alloc in the proposal inner loop.
/// <para>In y = logit x the density is log-concave at all shapes, ∝ e^(ay)/(1 + e^y)^(a+b), so <see cref="SampleOn"/> is
/// rejection under <see cref="LogConcaveEnvelope"/> and the quantiles are Newton on the log tail there. Both solve in y, so
/// x and 1 − x keep their digits: an x within 2⁻⁵³ of 1 reads as 1, and <c>Beta.Shapes(b, a).InverseLogCdf(logS)</c> is the
/// complement at full precision.</para>
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

    /// <inheritdoc/>
    public double SampleOn(Xoshiro256PlusPlus rng, double lower, double upper)
    {
        ArgumentNullException.ThrowIfNull(rng);
        (lower, upper) = Tails.Clip(lower, upper, 0.0, 1.0);
        double mode = A / (A + B), anchor = Math.Clamp(mode, lower, upper);
        double complement = anchor == mode ? B / (A + B) : 1.0 - anchor;   // the mode's holds its digits where b is tiny
        double left = LogConcaveEnvelope.LogRatio(lower, anchor) + LogConcaveEnvelope.LogRatio(complement, 1.0 - lower);
        double right = LogConcaveEnvelope.LogRatio(anchor, upper) + LogConcaveEnvelope.LogRatio(1.0 - upper, complement);
        double reach = LogConcaveEnvelope.Reach(Math.Abs(A - (A + B) * anchor), (A + B) * anchor * complement);
        double delta = LogConcaveEnvelope.Draw(rng, new LogitRatio(A, B, anchor, complement), left, right, reach, reach);
        return Math.Clamp(FromOffset(anchor, complement, delta), lower, upper);
    }

    // The x at logit x₀ + δ, with E = e^δ − 1: x₀ + x₀c₀E/(1 + x₀E) where δ is small, so x₀ keeps its digits, and the form
    // that does not cancel beyond that. Above the anchor these are written in the complement, whose exponential cannot overflow.
    private static double FromOffset(double x0, double c0, double delta)
    {
        if (delta < 0.0)
        {
            if (delta > -0.5)
            {
                double e = SpecialFunctions.Expm1(delta);
                return x0 + x0 * c0 * e / (1.0 + x0 * e);
            }
            double ratio = x0 * Math.Exp(delta);
            return ratio / (c0 + ratio);
        }
        if (delta < 0.5)
        {
            double e = SpecialFunctions.Expm1(-delta);
            return x0 - x0 * c0 * e / (1.0 + c0 * e);
        }
        double opposite = c0 * Math.Exp(-delta);
        return 1.0 - opposite / (x0 + opposite);
    }

    // h(δ) on logit x = logit x₀ + δ: aδ − (a + b) ln(c₀ + x₀e^δ), written below and above the anchor so neither overflows.
    internal readonly struct LogitRatio(double a, double b, double x0, double c0) : IConcaveLogRatio
    {
        public double Value(double delta)
        {
            if (delta < 0.0)
            {
                double below = delta > -0.5 ? SpecialFunctions.Log1p(x0 * SpecialFunctions.Expm1(delta)) : Math.Log(c0 + x0 * Math.Exp(delta));
                return a * delta - (a + b) * below;
            }
            double above = delta < 0.5 ? SpecialFunctions.Log1p(c0 * SpecialFunctions.Expm1(-delta)) : Math.Log(x0 + c0 * Math.Exp(-delta));
            return -b * delta - (a + b) * above;
        }

        public double Slope(double delta)
        {
            if (delta < 0.0)
            {
                double ratio = x0 * Math.Exp(delta);
                return a - (a + b) * (ratio / (c0 + ratio));
            }
            double opposite = c0 * Math.Exp(-delta);
            return (a + b) * (opposite / (x0 + opposite)) - b;
        }
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
        return Invert(logP, SpecialFunctions.Log1mExp(logP));
    }

    /// <inheritdoc/>
    public double InverseLogSurvival(double logS)
    {
        Tails.CheckLogProbability(logS);
        return Invert(SpecialFunctions.Log1mExp(logS), logS);
    }

    // Solved on the smaller tail. Each has an asymptote in y, I_x ≈ xᵃ/(aB) and its reflection, whose root is a bound on
    // the far side and the value itself where x or 1 − x underflows; the logit's normal approximation starts the bulk.
    private double Invert(double logLower, double logUpper)
    {
        if (double.IsNegativeInfinity(logLower)) return 0.0;
        if (double.IsNegativeInfinity(logUpper)) return 1.0;
        double logBeta = SpecialFunctions.LogBeta(A, B);
        double logX = (logLower + Math.Log(A) + logBeta) / A, logC = (logUpper + Math.Log(B) + logBeta) / B;
        if (logX < -700.0) return Math.Exp(logX);
        if (logC < -700.0) return 1.0;

        bool lower = logLower <= logUpper;
        double z = lower ? Normal.StandardInverseLogCdf(logLower) : -Normal.StandardInverseLogCdf(logUpper);
        double centre = SpecialFunctions.Digamma(A) - SpecialFunctions.Digamma(B) + z * Math.Sqrt(1.0 / A + 1.0 / B);
        double start = lower
            ? Math.Max(logX < 0.0 ? logX - SpecialFunctions.Log1mExp(logX) : double.NegativeInfinity, centre)
            : Math.Min(logC < 0.0 ? SpecialFunctions.Log1mExp(logC) - logC : double.PositiveInfinity, centre);

        double y = TailInversion.Solve(new LogitTail(A, B, logBeta, lower), lower ? logLower : logUpper, start,
                                       double.NegativeInfinity, double.PositiveInfinity, lower, scale: 1.0, out double correction);
        (double x, double c, _, _) = FromLogit(y);
        return x + x * c * correction;
    }

    // (x, 1 − x, ln x, ln(1 − x)) from y = logit x, for any y.
    private static (double X, double C, double LogX, double LogC) FromLogit(double y)
    {
        double r = Math.Exp(-Math.Abs(y)), logOnePlus = SpecialFunctions.Log1p(r);
        return y >= 0.0 ? (1.0 / (1.0 + r), r / (1.0 + r), -logOnePlus, -y - logOnePlus)
                        : (r / (1.0 + r), 1.0 / (1.0 + r), y - logOnePlus, -logOnePlus);
    }

    // ln I or its complement in y = logit x, whose slope is ±x(1 − x)f(x)/tail.
    private readonly struct LogitTail(double a, double b, double logBeta, bool lower) : ILogTail
    {
        public double At(double y, out double slope)
        {
            (double x, double c, double logX, double logC) = FromLogit(y);
            (double logI, double logComplement) = SpecialFunctions.LogRegularizedBeta(x, c, logX, logC, a, b);
            double tail = lower ? logI : logComplement, rate = Math.Exp(a * logX + b * logC - logBeta - tail);
            slope = lower ? rate : -rate;
            return tail;
        }
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
