using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The asymmetric Laplace ALD(μ, σ, q), density <c>q(1 − q)/σ · exp(−ρ_q((x − μ)/σ))</c> with the check loss
/// <c>ρ_q(z) = z(q − 1[z &lt; 0])</c>: the quantile regression likelihood, whose q-quantile is μ (Yu and Moyeed 2001). It is a normal
/// variance-mean mixture, <c>μ + σ(θv + τ√v·ζ)</c> with v ~ Exp(1), θ = (1 − 2q)/(q(1 − q)) and τ² = 2/(q(1 − q)) (Kozumi and
/// Kobayashi 2011), which is how BASS's quantile noise reaches the core (E1). Two exponential tails meet at μ, so every function is in
/// closed form, in logs. A readonly value type.
/// </summary>
public readonly struct AsymmetricLaplace : IContinuousDistribution
{
    private AsymmetricLaplace(double location, double scale, double level)
    {
        Location = location;
        Scale = scale;
        Level = level;
    }

    /// <summary>ALD(<paramref name="location"/>, <paramref name="scale"/>, <paramref name="level"/>): μ finite, σ positive and finite,
    /// q in (0, 1).</summary>
    public static AsymmetricLaplace Of(double location, double scale, double level)
    {
        if (!double.IsFinite(location)) throw new ArgumentOutOfRangeException(nameof(location), "μ must be finite.");
        if (!(scale > 0.0 && double.IsFinite(scale))) throw new ArgumentOutOfRangeException(nameof(scale), "σ must be positive and finite.");
        if (!(level > 0.0 && level < 1.0)) throw new ArgumentOutOfRangeException(nameof(level), "q lies in (0, 1).");
        return new AsymmetricLaplace(location, scale, level);
    }

    /// <summary>μ, the q-quantile.</summary>
    public double Location { get; }

    /// <summary>σ.</summary>
    public double Scale { get; }

    /// <summary>q, the level whose quantile μ is.</summary>
    public double Level { get; }

    /// <summary><c>μ + σ(1 − 2q)/(q(1 − q))</c>.</summary>
    public double Mean => Location + Scale * (1.0 - 2.0 * Level) / (Level * (1.0 - Level));

    /// <summary><c>σ²(1 − 2q + 2q²)/(q²(1 − q)²)</c>.</summary>
    public double Variance
    {
        get
        {
            double q = Level, d = q * (1.0 - q);
            return Scale * Scale * (1.0 - 2.0 * q + 2.0 * q * q) / (d * d);
        }
    }

    // ln q and ln(1 − q).
    private double LogLevel => Math.Log(Level);
    private double LogUpper => SpecialFunctions.Log1p(-Level);

    /// <summary>Log density: <c>ln(q(1 − q)/σ) − ρ_q(z)</c>, z = (x − μ)/σ.</summary>
    public double LogPdf(double x)
    {
        double z = (x - Location) / Scale;
        return LogLevel + LogUpper - Math.Log(Scale) - z * (Level - (z < 0.0 ? 1.0 : 0.0));
    }

    public double Cdf(double x) => Math.Exp(LogCdf(x));
    public double Survival(double x) => Math.Exp(LogSurvival(x));

    /// <summary>ln F: <c>ln q + (1 − q)z</c> below μ, <c>ln(1 − (1 − q)e^(−qz))</c> above.</summary>
    public double LogCdf(double x)
    {
        double z = (x - Location) / Scale;
        return z < 0.0 ? LogLevel + (1.0 - Level) * z : SpecialFunctions.Log1mExp(LogUpper - Level * z);
    }

    /// <summary>ln S: <c>ln(1 − q) − qz</c> above μ, <c>ln(1 − qe^((1 − q)z))</c> below.</summary>
    public double LogSurvival(double x)
    {
        double z = (x - Location) / Scale;
        return z >= 0.0 ? LogUpper - Level * z : SpecialFunctions.Log1mExp(LogLevel + (1.0 - Level) * z);
    }

    public double IntervalProbability(double lower, double upper) => Math.Exp(LogIntervalProbability(lower, upper));

    public double LogIntervalProbability(double lower, double upper)
    {
        if (Tails.IsEmpty(lower, upper)) return double.NegativeInfinity;
        return Tails.LogInterval(LogCdf(lower), LogSurvival(lower), LogCdf(upper), LogSurvival(upper));
    }

    /// <summary>One draw: the left tail with probability q, μ − σE/(1 − q), else μ + σE/q, E a unit exponential.</summary>
    public double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        double e = -Math.Log(1.0 - rng.NextDouble());
        return rng.NextDouble() < Level ? Location - Scale * e / (1.0 - Level) : Location + Scale * e / Level;
    }

    /// <summary>One draw on [<paramref name="lower"/>, <paramref name="upper"/>]: a side by its mass on the interval, then that side's
    /// exponential cut to it.</summary>
    public double SampleOn(Xoshiro256PlusPlus rng, double lower, double upper)
    {
        ArgumentNullException.ThrowIfNull(rng);
        (lower, upper) = Tails.Clip(lower, upper, double.NegativeInfinity, double.PositiveInfinity);
        double left = upper <= Location ? 1.0 : lower >= Location ? 0.0
            : Math.Exp(LogIntervalProbability(lower, Location) - LogIntervalProbability(lower, upper));
        if (rng.NextDouble() < left)
        {
            // Below μ, the distance below min(upper, μ) is exponential at rate (1 − q)/σ, cut where it reaches lower.
            double top = Math.Min(upper, Location), rate = (1.0 - Level) / Scale;
            double e = LogConcaveEnvelope.TruncatedExponential(rng.NextDouble(), 1.0, rate * (top - lower));
            return Math.Clamp(top - e / rate, lower, upper);
        }
        double bottom = Math.Max(lower, Location), up = Level / Scale;
        double f = LogConcaveEnvelope.TruncatedExponential(rng.NextDouble(), 1.0, up * (upper - bottom));
        return Math.Clamp(bottom + f / up, lower, upper);
    }

    /// <summary>The p-quantile: <c>μ + σ·ln(p/q)/(1 − q)</c> up to q, <c>μ − σ·ln((1 − p)/(1 − q))/q</c> above.</summary>
    public double Quantile(double p)
    {
        Tails.CheckProbability(p);
        return p <= Level ? InverseLogCdf(Math.Log(p)) : InverseLogSurvival(SpecialFunctions.Log1p(-p));
    }

    /// <inheritdoc/>
    public double InverseLogCdf(double logP)
    {
        Tails.CheckLogProbability(logP);
        return logP <= LogLevel ? Location + Scale * (logP - LogLevel) / (1.0 - Level)
                                : Location - Scale * (SpecialFunctions.Log1mExp(logP) - LogUpper) / Level;
    }

    /// <inheritdoc/>
    public double InverseLogSurvival(double logS)
    {
        Tails.CheckLogProbability(logS);
        return logS <= LogUpper ? Location - Scale * (logS - LogUpper) / Level
                                : Location + Scale * (SpecialFunctions.Log1mExp(logS) - LogLevel) / (1.0 - Level);
    }
}
