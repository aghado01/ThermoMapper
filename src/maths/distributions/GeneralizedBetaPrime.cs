using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The generalized beta prime GBP(p, a, b) on v &gt; 0: density p·v^(pa−1)·(1 + v^p)^(−(a+b))/B(a, b), p, a, b &gt; 0 (McDonald 1984's
/// beta of the second kind with a power; the mixing law of Armagan, Dunson and Clyde 2011). y = v^p/(1 + v^p) is Beta(a, b), so the tails are
/// the beta's (<see cref="Beta"/>), taken in logit y = p·ln v so that both tails keep their digits, its quantiles are the beta's solve in
/// that logit, and a draw is (G_a/G_b)^(1/p) from two unit gammas. In u = ln v the density e^(pau)·(1 + e^(pu))^(−(a+b)), times p/B(a, b),
/// is log-concave, with its mode where e^(pu) = a/b and slopes pa as u → −∞ and −pb as u → +∞. Its moments are E[v^k] = B(a + k/p, b − k/p)/B(a, b) for −pa &lt; k &lt; pb and infinite above.
/// The horseshoe's mixing law is p = 1, a = b = ½, under which √v is half-Cauchy (Carvalho, Polson and Scott 2010); the Strawderman–Berger's
/// is p = ½, a = 1, b = ½. A readonly value type.
/// </summary>
public readonly struct GeneralizedBetaPrime : IContinuousDistribution
{
    private readonly double _logNormalizer;   // ln B(a, b) − ln p

    private GeneralizedBetaPrime(double p, double a, double b)
    {
        P = p;
        A = a;
        B = b;
        _logNormalizer = SpecialFunctions.LogBeta(a, b) - Math.Log(p);
    }

    /// <summary>GBP(<paramref name="p"/>, <paramref name="a"/>, <paramref name="b"/>), each positive and finite.</summary>
    public static GeneralizedBetaPrime Of(double p, double a, double b)
    {
        if (!(p > 0.0 && double.IsFinite(p))) throw new ArgumentOutOfRangeException(nameof(p), "The power p must be positive and finite.");
        if (!(a > 0.0 && double.IsFinite(a))) throw new ArgumentOutOfRangeException(nameof(a), "Shape a must be positive and finite.");
        if (!(b > 0.0 && double.IsFinite(b))) throw new ArgumentOutOfRangeException(nameof(b), "Shape b must be positive and finite.");
        return new GeneralizedBetaPrime(p, a, b);
    }

    /// <summary>The power p.</summary>
    public double P { get; }

    /// <summary>The shape a, which sets the density's power pa − 1 at 0.</summary>
    public double A { get; }

    /// <summary>The shape b, which sets the tail's power −pb − 1.</summary>
    public double B { get; }

    /// <summary>The mode of u = ln v, ln(a/b)/p.</summary>
    public double ModeInLog => Math.Log(A / B) / P;

    /// <summary>ln of u = ln v's density at <paramref name="u"/>: pa·u − (a + b)·ln(1 + e^(pu)) − ln B(a, b) + ln p.</summary>
    public double LogDensityInLog(double u) => P * A * u - (A + B) * Softplus(P * u) - _logNormalizer;

    /// <summary>ln E[v^<paramref name="k"/>], ln B(a + k/p, b − k/p) − ln B(a, b) for −pa &lt; k &lt; pb; +∞ outside.</summary>
    public double LogMoment(double k)
    {
        double shift = k / P;
        return A + shift > 0.0 && B - shift > 0.0
            ? SpecialFunctions.LogBeta(A + shift, B - shift) - SpecialFunctions.LogBeta(A, B)
            : double.PositiveInfinity;
    }

    /// <summary>E[v], finite for b &gt; 1/p; +∞ otherwise.</summary>
    public double Mean => Math.Exp(LogMoment(1.0));

    /// <summary>Var[v], finite for b &gt; 2/p; +∞ otherwise.</summary>
    public double Variance
    {
        get
        {
            double second = LogMoment(2.0);
            if (double.IsPositiveInfinity(second)) return double.PositiveInfinity;
            double mean = Mean;
            return Math.Exp(second) - mean * mean;
        }
    }

    /// <inheritdoc/>
    public double LogPdf(double x)
    {
        if (!(x > 0.0) || double.IsPositiveInfinity(x)) return double.NegativeInfinity;
        double u = Math.Log(x);
        return LogDensityInLog(u) - u;
    }

    /// <inheritdoc/>
    public double Cdf(double x) => Math.Exp(LogCdf(x));

    /// <inheritdoc/>
    public double Survival(double x) => Math.Exp(LogSurvival(x));

    /// <inheritdoc/>
    public double LogCdf(double x) => LogTails(x).LogCdf;

    /// <inheritdoc/>
    public double LogSurvival(double x) => LogTails(x).LogSurvival;

    // (ln F, ln S) at x through Beta(a, b) at logit y = p·ln x.
    private (double LogCdf, double LogSurvival) LogTails(double x)
    {
        if (double.IsNaN(x)) return (double.NaN, double.NaN);
        if (!(x > 0.0)) return (double.NegativeInfinity, 0.0);
        if (double.IsPositiveInfinity(x)) return (0.0, double.NegativeInfinity);
        (double y, double c, double logY, double logC) = Beta.FromLogit(P * Math.Log(x));
        return SpecialFunctions.LogRegularizedBeta(y, c, logY, logC, A, B);
    }

    /// <inheritdoc/>
    public double IntervalProbability(double lower, double upper) => Math.Exp(LogIntervalProbability(lower, upper));

    /// <inheritdoc/>
    public double LogIntervalProbability(double lower, double upper)
    {
        if (Tails.IsEmpty(lower, upper) || upper <= 0.0) return double.NegativeInfinity;
        (double fl, double sl) = LogTails(lower);
        (double fu, double su) = LogTails(upper);
        return Tails.LogInterval(fl, sl, fu, su);
    }

    /// <summary>One draw, (G_a/G_b)^(1/p) taken in logs, so neither a tiny nor a huge ratio over- or underflows before the power.</summary>
    public double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        double ga = Gamma.SampleCore(rng, A, 1.0), gb = Gamma.SampleCore(rng, B, 1.0);
        return Math.Exp((Math.Log(ga) - Math.Log(gb)) / P);
    }

    /// <summary>One draw on [<paramref name="lower"/>, <paramref name="upper"/>] by inversion: the quantile of a uniform between the ends'
    /// masses, on the side whose tails are small.</summary>
    public double SampleOn(Xoshiro256PlusPlus rng, double lower, double upper)
    {
        ArgumentNullException.ThrowIfNull(rng);
        (lower, upper) = Tails.Clip(lower, upper, 0.0, double.PositiveInfinity);
        (double fl, double sl) = LogTails(lower);
        (double fu, double su) = LogTails(upper);
        double w = rng.NextDouble();
        if (fu <= -SpecialFunctions.Ln2)
        {
            double ratio = Math.Exp(fl - fu);
            return Math.Clamp(InverseLogCdf(fu + Math.Log(ratio + w * (1.0 - ratio))), lower, upper);
        }
        double back = Math.Exp(su - sl);
        return Math.Clamp(InverseLogSurvival(sl + Math.Log(back + w * (1.0 - back))), lower, upper);
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
        return FromTails(logP, SpecialFunctions.Log1mExp(logP));
    }

    /// <inheritdoc/>
    public double InverseLogSurvival(double logS)
    {
        Tails.CheckLogProbability(logS);
        return FromTails(SpecialFunctions.Log1mExp(logS), logS);
    }

    // v at its two log tails: the beta solves in its own logit, p·ln v, so v and 1/v each keep their digits.
    private double FromTails(double logLower, double logUpper) => Math.Exp(Beta.Shapes(A, B).InverseLogit(logLower, logUpper) / P);

    // ln(1 + eˣ) without overflow.
    internal static double Softplus(double x) => x > 0.0 ? x + SpecialFunctions.Log1p(Math.Exp(-x)) : SpecialFunctions.Log1p(Math.Exp(x));
}
