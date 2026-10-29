using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The location-scale Student-t distribution: <c>X = location + scale·T</c> with <c>T</c> a standard t on
/// <see cref="DegreesOfFreedom"/> ν. Built with <see cref="Standard"/> or <see cref="LocationScale"/>; the scale
/// is the t's scale, not its standard deviation (which is <c>scale·√(ν/(ν−2))</c>, finite only for ν &gt; 2).
/// The tails are <c>P(T &gt; t) = ½·I_x(ν/2, ½)</c> with <c>x = ν/(ν + t²)</c>, evaluated from <c>x</c> and
/// <c>1 − x</c> each formed without cancellation and in log space, so <c>t = 10²⁰⁰</c> under ν = 1 keeps its
/// digits. Accuracy degrades for very large ν, where <c>ln B(ν/2, ½)</c> cancels; use <see cref="Normal"/> there.
/// <para>Interval moments are by <see cref="GaussLegendre"/> at every ν, ν = 1 included, over a finite interval: panels
/// split at the location and graded geometrically away from it, so each panel stays three half-widths from the density's
/// poles at location ± i·scale·√ν, with 20 nodes per panel, and the error estimated by the 40-node rule's difference.</para>
/// </summary>
public readonly struct StudentT : IContinuousDistribution, IIntervalMoments
{
    private const int QuadratureNodes = 20;

    private readonly double _logNormalizer;   // ln Γ((ν+1)/2) − ln Γ(ν/2) − ½ ln(νπ) − ln scale

    /// <summary>Degrees of freedom ν.</summary>
    public double DegreesOfFreedom { get; }

    /// <summary>Location (the median; the mean when ν &gt; 1).</summary>
    public double Location { get; }

    /// <summary>Scale of the t (not its standard deviation).</summary>
    public double Scale { get; }

    private StudentT(double nu, double location, double scale)
    {
        DegreesOfFreedom = nu;
        Location = location;
        Scale = scale;
        _logNormalizer = SpecialFunctions.LogGamma(0.5 * (nu + 1.0)) - SpecialFunctions.LogGamma(0.5 * nu)
                       - 0.5 * Math.Log(nu * Math.PI) - Math.Log(scale);
    }

    /// <summary>The standard t on ν degrees of freedom.</summary>
    public static StudentT Standard(double degreesOfFreedom) => LocationScale(degreesOfFreedom, 0.0, 1.0);

    /// <summary><c>location + scale·T_ν</c>.</summary>
    public static StudentT LocationScale(double degreesOfFreedom, double location, double scale)
    {
        if (!(degreesOfFreedom > 0.0) || double.IsPositiveInfinity(degreesOfFreedom))
            throw new ArgumentOutOfRangeException(nameof(degreesOfFreedom), "Degrees of freedom must be positive and finite.");
        if (!double.IsFinite(location)) throw new ArgumentOutOfRangeException(nameof(location), "Location must be finite.");
        if (!(scale > 0.0) || double.IsPositiveInfinity(scale))
            throw new ArgumentOutOfRangeException(nameof(scale), "Scale must be positive and finite.");
        return new StudentT(degreesOfFreedom, location, scale);
    }

    /// <summary>The mean, <see cref="Location"/> for ν &gt; 1; NaN (undefined) otherwise.</summary>
    public double Mean => DegreesOfFreedom > 1.0 ? Location : double.NaN;

    /// <summary><c>scale²·ν/(ν − 2)</c> for ν &gt; 2; +∞ for 1 &lt; ν ≤ 2; NaN (undefined) for ν ≤ 1.</summary>
    public double Variance => DegreesOfFreedom > 2.0 ? Scale * Scale * DegreesOfFreedom / (DegreesOfFreedom - 2.0)
                            : DegreesOfFreedom > 1.0 ? double.PositiveInfinity : double.NaN;

    /// <summary>One draw: <c>location + scale·Z·√(ν/V)</c> with <c>Z ~ N(0,1)</c> and <c>V ~ χ²_ν</c>.</summary>
    public double Sample(Xoshiro256PlusPlus rng)
    {
        double z = Normal.SampleStandard(rng);
        double v = 2.0 * Gamma.SampleStandard(rng, 0.5 * DegreesOfFreedom);
        return Location + Scale * z * Math.Sqrt(DegreesOfFreedom / v);
    }

    public double LogPdf(double x)
    {
        double t = (x - Location) / Scale;
        return _logNormalizer - 0.5 * (DegreesOfFreedom + 1.0) * LogOnePlusTSquaredOverNu(t);
    }

    public double Cdf(double x) => Math.Exp(LogCdf(x));
    public double Survival(double x) => Math.Exp(LogSurvival(x));

    public double LogCdf(double x)
    {
        double t = (x - Location) / Scale;
        return t < 0.0 ? TailLogs(-t).LogUpper : SpecialFunctions.Log1mExp(TailLogs(t).LogUpper);
    }

    public double LogSurvival(double x)
    {
        double t = (x - Location) / Scale;
        return t > 0.0 ? TailLogs(t).LogUpper : SpecialFunctions.Log1mExp(TailLogs(-t).LogUpper);
    }

    public double IntervalProbability(double lower, double upper) => Math.Exp(LogIntervalProbability(lower, upper));

    public double LogIntervalProbability(double lower, double upper)
    {
        if (Tails.IsEmpty(lower, upper)) return double.NegativeInfinity;
        double a = (lower - Location) / Scale, b = (upper - Location) / Scale;
        if (a >= 0.0) return SameTail(TailLogs(a).LogUpper, TailLogs(b).LogUpper);     // S(a) − S(b)
        if (b <= 0.0) return SameTail(TailLogs(-b).LogUpper, TailLogs(-a).LogUpper);   // F(b) − F(a)
        return Math.Log(0.5 * (Math.Exp(TailLogs(-a).LogCentral) + Math.Exp(TailLogs(b).LogCentral)));
    }

    /// <inheritdoc/>
    /// <remarks>The interval must be finite: truncate the distribution to one. The moments are the 40-node rule's.</remarks>
    public double IntervalMoments(double lower, double upper, double centre, double logScale, Span<double> moments)
    {
        moments.Clear();
        if (Tails.IsEmpty(lower, upper) || moments.IsEmpty) return 0.0;
        if (!double.IsFinite(lower) || !double.IsFinite(upper))
            throw new ArgumentException($"Student-t interval moments integrate over a finite interval, got [{lower}, {upper}]; truncate the distribution.");

        // Breakpoints in t: the interval's ends, 0, and ±g·2ᵏ with g = ½·min(√ν, 1), so a panel [g·2ᵏ, g·2ᵏ⁺¹] is three
        // half-widths from the poles at ±i√ν and the central panels four.
        double a = (lower - Location) / Scale, b = (upper - Location) / Scale;
        double g = 0.5 * Math.Min(Math.Sqrt(DegreesOfFreedom), 1.0), far = Math.Max(Math.Abs(a), Math.Abs(b));
        var breaks = new System.Collections.Generic.List<double> { a, b };
        if (a < 0.0 && b > 0.0) breaks.Add(0.0);
        for (double r = g; r < far; r *= 2.0)
        {
            if (r > a && r < b) breaks.Add(r);
            if (-r > a && -r < b) breaks.Add(-r);
        }
        breaks.Sort();

        int n = moments.Length;
        Span<double> coarse = n <= 32 ? stackalloc double[n] : new double[n];
        coarse.Clear();
        for (int p = 0; p + 1 < breaks.Count; p++)
        {
            double x0 = Location + Scale * breaks[p], x1 = Location + Scale * breaks[p + 1];
            Panel(GaussLegendre.Of(QuadratureNodes), x0, x1, centre, logScale, coarse);
            Panel(GaussLegendre.Of(2 * QuadratureNodes), x0, x1, centre, logScale, moments);
        }
        double error = 0.0;
        for (int k = 0; k < n; k++) error = Math.Max(error, Math.Abs(moments[k] - coarse[k]));
        return error;
    }

    // Adds the rule's ∫ (x − c)ᵏ f(x)·e^(−logScale) dx over [x0, x1] into sums.
    private void Panel(GaussLegendre rule, double x0, double x1, double centre, double logScale, Span<double> sums)
    {
        double mid = 0.5 * (x0 + x1), half = 0.5 * (x1 - x0);
        ReadOnlySpan<double> nodes = rule.Nodes, weights = rule.Weights;
        for (int i = 0; i < nodes.Length; i++)
        {
            double x = mid + half * nodes[i];
            double w = half * weights[i] * Math.Exp(LogPdf(x) - logScale), d = x - centre;
            for (int k = 0; k < sums.Length; k++)
            {
                sums[k] += w;
                w *= d;
            }
        }
    }

    // For t ≥ 0: ln P(T > t) = ln ½ + ln I_x(ν/2, ½), and ln P(|T| < t) = ln I_y(½, ν/2), with x = ν/(ν + t²) and
    // y = t²/(ν + t²). Both are formed from whichever of t²/ν and ν/t² is below 1, so neither overflows.
    private (double LogUpper, double LogCentral) TailLogs(double t)
    {
        double nu = DegreesOfFreedom;
        if (double.IsPositiveInfinity(t)) return (double.NegativeInfinity, 0.0);
        if (t == 0.0) return (-SpecialFunctions.Ln2, double.NegativeInfinity);

        double x, y, logX, logY;
        if (t * t < nu)
        {
            double q = t * t / nu;
            double l1q = SpecialFunctions.Log1p(q);
            x = 1.0 / (1.0 + q);
            y = q / (1.0 + q);
            logX = -l1q;
            logY = Math.Log(q) - l1q;
        }
        else
        {
            double r = nu / t / t;
            double l1r = SpecialFunctions.Log1p(r);
            x = r / (1.0 + r);
            y = 1.0 / (1.0 + r);
            logX = Math.Log(nu) - 2.0 * Math.Log(t) - l1r;
            logY = -l1r;
        }
        (double logI, double logComplement) = SpecialFunctions.LogRegularizedBeta(x, y, logX, logY, 0.5 * nu, 0.5);
        return (logI - SpecialFunctions.Ln2, logComplement);
    }

    // ln(1 + t²/ν) without overflowing t².
    private double LogOnePlusTSquaredOverNu(double t)
    {
        double at = Math.Abs(t), nu = DegreesOfFreedom;
        return at * at < nu
            ? SpecialFunctions.Log1p(at * at / nu)
            : 2.0 * Math.Log(at) - Math.Log(nu) + SpecialFunctions.Log1p(nu / at / at);
    }

    private static double SameTail(double logHi, double logLo)
        => double.IsNegativeInfinity(logHi) ? double.NegativeInfinity
         : logHi + SpecialFunctions.Log1mExp(logLo - logHi);
}
