using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The normal distribution N(mean, sd²), built with <see cref="MeanSd"/>, <see cref="MeanVariance"/> or
/// <see cref="Standard"/>. <see cref="SampleStandard"/> is the standard-normal draw that also backs the
/// Marsaglia–Tsang <see cref="Gamma"/> sampler. The tails come from <c>Φ(−z) = ½·Q(½, z²/2)</c> in log space, so
/// <c>ln Φ(−40)</c> is exact where <c>Φ(−40)</c> underflows, and an interval straddling the mean is the sum of
/// two error functions rather than a difference of CDFs near ½. A readonly value type: zero-alloc in inner loops.
/// <para>Interval moments are analytic at every power: the density's Taylor series about the centre of each short
/// stretch of the interval, φ(z + u) = φ(z)·Σₙ (−1)ⁿHeₙ(z)uⁿ/n!, integrated term by term and summed to convergence. A
/// stretch's half-width h keeps h·(1 + |z|) ≤ 1, so the series never cancels by more than e³, where the moments' recurrence
/// in Φ and φ loses every digit on a narrow interval.</para>
/// </summary>
public readonly struct Normal : IContinuousDistribution, IIntervalMoments
{
    private static readonly double LogSqrt2Pi = 0.5 * Math.Log(2.0 * Math.PI);

    // Where the density falls below e^−Reach of its largest value on an interval, the rest of the interval is left out
    // of the moments.
    private const double Reach = 92.0;
    private const int MaxSeriesTerms = 400;

    /// <summary>The mean.</summary>
    public double Mean { get; }

    /// <summary>The standard deviation.</summary>
    public double Sd { get; }

    private Normal(double mean, double sd)
    {
        Mean = mean;
        Sd = sd;
    }

    /// <summary>N(0, 1).</summary>
    public static Normal Standard => new(0.0, 1.0);

    /// <summary>N(mean, sd²) from the standard deviation.</summary>
    public static Normal MeanSd(double mean, double sd)
    {
        if (!double.IsFinite(mean)) throw new ArgumentOutOfRangeException(nameof(mean), "Mean must be finite.");
        if (!(sd > 0.0) || double.IsPositiveInfinity(sd))
            throw new ArgumentOutOfRangeException(nameof(sd), "Standard deviation must be positive and finite.");
        return new Normal(mean, sd);
    }

    /// <summary>N(mean, variance) from the variance.</summary>
    public static Normal MeanVariance(double mean, double variance)
    {
        if (!(variance > 0.0) || double.IsPositiveInfinity(variance))
            throw new ArgumentOutOfRangeException(nameof(variance), "Variance must be positive and finite.");
        return MeanSd(mean, Math.Sqrt(variance));
    }

    /// <summary>The variance <c>sd²</c>.</summary>
    public double Variance => Sd * Sd;

    /// <summary>One standard-normal draw via the Box–Muller transform (cosine branch; two uniforms).</summary>
    public static double SampleStandard(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        double u1 = 1.0 - rng.NextDouble();   // (0,1] avoids log(0)
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    /// <summary>One draw from this distribution.</summary>
    public double Sample(Xoshiro256PlusPlus rng) => Mean + Sd * SampleStandard(rng);

    public double LogPdf(double x)
    {
        double z = (x - Mean) / Sd;
        return -0.5 * z * z - Math.Log(Sd) - LogSqrt2Pi;
    }

    public double Cdf(double x) => Math.Exp(LogCdf(x));
    public double Survival(double x) => Math.Exp(LogSurvival(x));
    public double LogCdf(double x) => StandardLogCdf((x - Mean) / Sd);
    public double LogSurvival(double x) => StandardLogCdf((Mean - x) / Sd);
    public double IntervalProbability(double lower, double upper) => Math.Exp(LogIntervalProbability(lower, upper));

    public double LogIntervalProbability(double lower, double upper)
    {
        if (Tails.IsEmpty(lower, upper)) return double.NegativeInfinity;
        double a = (lower - Mean) / Sd, b = (upper - Mean) / Sd;
        if (a >= 0.0) return SameTail(LogLowerTail(a), LogLowerTail(b));      // Φ(−a) − Φ(−b)
        if (b <= 0.0) return SameTail(LogLowerTail(-b), LogLowerTail(-a));    // Φ(b) − Φ(a)
        return Math.Log(0.5 * (Erf(-a) + Erf(b)));                             // straddles the mean
    }

    /// <inheritdoc/>
    public double IntervalMoments(double lower, double upper, double centre, double logScale, Span<double> moments)
    {
        moments.Clear();
        if (Tails.IsEmpty(lower, upper) || moments.IsEmpty) return 0.0;
        double a = (lower - Mean) / Sd, b = (upper - Mean) / Sd;
        double near = a > 0.0 ? a : b < 0.0 ? b : 0.0;
        double reach = Math.Sqrt(near * near + 2.0 * Reach);
        a = Math.Max(a, -reach);
        b = Math.Min(b, reach);
        if (!(b > a)) return 0.0;

        int n = moments.Length;
        Span<double> piece = n <= 32 ? stackalloc double[n] : new double[n];
        Span<double> powers = n <= 32 ? stackalloc double[n] : new double[n];
        double zc = (centre - Mean) / Sd;
        for (double z0 = a; z0 < b;)
        {
            double h = Math.Min(1.0, 1.0 / (1.0 + Math.Max(Math.Abs(z0), Math.Abs(z0 + 2.0))));
            double z1 = z0 + 2.0 * h;
            if (!(z1 > z0)) throw new ArithmeticException($"The interval [{lower}, {upper}] lies too far in the tail of N({Mean}, {Sd}²) for its moments.");
            if (z1 > b) z1 = b;
            double zm = 0.5 * (z0 + z1), half = 0.5 * (z1 - z0);
            double weight = Math.Exp(-0.5 * zm * zm - LogSqrt2Pi - logScale);
            z0 = z1;
            if (weight == 0.0) continue;
            SeriesPiece(zm, half, piece);
            // ∫ (z − zc)ᵏ φ over the stretch = Σₘ C(k, m)·dᵏ⁻ᵐ·∫ (z − zm)ᵐ φ, with d = zm − zc.
            double d = zm - zc;
            powers[0] = 1.0;
            for (int k = 1; k < n; k++) powers[k] = powers[k - 1] * d;
            for (int k = 0; k < n; k++)
            {
                double sum = 0.0, binomial = 1.0;
                for (int m = 0; m <= k; m++)
                {
                    sum += binomial * powers[k - m] * piece[m];
                    binomial = binomial * (k - m) / (m + 1);
                }
                moments[k] += weight * sum;
            }
        }
        double scale = 1.0;
        for (int k = 0; k < n; k++)
        {
            moments[k] *= scale;
            scale *= Sd;
        }
        return 0.0;
    }

    // piece[m] = ∫ uᵐ φ(zm + u)/φ(zm) du over [−h, h] = (−1)ᵐ hᵐ⁺¹ Σⱼ Hⱼ·2/(m + j + 1), over j of m's parity, with
    // Hⱼ = Heⱼ(zm)hʲ/j! by Heⱼ₊₁ = z·Heⱼ − j·Heⱼ₋₁.
    private static void SeriesPiece(double zm, double h, Span<double> piece)
    {
        int n = piece.Length;
        piece.Clear();
        double previous = 0.0, current = 1.0;
        for (int j = 0; ; j++)
        {
            for (int m = j & 1; m < n; m += 2) piece[m] += current * 2.0 / (m + j + 1);
            double next = (zm * h * current - h * h * previous) / (j + 1);
            previous = current;
            current = next;
            if (j >= n && Math.Abs(previous) + Math.Abs(current) <= 1e-20) break;
            if (j >= MaxSeriesTerms) throw new ArithmeticException($"The normal moment series did not converge at z = {zm}, h = {h}.");
        }
        double power = h;
        for (int m = 0; m < n; m++)
        {
            piece[m] *= (m & 1) == 0 ? power : -power;
            power *= h;
        }
    }

    /// <summary>The quantile function: the x with <c>P(X ≤ x) = p</c>; −∞ at 0, +∞ at 1, NaN outside [0, 1].</summary>
    public double Quantile(double p) => Mean + Sd * StandardQuantile(p);

    /// <summary>
    /// The standard normal quantile Φ⁻¹(p), by Wichura's algorithm AS 241 (PPND16, 1988): rational approximations
    /// in three regions, accurate to about 1 part in 10¹⁶. −∞ at 0, +∞ at 1, NaN outside [0, 1].
    /// </summary>
    public static double StandardQuantile(double p)
    {
        if (!(p >= 0.0 && p <= 1.0)) return double.NaN;
        if (p == 0.0) return double.NegativeInfinity;
        if (p == 1.0) return double.PositiveInfinity;

        double q = p - 0.5;
        if (Math.Abs(q) <= 0.425)
        {
            double r = 0.180625 - q * q;
            return q * (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r + 67265.770927008700853) * r
                              + 45921.953931549871457) * r + 13731.693765509461125) * r + 1971.5909503065514427) * r
                              + 133.14166789178437745) * r + 3.387132872796366608)
                     / (((((((r * 5226.495278852545925 + 28729.085735721942674) * r + 39307.89580009271061) * r
                              + 21213.794301586595867) * r + 5394.1960214247511077) * r + 687.1870074920579083) * r
                              + 42.313330701600911252) * r + 1.0);
        }

        double t = Math.Sqrt(-Math.Log(q < 0.0 ? p : 1.0 - p));
        double value;
        if (t <= 5.0)
        {
            t -= 1.6;
            value = (((((((t * 7.7454501427834140764e-4 + 0.0227238449892691845833) * t + 0.24178072517745061177) * t
                          + 1.27045825245236838258) * t + 3.64784832476320460504) * t + 5.7694972214606914055) * t
                          + 4.6303378461565452959) * t + 1.42343711074968357734)
                  / (((((((t * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * t + 0.0151986665636164571966) * t
                          + 0.14810397642748007459) * t + 0.68976733498510000455) * t + 1.6763848301838038494) * t
                          + 2.05319162663775882187) * t + 1.0);
        }
        else
        {
            t -= 5.0;
            value = (((((((t * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * t + 0.0012426609473880784386) * t
                          + 0.026532189526576123093) * t + 0.29656057182850489123) * t + 1.7848265399172913358) * t
                          + 5.4637849111641143699) * t + 6.6579046435011037772)
                  / (((((((t * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * t + 1.8463183175100546818e-5) * t
                          + 7.868691311456132591e-4) * t + 0.0148753612908506148525) * t + 0.13692988092273580531) * t
                          + 0.59983220655588793769) * t + 1.0);
        }
        return q < 0.0 ? -value : value;
    }

    /// <summary><c>ln Φ(z)</c> of the standard normal.</summary>
    public static double StandardLogCdf(double z)
        => z < 0.0 ? LogLowerTail(-z) : SpecialFunctions.Log1mExp(LogLowerTail(z));

    // ln Φ(−t) for t ≥ 0: Φ(−t) = ½·erfc(t/√2) = ½·Q(½, t²/2).
    private static double LogLowerTail(double t)
        => double.IsPositiveInfinity(t) ? double.NegativeInfinity
         : SpecialFunctions.LogRegularizedGamma(0.5, 0.5 * t * t).LogQ - SpecialFunctions.Ln2;

    // erf(t/√2) = P(½, t²/2) for t ≥ 0, accurate near 0 where 1 − 2Φ(−t) is not.
    private static double Erf(double t)
        => double.IsPositiveInfinity(t) ? 1.0 : Math.Exp(SpecialFunctions.LogRegularizedGamma(0.5, 0.5 * t * t).LogP);

    // ln(e^hi − e^lo) for the two log tail masses hi ≥ lo.
    private static double SameTail(double logHi, double logLo)
        => double.IsNegativeInfinity(logHi) ? double.NegativeInfinity
         : logHi + SpecialFunctions.Log1mExp(logLo - logHi);
}
