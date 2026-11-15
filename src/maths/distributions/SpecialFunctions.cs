using System;

namespace Maths.Distributions;

/// <summary>
/// Special functions backing the continuous distributions: the log-gamma function, and the regularized incomplete
/// gamma and beta functions in log space, from which the Gamma, inverse-gamma, normal, Student-t and Beta CDFs
/// follow. LogGamma is the Lanczos approximation (g = 7) with the reflection formula for <c>x &lt; 0.5</c>. The
/// incomplete functions use the power series and the Lentz continued fraction (Numerical Recipes §6.2, §6.4),
/// each evaluated for the tail it represents directly, so a tail probability far below 1e-308 keeps its digits as
/// a log; the other tail is its accurate complement.
/// </summary>
public static partial class SpecialFunctions
{
    private static readonly double[] LanczosG7 =
    {
        0.99999999999980993, 676.5203681218851, -1259.1392167224028,
        771.32342877765313, -176.61502916214059, 12.507343278686905,
        -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7,
    };

    private const double Eps = 2.220446049250313e-16;
    private const double Tiny = 1e-300;
    private const int MaxIterations = 100_000;
    internal const double Ln2 = 0.69314718055994530942;

    /// <summary>Natural log of the gamma function, <c>ln Γ(x)</c>, for <c>x &gt; 0</c>.</summary>
    public static double LogGamma(double x)
    {
        if (x < 0.5)
        {
            // Below about 5.6·10⁻³⁰⁹, π/sin(πx) ≈ 1/x overflows; there ln Γ(x) = −ln x − γx + O(x²) is −ln x to the last bit.
            double reflection = Math.PI / Math.Sin(Math.PI * x);
            if (double.IsPositiveInfinity(reflection) && x > 0.0) return -Math.Log(x);
            return Math.Log(reflection) - LogGamma(1.0 - x);
        }

        x -= 1.0;
        double a = LanczosG7[0];
        double t = x + 7.5;
        for (int i = 1; i < LanczosG7.Length; i++)
            a += LanczosG7[i] / (x + i);

        return 0.5 * Math.Log(2.0 * Math.PI) + (x + 0.5) * Math.Log(t) - t + Math.Log(a);
    }

    /// <summary>
    /// The log of the multivariate gamma function, <c>ln Γ_k(a) = k(k − 1)/4 · ln π + Σⱼ ln Γ(a − j/2)</c> over j = 0 … k − 1, for
    /// <c>a &gt; (k − 1)/2</c>: the inverse Wishart's normalizer. Γ₁ is Γ.
    /// </summary>
    public static double LogMultivariateGamma(int k, double a)
    {
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k), "The dimension is at least 1.");
        if (!(a > 0.5 * (k - 1))) throw new ArgumentOutOfRangeException(nameof(a), "Γ_k(a) needs a > (k − 1)/2.");
        double sum = 0.25 * k * (k - 1) * Math.Log(Math.PI);
        for (int j = 0; j < k; j++) sum += LogGamma(a - 0.5 * j);
        return sum;
    }

    /// <summary>
    /// The digamma function ψ(x) = d ln Γ(x)/dx. For x &gt; 0 it climbs by ψ(x) = ψ(x + 1) − 1/x to x ≥ 15, where the
    /// asymptotic series ln x − 1/(2x) − Σ B₂ₖ/(2k x²ᵏ), to k = 5, is within 2·10⁻¹⁶ of ψ, so ψ keeps about 10⁻¹⁵ absolute
    /// accuracy, the recurrence's rounding. Below zero it reflects, ψ(x) = ψ(1 − x) − π/tan(πx). NaN at the poles 0, −1,
    /// −2, ….
    /// </summary>
    public static double Digamma(double x)
    {
        if (double.IsNaN(x) || double.IsNegativeInfinity(x)) return double.NaN;
        if (double.IsPositiveInfinity(x)) return double.PositiveInfinity;
        if (x <= 0.0)
            return x == Math.Floor(x) ? double.NaN : Digamma(1.0 - x) - Math.PI / Math.Tan(Math.PI * x);

        double result = 0.0;
        while (x < 15.0)
        {
            result -= 1.0 / x;
            x += 1.0;
        }
        double f = 1.0 / (x * x);
        return result + Math.Log(x) - 0.5 / x
               - f * (1.0 / 12 - f * (1.0 / 120 - f * (1.0 / 252 - f * (1.0 / 240 - f * (1.0 / 132)))));
    }

    /// <summary>
    /// <c>ln B(a, b) = ln Γ(a) + ln Γ(b) − ln Γ(a + b)</c>, without the cancellation of those three terms. With s the smaller shape and
    /// l the larger, it is ln Γ(s) − <see cref="LogGammaRatio"/>(l, s), which keeps its digits however far apart the shapes are (the
    /// three-term form returned 0 for ln B(10¹⁶, 1), whose value is −36.84). Where both are at least 10 it is R's <c>lbeta</c>
    /// arrangement, −½·ln l + ½·ln 2π + μ(s) + μ(l) − μ(s + l) + (s − ½)·ln(s/(s + l)) + l·ln(1 − s/(s + l)), with μ Stirling's
    /// correction: no two terms of size s·ln s cancel there, where the ratio form would lose up to three digits (1.2·10⁻¹³ near
    /// 10²⁰⁰). Against ln B to 400 digits at 2428 points from 10⁻¹² to 10³⁰⁰ it is within 1.7·10⁻¹⁵ of max(1, |ln B|). Shapes that are
    /// not finite and positive take the three-term form, as before.
    /// </summary>
    public static double LogBeta(double a, double b)
    {
        if (!(a > 0.0 && b > 0.0 && double.IsFinite(a) && double.IsFinite(b))) return LogGamma(a) + LogGamma(b) - LogGamma(a + b);
        double s = Math.Min(a, b), l = Math.Max(a, b);
        if (s < StirlingBetaFloor) return LogGamma(s) - LogGammaRatio(l, s);
        // s + l overflows only where both are near the double's end; there the weight is taken from l/s.
        double total = s + l, w = double.IsPositiveInfinity(total) ? 1.0 / (1.0 + l / s) : s / total;
        return -0.5 * Math.Log(l) + 0.5 * Math.Log(2.0 * Math.PI) + StirlingCorrection(s) + StirlingCorrection(l) - StirlingCorrection(total)
               + (s - 0.5) * Math.Log(w) + l * Log1p(-w);
    }

    // At and above this smaller shape LogBeta takes R's lbeta arrangement; there Stirling's correction to eight terms is within 10⁻¹⁷.
    private const double StirlingBetaFloor = 10.0;

    // μ(x) = ln Γ(x) − (x − ½)·ln x + x − ½·ln 2π, x ≥ 10: Stirling's series Σ B₂ₖ/(2k(2k − 1)·x^(2k−1)) to eight terms; 0 at +∞.
    private static double StirlingCorrection(double x)
    {
        double f = 1.0 / (x * x);
        return (1.0 / 12 - f * (1.0 / 360 - f * (1.0 / 1260 - f * (1.0 / 1680 - f * (1.0 / 1188 - f * (691.0 / 360360 - f * (1.0 / 156
                - f * 3617.0 / 122400))))))) / x;
    }

    /// <summary>Regularized lower incomplete gamma <c>P(a, x) = γ(a, x)/Γ(a)</c>, the Gamma(a, rate 1) CDF.</summary>
    public static double RegularizedGammaP(double a, double x) => Math.Exp(LogRegularizedGamma(a, x).LogP);

    /// <summary>Regularized upper incomplete gamma <c>Q(a, x) = 1 − P(a, x)</c>, computed directly in its tail.</summary>
    public static double RegularizedGammaQ(double a, double x) => Math.Exp(LogRegularizedGamma(a, x).LogQ);

    /// <summary><c>ln P(a, x)</c>.</summary>
    public static double LogRegularizedGammaP(double a, double x) => LogRegularizedGamma(a, x).LogP;

    /// <summary><c>ln Q(a, x)</c>; finite where <c>Q</c> itself underflows.</summary>
    public static double LogRegularizedGammaQ(double a, double x) => LogRegularizedGamma(a, x).LogQ;

    /// <summary>Regularized incomplete beta <c>I_x(a, b)</c>, the Beta(a, b) CDF.</summary>
    public static double RegularizedBeta(double x, double a, double b) => Math.Exp(LogRegularizedBetaPair(x, a, b).LogI);

    /// <summary><c>ln I_x(a, b)</c>.</summary>
    public static double LogRegularizedBeta(double x, double a, double b) => LogRegularizedBetaPair(x, a, b).LogI;

    /// <summary><c>ln(1 − I_x(a, b)) = ln I_{1−x}(b, a)</c>, computed directly where it is the small tail.</summary>
    public static double LogRegularizedBetaComplement(double x, double a, double b)
        => LogRegularizedBetaPair(x, a, b).LogComplement;

    // (ln P, ln Q). The series represents P (x < a + 1), the continued fraction Q (x ≥ a + 1).
    internal static (double LogP, double LogQ) LogRegularizedGamma(double a, double x)
    {
        if (!(a > 0.0)) throw new ArgumentOutOfRangeException(nameof(a), "Shape must be positive.");
        if (!(x >= 0.0)) throw new ArgumentOutOfRangeException(nameof(x), "x must be non-negative.");
        if (x == 0.0) return (double.NegativeInfinity, 0.0);
        if (double.IsPositiveInfinity(x)) return (0.0, double.NegativeInfinity);

        double logPrefix = a * Math.Log(x) - x - LogGamma(a);   // ln(xᵃ e⁻ˣ / Γ(a))
        if (x < a + 1.0)
        {
            double ap = a, del = 1.0 / a, sum = del;
            for (int n = 0; n < MaxIterations; n++)
            {
                ap += 1.0;
                del *= x / ap;
                sum += del;
                if (Math.Abs(del) < Math.Abs(sum) * Eps)
                {
                    double logP = logPrefix + Math.Log(sum);
                    return (logP, Log1mExp(logP));
                }
            }
            throw new ArithmeticException($"Incomplete gamma series did not converge for a = {a}, x = {x}.");
        }

        double b = x + 1.0 - a, c = 1.0 / Tiny, d = 1.0 / b, h = d;
        for (int i = 1; i <= MaxIterations; i++)
        {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (Math.Abs(d) < Tiny) d = Tiny;
            c = b + an / c;
            if (Math.Abs(c) < Tiny) c = Tiny;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (Math.Abs(del - 1.0) <= Eps)
            {
                double logQ = logPrefix + Math.Log(h);
                return (Log1mExp(logQ), logQ);
            }
        }
        throw new ArithmeticException($"Incomplete gamma continued fraction did not converge for a = {a}, x = {x}.");
    }

    private static (double LogI, double LogComplement) LogRegularizedBetaPair(double x, double a, double b)
    {
        if (!(x >= 0.0 && x <= 1.0)) throw new ArgumentOutOfRangeException(nameof(x), "x must be in [0, 1].");
        double y = 1.0 - x;
        return LogRegularizedBeta(x, y, Math.Log(x), Math.Log(y), a, b);
    }

    // (ln I_x(a, b), ln I_y(b, a)) with y = 1 − x supplied separately, so a caller that knows both x and 1 − x
    // without cancellation (the Student-t CDF) keeps them; the logs are supplied too, for x below double range.
    internal static (double LogI, double LogComplement) LogRegularizedBeta(
        double x, double y, double logX, double logY, double a, double b)
    {
        if (!(a > 0.0)) throw new ArgumentOutOfRangeException(nameof(a), "Shape a must be positive.");
        if (!(b > 0.0)) throw new ArgumentOutOfRangeException(nameof(b), "Shape b must be positive.");
        if (double.IsNegativeInfinity(logX)) return (double.NegativeInfinity, 0.0);
        if (double.IsNegativeInfinity(logY)) return (0.0, double.NegativeInfinity);

        double logFront = a * logX + b * logY - LogBeta(a, b);
        if (x < (a + 1.0) / (a + b + 2.0))
        {
            double logI = logFront + Math.Log(BetaContinuedFraction(x, a, b)) - Math.Log(a);
            return (logI, Log1mExp(logI));
        }
        double logC = logFront + Math.Log(BetaContinuedFraction(y, b, a)) - Math.Log(b);
        return (Log1mExp(logC), logC);
    }

    private static double BetaContinuedFraction(double x, double a, double b)
    {
        double qab = a + b, qap = a + 1.0, qam = a - 1.0;
        double c = 1.0, d = 1.0 - qab * x / qap;
        if (Math.Abs(d) < Tiny) d = Tiny;
        d = 1.0 / d;
        double h = d;
        for (int m = 1; m <= MaxIterations; m++)
        {
            int m2 = 2 * m;
            double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = 1.0 + aa * d;
            if (Math.Abs(d) < Tiny) d = Tiny;
            c = 1.0 + aa / c;
            if (Math.Abs(c) < Tiny) c = Tiny;
            d = 1.0 / d;
            h *= d * c;
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
            d = 1.0 + aa * d;
            if (Math.Abs(d) < Tiny) d = Tiny;
            c = 1.0 + aa / c;
            if (Math.Abs(c) < Tiny) c = Tiny;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (Math.Abs(del - 1.0) <= Eps) return h;
        }
        throw new ArithmeticException($"Incomplete beta continued fraction did not converge for a = {a}, b = {b}, x = {x}.");
    }

    /// <summary><c>ln(1 + x)</c>, accurate for small <c>x</c> (Goldberg's correction). Past x ≈ 2.5·10³⁰⁵ the correction's product
    /// ln(u)·x overflows; there u is x, the correction is 1, and it returns ln u.</summary>
    public static double Log1p(double x)
    {
        if (double.IsPositiveInfinity(x)) return x;
        double u = 1.0 + x;
        if (u == 1.0) return x;
        double corrected = Math.Log(u) * x / (u - 1.0);
        return double.IsPositiveInfinity(corrected) ? Math.Log(u) : corrected;
    }

    /// <summary><c>eˣ − 1</c>, accurate for small <c>x</c> (Kahan's correction).</summary>
    internal static double Expm1(double x)
    {
        double u = Math.Exp(x);
        if (u == 1.0) return x;
        double um1 = u - 1.0;
        if (um1 == -1.0) return -1.0;
        if (double.IsPositiveInfinity(u)) return u;
        return um1 * x / Math.Log(u);
    }

    /// <summary><c>ln(1 − eᵘ)</c> for <c>u ≤ 0</c> (Mächler 2012): the log of a complement.</summary>
    internal static double Log1mExp(double u)
        => u > -Ln2 ? Math.Log(-Expm1(u)) : Log1p(-Math.Exp(u));
}
