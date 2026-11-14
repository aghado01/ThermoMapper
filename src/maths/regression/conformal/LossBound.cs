using System;
using System.Collections.Generic;
using System.Numerics;

namespace Maths.Regression.Conformal;

/// <summary>
/// The bound a calibration set's losses imply on the risk (the conformal layer's slice 18): with n losses ℓᵢ in [0, B], sorted, and
/// ℓ₍ₙ₊₁₎ = B, the random variable L⁺ = Σ Uᵢ ℓ₍ᵢ₎ with U ~ Dir(1, …, 1) (Snell and Griffiths, arXiv 2502.13228).
/// <list type="bullet">
/// <item><see cref="Gaffke(double)"/>, b*_β, its β-quantile, is Gaffke's (2005) upper confidence bound on the risk J = E ℓ: over replicated
/// calibration sets of independent losses at most B, P(J ≤ b*_β) ≥ β (Vlassis and Thomas, arXiv 2607.08415, Theorem 1), with neither
/// identical laws nor continuity needed. At one loss it is ℓ + β(B − ℓ) (<c>gaffke_one</c>). It is a confidence bound, not a posterior:
/// Snell and Griffiths's claim that it bounds the posterior under every prior (their Theorem 4.3) fails as stated (the Lean's finding).</item>
/// <item><see cref="Mean"/>, E L⁺ = (Σ ℓ + B)/(n + 1), is conformal risk control's reading (<c>crc_orderings</c>).</item>
/// </list>
/// P(L⁺ &gt; x) is the divided difference of (t − x)₊ⁿ on the n + 1 knots, confluent where they tie (the B-spline's tail), read here in exact
/// rational arithmetic on the losses' binary values, so ties and discrete losses read exactly: at 0/1 losses with B = 1 the law is
/// Beta(k + 1, n − k), k the losses at 1. Each reading costs O(n²) big-number steps; the bound bisects on it.
/// </summary>
public sealed class LossBound
{
    private readonly double[] _knots;   // the losses sorted, then B
    private readonly double _bound;

    /// <summary>The law of L⁺ for <paramref name="losses"/>, each in [0, <paramref name="bound"/>].</summary>
    public LossBound(IReadOnlyList<double> losses, double bound)
    {
        ArgumentNullException.ThrowIfNull(losses);
        if (!(bound >= 0.0 && double.IsFinite(bound))) throw new ArgumentOutOfRangeException(nameof(bound), "B is finite and nonnegative.");
        _knots = new double[losses.Count + 1];
        for (int i = 0; i < losses.Count; i++)
        {
            double l = losses[i];
            if (!(l >= 0.0 && l <= bound)) throw new ArgumentOutOfRangeException(nameof(losses), $"Loss {i} is {l}: each loss lies in [0, B].");
            _knots[i] = l;
        }
        _knots[^1] = bound;
        Array.Sort(_knots);
        _bound = bound;
    }

    /// <summary>n, the losses.</summary>
    public int Count => _knots.Length - 1;

    /// <summary>E L⁺ = (Σ ℓ + B)/(n + 1), conformal risk control's reading.</summary>
    public double Mean
    {
        get
        {
            double sum = 0.0;
            foreach (double t in _knots) sum += t;
            return sum / _knots.Length;
        }
    }

    /// <summary>P(L⁺ &gt; x), exactly, rounded once to a double.</summary>
    public double Survival(double x)
    {
        (BigInteger num, BigInteger den) = SurvivalExact(x);
        return Ratio(num, den);
    }

    /// <summary>P(L⁺ ≤ x).</summary>
    public double Cdf(double x)
    {
        (BigInteger num, BigInteger den) = SurvivalExact(x);
        return Ratio(den - num, den);
    }

    /// <summary>
    /// b*_β: the least x with P(L⁺ ≤ x) ≥ β, read exactly and found by bisection to the doubles' resolution, so the double returned is at
    /// least the exact quantile. Gaffke's upper confidence bound on the risk at confidence β.
    /// </summary>
    public double Gaffke(double beta)
    {
        if (!(beta > 0.0 && beta < 1.0)) throw new ArgumentOutOfRangeException(nameof(beta), "β lies in (0, 1).");
        double lo = _knots[0], hi = _bound;
        if (Reaches(lo, beta)) return lo;
        // Invariant: the CDF falls short of β at lo and reaches it at hi.
        while (true)
        {
            double mid = lo + 0.5 * (hi - lo);
            if (mid <= lo || mid >= hi) return hi;
            if (Reaches(mid, beta)) hi = mid; else lo = mid;
        }
    }

    /// <summary>Gaffke's bound for <paramref name="losses"/> at confidence <paramref name="beta"/> (<see cref="Gaffke(double)"/>).</summary>
    public static double Gaffke(IReadOnlyList<double> losses, double bound, double beta) => new LossBound(losses, bound).Gaffke(beta);

    // P(L⁺ ≤ x) ≥ β, exactly: den − num ≥ β·den.
    private bool Reaches(double x, double beta)
    {
        (BigInteger num, BigInteger den) = SurvivalExact(x);
        return AtLeast(den - num, den, beta);
    }

    // P(L⁺ > x) = DD[t₀, …, tₙ] of (t − x)₊ⁿ, on integer coordinates T = t·2^s, X = x·2^s, where the order-n divided difference is the same:
    // the confluent table, the derivative f^(k)/k! = C(n, k)(T − X)₊^(n − k) where k + 1 knots tie, (T − X)₊⁰ = 1{T > X}.
    private (BigInteger Num, BigInteger Den) SurvivalExact(double x)
    {
        if (double.IsNaN(x)) throw new ArgumentException("x is NaN.", nameof(x));
        int n = _knots.Length - 1;
        if (x < _knots[0]) return (BigInteger.One, BigInteger.One);
        if (x >= _knots[^1]) return (BigInteger.Zero, BigInteger.One);
        int scale = 0;
        foreach (double t in _knots) scale = Math.Max(scale, Fraction(t));
        scale = Math.Max(scale, Fraction(x));
        var T = new BigInteger[n + 1];
        for (int i = 0; i <= n; i++) T[i] = Integer(_knots[i], scale);
        BigInteger X = Integer(x, scale);
        int low = 0;
        while (low <= n && T[low] == T[0]) low++;
        if (T[n] == T[low < n + 1 ? low : n])
        {
            // Two values a < b, with c_a = low knots at a and c_b = n + 1 − low at b: L⁺ = a + (b − a)·G with G ~ Beta(c_b, c_a), so
            // P(L⁺ > x) = P(Bin(n, u) ≤ c_b − 1), u = (x − a)/(b − a), summed exactly over the integer coordinates.
            BigInteger a = T[0], b = T[n], span = b - a, up = X - a, down = b - X, sum = BigInteger.Zero;
            int cb = n + 1 - low;
            for (int j = 0; j <= cb - 1; j++) sum += Binomial(n, j) * BigInteger.Pow(up, j) * BigInteger.Pow(down, n - j);
            return (sum, BigInteger.Pow(span, n));
        }
        var num = new BigInteger[n + 1];
        var den = new BigInteger[n + 1];
        for (int i = 0; i <= n; i++)
        {
            num[i] = T[i] > X ? BigInteger.Pow(T[i] - X, n) : BigInteger.Zero;
            den[i] = BigInteger.One;
        }
        for (int k = 1; k <= n; k++)
        {
            for (int i = n; i >= k; i--)
            {
                if (T[i] == T[i - k])
                {
                    num[i] = T[i] > X ? Binomial(n, k) * BigInteger.Pow(T[i] - X, n - k) : BigInteger.Zero;
                    den[i] = BigInteger.One;
                }
                else
                {
                    BigInteger a = num[i] * den[i - 1] - num[i - 1] * den[i], b = den[i] * den[i - 1] * (T[i] - T[i - k]);
                    BigInteger g = BigInteger.GreatestCommonDivisor(a, b);
                    if (!g.IsZero && !g.IsOne) { a /= g; b /= g; }
                    num[i] = a;
                    den[i] = b;
                }
            }
        }
        BigInteger rn = num[n], rd = den[n];
        if (rd.Sign < 0) { rn = -rn; rd = -rd; }
        if (rn.IsZero) rd = BigInteger.One;
        return (rn, rd);
    }

    // The binary places t needs to be a whole number at scale 2^s.
    private static int Fraction(double t)
    {
        if (t == 0.0) return 0;
        long bits = BitConverter.DoubleToInt64Bits(t);
        int exponent = (int)((bits >> 52) & 0x7FF);
        long mantissa = exponent == 0 ? bits & ((1L << 52) - 1) : (bits & ((1L << 52) - 1)) | (1L << 52);
        int shift = exponent == 0 ? 1074 : 1075 - exponent;
        while (shift > 0 && (mantissa & 1) == 0) { mantissa >>= 1; shift--; }
        return Math.Max(0, shift);
    }

    private static BigInteger Integer(double t, int scale)
    {
        if (t == 0.0) return BigInteger.Zero;
        long bits = BitConverter.DoubleToInt64Bits(t);
        int exponent = (int)((bits >> 52) & 0x7FF);
        long mantissa = exponent == 0 ? bits & ((1L << 52) - 1) : (bits & ((1L << 52) - 1)) | (1L << 52);
        int shift = exponent == 0 ? 1074 : 1075 - exponent;   // t = mantissa·2^(−shift)
        BigInteger value = mantissa;
        int up = scale - shift;
        return up >= 0 ? value << up : value >> -up;   // exact: scale is at least t's binary places
    }

    private static BigInteger Binomial(int n, int k)
    {
        BigInteger result = BigInteger.One;
        for (int i = 1; i <= k; i++) result = result * (n - k + i) / i;
        return result;
    }

    internal static double Ratio(BigInteger num, BigInteger den)
    {
        if (num.IsZero) return 0.0;
        if (num.Sign < 0) return -Ratio(-num, den);
        return Math.Exp(BigInteger.Log(num) - BigInteger.Log(den));
    }

    // count ≥ level·total exactly, the level read as the double it is.
    internal static bool AtLeast(BigInteger count, BigInteger total, double level)
    {
        if (level <= 0.0) return true;
        if (level >= 1.0) return count >= total;
        long bits = BitConverter.DoubleToInt64Bits(level);
        int exponent = (int)((bits >> 52) & 0x7FF);
        long fraction = bits & ((1L << 52) - 1);
        long mantissa = exponent == 0 ? fraction : fraction | (1L << 52);
        int shift = exponent == 0 ? 1074 : 1075 - exponent;
        return (count << shift) >= mantissa * total;
    }
}

/// <summary>A threshold a risk rule chose (<see cref="RiskControl"/>): its index on the grid, the reading that chose it, and whether any
/// threshold met α (the last threshold where none did).</summary>
public sealed record RiskThreshold(int Index, double Reading, bool Met);

/// <summary>
/// Threshold rules that control a loss (the conformal layer's slice 18), on a loss matrix: unit i's loss at threshold j, nonincreasing in
/// j, each in [0, B] (Angelopoulos et al., conformal risk control; Snell and Griffiths). A set family indexed by a threshold, with
/// ℓᵢ(λ) = 1{yᵢ ∉ C_λ(xᵢ)}, gives miscoverage; any bounded monotone loss serves.
/// <list type="bullet">
/// <item><see cref="Mean"/>: the first threshold with (Σ ℓ + B)/(n + 1) ≤ α. Over the orderings of any data set, the test unit's loss at it
/// averages at most α, so the expected test loss is at most α (<c>crc_orderings</c>; Angelopoulos et al.'s Theorem 1), given the losses at
/// the last threshold at most α.</item>
/// <item><see cref="HighProbability"/>: the first threshold with b*_β ≤ α (<see cref="LossBound.Gaffke(double)"/>). Its risk exceeds α only
/// where the last threshold with risk above α has b*_β below its risk (<c>fail_le_fixed</c>), with chance at most 1 − β by Gaffke's
/// calibration (Vlassis and Thomas's Theorem 1), for independent units.</item>
/// </list>
/// </summary>
public static class RiskControl
{
    /// <summary>Conformal risk control: the first threshold whose (Σᵢ ℓᵢⱼ + B)/(n + 1) is at most α.</summary>
    public static RiskThreshold Mean(double[,] losses, double bound, double alpha)
    {
        Check(losses, bound, alpha);
        int n = losses.GetLength(0), thresholds = losses.GetLength(1);
        double reading = double.NaN;
        for (int j = 0; j < thresholds; j++)
        {
            double sum = bound;
            for (int i = 0; i < n; i++) sum += losses[i, j];
            reading = sum / (n + 1);
            if (reading <= alpha) return new RiskThreshold(j, reading, true);
        }
        return new RiskThreshold(thresholds - 1, reading, false);
    }

    /// <summary>The high-probability rule: the first threshold whose Gaffke bound b*_β is at most α, found by bisection, the bound falling
    /// as the losses do.</summary>
    public static RiskThreshold HighProbability(double[,] losses, double bound, double alpha, double beta)
    {
        Check(losses, bound, alpha);
        if (!(beta > 0.0 && beta < 1.0)) throw new ArgumentOutOfRangeException(nameof(beta), "β lies in (0, 1).");
        int n = losses.GetLength(0), thresholds = losses.GetLength(1);
        var column = new double[n];
        double BoundAt(int j)
        {
            for (int i = 0; i < n; i++) column[i] = losses[i, j];
            return LossBound.Gaffke(column, bound, beta);
        }
        double last = BoundAt(thresholds - 1);
        if (!(last <= alpha)) return new RiskThreshold(thresholds - 1, last, false);
        int lo = 0, hi = thresholds - 1;
        double reading = last;
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            double b = BoundAt(mid);
            if (b <= alpha) { hi = mid; reading = b; } else lo = mid + 1;
        }
        return new RiskThreshold(lo, lo == thresholds - 1 ? last : reading, true);
    }

    private static void Check(double[,] losses, double bound, double alpha)
    {
        ArgumentNullException.ThrowIfNull(losses);
        if (!(bound >= 0.0 && double.IsFinite(bound))) throw new ArgumentOutOfRangeException(nameof(bound), "B is finite and nonnegative.");
        if (!(alpha >= 0.0 && double.IsFinite(alpha))) throw new ArgumentOutOfRangeException(nameof(alpha), "α is finite and nonnegative.");
        int n = losses.GetLength(0), thresholds = losses.GetLength(1);
        if (thresholds == 0) throw new ArgumentException("At least one threshold is required.", nameof(losses));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < thresholds; j++)
            {
                double l = losses[i, j];
                if (!(l >= 0.0 && l <= bound)) throw new ArgumentOutOfRangeException(nameof(losses), $"Unit {i}'s loss at threshold {j} is {l}: each lies in [0, B].");
                if (j > 0 && l > losses[i, j - 1]) throw new ArgumentException($"Unit {i}'s losses rise at threshold {j}: they are nonincreasing in the threshold.", nameof(losses));
            }
    }
}
