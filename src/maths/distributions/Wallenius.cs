using System;
using System.Numerics;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>How <see cref="Wallenius"/> evaluates a set's probability.</summary>
public enum WalleniusMethod
{
    /// <summary><see cref="DrawOrders"/> up to K = 10, <see cref="Quadrature"/> above.</summary>
    Auto,

    /// <summary>The sum over the K! draw orders, accumulated over the 2^K subsets of the set: every term and every
    /// denominator a sum of positives, so exact to a few ulps at any weights, in O(K·2^K).</summary>
    DrawOrders,

    /// <summary>Adaptive Gauss–Kronrod quadrature of Chesson's positive integrand, in O(K) per node.</summary>
    Quadrature,
}

/// <summary>
/// Wallenius' multivariate noncentral hypergeometric distribution with one ball per variable (Wallenius 1963, in
/// Chesson's 1976 integral form): the law of the set S of K distinct indices drawn by weighted sampling without
/// replacement, each draw choosing among the indices not yet drawn in proportion to its weight z, so that
/// P(S) = ∫₀¹ Πᵢ∈S (1 − t^αᵢ) dt with αᵢ = zᵢ / Σ_{k∉S} zₖ. Its inclusion–exclusion expansion (R BASS's
/// <c>dmwnchBass</c>) alternates in sign and loses digits under skewed weights, so neither evaluation here uses it:
/// the draw-order recursion and the quadrature both sum positive terms, and both work in log space, so log P(S) keeps
/// its digits where P(S) underflows. K = 1 is zᵢ/Σz, and the set of every index has probability 1.
/// </summary>
public static class Wallenius
{
    private const int AutoDrawOrdersMaxOrder = 10;
    private const int DrawOrdersMaxOrder = 20;

    // The quadrature runs on v = ln u, t = e^(−u). Below u = LowerCut/(1 + α_max) the integrand's share of P(S) is under
    // e²·LowerCut^(K+1); above u = 60 + 10K, under e³·Γ(K + 1, 61 + 10K). Both are below 1e-17 for K ≥ 2.
    private const double LowerCut = 1e-6;
    private const double QuadratureTolerance = 1e-14;
    private const int MaxPanels = 512;

    /// <summary>ln P(S): the log probability that weighted sampling without replacement draws exactly the indices in
    /// <paramref name="set"/> (in any order) from <paramref name="weights"/>.</summary>
    public static double LogProbability(
        ReadOnlySpan<double> weights, ReadOnlySpan<int> set, WalleniusMethod method = WalleniusMethod.Auto)
    {
        Validate(weights, set);
        int k = set.Length;
        if (k == weights.Length) return 0.0;

        // The weight left in the urn, summed directly rather than as a difference from the total, which cancels when
        // the set holds the large weights.
        double rest = 0.0;
        for (int u = 0; u < weights.Length; u++)
            if (!Contains(set, u)) rest += weights[u];

        if (k == 1)
        {
            double z = weights[set[0]], ratio = rest / z;
            return double.IsPositiveInfinity(ratio) ? Math.Log(z) - Math.Log(rest) : -SpecialFunctions.Log1p(ratio);
        }

        switch (method)
        {
            case WalleniusMethod.Auto:
                return k <= AutoDrawOrdersMaxOrder ? LogDrawOrders(weights, set, rest) : LogQuadrature(weights, set, rest);
            case WalleniusMethod.DrawOrders:
                if (k > DrawOrdersMaxOrder)
                    throw new ArgumentOutOfRangeException(nameof(set), $"The draw-order recursion runs to K = {DrawOrdersMaxOrder}.");
                return LogDrawOrders(weights, set, rest);
            case WalleniusMethod.Quadrature:
                return LogQuadrature(weights, set, rest);
            default:
                throw new ArgumentOutOfRangeException(nameof(method));
        }
    }

    /// <summary>P(S); see <see cref="LogProbability"/>.</summary>
    public static double Probability(
        ReadOnlySpan<double> weights, ReadOnlySpan<int> set, WalleniusMethod method = WalleniusMethod.Auto)
        => Math.Exp(LogProbability(weights, set, method));

    /// <summary>One draw: fills <paramref name="set"/> with K = <c>set.Length</c> distinct indices in draw order, each
    /// chosen among the indices not yet drawn in proportion to its weight, one uniform per draw. The drawn set's law
    /// is <see cref="LogProbability"/>'s.</summary>
    public static void Sample(Xoshiro256PlusPlus rng, ReadOnlySpan<double> weights, Span<int> set)
    {
        ArgumentNullException.ThrowIfNull(rng);
        ValidateWeights(weights);
        if (set.Length < 1 || set.Length > weights.Length)
            throw new ArgumentOutOfRangeException(nameof(set), "A draw takes between 1 and p indices.");

        for (int d = 0; d < set.Length; d++)
        {
            ReadOnlySpan<int> drawn = set[..d];
            double remaining = 0.0;
            for (int u = 0; u < weights.Length; u++)
                if (!Contains(drawn, u)) remaining += weights[u];

            // The last index still in the urn takes a target that rounding leaves past the running sum.
            double target = rng.NextDouble() * remaining, cumulative = 0.0;
            int pick = -1;
            for (int u = 0; u < weights.Length; u++)
            {
                if (Contains(drawn, u)) continue;
                pick = u;
                cumulative += weights[u];
                if (target < cumulative) break;
            }
            set[d] = pick;
        }
    }

    // G(U) = Σ_{i∈U} zᵢ·G(U∖i) / (R + Z_U), G(∅) = 1, is the probability that the first |U| draws from U and the rest R
    // take exactly U: the sum over U's draw orders. Scaled by Π_{i∈U} min(zᵢ/R, 1) it is
    // H(U) = Σ_{i∈U} max(R, zᵢ)·H(U∖i) / (R + Z_U), which lies in [1/(K+1)!, K!], so it neither overflows nor underflows.
    private static double LogDrawOrders(ReadOnlySpan<double> weights, ReadOnlySpan<int> set, double rest)
    {
        int k = set.Length, n = 1 << k;
        Span<double> h = k <= AutoDrawOrdersMaxOrder ? stackalloc double[n] : new double[n];
        Span<double> mass = k <= AutoDrawOrdersMaxOrder ? stackalloc double[n] : new double[n];
        h[0] = 1.0;
        mass[0] = 0.0;
        for (int mask = 1; mask < n; mask++)
        {
            mass[mask] = mass[mask & (mask - 1)] + weights[set[BitOperations.TrailingZeroCount(mask)]];
            double sum = 0.0;
            for (int left = mask; left != 0; left &= left - 1)
            {
                int i = BitOperations.TrailingZeroCount(left);
                sum += Math.Max(rest, weights[set[i]]) * h[mask & ~(1 << i)];
            }
            h[mask] = sum / (rest + mass[mask]);
        }

        double logRest = Math.Log(rest), logScale = 0.0;
        foreach (int s in set) logScale += Math.Min(Math.Log(weights[s]) - logRest, 0.0);
        return logScale + Math.Log(h[n - 1]);
    }

    // ln P(S) = Σ ln min(αᵢ, 1) + ln ∫ e^(v − e^v) Πᵢ fᵢ(e^v) dv, with fᵢ(u) = (1 − e^(−αᵢu)) / min(αᵢ, 1) ≤ max(u, 1),
    // so the integral neither underflows nor overflows where P(S) does.
    private static double LogQuadrature(ReadOnlySpan<double> weights, ReadOnlySpan<int> set, double rest)
    {
        int k = set.Length;
        Span<double> alpha = k <= 64 ? stackalloc double[k] : new double[k];
        double logScale = 0.0, alphaMax = 0.0;
        for (int i = 0; i < k; i++)
        {
            double a = alpha[i] = weights[set[i]] / rest;
            logScale += Math.Log(Math.Min(a, 1.0));
            if (!double.IsPositiveInfinity(a)) alphaMax = Math.Max(alphaMax, a);
        }
        double lower = Math.Log(LowerCut) - SpecialFunctions.Log1p(alphaMax);
        double upper = Math.Log(60.0 + 10.0 * k);
        return logScale + Math.Log(AdaptiveKronrod(alpha, lower, upper));
    }

    private static double Integrand(ReadOnlySpan<double> alpha, double v)
    {
        double u = Math.Exp(v), value = Math.Exp(v - u);
        foreach (double a in alpha)
        {
            double f = -SpecialFunctions.Expm1(-a * u);
            value *= a < 1.0 ? f / a : f;
        }
        return value;
    }

    // Global adaptive Gauss–Kronrod (G7–K15): panels of width at most 1 in v, the panel with the largest |K15 − G7|
    // bisected until the summed estimate is within QuadratureTolerance of the integral.
    private static double AdaptiveKronrod(ReadOnlySpan<double> alpha, double lower, double upper)
    {
        Span<double> left = stackalloc double[MaxPanels], right = stackalloc double[MaxPanels];
        Span<double> value = stackalloc double[MaxPanels], error = stackalloc double[MaxPanels];
        int panels = Math.Clamp((int)Math.Ceiling(upper - lower), 1, MaxPanels / 2);
        double width = (upper - lower) / panels;
        for (int j = 0; j < panels; j++)
        {
            left[j] = lower + j * width;
            right[j] = j == panels - 1 ? upper : lower + (j + 1) * width;
            (value[j], error[j]) = Kronrod15(alpha, left[j], right[j]);
        }

        while (true)
        {
            double total = 0.0, totalError = 0.0;
            int worst = 0;
            for (int j = 0; j < panels; j++)
            {
                total += value[j];
                totalError += error[j];
                if (error[j] > error[worst]) worst = j;
            }
            if (totalError <= QuadratureTolerance * total) return total;
            if (panels == MaxPanels)
                throw new ArithmeticException($"Wallenius quadrature did not converge in {MaxPanels} panels.");

            double mid = 0.5 * (left[worst] + right[worst]);
            left[panels] = mid;
            right[panels] = right[worst];
            right[worst] = mid;
            (value[worst], error[worst]) = Kronrod15(alpha, left[worst], mid);
            (value[panels], error[panels]) = Kronrod15(alpha, mid, right[panels]);
            panels++;
        }
    }

    // QUADPACK's qk15 nodes and weights: Kronrod abscissae, the Gauss points at odd indices and the centre.
    private static readonly double[] KronrodNodes =
    {
        0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
        0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
        0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
        0.207784955007898467600689403773245, 0.0,
    };

    private static readonly double[] KronrodWeights =
    {
        0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
        0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
        0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
        0.204432940075298892414161999234649, 0.209482141084727828012999174891714,
    };

    private static readonly double[] GaussWeights =
    {
        0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
        0.381830050505118944950369775488975, 0.417959183673469387755102040816327,
    };

    private static (double Value, double Error) Kronrod15(ReadOnlySpan<double> alpha, double a, double b)
    {
        double centre = 0.5 * (a + b), half = 0.5 * (b - a);
        double fc = Integrand(alpha, centre);
        double kronrod = KronrodWeights[7] * fc, gauss = GaussWeights[3] * fc;
        for (int j = 0; j < 7; j++)
        {
            double dx = half * KronrodNodes[j];
            double pair = Integrand(alpha, centre - dx) + Integrand(alpha, centre + dx);
            kronrod += KronrodWeights[j] * pair;
            if ((j & 1) == 1) gauss += GaussWeights[j / 2] * pair;
        }
        return (kronrod * half, Math.Abs(kronrod - gauss) * half);
    }

    private static bool Contains(ReadOnlySpan<int> set, int u)
    {
        foreach (int s in set)
            if (s == u) return true;
        return false;
    }

    private static void ValidateWeights(ReadOnlySpan<double> weights)
    {
        if (weights.Length < 1) throw new ArgumentException("At least one weight is required.", nameof(weights));
        foreach (double z in weights)
            if (!(z > 0.0) || double.IsPositiveInfinity(z))
                throw new ArgumentOutOfRangeException(nameof(weights), "Weights must be positive and finite.");
    }

    private static void Validate(ReadOnlySpan<double> weights, ReadOnlySpan<int> set)
    {
        ValidateWeights(weights);
        if (set.Length < 1 || set.Length > weights.Length)
            throw new ArgumentOutOfRangeException(nameof(set), "A set holds between 1 and p indices.");
        for (int i = 0; i < set.Length; i++)
        {
            if ((uint)set[i] >= (uint)weights.Length)
                throw new ArgumentOutOfRangeException(nameof(set), $"Index {set[i]} lies outside the {weights.Length} weights.");
            for (int j = 0; j < i; j++)
                if (set[j] == set[i]) throw new ArgumentException($"Index {set[i]} appears twice in the set.", nameof(set));
        }
    }
}
