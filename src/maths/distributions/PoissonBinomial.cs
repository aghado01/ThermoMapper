using System;

namespace Maths.Distributions;

/// <summary>
/// The Poisson-binomial law: the number K of successes among independent coins with probabilities p₁, …, p_p. Its pmf
/// comes by dynamic programming over the coins, each coin mixing the counts so far as fᵏ ← fᵏ(1 − pⱼ) + fᵏ⁻¹pⱼ, truncated at
/// the largest count asked for: O(p·k_max), every entry a sum of nonnegative products, so each keeps its relative precision
/// with no cancellation. KHAOS's birth reads it for the validity of an inclusion draw, P(1 ≤ K ≤ q_max) (E4).
/// </summary>
public static class PoissonBinomial
{
    /// <summary>
    /// Fills <paramref name="pmf"/>[k] with P(K = k) for k below its length, for coins with success probabilities
    /// <paramref name="probabilities"/>, each in [0, 1]; counts above the coins' number have probability 0.
    /// </summary>
    public static void Pmf(ReadOnlySpan<double> probabilities, Span<double> pmf)
    {
        if (pmf.Length == 0) throw new ArgumentException("The pmf needs an entry.", nameof(pmf));
        pmf.Clear();
        pmf[0] = 1.0;
        int top = pmf.Length - 1;
        for (int j = 0; j < probabilities.Length; j++)
        {
            double p = probabilities[j];
            if (!(p >= 0.0 && p <= 1.0)) throw new ArgumentOutOfRangeException(nameof(probabilities), $"Coin {j}'s probability {p} lies outside [0, 1].");
            double q = 1.0 - p;
            for (int k = Math.Min(j + 1, top); k > 0; k--) pmf[k] = pmf[k] * q + pmf[k - 1] * p;
            pmf[0] *= q;
        }
    }

    /// <summary>P(<paramref name="lower"/> ≤ K ≤ <paramref name="upper"/>), summed over the counts from the pmf.</summary>
    public static double Probability(ReadOnlySpan<double> probabilities, int lower, int upper)
    {
        if (lower < 0) throw new ArgumentOutOfRangeException(nameof(lower), "A count is not negative.");
        if (upper < lower) return 0.0;
        int top = Math.Min(upper, probabilities.Length);
        if (lower > top)
        {
            Pmf(probabilities, stackalloc double[1]);   // validates the coins
            return 0.0;
        }
        Span<double> pmf = top < 256 ? stackalloc double[top + 1] : new double[top + 1];
        Pmf(probabilities, pmf);
        double sum = 0.0;
        for (int k = lower; k <= top; k++) sum += pmf[k];
        return sum;
    }

    /// <summary>P(K ≤ <paramref name="k"/>).</summary>
    public static double Cdf(ReadOnlySpan<double> probabilities, int k) => k < 0 ? 0.0 : Probability(probabilities, 0, k);
}
