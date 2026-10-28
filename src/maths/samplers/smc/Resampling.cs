using System;
using Maths.Rng;

namespace Maths.Samplers.Smc;

/// <summary>How a weighted population is resampled to equal weights (Douc, Cappé and Moulines 2005).</summary>
public enum ResamplingScheme
{
    /// <summary>N independent draws from the weights: unbiased, the highest variance. <see cref="SmcOptions"/>' default,
    /// since it is the scheme the SMC limit theory and the normalizer's genealogy error describe.</summary>
    Multinomial,

    /// <summary>One uniform in each of N equal strata of [0, 1).</summary>
    Stratified,

    /// <summary>One uniform, shifted across the N strata: the lowest variance in practice, and each particle gets
    /// ⌊N·W_i⌋ or ⌈N·W_i⌉ offspring, but without the general convergence theory the other schemes have (Gerber,
    /// Chopin and Whiteley 2019).</summary>
    Systematic,

    /// <summary>⌊N·W_i⌋ offspring each, then the remainder drawn multinomially from the residual weights.</summary>
    Residual,
}

/// <summary>
/// Resampling: maps normalized weights to N ancestor indices, each particle's expected offspring count N·W_i under
/// every scheme. Ancestors come out in ascending order, so offspring of one particle are adjacent.
/// </summary>
public static class Resampling
{
    /// <summary>Writes one ancestor index per offspring into <paramref name="ancestors"/>, whose length is the population size.</summary>
    public static void Ancestors(ResamplingScheme scheme, ReadOnlySpan<double> weights, Xoshiro256PlusPlus rng, Span<int> ancestors)
    {
        ArgumentNullException.ThrowIfNull(rng);
        int n = ancestors.Length;
        if (weights.Length == 0) throw new ArgumentException("There are no particles.", nameof(weights));
        if (n == 0) return;
        switch (scheme)
        {
            case ResamplingScheme.Multinomial:
                Merge(weights, SortedUniforms(n, rng), ancestors);
                return;
            case ResamplingScheme.Stratified:
            {
                var u = new double[n];
                for (int k = 0; k < n; k++) u[k] = (k + rng.NextDouble()) / n;
                Merge(weights, u, ancestors);
                return;
            }
            case ResamplingScheme.Systematic:
            {
                double shift = rng.NextDouble();
                var u = new double[n];
                for (int k = 0; k < n; k++) u[k] = (k + shift) / n;
                Merge(weights, u, ancestors);
                return;
            }
            case ResamplingScheme.Residual:
                Residual(weights, rng, ancestors);
                return;
            default:
                throw new ArgumentOutOfRangeException(nameof(scheme));
        }
    }

    // ⌊N·W_i⌋ copies of each, then the remaining R offspring drawn from the residual weights (N·W_i − ⌊N·W_i⌋)/R.
    private static void Residual(ReadOnlySpan<double> weights, Xoshiro256PlusPlus rng, Span<int> ancestors)
    {
        int n = ancestors.Length, m = weights.Length;
        var counts = new int[m];
        var residual = new double[m];
        int assigned = 0;
        for (int i = 0; i < m; i++)
        {
            double expected = n * weights[i];
            counts[i] = (int)Math.Floor(expected);
            residual[i] = expected - counts[i];
            assigned += counts[i];
        }
        int remaining = n - assigned;
        if (remaining > 0)
        {
            double total = 0.0;
            foreach (double r in residual) total += r;
            for (int i = 0; i < m; i++) residual[i] /= total;
            var extra = new int[remaining];
            Merge(residual, SortedUniforms(remaining, rng), extra);
            foreach (int a in extra) counts[a]++;
        }
        int k = 0;
        for (int i = 0; i < m && k < n; i++)
            for (int c = 0; c < counts[i] && k < n; c++) ancestors[k++] = i;
    }

    // n sorted uniforms on [0, 1) from normalized cumulative exponential spacings, without a sort.
    private static double[] SortedUniforms(int n, Xoshiro256PlusPlus rng)
    {
        var u = new double[n];
        double sum = 0.0;
        for (int k = 0; k < n; k++)
        {
            sum += -Math.Log(1.0 - rng.NextDouble());
            u[k] = sum;
        }
        double total = sum - Math.Log(1.0 - rng.NextDouble());
        for (int k = 0; k < n; k++) u[k] /= total;
        return u;
    }

    // Walks ascending uniforms through the cumulative weights. A zero-weight particle is never chosen, not even the
    // last when rounding leaves the cumulative sum just short of 1.
    private static void Merge(ReadOnlySpan<double> weights, double[] sortedUniforms, Span<int> ancestors)
    {
        int last = weights.Length - 1;
        while (last > 0 && weights[last] == 0.0) last--;
        int i = 0;
        double cumulative = weights[0];
        for (int k = 0; k < sortedUniforms.Length; k++)
        {
            double u = sortedUniforms[k];
            while (i < last && (u >= cumulative || weights[i] == 0.0))
            {
                i++;
                cumulative += weights[i];
            }
            ancestors[k] = i;
        }
    }
}
