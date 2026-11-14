using System;
using Maths.Rng;

namespace Clustering.Geometric.KMeans;

/// <summary>
/// k-means++ seeding (Arthur &amp; Vassilvitskii 2007): the first seed is a uniform draw from the data, each later one
/// a draw with probability proportional to the squared distance to the nearest seed already chosen. Draws come from a
/// pinned <see cref="Xoshiro256PlusPlus"/> stream, so a seed set is a function of the data and the generator state.
/// </summary>
public static class KMeansPlusPlusSeeding
{
    /// <summary>
    /// Returns <paramref name="k"/> seeds, each a copy of a data row. A row already chosen has squared distance zero
    /// and is drawn again only when every row coincides with a seed, in which case the draw is uniform.
    /// Costs one pass over the data per seed.
    /// </summary>
    public static double[][] Seed(double[][] data, int k, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(rng);
        int n = data.Length;
        if (n == 0) throw new ArgumentException("Data cannot be empty.", nameof(data));
        if (k < 1 || k > n) throw new ArgumentOutOfRangeException(nameof(k), "k must be between 1 and n.");

        var seeds = new double[k][];
        var nearest = new double[n];
        Array.Fill(nearest, double.PositiveInfinity);

        int pick = UniformIndex(rng, n);
        for (int c = 0; c < k; c++)
        {
            double[] seed = (double[])data[pick].Clone();
            seeds[c] = seed;
            if (c == k - 1) break;

            double total = 0.0;
            for (int i = 0; i < n; i++)
            {
                double dist = SquaredDistance(data[i], seed);
                if (dist < nearest[i]) nearest[i] = dist;
                total += nearest[i];
            }

            if (!(total > 0.0))
            {
                pick = UniformIndex(rng, n);
                continue;
            }

            double target = rng.NextDouble() * total;
            double cumulative = 0.0;
            int lastPositive = -1;
            pick = -1;
            for (int i = 0; i < n; i++)
            {
                if (!(nearest[i] > 0.0)) continue;
                lastPositive = i;
                cumulative += nearest[i];
                if (cumulative > target)
                {
                    pick = i;
                    break;
                }
            }

            // The running sum can end a rounding error short of the total; the last eligible row takes the remainder.
            if (pick < 0) pick = lastPositive;
        }

        return seeds;
    }

    private static int UniformIndex(Xoshiro256PlusPlus rng, int bound)
        => (int)Math.BigMul(rng.NextUInt64(), (ulong)bound, out _);

    private static double SquaredDistance(double[] a, double[] b)
    {
        double sum = 0.0;
        for (int i = 0; i < a.Length; i++)
        {
            double diff = a[i] - b[i];
            sum += diff * diff;
        }
        return sum;
    }
}
