// ============================================================================
// TDA.Mapper — MapperKMeans.cs
// ============================================================================
// IClusterer adapter for KMeans++ (Clustering.KMeans.KMeansPlusPlus).
// K-means is the default clusterer for Mapper preimages — each preimage is
// small (typically 10–200 points) so Lloyd iteration cost is negligible.
// Pass k=0 for auto-k = clamp(√(N/3), 2, 10); pass k>0 to fix K.
// ============================================================================

#nullable enable
using System;
using Clustering.Geometric.KMeans;
using TDA.Mapper;

namespace TDA.Mapper.Clusterers;

public sealed class KMeansClusterer : IClusterer
{
    /// <summary>Fixed k, or 0 for auto-k = clamp(√(N/3), 2, 10).</summary>
    public int K { get; }
    public int MaxIterations { get; }
    public int Seed { get; }

    public string Name => K > 0 ? $"KMeans (k={K})" : "KMeans (auto-k)";

    public KMeansClusterer(int k = 0, int maxIterations = 100, int seed = 42)
    {
        if (k < 0) throw new ArgumentOutOfRangeException(nameof(k), "k must be >= 0 (0 = auto)");
        if (maxIterations < 1) throw new ArgumentOutOfRangeException(nameof(maxIterations));

        K = k;
        MaxIterations = maxIterations;
        Seed = seed;
    }

    public ClusterResult Cluster(double[][] subset)
    {
        if (subset is null || subset.Length == 0)
            return new ClusterResult(Array.Empty<int>(), 0);

        // Trivial subsets: single point or pair → no clustering needed.
        if (subset.Length == 1)
            return new ClusterResult(new[] { 0 }, 1);

        int actualK = K > 0
            ? Math.Min(K, subset.Length)
            : Math.Clamp((int)Math.Sqrt(subset.Length / 3.0), 2, Math.Min(10, subset.Length));

        var result = KMeansPlusPlus.Cluster(
            data: subset,
            k: actualK,
            maxIterations: MaxIterations,
            seed: Seed);

        return new ClusterResult(result.Labels, actualK);
    }
}
