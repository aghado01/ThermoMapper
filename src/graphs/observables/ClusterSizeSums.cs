using System;

namespace Graphs.Observables;

/// <summary>
/// Means of sums over clusters, ⟨Σ_k g(c_k)⟩, reduced from an exact cluster-size histogram — a degree-0 reduction
/// over counts n_s of clusters at each size s, summed over draws. Each such sum is linear in the counts, so its mean
/// is exact: ⟨Σ_k g(c_k)⟩ = Σ_s g(s)·n_s / draws. Graph-intrinsic and model-agnostic, beside
/// <see cref="ClusterSizeEntropy"/>; an SW run's <c>ExactClusterSizeHistogram</c> is the feed.
/// </summary>
/// <remarks>
/// The histogram is indexed by size: entry s counts clusters of size s, and entry 0 is ignored. Returns 0 when there
/// are no draws, as the Potts reductions do.
/// </remarks>
public static class ClusterSizeSums
{
    /// <summary>⟨Σ_k g(c_k)⟩ for any per-cluster function <paramref name="g"/> of the size.</summary>
    /// <param name="countsBySize">n_s, the clusters of size s over the draws.</param>
    /// <param name="draws">The number of draws the counts were summed over (0 ⇒ 0).</param>
    /// <param name="g">The per-cluster function of the size.</param>
    public static double Mean(ReadOnlySpan<long> countsBySize, long draws, Func<int, double> g)
    {
        ArgumentNullException.ThrowIfNull(g);
        if (draws <= 0)
            return 0.0;

        double sum = 0.0;
        for (int s = 1; s < countsBySize.Length; s++)
            if (countsBySize[s] != 0)
                sum += g(s) * countsBySize[s];
        return sum / draws;
    }

    /// <summary>S₂ = ⟨Σ_k c_k²⟩/N², the chance that two nodes drawn with replacement share a cluster.</summary>
    /// <param name="countsBySize">n_s, the clusters of size s over the draws.</param>
    /// <param name="draws">The number of draws (0 ⇒ 0).</param>
    /// <param name="nodeCount">N, the number of nodes.</param>
    public static double SecondMoment(ReadOnlySpan<long> countsBySize, long draws, int nodeCount)
    {
        if (nodeCount <= 0)
            throw new ArgumentOutOfRangeException(nameof(nodeCount), "The node count must be positive.");
        return Mean(countsBySize, draws, s => (double)s * s) / ((double)nodeCount * nodeCount);
    }

    /// <summary>
    /// The size-floored count F_n = ⟨Σ_k (c_k − n)₊⟩/n: the clusters' mass above the floor <paramref name="floor"/>,
    /// in units of the floor, at any floor chosen after the run.
    /// </summary>
    /// <param name="countsBySize">n_s, the clusters of size s over the draws.</param>
    /// <param name="draws">The number of draws (0 ⇒ 0).</param>
    /// <param name="floor">The floor n &gt; 0.</param>
    public static double FlooredCount(ReadOnlySpan<long> countsBySize, long draws, double floor)
    {
        if (!(floor > 0.0))
            throw new ArgumentOutOfRangeException(nameof(floor), "The floor must be positive.");
        return Mean(countsBySize, draws, s => Math.Max(0.0, s - floor)) / floor;
    }

    /// <summary>
    /// ⟨Σ_k c_k ln c_k⟩, which carries the expected partition entropy: ln N − ⟨Σ_k c_k ln c_k⟩/N nats.
    /// </summary>
    /// <param name="countsBySize">n_s, the clusters of size s over the draws.</param>
    /// <param name="draws">The number of draws (0 ⇒ 0).</param>
    public static double SizeLogSize(ReadOnlySpan<long> countsBySize, long draws)
        => Mean(countsBySize, draws, s => s * Math.Log(s));
}
