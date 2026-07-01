using System;
using Maths.Distance;

namespace Graphs.Primitives.Mst;

/// <summary>
/// Core-distance computation shared by density clustering (HDBSCAN) and
/// mutual-reachability graph construction.
/// </summary>
public static class CoreDistances
{
    /// <summary>
    /// Computes the core distance for each point: the distance to its
    /// <paramref name="minPts"/>-th nearest neighbor. Results are written to
    /// <paramref name="result"/> (length ≥ <paramref name="n"/>).
    /// </summary>
    public static void Compute<TMetric>(
        ReadOnlySpan<double> data,
        int                  n,
        int                  dim,
        int                  minPts,
        TMetric              metric,
        Span<double>         result)
        where TMetric : struct, IDistanceMetric
    {
        if (n < 1)
            throw new ArgumentOutOfRangeException(nameof(n), "n must be positive.");
        if (minPts < 2)
            throw new ArgumentOutOfRangeException(nameof(minPts), "minPts must be >= 2.");
        if (minPts > n)
            throw new ArgumentOutOfRangeException(nameof(minPts), "minPts must be <= n.");
        if (dim <= 0)
            throw new ArgumentOutOfRangeException(nameof(dim), "dim must be positive.");
        if (data.Length < n * dim)
            throw new ArgumentException("Data length is too small for the supplied shape.", nameof(data));
        if (result.Length < n)
            throw new ArgumentException("result must have length at least n.", nameof(result));

        var kBuf = new double[minPts];
        for (int i = 0; i < n; i++)
        {
            ReadOnlySpan<double> rowI = data.Slice(i * dim, dim);
            for (int j = 0; j < minPts; j++)
                kBuf[j] = double.PositiveInfinity;
            double kthLargest = double.PositiveInfinity;

            for (int j = 0; j < n; j++)
            {
                if (i == j) continue;
                ReadOnlySpan<double> rowJ = data.Slice(j * dim, dim);
                double d = metric.Distance(rowI, rowJ);
                if (d < kthLargest)
                {
                    int pos = minPts - 1;
                    kBuf[pos] = d;
                    while (pos > 0 && kBuf[pos] < kBuf[pos - 1])
                    {
                        (kBuf[pos], kBuf[pos - 1]) = (kBuf[pos - 1], kBuf[pos]);
                        pos--;
                    }
                    kthLargest = kBuf[minPts - 1];
                }
            }

            result[i] = kBuf[minPts - 1];
        }
    }
}
