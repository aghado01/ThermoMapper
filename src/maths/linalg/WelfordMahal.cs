// ============================================================================
// Maths/Linalg/WelfordMahal.cs
// ============================================================================
// Incremental diagonal Mahalanobis distance via Welford's online algorithm.
//
// Maintains per-dimension running mean and M2 (sum of squared deviations)
// using Welford's method — numerically stable, single-pass, O(D) memory.
//
// This is NOT a substitute for WeiszfeldScatter or metrics/Mahalanobis.cs:
//   - WeiszfeldScatter: full D×D robust batch scatter (requires all data)
//   - metrics/Mahalanobis.cs: full covariance-inverse distance (batch)
//   - This type: diagonal approximation, streaming/memory-constrained contexts
//
// The resulting distance is standardized Euclidean — each dimension is
// normalized by its sample standard deviation. This satisfies the triangle
// inequality and is a proper metric, but does not capture inter-dimensional
// correlations. Use when: (a) data arrives as a stream, (b) D is large enough
// that a full D×D covariance matrix is impractical, or (c) dimensions are
// approximately uncorrelated by construction (e.g. PCA-reduced features).
//
// Thread safety: none. Wrap in a lock or use separate instances per thread.
// ============================================================================

#nullable enable
using System;
using System.Runtime.CompilerServices;

namespace Maths.LinAlg;

/// <summary>
/// Incremental diagonal Mahalanobis distance estimator using Welford's
/// online algorithm. Update with observations one at a time; query distance
/// at any point with at least two observations accumulated.
/// </summary>
public sealed class OnlineMahalanobis
{
    private readonly double[] _mean;
    private readonly double[] _m2;       // sum of squared deviations from mean
    private readonly double[] _scratch;  // reused per Distance() call
    private int _count;

    /// <summary>Dimensionality of the observation space.</summary>
    public int Dimension { get; }

    /// <summary>Number of observations accumulated so far.</summary>
    public int Count => _count;

    /// <summary>
    /// Whether at least two observations have been accumulated.
    /// Distance() returns 0 for any query until this is true.
    /// </summary>
    public bool IsReady => _count >= 2;

    /// <param name="dimension">
    /// Number of dimensions per observation. Must be &gt; 0.
    /// </param>
    public OnlineMahalanobis(int dimension)
    {
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(dimension, 0, nameof(dimension));
        Dimension = dimension;
        _mean = new double[dimension];
        _m2 = new double[dimension];
        _scratch = new double[dimension];
    }

    /// <summary>
    /// Incorporates a new observation into the running statistics.
    /// Uses Welford's numerically stable update.
    /// </summary>
    /// <param name="observation">
    /// Array of length <see cref="Dimension"/>. Not mutated.
    /// </param>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public void Update(ReadOnlySpan<double> observation)
    {
        if (observation.Length != Dimension)
            throw new ArgumentException(
                $"Observation length {observation.Length} does not match Dimension {Dimension}.",
                nameof(observation));

        _count++;
        for (int i = 0; i < Dimension; i++)
        {
            double delta = observation[i] - _mean[i];
            _mean[i] += delta / _count;
            double delta2 = observation[i] - _mean[i];
            _m2[i] += delta * delta2;
        }
    }

    /// <summary>Array overload of <see cref="Update(ReadOnlySpan{double})"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Update(double[] observation)
        => Update(observation.AsSpan());

    /// <summary>
    /// Computes the diagonal Mahalanobis distance from
    /// <paramref name="observation"/> to the current running mean.
    /// Returns 0 if fewer than two observations have been accumulated.
    /// </summary>
    /// <param name="observation">
    /// Array of length <see cref="Dimension"/>. Not mutated.
    /// </param>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public double Distance(ReadOnlySpan<double> observation)
    {
        if (observation.Length != Dimension)
            throw new ArgumentException(
                $"Observation length {observation.Length} does not match Dimension {Dimension}.",
                nameof(observation));

        if (_count < 2) return 0.0;

        double sumSq = 0.0;
        double invCountMinus1 = 1.0 / (_count - 1);

        for (int i = 0; i < Dimension; i++)
        {
            double variance = _m2[i] * invCountMinus1;
            if (variance <= 0.0) continue;  // constant dimension — skip
            double diff = observation[i] - _mean[i];
            sumSq += (diff * diff) / variance;
        }

        return Math.Sqrt(sumSq);
    }

    /// <summary>Array overload of <see cref="Distance(ReadOnlySpan{double})"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public double Distance(double[] observation)
        => Distance(observation.AsSpan());

    /// <summary>
    /// Computes the diagonal Mahalanobis distance between two arbitrary
    /// observations using the current running statistics as the normalizer.
    /// Useful for pairwise distance queries against a reference distribution.
    /// Returns 0 if fewer than two observations have been accumulated.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        if (a.Length != Dimension)
            throw new ArgumentException(
                $"Length {a.Length} does not match Dimension {Dimension}.", nameof(a));
        if (b.Length != Dimension)
            throw new ArgumentException(
                $"Length {b.Length} does not match Dimension {Dimension}.", nameof(b));

        if (_count < 2) return 0.0;

        double sumSq = 0.0;
        double invCountMinus1 = 1.0 / (_count - 1);

        for (int i = 0; i < Dimension; i++)
        {
            double variance = _m2[i] * invCountMinus1;
            if (variance <= 0.0) continue;
            double diff = a[i] - b[i];
            sumSq += (diff * diff) / variance;
        }

        return Math.Sqrt(sumSq);
    }

    /// <summary>
    /// Returns a snapshot of the current per-dimension sample variances.
    /// Variance[i] = M2[i] / (Count - 1). Returns zeros if Count &lt; 2.
    /// </summary>
    public double[] GetVariances()
    {
        var variances = new double[Dimension];
        if (_count < 2) return variances;

        double invCountMinus1 = 1.0 / (_count - 1);
        for (int i = 0; i < Dimension; i++)
            variances[i] = _m2[i] * invCountMinus1;

        return variances;
    }

    /// <summary>
    /// Returns a snapshot of the current running mean.
    /// </summary>
    public double[] GetMean()
    {
        var mean = new double[Dimension];
        _mean.AsSpan().CopyTo(mean);
        return mean;
    }

    /// <summary>
    /// Resets all accumulated statistics. Dimension is preserved.
    /// </summary>
    public void Reset()
    {
        _count = 0;
        Array.Clear(_mean, 0, Dimension);
        Array.Clear(_m2, 0, Dimension);
        Array.Clear(_scratch, 0, Dimension);
    }
}
