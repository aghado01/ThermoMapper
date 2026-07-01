// ============================================================================
// Hashish/CosineVectors.cs
// ============================================================================
// Cosine similarity and distance utilities for dense distributional vectors
// produced by CooccurrenceStats.PpmiMatrix() or TfIdf.
//
// Distinct from metrics/Cosine.cs which is scoped to SPC graph initialization
// (angular distance, double[] only). This variant:
//   - operates on ReadOnlySpan<double> throughout
//   - exposes raw similarity for ranked neighbor queries
//   - provides an L2-normalize-in-place helper for pre-normalizing PPMI rows
//   - provides a pairwise distance matrix builder for direct SPC input
// ============================================================================

#nullable enable
using System;
using System.Buffers;
using System.Numerics.Tensors;
using System.Runtime.CompilerServices;

namespace Hashish;

public static class CosineVectors
{
    private const double Epsilon = 1e-12;

    /// <summary>
    /// Cosine similarity in [−1, 1].
    /// Returns 0 if either vector has zero norm.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Similarity(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        double sim = TensorPrimitives.CosineSimilarity<double>(a, b);
        return double.IsNaN(sim) ? 0.0 : Math.Clamp(sim, -1.0, 1.0);
    }

    /// <summary>
    /// Angular distance in [0, 1]: arccos(similarity) / π.
    /// Suitable as a drop-in edge weight for SPC graph construction.
    /// Returns 1.0 (max distance) if either vector has zero norm.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
        => Math.Acos(Similarity(a, b)) / Math.PI;

    /// <summary>
    /// Fast dot-product distance for pre-normalized (unit) vectors.
    /// Skips norm computation entirely — only valid if both inputs are
    /// already L2-normalized. Use after <see cref="NormalizeInPlace"/>.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double DistanceNormalized(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        double dot = TensorPrimitives.Dot<double>(a, b);
        return Math.Acos(Math.Clamp(dot, -1.0, 1.0)) / Math.PI;
    }

    /// <summary>
    /// Normalizes <paramref name="vector"/> to unit L2 norm in place.
    /// No-op if the norm is below epsilon (zero vector).
    /// Call on each PPMI row before building a distance matrix to enable
    /// the faster <see cref="DistanceNormalized"/> path.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static void NormalizeInPlace(Span<double> vector)
    {
        double norm = Math.Sqrt(TensorPrimitives.Dot<double>(vector, vector));
        if (norm < Epsilon) return;
        TensorPrimitives.Multiply<double>(vector, 1.0 / norm, vector);
    }

    /// <summary>
    /// Builds a symmetric pairwise angular distance matrix from a set of
    /// row vectors. Returns a flat row-major double[] of length n×n suitable
    /// for direct input to SPC's distance matrix constructor.
    ///
    /// Vectors are L2-normalized internally; originals are not mutated.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static double[] BuildDistanceMatrix(double[][] vectors)
    {
        ArgumentNullException.ThrowIfNull(vectors);
        int n = vectors.Length;
        if (n == 0) return Array.Empty<double>();

        int d = vectors[0].Length;

        // Normalize a copy of each row into a rented flat buffer.
        double[] normBuf = ArrayPool<double>.Shared.Rent(n * d);
        try
        {
            // Copy and normalize each row.
            for (int i = 0; i < n; i++)
            {
                var src = vectors[i].AsSpan(0, d);
                var dst = normBuf.AsSpan(i * d, d);
                src.CopyTo(dst);
                NormalizeInPlace(dst);
            }

            // Build upper triangle, mirror to lower.
            var matrix = new double[n * n];
            for (int i = 0; i < n; i++)
            {
                var rowI = normBuf.AsSpan(i * d, d);
                for (int j = i + 1; j < n; j++)
                {
                    double dist = DistanceNormalized(
                        rowI, normBuf.AsSpan(j * d, d));
                    matrix[i * n + j] = dist;
                    matrix[j * n + i] = dist;
                }
            }

            return matrix;
        }
        finally
        {
            ArrayPool<double>.Shared.Return(normBuf);
        }
    }
}
