using System;
using System.Buffers;

namespace Maths.Geometry.Estimators.Ambient
// ============================================================================
// Estimators/Euclidean/EuclideanMean.cs
// ============================================================================
// Flat-space specialisation of GeometricMean — closed-form, single-pass.
// Drop-in for GeometricMean when you don't need the manifold abstraction.
//
// Output shapes mirror GeometricMean exactly:
//   Compute(...)             -> void, destination mutated in place
//   ComputeWithScatter(...)  -> void, destination + scatterDestination mutated
//
// scatterDestination is a row-major D*D flat buffer (length D²),
// identical layout to GeometricMean.ComputeWithScatter.
// ============================================================================

{
    public static class EuclideanMean
    {
        // ── Location only ─────────────────────────────────────────────────────

        /// <summary>
        /// Computes the weighted Euclidean mean into <paramref name="destination"/>.
        /// Single pass, closed-form exact.
        /// </summary>
        public static void Compute(
            double[][] points,
            ReadOnlySpan<int> indices,
            ReadOnlySpan<double> weights,
            Span<double> destination)
        {
            ValidateInputs(points, indices, weights, destination, out _);
            if (ComputeCore(points, indices, weights, destination) <= 0)
                throw new ArgumentException(
                    "At least one weight must be positive.", nameof(weights));
        }

        // ── (Location, Scatter) ───────────────────────────────────────────────

        /// <summary>
        /// Computes the weighted Euclidean mean and sample covariance in two passes.
        /// <para>
        /// <paramref name="scatterDestination"/> is a row-major D×D flat buffer
        /// of length D². Layout identical to <c>GeometricMean.ComputeWithScatter</c>.
        /// </para>
        /// </summary>
        public static void ComputeWithScatter(
            double[][] points,
            ReadOnlySpan<int> indices,
            ReadOnlySpan<double> weights,
            Span<double> destination,
            Span<double> scatterDestination,
            double consistencyFactor = 1.0)
        {
            ValidateInputs(points, indices, weights, destination, out int dim);
            double totalWeight = ComputeCore(points, indices, weights, destination);
            if (totalWeight <= 0)
                throw new ArgumentException(
                    "At least one weight must be positive.", nameof(weights));

            scatterDestination.Clear();

            double[] vArr = ArrayPool<double>.Shared.Rent(dim);
            try
            {
                Span<double> v = vArr.AsSpan(0, dim);
                for (int i = 0; i < indices.Length; i++)
                {
                    double w = weights[i];
                    if (w == 0) continue;
                    double[] pt = points[indices[i]];
                    for (int d = 0; d < dim; d++) v[d] = pt[d] - destination[d];
                    for (int r = 0; r < dim; r++)
                        for (int c = 0; c < dim; c++)
                            scatterDestination[r * dim + c] += w * v[r] * v[c];
                }

                double scale = consistencyFactor / totalWeight;
                for (int k = 0; k < scatterDestination.Length; k++)
                    scatterDestination[k] *= scale;
            }
            finally
            {
                ArrayPool<double>.Shared.Return(vArr);
            }
        }

        // ── Shared validation + core computation ──────────────────────────────

        // Computes the weighted mean into destination and returns the total weight.
        // Caller is responsible for throwing if the returned value is <= 0.
        private static double ComputeCore(
            double[][] points,
            ReadOnlySpan<int> indices,
            ReadOnlySpan<double> weights,
            Span<double> destination)
        {
            destination.Clear();
            double totalWeight = 0;
            int dim = destination.Length;
            for (int i = 0; i < indices.Length; i++)
            {
                double w = weights[i];
                if (w == 0) continue;
                totalWeight += w;
                double[] pt = points[indices[i]];
                for (int d = 0; d < dim; d++) destination[d] += w * pt[d];
            }
            if (totalWeight > 0)
            {
                double inv = 1.0 / totalWeight;
                for (int d = 0; d < dim; d++) destination[d] *= inv;
            }
            return totalWeight;
        }

        internal static void ValidateInputs(
            double[][] points,
            ReadOnlySpan<int> indices,
            ReadOnlySpan<double> weights,
            Span<double> destination,
            out int dim)
        {
            if (points == null)
                throw new ArgumentNullException(nameof(points));
            if (indices.Length == 0)
                throw new ArgumentException(
                    "At least one point index is required.", nameof(indices));
            if (indices.Length != weights.Length)
                throw new ArgumentException(
                    "Indices and weights must have the same length.", nameof(weights));

            int first = indices[0];
            if ((uint)first >= (uint)points.Length || points[first] == null)
                throw new ArgumentOutOfRangeException(nameof(indices));

            dim = points[first].Length;
            if (dim == 0)
                throw new ArgumentException(
                    "Feature vectors must not be empty.", nameof(points));
            if (destination.Length < dim)
                throw new ArgumentException(
                    "Destination span is shorter than the feature dimension.",
                    nameof(destination));

            for (int i = 0; i < indices.Length; i++)
            {
                int idx = indices[i];
                if ((uint)idx >= (uint)points.Length || points[idx] == null)
                    throw new ArgumentOutOfRangeException(nameof(indices));
                if (points[idx].Length != dim)
                    throw new ArgumentException(
                        "All selected feature vectors must have the same dimension.",
                        nameof(points));
                double w = weights[i];
                if (double.IsNaN(w) || double.IsInfinity(w) || w < 0)
                    throw new ArgumentOutOfRangeException(nameof(weights),
                        "Weights must be finite and non-negative.");
            }
        }
    }
}
