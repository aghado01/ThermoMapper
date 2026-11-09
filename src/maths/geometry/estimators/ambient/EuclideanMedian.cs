using System;
using System.Buffers;
using Maths.Geometry.Solver;

namespace Maths.Geometry.Estimators.Ambient
// ============================================================================
// Estimators/Euclidean/EuclideanMedian.cs
// ============================================================================
// Flat-space specialisation of GeometricMedian — regularized Weiszfeld IRLS.
// Drop-in for GeometricMedian when you don't need the manifold abstraction.
//
// Output shapes mirror GeometricMedian exactly:
//   Compute(...)             -> void, destination mutated in place
//   ComputeWithScatter(...)  -> void, destination + scatterDestination mutated
//
// Warm-start: caller pre-populates destination (typically from EuclideanMean).
// tolerance and epsilon are read in the unit that `scale` picks: by default the
// spread σ, the weighted median distance from the iterate to the data, recomputed
// each iteration (IrlsScale.Data keeps them in the data's units).
// scatterDestination is a row-major D*D flat buffer (length D²),
// identical layout to GeometricMedian.ComputeWithScatter.
// ============================================================================

{
    public static class EuclideanMedian
    {
        // ── Location only ─────────────────────────────────────────────────────

        /// <summary>
        /// Computes the weighted geometric median into <paramref name="destination"/>.
        /// <paramref name="destination"/> is the warm-start initial value and the output.
        /// </summary>
        public static void Compute(
            double[][] points,
            ReadOnlySpan<int> indices,
            ReadOnlySpan<double> weights,
            Span<double> destination,
            int maxIterations = 200,
            double tolerance = 1e-8,
            double epsilon = 1e-10,
            IrlsScale scale = IrlsScale.Spread)
        {
            EuclideanMean.ValidateInputs(
                points, indices, weights, destination, out int dim);
            ValidateIterationParams(maxIterations, tolerance, epsilon, scale);

            double[] scratchArr = ArrayPool<double>.Shared.Rent(dim);
            double[] spreadArr = ArrayPool<double>.Shared.Rent(3 * indices.Length);
            try
            {
                IrlsLoop(points, indices, weights, destination,
                         scratchArr.AsSpan(0, dim), spreadArr.AsSpan(0, 3 * indices.Length),
                         dim, maxIterations, tolerance, epsilon, scale);
            }
            finally
            {
                ArrayPool<double>.Shared.Return(scratchArr);
                ArrayPool<double>.Shared.Return(spreadArr);
            }
        }

        // ── (Location, Scatter) ───────────────────────────────────────────────

        /// <summary>
        /// Computes the weighted geometric median and Weiszfeld scatter in one pass.
        /// <para>
        /// <paramref name="destination"/> is the warm-start and the location output.
        /// <paramref name="scatterDestination"/> is a row-major D×D flat buffer
        /// of length D². Layout identical to <c>GeometricMedian.ComputeWithScatter</c>.
        /// </para>
        /// </summary>
        public static void ComputeWithScatter(
            double[][] points,
            ReadOnlySpan<int> indices,
            ReadOnlySpan<double> weights,
            Span<double> destination,
            Span<double> scatterDestination,
            int maxIterations = 200,
            double tolerance = 1e-8,
            double epsilon = 1e-10,
            double consistencyFactor = 1.0,
            IrlsScale scale = IrlsScale.Spread)
        {
            EuclideanMean.ValidateInputs(
                points, indices, weights, destination, out int dim);
            ValidateIterationParams(maxIterations, tolerance, epsilon, scale);

            double[] scratchArr = ArrayPool<double>.Shared.Rent(dim);
            double[] irlsWArr = ArrayPool<double>.Shared.Rent(indices.Length);
            double[] spreadArr = ArrayPool<double>.Shared.Rent(3 * indices.Length);
            try
            {
                Span<double> scratch = scratchArr.AsSpan(0, dim);
                Span<double> irlsW = irlsWArr.AsSpan(0, indices.Length);

                IrlsLoop(points, indices, weights, destination,
                         scratch, spreadArr.AsSpan(0, 3 * indices.Length),
                         dim, maxIterations, tolerance, epsilon, scale,
                         finalIrlsWeights: irlsW);

                AccumulateScatter(points, indices, irlsW, destination,
                                  scatterDestination, dim, consistencyFactor);
            }
            finally
            {
                ArrayPool<double>.Shared.Return(scratchArr);
                ArrayPool<double>.Shared.Return(irlsWArr);
                ArrayPool<double>.Shared.Return(spreadArr);
            }
        }

        // ── IRLS core ─────────────────────────────────────────────────────────

        private static void IrlsLoop(
            double[][] points,
            ReadOnlySpan<int> indices,
            ReadOnlySpan<double> weights,
            Span<double> destination,
            Span<double> scratch,
            Span<double> spreadScratch,
            int dim,
            int maxIterations,
            double tolerance,
            double epsilon,
            IrlsScale scale,
            Span<double> finalIrlsWeights = default)
        {
            if (scale == IrlsScale.Spread)
            {
                SpreadLoop(points, indices, weights, destination, scratch, spreadScratch,
                           dim, maxIterations, tolerance, epsilon, finalIrlsWeights);
                return;
            }

            for (int iter = 0; iter < maxIterations; iter++)
            {
                scratch.Clear();
                double denominator = 0;

                for (int i = 0; i < indices.Length; i++)
                {
                    double w = weights[i];
                    if (w == 0) continue;

                    double[] pt = points[indices[i]];
                    double r = EuclideanDistance(pt, destination, dim);
                    double wi = w / Math.Max(r, epsilon);

                    denominator += wi;
                    for (int d = 0; d < dim; d++) scratch[d] += wi * pt[d];
                }

                if (denominator <= 0) break;

                double inv = 1.0 / denominator;
                double shift = 0;
                double pNormSq = 0;

                for (int d = 0; d < dim; d++)
                {
                    double next = scratch[d] * inv;
                    double diff = next - destination[d];
                    shift += diff * diff;
                    pNormSq += destination[d] * destination[d];
                    scratch[d] = next;
                }

                scratch.CopyTo(destination);

                double normScale = 1.0 + Math.Sqrt(pNormSq);
                if (Math.Sqrt(shift) <= tolerance * normScale) break;
            }

            // Recompute weights at the converged location so scatter reflects
            // the final position, not the pre-update position of the last step.
            if (!finalIrlsWeights.IsEmpty)
            {
                for (int i = 0; i < indices.Length; i++)
                {
                    double w = weights[i];
                    if (w == 0) { finalIrlsWeights[i] = 0.0; continue; }
                    double r = EuclideanDistance(points[indices[i]], destination, dim);
                    finalIrlsWeights[i] = w / Math.Max(r, epsilon);
                }
            }
        }

        // IrlsScale.Spread: the floor is ε·σ and the stop τ·σ, with σ the spread
        // (Irls.SpreadOf) at the current iterate. The step accumulates the differences
        // xᵢ − p, so a centre far from the origin costs only the rounding of p, and the
        // stopping radius never falls below Irls.RoundingFloor·‖p‖.
        private static void SpreadLoop(
            double[][] points,
            ReadOnlySpan<int> indices,
            ReadOnlySpan<double> weights,
            Span<double> destination,
            Span<double> step,
            Span<double> spreadScratch,
            int dim,
            int maxIterations,
            double tolerance,
            double epsilon,
            Span<double> finalIrlsWeights)
        {
            int n = indices.Length;
            Span<double> radii = spreadScratch.Slice(0, n);
            Span<double> keys = spreadScratch.Slice(n, n);
            Span<double> keyW = spreadScratch.Slice(2 * n, n);

            for (int iter = 0; iter < maxIterations; iter++)
            {
                double spread = Radii(points, indices, weights, destination, radii, keys, keyW, dim);
                if (spread == 0) break;   // every weighted point sits on the iterate
                double floor = epsilon * spread;

                step.Clear();
                double denominator = 0;
                for (int i = 0; i < n; i++)
                {
                    double w = weights[i];
                    if (w == 0) continue;

                    double[] pt = points[indices[i]];
                    double wi = w / Math.Max(radii[i], floor);
                    denominator += wi;
                    for (int d = 0; d < dim; d++) step[d] += wi * (pt[d] - destination[d]);
                }

                if (denominator <= 0) break;

                double inv = 1.0 / denominator;
                double shift = 0;
                double pNormSq = 0;
                for (int d = 0; d < dim; d++)
                {
                    double s = step[d] * inv;
                    shift += s * s;
                    destination[d] += s;
                    pNormSq += destination[d] * destination[d];
                }

                double radius = Math.Max(tolerance * spread, Irls.RoundingFloor * Math.Sqrt(pNormSq));
                if (Math.Sqrt(shift) <= radius) break;
            }

            if (!finalIrlsWeights.IsEmpty)
            {
                // With σ = 0 every weighted point sits on the location and any finite weights give a zero scatter.
                double spread = Radii(points, indices, weights, destination, radii, keys, keyW, dim);
                double floor = epsilon * spread;
                for (int i = 0; i < n; i++)
                {
                    double w = weights[i];
                    finalIrlsWeights[i] = w == 0 ? 0.0 : spread == 0 ? w : w / Math.Max(radii[i], floor);
                }
            }
        }

        // Distances from the iterate to the weighted points into radii; returns the spread.
        private static double Radii(
            double[][] points,
            ReadOnlySpan<int> indices,
            ReadOnlySpan<double> weights,
            ReadOnlySpan<double> center,
            Span<double> radii,
            Span<double> keys,
            Span<double> keyW,
            int dim)
        {
            for (int i = 0; i < indices.Length; i++)
                radii[i] = weights[i] == 0 ? 0.0 : EuclideanDistance(points[indices[i]], center, dim);
            return Irls.SpreadOf(radii, weights, keys, keyW);
        }

        // ── Weiszfeld scatter accumulation ────────────────────────────────────
        // Sigma = (cD / sum_i w_i) * sum_i w_i * (x_i - mu)(x_i - mu)^T
        private static void AccumulateScatter(
            double[][] points,
            ReadOnlySpan<int> indices,
            ReadOnlySpan<double> irlsWeights,
            ReadOnlySpan<double> location,
            Span<double> scatterDst,
            int dim,
            double consistencyFactor)
        {
            scatterDst.Clear();
            double sumW = 0;

            double[] vArr = ArrayPool<double>.Shared.Rent(dim);
            try
            {
                Span<double> v = vArr.AsSpan(0, dim);
                for (int i = 0; i < indices.Length; i++)
                {
                    double w = irlsWeights[i];
                    if (w == 0) continue;
                    sumW += w;
                    double[] pt = points[indices[i]];
                    for (int d = 0; d < dim; d++) v[d] = pt[d] - location[d];
                    for (int r = 0; r < dim; r++)
                        for (int c = 0; c < dim; c++)
                            scatterDst[r * dim + c] += w * v[r] * v[c];
                }

                if (sumW > 0)
                {
                    double scale = consistencyFactor / sumW;
                    for (int k = 0; k < scatterDst.Length; k++)
                        scatterDst[k] *= scale;
                }
            }
            finally
            {
                ArrayPool<double>.Shared.Return(vArr);
            }
        }

        // ── Helpers ───────────────────────────────────────────────────────────

        private static double EuclideanDistance(
            double[] point,
            ReadOnlySpan<double> center,
            int dim)
        {
            double sq = 0;
            for (int d = 0; d < dim; d++)
            {
                double diff = point[d] - center[d];
                sq += diff * diff;
            }
            return Math.Sqrt(sq);
        }

        private static void ValidateIterationParams(
            int maxIterations, double tolerance, double epsilon, IrlsScale scale)
        {
            if (scale != IrlsScale.Spread && scale != IrlsScale.Data)
                throw new ArgumentOutOfRangeException(nameof(scale));
            if (maxIterations < 1)
                throw new ArgumentOutOfRangeException(nameof(maxIterations));
            if (!(tolerance > 0) || double.IsNaN(tolerance) || double.IsInfinity(tolerance))
                throw new ArgumentOutOfRangeException(nameof(tolerance));
            if (!(epsilon > 0) || double.IsNaN(epsilon) || double.IsInfinity(epsilon))
                throw new ArgumentOutOfRangeException(nameof(epsilon));
        }
    }
}
