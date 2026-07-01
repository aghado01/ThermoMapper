using System;
using System.Buffers;
using Maths.Geometry;
using Maths.Geometry.Solver;

namespace Maths.Geometry.Estimators.Intrinsic
// ============================================================================
// Estimators/Manifolds/ManifoldMean.cs
// ============================================================================
// Hot-path Riemannian (Fréchet) mean via Karcher flow.
// For flat manifolds: single-pass weighted average (no iteration).
// For curved manifolds: iterative tangent-space round-trip.
//
// Downstream callers (GMM init, SPC handoff, etc.) call this directly.
// Zero heap allocation on the hot path — all scratch via ArrayPool.
// ============================================================================

{
    public static class GeometricMean
    {
        // ── Location only ─────────────────────────────────────────────────────

        /// <summary>
        /// Computes the Riemannian (Fréchet) mean in place.
        /// <paramref name="destination"/> is both the warm-start and the output.
        /// </summary>
        public static void Compute<TManifold>(
            TManifold              manifold,
            ReadOnlySpan<double[]> data,
            ReadOnlySpan<double>   weights,
            double[]               destination,
            IrlsOptions            opts = default)
            where TManifold : struct, IRiemannianManifold
        {
            // L2Loss.IsClosedForm == true → Irls takes the Karcher short-circuit path.
            Irls.Solve<TManifold, L2Loss>(
                manifold, data, weights, destination, opts,
                finalIrlsWeights: default);
        }

        // ── (Location, Scatter) ───────────────────────────────────────────────

        /// <summary>
        /// Computes the Riemannian mean and the Karcher (L2-weighted) scatter
        /// in a single pass.
        /// <para>
        /// <paramref name="scatterDestination"/> must be a caller-supplied
        /// row-major D×D flat buffer of length D².
        /// </para>
        /// </summary>
        public static void ComputeWithScatter<TManifold>(
            TManifold              manifold,
            ReadOnlySpan<double[]> data,
            ReadOnlySpan<double>   weights,
            double[]               locationDestination,
            Span<double>           scatterDestination,
            IrlsOptions            opts = default,
            double                 consistencyFactor = 1.0)
            where TManifold : struct, IRiemannianManifold
        {
            int      n       = data.Length;
            double[] rentedW = ArrayPool<double>.Shared.Rent(n);
            try
            {
                Span<double> finalW = rentedW.AsSpan(0, n);
                Irls.Solve<TManifold, L2Loss>(
                    manifold, data, weights, locationDestination, opts, finalW);
                KarcherScatter.Accumulate(
                    manifold, data, locationDestination, finalW,
                    scatterDestination, consistencyFactor);
            }
            finally
            {
                ArrayPool<double>.Shared.Return(rentedW);
            }
        }
    }
}
