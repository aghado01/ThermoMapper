using System;
using System.Buffers;
using Maths.Geometry;
using Maths.Geometry.Solver;

namespace Maths.Geometry.Estimators.Intrinsic
// ============================================================================
// Estimators/Manifolds/ManifoldMedian.cs
// ============================================================================
// Hot-path geometric median via hybrid Weiszfeld / projected subgradient IRLS.
// Convergence guaranteed for non-degenerate configurations; singularity at
// coincident data points handled via OptimalityCheck or Regularise policy
// (configured in IrlsOptions).
//
// Downstream callers (GMM init, SPC handoff, etc.) call this directly.
// Zero heap allocation on the hot path — all scratch via ArrayPool.
// ============================================================================

{
    public static class GeometricMedian
    {
        // ── Location only ─────────────────────────────────────────────────────

        /// <summary>
        /// Computes the geometric median in place.
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
            Irls.Solve<TManifold, L1Loss>(
                manifold, data, weights, destination, opts,
                finalIrlsWeights: default);
        }

        // ── Warm-start selection ──────────────────────────────────────────────

        /// <summary>
        /// Index of the medoid: the data point minimizing the total distance to all others.
        /// The natural warm-start for <see cref="Compute{TManifold}"/> — a sample inside
        /// the majority cluster keeps the IRLS iteration within the manifold's injectivity radius,
        /// whereas an arbitrary point can sit on the cut locus of the rest, where the log map
        /// degenerates and Weiszfeld stalls. O(n²) distance evaluations; ties resolve to the
        /// lowest index.
        /// </summary>
        public static int MedoidIndex<TManifold>(
            TManifold              manifold,
            ReadOnlySpan<double[]> data)
            where TManifold : struct, IRiemannianManifold
        {
            int      n      = data.Length;
            double[] rented = ArrayPool<double>.Shared.Rent(n);
            try
            {
                Span<double> totals = rented.AsSpan(0, n);
                totals.Clear();
                for (int i = 0; i < n; i++)
                {
                    for (int j = i + 1; j < n; j++)
                    {
                        double distance = manifold.Distance(data[i], data[j]);
                        totals[i] += distance;
                        totals[j] += distance;
                    }
                }

                int medoid = 0;
                for (int i = 1; i < n; i++)
                {
                    if (totals[i] < totals[medoid])
                        medoid = i;
                }
                return medoid;
            }
            finally
            {
                ArrayPool<double>.Shared.Return(rented);
            }
        }

        // ── (Location, Scatter) ───────────────────────────────────────────────

        /// <summary>
        /// Computes the geometric median and the Weiszfeld (L1-weighted) scatter
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
                Irls.Solve<TManifold, L1Loss>(
                    manifold, data, weights, locationDestination, opts, finalW);
                WeiszfeldScatter.Accumulate(
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
