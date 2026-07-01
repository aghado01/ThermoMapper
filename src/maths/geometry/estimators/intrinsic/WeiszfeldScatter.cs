// ============================================================================
// Estimators/Manifolds/WeiszfeldScatter.cs
// ============================================================================
// Robust scatter companion to GeometricMedian.
// Weights are the converged L1 IRLS weights w_i = 1 / d(mu, x_i), so points
// far from the median contribute less to the scatter estimate — consistent with
// the median's breakdown-point guarantees.
//
// Sigma = (cD / sum_i w_i) * sum_i w_i * v_i * v_i^T
// where v_i = log_mu(x_i) and w_i are the L1 IRLS weights at convergence.
//
// Can also be called standalone with a fixed location estimate.
// ============================================================================
using System;
using System.Buffers;
using Maths.Geometry;
namespace Maths.Geometry.Estimators.Intrinsic
{
    public static class WeiszfeldScatter
    {
        /// <summary>
        /// Accumulates scatter from converged IRLS weights into
        /// <paramref name="scatterDst"/> (row-major D×D, length D²).
        /// Called internally by <see cref="GeometricMedian.ComputeWithScatter"/>.
        /// </summary>
        internal static void Accumulate<TManifold>(
            TManifold              manifold,
            ReadOnlySpan<double[]> data,
            double[]               location,
            ReadOnlySpan<double>   finalIrlsWeights,
            Span<double>           scatterDst,
            double                 consistencyFactor)
            where TManifold : struct, IRiemannianManifold
            => ScatterAccumulator.Accumulate(
                manifold, data, location,
                finalIrlsWeights, scatterDst, consistencyFactor);

        /// <summary>
        /// Standalone scatter: computes robust Weiszfeld scatter around a fixed
        /// <paramref name="location"/> without re-running median estimation.
        /// Runs one pass computing L1 weights from distances, then accumulates.
        /// </summary>
        public static void Compute<TManifold>(
            TManifold              manifold,
            ReadOnlySpan<double[]> data,
            double[]               location,
            Span<double>           scatterDestination,
            ReadOnlySpan<double>   outerWeights = default,
            double                 epsilon = 1e-10,
            double                 consistencyFactor = 1.0)
            where TManifold : struct, IRiemannianManifold
        {
            int      n       = data.Length;
            double[] rentedW = ArrayPool<double>.Shared.Rent(n);
            try
            {
                Span<double>         w       = rentedW.AsSpan(0, n);
                ReadOnlySpan<double> locSpan = location;

                for (int i = 0; i < n; i++)
                {
                    double r  = manifold.Distance(locSpan, data[i]);
                    double ow = outerWeights.IsEmpty ? 1.0 : outerWeights[i];
                    w[i] = r > epsilon ? ow / r : 0.0;
                }

                ScatterAccumulator.Accumulate(
                    manifold, data, location, w,
                    scatterDestination, consistencyFactor);
            }
            finally
            {
                ArrayPool<double>.Shared.Return(rentedW);
            }
        }
    }
}
