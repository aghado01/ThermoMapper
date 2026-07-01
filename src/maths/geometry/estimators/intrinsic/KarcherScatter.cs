// ============================================================================
// Estimators/Manifolds/KarcherScatter.cs
// ============================================================================
// Riemannian sample covariance — the scatter companion to GeometricMean.
// Weights are the converged L2 IRLS weights (uniform for the mean, non-uniform
// when outer sample weights are supplied).
//
// Sigma = (cD / sum_i w_i) * sum_i w_i * v_i * v_i^T
// where v_i = log_mu(x_i).
//
// Can also be called standalone with a fixed location estimate (e.g. when
// the mean was computed externally and only scatter is needed).
// ============================================================================
using System;
using System.Buffers;
using Maths.Geometry;
namespace Maths.Geometry.Estimators.Intrinsic
{
    public static class KarcherScatter
    {
        /// <summary>
        /// Accumulates scatter from converged IRLS weights into
        /// <paramref name="scatterDst"/> (row-major D×D, length D²).
        /// Called internally by <see cref="GeometricMean.ComputeWithScatter"/>.
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
        /// Standalone scatter: computes the Riemannian covariance around a
        /// fixed <paramref name="location"/> without re-running mean estimation.
        /// Outer <paramref name="weights"/> are used directly (uniform if default).
        /// </summary>
        public static void Compute<TManifold>(
            TManifold              manifold,
            ReadOnlySpan<double[]> data,
            double[]               location,
            Span<double>           scatterDestination,
            ReadOnlySpan<double>   weights = default,
            double                 consistencyFactor = 1.0)
            where TManifold : struct, IRiemannianManifold
        {
            int      n       = data.Length;
            double[] rentedW = ArrayPool<double>.Shared.Rent(n);
            try
            {
                Span<double> w = rentedW.AsSpan(0, n);
                if (weights.IsEmpty)
                    w.Fill(1.0);
                else
                    weights.CopyTo(w);

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
