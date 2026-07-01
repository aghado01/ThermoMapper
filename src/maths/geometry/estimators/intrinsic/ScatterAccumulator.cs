// ============================================================================
// Estimators/Shared/ScatterAccumulator.cs   (internal primitive)
// ============================================================================
// Single implementation of the weighted outer-product scatter accumulation
// used by both KarcherScatter and WeiszfeldScatter.
//
// Sigma = (consistencyFactor / sum_i w_i) * sum_i w_i * v_i * v_i^T
// where v_i = log_mu(x_i).
//
// scatterDst is a caller-supplied row-major D×D flat buffer (length D²).
// ============================================================================
using System;
using System.Buffers;
using Maths.Geometry;

namespace Maths.Geometry.Estimators.Intrinsic
{
    internal static class ScatterAccumulator
    {
        public static void Accumulate<TManifold>(
            TManifold              manifold,
            ReadOnlySpan<double[]> data,
            double[]               location,
            ReadOnlySpan<double>   weights,
            Span<double>           scatterDst,
            double                 consistencyFactor)
            where TManifold : struct, IRiemannianManifold
        {
            int      dim    = manifold.Dimension;
            double[] logArr = ArrayPool<double>.Shared.Rent(dim);
            try
            {
                Span<double>         v       = logArr.AsSpan(0, dim);
                ReadOnlySpan<double> locSpan = location;
                scatterDst.Clear();
                double sumW = 0;

                for (int i = 0; i < data.Length; i++)
                {
                    double w = weights[i];
                    if (w == 0) continue;
                    sumW += w;
                    manifold.LogMap(locSpan, data[i], v);
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
                ArrayPool<double>.Shared.Return(logArr);
            }
        }
    }
}
