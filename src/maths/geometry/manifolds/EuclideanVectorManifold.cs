// ============================================================================
// Manifolds/EuclideanVectorManifold.cs
// ============================================================================
// Flat Euclidean R^d.  Points are passed as ReadOnlySpan<double> / Span<double>
// at all call sites.  The 1D case is EuclideanVectorManifold(dimension: 1).
// ============================================================================
using System;
using System.Numerics.Tensors;

namespace Maths.Geometry
{
    public readonly struct EuclideanVectorManifold : IRiemannianManifold
    {
        public static bool IsFlat => true;
        public int Dimension { get; }

        public EuclideanVectorManifold(int dimension) => Dimension = dimension;

        public double Distance(ReadOnlySpan<double> p, ReadOnlySpan<double> q)
            => TensorPrimitives.Distance<double>(p, q);

        public void LogMap(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst)
            => TensorPrimitives.Subtract<double>(q, p, dst);

        public void ExpMap(ReadOnlySpan<double> p, ReadOnlySpan<double> v, Span<double> dst)
            => TensorPrimitives.Add<double>(p, v, dst);

        public void AddScaled(Span<double> dst, ReadOnlySpan<double> v, double scalar)
        {
            for (int i = 0; i < Dimension; i++) dst[i] += v[i] * scalar;
        }

        public double Norm(ReadOnlySpan<double> p, ReadOnlySpan<double> v)
        {
            double sumSq = 0;
            for (int i = 0; i < Dimension; i++) sumSq += v[i] * v[i];
            return Math.Sqrt(sumSq);
        }
    }
}
