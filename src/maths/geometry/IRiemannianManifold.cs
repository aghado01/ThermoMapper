// ============================================================================
// Manifolds/IRiemannianManifold.cs
// ============================================================================
// Decision: non-generic, span-based throughout.
// ProductManifold already implements this shape. EuclideanVectorManifold is
// updated to match (see EuclideanVectorManifold.cs).
// Points are represented as double[] (or ReadOnlySpan<double> at call sites).
// ============================================================================
using System;

namespace Maths.Geometry
{
    public interface IRiemannianManifold
    {
        /// <summary>
        /// Compile-time flat-space flag. JIT erases curved-path dead branches
        /// for EuclideanVectorManifold specialisations.
        /// </summary>
        static abstract bool IsFlat { get; }

        /// <summary>
        /// Tangent-space dimension D. For product manifolds, the sum of
        /// per-factor dimensions.
        /// </summary>
        int Dimension { get; }

        /// <summary>Geodesic distance d(p, q).</summary>
        double Distance(ReadOnlySpan<double> p, ReadOnlySpan<double> q);

        /// <summary>
        /// Writes log_p(q) into <paramref name="dst"/> (length = Dimension).
        /// </summary>
        void LogMap(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst);

        /// <summary>
        /// Writes exp_p(v) into <paramref name="dst"/> (length = Dimension).
        /// </summary>
        void ExpMap(ReadOnlySpan<double> p, ReadOnlySpan<double> v, Span<double> dst);

        /// <summary>dst += scalar * v  (in-place tangent accumulation).</summary>
        void AddScaled(Span<double> dst, ReadOnlySpan<double> v, double scalar);

        /// <summary>
        /// Riemannian norm ||v||_p. For flat manifolds this is just the
        /// Euclidean norm of the tangent vector; for curved manifolds it equals
        /// the geodesic length of the corresponding exp_p(v) path.
        /// </summary>
        double Norm(ReadOnlySpan<double> p, ReadOnlySpan<double> v);
    }
}
