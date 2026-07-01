// ============================================================================
// Geometry/ProductManifold.cs
// ============================================================================
using System;

namespace Maths.Geometry
{
    public readonly struct ProductManifold<TA, TB> : IRiemannianManifold
        where TA : struct, IRiemannianManifold
        where TB : struct, IRiemannianManifold
    {
        public static bool IsFlat => TA.IsFlat && TB.IsFlat;

        private readonly TA _a;
        private readonly TB _b;
        private readonly int _dimA;
        public int Dimension { get; }

        public ProductManifold(TA a, TB b)
        {
            _a = a;
            _b = b;
            _dimA = a.Dimension;
            Dimension = a.Dimension + b.Dimension;
        }

        public double Distance(ReadOnlySpan<double> p, ReadOnlySpan<double> q)
        {
            double dA = _a.Distance(p.Slice(0, _dimA), q.Slice(0, _dimA));
            double dB = _b.Distance(p.Slice(_dimA), q.Slice(_dimA));
            return Math.Sqrt(dA * dA + dB * dB);
        }

        public void LogMap(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst)
        {
            _a.LogMap(p.Slice(0, _dimA), q.Slice(0, _dimA), dst.Slice(0, _dimA));
            _b.LogMap(p.Slice(_dimA), q.Slice(_dimA), dst.Slice(_dimA));
        }

        public void ExpMap(ReadOnlySpan<double> p, ReadOnlySpan<double> v, Span<double> dst)
        {
            _a.ExpMap(p.Slice(0, _dimA), v.Slice(0, _dimA), dst.Slice(0, _dimA));
            _b.ExpMap(p.Slice(_dimA), v.Slice(_dimA), dst.Slice(_dimA));
        }

        public void AddScaled(Span<double> dst, ReadOnlySpan<double> v, double scalar)
        {
            for (int i = 0; i < Dimension; i++) dst[i] += v[i] * scalar;
        }

        public double Norm(ReadOnlySpan<double> p, ReadOnlySpan<double> v)
        {
            double a = _a.Norm(p.Slice(0, _dimA), v.Slice(0, _dimA));
            double b = _b.Norm(p.Slice(_dimA), v.Slice(_dimA));
            return Math.Sqrt(a * a + b * b);
        }
    }
}
