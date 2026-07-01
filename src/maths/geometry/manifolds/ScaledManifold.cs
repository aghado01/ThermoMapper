using System;

namespace Maths.Geometry
{
    /// <summary>
    /// Metric-scaling wrapper: equips an inner manifold with the conformally rescaled metric
    /// g' = scale·g. Scaling a Riemannian metric by a positive constant leaves the geodesics — and
    /// hence Exp and Log — unchanged (Kisung You, "Constant Metric Scaling in Riemannian Computation",
    /// arXiv 2601.10992); only lengths rescale, so Distance and Norm pick up a √scale factor.
    ///
    /// <para>Its purpose is to weight the factors of a <see cref="ProductManifold{TA,TB}"/>: with
    /// factor scales α and 2−α the product distance becomes √(α·d_A² + (2−α)·d_B²), exactly the
    /// scale-calibrated metric of robust distributed PCA (You 2605.20681 §5).</para>
    /// </summary>
    public readonly struct ScaledManifold<TInner> : IRiemannianManifold
        where TInner : struct, IRiemannianManifold
    {
        private readonly TInner _inner;
        private readonly double _sqrtScale;

        public ScaledManifold(TInner inner, double scale)
        {
            if (!double.IsFinite(scale) || scale <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(scale), scale, "Metric scale must be finite and positive.");
            _inner = inner;
            _sqrtScale = Math.Sqrt(scale);
        }

        public static bool IsFlat => TInner.IsFlat;
        public int Dimension => _inner.Dimension;

        public double Distance(ReadOnlySpan<double> p, ReadOnlySpan<double> q)
            => _sqrtScale * _inner.Distance(p, q);

        public double Norm(ReadOnlySpan<double> p, ReadOnlySpan<double> v)
            => _sqrtScale * _inner.Norm(p, v);

        // Geodesics (hence Log/Exp) are invariant under constant metric scaling — delegate verbatim.
        public void LogMap(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst)
            => _inner.LogMap(p, q, dst);

        public void ExpMap(ReadOnlySpan<double> p, ReadOnlySpan<double> v, Span<double> dst)
            => _inner.ExpMap(p, v, dst);

        public void AddScaled(Span<double> dst, ReadOnlySpan<double> v, double scalar)
            => _inner.AddScaled(dst, v, scalar);
    }
}
