using System;
using Maths.Distance.Geodesic;

namespace Maths.Geometry
{
    /// <summary>
    /// Unit sphere S^(n-1) embedded in R^n. Points are unit vectors (normalized
    /// internally, matching <see cref="SphericalGeodesic"/> / SphericalGeodesicMetric);
    /// tangent vectors at p live in the ambient R^n, orthogonal to p.
    ///
    /// <para><b>Dimension vs intrinsic dimension — the off-by-one ward.</b>
    /// <see cref="Dimension"/> = n = the ambient coordinate count, so it matches the
    /// IRiemannianManifold buffer-length contract (points/tangents/dst are length n).
    /// The INTRINSIC manifold dimension is m = n - 1 = <see cref="IntrinsicDimension"/>;
    /// that is the value the heat-kernel Van Vleck exponent (m-1)/2 and the spherical
    /// calibration consume — NEVER <see cref="Dimension"/>. On the Poincaré ball the two
    /// coincide (a d-coordinate ball is a d-manifold); on the embedded sphere they differ
    /// by one, and passing n where m is expected silently over-counts the curvature.</para>
    ///
    /// <para>On the embedded sphere the induced metric is the ambient Euclidean metric
    /// restricted to T_p, so Norm(p,v) = ||v||_2 directly and ||log_p(q)||_2 = d(p,q)
    /// with no conformal factor — the Poincaré ||log|| = d/λ bug does not recur here.</para>
    /// </summary>
    public readonly struct SphericalManifold : IRiemannianManifold
    {
        public static bool IsFlat => false;

        /// <summary>Ambient coordinate count n (the sphere is S^(n-1)).</summary>
        public int Dimension { get; }

        /// <summary>Intrinsic manifold dimension m = n - 1. Feed THIS to the spherical
        /// Van Vleck / calibration, never <see cref="Dimension"/>.</summary>
        public int IntrinsicDimension => Dimension - 1;

        public SphericalManifold(int ambientDimension)
        {
            if (ambientDimension < 2)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(ambientDimension), ambientDimension,
                    "Spherical manifold requires at least 2 ambient coordinates (S^1).");
            }

            Dimension = ambientDimension;
        }

        public double Distance(ReadOnlySpan<double> p, ReadOnlySpan<double> q)
            => SphericalGeodesic.Distance(p, q);

        public void LogMap(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst)
        {
            ValidateTriple(p, q, dst);

            Span<double> pu = Dimension <= 256 ? stackalloc double[Dimension] : new double[Dimension];
            Span<double> qu = Dimension <= 256 ? stackalloc double[Dimension] : new double[Dimension];
            Normalize(p, pu);
            Normalize(q, qu);

            double dot = Math.Clamp(Dot(pu, qu), -1.0, 1.0);
            double theta = Math.Acos(dot);

            if (theta < 1e-12)
            {
                dst.Clear();                 // p == q
                return;
            }

            double sinTheta = Math.Sin(theta);
            if (sinTheta < 1e-12)
            {
                // p ~= -q : antipode / first conjugate point (the caustic). The geodesic
                // direction is undefined here; refuse to fabricate one. Callers flag
                // near-antipodal pairs as outside the short-time parametrix.
                dst.Clear();
                return;
            }

            // log_p(q) = (theta / sin theta) * (q - <p,q> p); ||result||_2 = theta.
            double scale = theta / sinTheta;
            for (int i = 0; i < Dimension; i++)
                dst[i] = (qu[i] - dot * pu[i]) * scale;
        }

        public void ExpMap(ReadOnlySpan<double> p, ReadOnlySpan<double> v, Span<double> dst)
        {
            ValidateTriple(p, v, dst);

            Span<double> pu = Dimension <= 256 ? stackalloc double[Dimension] : new double[Dimension];
            Normalize(p, pu);

            double theta = EuclideanNorm(v);
            if (theta < 1e-12)
            {
                pu.CopyTo(dst);              // exp_p(0) = p
                return;
            }

            // exp_p(v) = cos(theta) p + sin(theta) v/||v||, with theta = ||v||.
            double cos = Math.Cos(theta);
            double sinOverTheta = Math.Sin(theta) / theta;
            for (int i = 0; i < Dimension; i++)
                dst[i] = (cos * pu[i]) + (sinOverTheta * v[i]);

            NormalizeInPlace(dst);           // guard against numerical drift off the sphere
        }

        public void AddScaled(Span<double> dst, ReadOnlySpan<double> v, double scalar)
        {
            if (dst.Length < Dimension || v.Length < Dimension)
                throw new ArgumentException("Tangent buffers must match the manifold dimension.");

            for (int i = 0; i < Dimension; i++)
                dst[i] += v[i] * scalar;
        }

        public double Norm(ReadOnlySpan<double> p, ReadOnlySpan<double> v)
        {
            if (p.Length < Dimension || v.Length < Dimension)
                throw new ArgumentException("Point and tangent vector must match the manifold dimension.");

            // Induced metric on the embedded sphere is the ambient Euclidean metric —
            // no conformal factor (cf. PoincareBallManifold.Norm).
            return EuclideanNorm(v);
        }

        private void ValidateTriple(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst)
        {
            if (p.Length < Dimension || q.Length < Dimension || dst.Length < Dimension)
                throw new ArgumentException("Point and tangent buffers must match the manifold dimension.");
        }

        private void Normalize(ReadOnlySpan<double> src, Span<double> dst)
        {
            double norm = EuclideanNorm(src);
            if (norm < 1e-12)
            {
                dst.Clear();                 // degenerate zero vector -> first-axis pole
                dst[0] = 1.0;
                return;
            }

            for (int i = 0; i < Dimension; i++)
                dst[i] = src[i] / norm;
        }

        private void NormalizeInPlace(Span<double> v)
        {
            double norm = EuclideanNorm(v);
            if (norm < 1e-12)
            {
                v.Clear();
                v[0] = 1.0;
                return;
            }

            for (int i = 0; i < Dimension; i++)
                v[i] /= norm;
        }

        private static double EuclideanNorm(ReadOnlySpan<double> v)
            => Math.Sqrt(Dot(v, v));

        private static double Dot(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
        {
            double sum = 0.0;
            for (int i = 0; i < a.Length; i++)
                sum += a[i] * b[i];
            return sum;
        }
    }
}
