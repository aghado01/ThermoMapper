using System;
using Maths.Distance.Geodesic;

namespace Maths.Geometry
{
    /// <summary>
    /// Hyperboloid (Lorentz) model of hyperbolic space H^d embedded in R^(d+1).
    /// Points are vectors on the upper sheet of the hyperboloid satisfying
    /// &lt;p, p&gt;_L = -1 and p_{N-1} > 0.
    /// Tangent vectors at p satisfy &lt;p, v&gt;_L = 0.
    ///
    /// <para><b>Dimension vs intrinsic dimension:</b>
    /// <see cref="Dimension"/> = n = the ambient coordinate count, matching the
    /// IRiemannianManifold buffer-length contract.
    /// The INTRINSIC manifold dimension is d = n - 1 = <see cref="IntrinsicDimension"/>.</para>
    /// </summary>
    public readonly struct LorentzHyperboloidManifold : IRiemannianManifold
    {
        public static bool IsFlat => false;

        /// <summary>Ambient coordinate count n (the hyperboloid is in R^n).</summary>
        public int Dimension { get; }

        /// <summary>Intrinsic manifold dimension d = n - 1.</summary>
        public int IntrinsicDimension => Dimension - 1;

        public LorentzHyperboloidManifold(int ambientDimension)
        {
            if (ambientDimension < 2)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(ambientDimension), ambientDimension,
                    "Lorentz hyperboloid manifold requires at least 2 ambient coordinates (H^1).");
            }

            Dimension = ambientDimension;
        }

        public double Distance(ReadOnlySpan<double> p, ReadOnlySpan<double> q)
            => Lorentz.Distance(p, q);

        public void LogMap(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst)
        {
            ValidateTriple(p, q, dst);

            double alpha = -Lorentz.InnerProduct(p, q);
            if (alpha <= 1.0 + 1e-12)
            {
                dst.Clear(); // p == q
                return;
            }

            double dist = Math.Log(alpha + Math.Sqrt(alpha * alpha - 1.0));
            double denom = Math.Sqrt(alpha * alpha - 1.0);
            double scale = dist / denom;

            for (int i = 0; i < Dimension; i++)
            {
                dst[i] = scale * (q[i] - alpha * p[i]);
            }
        }

        public void ExpMap(ReadOnlySpan<double> p, ReadOnlySpan<double> v, Span<double> dst)
        {
            ValidateTriple(p, v, dst);

            double tangentNorm = Norm(p, v);
            if (tangentNorm < 1e-12)
            {
                p.CopyTo(dst); // exp_p(0) = p
                return;
            }

            double cosh = Math.Cosh(tangentNorm);
            double sinhOverNorm = Math.Sinh(tangentNorm) / tangentNorm;

            for (int i = 0; i < Dimension; i++)
            {
                dst[i] = (cosh * p[i]) + (sinhOverNorm * v[i]);
            }

            // Project/normalize back to the hyperboloid to prevent numerical drift
            Lorentz.ProjectToHyperboloid(dst);
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

            double ip = Lorentz.InnerProduct(v, v);
            return Math.Sqrt(Math.Max(0.0, ip));
        }

        /// <summary>
        /// Projects an ambient vector u in R^n onto the tangent space T_p H^d:
        /// Proj_p(u) = u + &lt;p, u&gt;_L * p.
        /// </summary>
        public void ProjectTangent(ReadOnlySpan<double> p, ReadOnlySpan<double> u, Span<double> dst)
        {
            if (p.Length < Dimension || u.Length < Dimension || dst.Length < Dimension)
                throw new ArgumentException("Buffers must match the manifold dimension.");

            double ip = Lorentz.InnerProduct(p, u);
            for (int i = 0; i < Dimension; i++)
            {
                dst[i] = u[i] + ip * p[i];
            }
        }

        private void ValidateTriple(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst)
        {
            if (p.Length < Dimension || q.Length < Dimension || dst.Length < Dimension)
                throw new ArgumentException("Point and tangent buffers must match the manifold dimension.");
        }
    }
}
