using System;
using Maths.Distance.Geodesic;

namespace Maths.Geometry
{
    public readonly struct PoincareBallManifold : IRiemannianManifold
    {
        public static bool IsFlat => false;

        public int Dimension { get; }

        public PoincareBallManifold(int dimension)
        {
            if (dimension < 1)
                throw new ArgumentOutOfRangeException(nameof(dimension), dimension, "Dimension must be positive.");

            Dimension = dimension;
        }

        public double Distance(ReadOnlySpan<double> p, ReadOnlySpan<double> q)
            => Poincare.Distance(p, q);

        public void LogMap(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst)
        {
            ValidatePointPair(p, q, dst);

            Span<double> delta = Dimension <= 256 ? stackalloc double[Dimension] : new double[Dimension];
            NegatedMobiusAdd(p, q, delta);

            double deltaNorm = EuclideanNorm(delta);
            if (deltaNorm < 1e-12)
            {
                dst.Clear();
                return;
            }

            double scale = (2.0 / ConformalFactor(p)) * Atanh(ClampUnitInterval(deltaNorm)) / deltaNorm;
            for (int i = 0; i < Dimension; i++)
                dst[i] = delta[i] * scale;
        }

        public void ExpMap(ReadOnlySpan<double> p, ReadOnlySpan<double> v, Span<double> dst)
        {
            ValidatePointPair(p, v, dst);

            double tangentNorm = EuclideanNorm(v);
            if (tangentNorm < 1e-12)
            {
                p.CopyTo(dst);
                return;
            }

            double scale = Math.Tanh(ConformalFactor(p) * tangentNorm / 2.0) / tangentNorm;
            Span<double> step = Dimension <= 256 ? stackalloc double[Dimension] : new double[Dimension];
            for (int i = 0; i < Dimension; i++)
                step[i] = v[i] * scale;

            MobiusAdd(p, step, dst);
            Poincare.ClampToBall(dst);
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

            return ConformalFactor(p) * EuclideanNorm(v);
        }

        private void ValidatePointPair(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst)
        {
            if (p.Length < Dimension || q.Length < Dimension || dst.Length < Dimension)
            {
                throw new ArgumentException(
                    "Point and tangent buffers must match the manifold dimension.");
            }
        }

        private static void NegatedMobiusAdd(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst)
        {
            Span<double> negP = p.Length <= 256 ? stackalloc double[p.Length] : new double[p.Length];
            for (int i = 0; i < p.Length; i++)
                negP[i] = -p[i];

            MobiusAdd(negP, q, dst);
        }

        private static void MobiusAdd(ReadOnlySpan<double> x, ReadOnlySpan<double> y, Span<double> dst)
        {
            if (x.Length != y.Length || dst.Length < x.Length)
                throw new ArgumentException("Mobius addition buffers must have the same dimension.");

            double xNormSq = Dot(x, x);
            double yNormSq = Dot(y, y);
            double xy = Dot(x, y);
            double denominator = 1.0 + (2.0 * xy) + (xNormSq * yNormSq);
            if (Math.Abs(denominator) < 1e-12)
            {
                dst.Clear();
                return;
            }

            double xScale = 1.0 + (2.0 * xy) + yNormSq;
            double yScale = 1.0 - xNormSq;
            for (int i = 0; i < x.Length; i++)
                dst[i] = ((xScale * x[i]) + (yScale * y[i])) / denominator;
        }

        private static double ConformalFactor(ReadOnlySpan<double> p)
        {
            double limit = 1.0 - Poincare.BoundaryMargin;
            double limitSq = limit * limit;
            double normSq = Math.Min(Dot(p, p), limitSq);
            return 2.0 / Math.Max(1e-12, 1.0 - normSq);
        }

        private static double EuclideanNorm(ReadOnlySpan<double> vector)
            => Math.Sqrt(Dot(vector, vector));

        private static double Dot(ReadOnlySpan<double> left, ReadOnlySpan<double> right)
        {
            double sum = 0.0;
            for (int i = 0; i < left.Length; i++)
                sum += left[i] * right[i];
            return sum;
        }

        private static double Atanh(double x)
            => 0.5 * Math.Log((1.0 + x) / Math.Max(1e-12, 1.0 - x));

        private static double ClampUnitInterval(double x)
            => Math.Min(Math.Max(x, 0.0), 1.0 - 1e-12);
    }
}
