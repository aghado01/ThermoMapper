using System;
using System.Buffers;
using System.Runtime.CompilerServices;
using System.Threading;

namespace Maths.Distance
{
    public static class Mahalanobis
    {
        // Three-tier scratch buffer allocation — avoids O(N²) heap allocations
        // during graph construction while remaining thread-safe under Parallel.For.
        //
        // Tier 1 (dim ≤ 64):  stackalloc — zero heap cost, JIT-inlinable.
        // Tier 2 (dim ≤ 512): ThreadLocal fixed arrays — one allocation per thread, ever.
        // Tier 3 (dim > 512): ArrayPool rent/return per call.
        //
        // ThreadLocal note: _mahalScratch allocates once per thread that touches it.
        // Primary caller is parallel graph construction (thread count ≤ ProcessorCount).
        // Calling from many short-lived threads (e.g., test isolation) is correct but
        // wastes memory. Write-before-read invariant: MahalanobisCore always writes
        // diff[] before reading it and temp[i] before the final dot — scratch is safe.
        private static readonly ThreadLocal<(double[] diff, double[] temp)> _mahalScratch =
            new(() => (new double[512], new double[512]));

        /// <summary>Explicit length guard: inner loop writes <c>diff[i] = a[i] - b[i]</c> over <c>a.Length</c> — a shorter <c>b</c> would throw inside the loop with no context on the mismatch.</summary>
        public static double Distance(double[] a, double[] b, double[,] invCov)
        {
            if (a == null || b == null || invCov == null)
                throw new ArgumentNullException(a == null ? nameof(a) : b == null ? nameof(b) : nameof(invCov));
            if (a.Length != b.Length)
                throw new ArgumentException("Feature vectors must have the same dimensionality.");
            return Distance(a, b, invCov, a.Length);
        }

        public static double Distance(double[] a, double[] b, double[,] invCov, int dim)
            => Math.Sqrt(DispatchQuadraticForm(a, b, invCov, dim));

        /// <summary>
        /// Returns the squared Mahalanobis distance: (a−b)ᵀ Σ⁻¹ (a−b).
        /// Avoids the sqrt of <see cref="Distance(double[],double[],double[,])"/>;
        /// prefer this when the squared form is all that is needed (e.g. log-pdf evaluation).
        /// </summary>
        public static double DistanceSquared(double[] a, double[] b, double[,] invCov)
        {
            if (a == null || b == null || invCov == null)
                throw new ArgumentNullException(a == null ? nameof(a) : b == null ? nameof(b) : nameof(invCov));
            if (a.Length != b.Length)
                throw new ArgumentException("Feature vectors must have the same dimensionality.");
            return DispatchQuadraticForm(a, b, invCov, a.Length);
        }

        /// <summary>
        /// Dimension-dispatched entry point. Selects the scratch-allocation tier and
        /// delegates to <see cref="QuadraticFormCore"/>.
        /// </summary>
        private static double DispatchQuadraticForm(double[] a, double[] b, double[,] invCov, int dim)
        {
            if (dim <= 64)
            {
                Span<double> diff = stackalloc double[dim];
                Span<double> temp = stackalloc double[dim];
                return QuadraticFormCore(a, b, invCov, diff, temp, dim);
            }

            if (dim <= 512)
            {
                var (diffArr, tempArr) = _mahalScratch.Value;
                return QuadraticFormCore(a, b, invCov,
                    diffArr.AsSpan(0, dim),
                    tempArr.AsSpan(0, dim), dim);
            }

            // dim > 512 — pool rent, guaranteed return
            double[] rentedDiff = ArrayPool<double>.Shared.Rent(dim);
            double[] rentedTemp = ArrayPool<double>.Shared.Rent(dim);
            try
            {
                return QuadraticFormCore(a, b, invCov,
                    rentedDiff.AsSpan(0, dim),
                    rentedTemp.AsSpan(0, dim), dim);
            }
            finally
            {
                ArrayPool<double>.Shared.Return(rentedDiff);
                ArrayPool<double>.Shared.Return(rentedTemp);
            }
        }

        /// <summary>
        /// Inner kernel. Returns the raw quadratic form max(0, (a−b)ᵀ Σ⁻¹ (a−b)).
        /// Callers apply sqrt for Euclidean-equivalent distance, or use directly for D².
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double QuadraticFormCore(
            double[] a, double[] b, double[,] invCov,
            Span<double> diff, Span<double> temp, int dim)
        {
            for (int i = 0; i < dim; i++)
                diff[i] = a[i] - b[i];

            for (int i = 0; i < dim; i++)
            {
                double s = 0;
                for (int j = 0; j < dim; j++)
                    s += invCov[i, j] * diff[j];
                temp[i] = s;
            }

            double quad = 0;
            for (int i = 0; i < dim; i++)
                quad += diff[i] * temp[i];

            return Math.Max(0, quad);
        }
    }
}
