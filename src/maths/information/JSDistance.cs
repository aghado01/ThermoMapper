using System;
using System.Buffers;
using System.Numerics.Tensors;
using System.Runtime.CompilerServices;
using System.Threading;

namespace Maths.Information
{
    public static class JSDistance
    {
        // Three-tier scratch buffers — same pattern as Mahalanobis. JSD needs
        // two len-sized buffers: one for m = (p+q)/2 (later repurposed for the
        // q-side term) and one for the p-side term.
        private static readonly ThreadLocal<(double[] a, double[] b)> _jsdScratch =
            new(() => (new double[512], new double[512]));

        /// <summary>Spans are sized from <c>p.Length</c>; a shorter <c>q</c> throws at <c>q.AsSpan(0, len)</c> with an unhelpful bounds error rather than at the call boundary.</summary>
        public static double Distance(double[] p, double[] q)
        {
            int len = p.Length;
            if (len == 0) return 0.0;

            if (len <= 64)
            {
                Span<double> a = stackalloc double[len];
                Span<double> b = stackalloc double[len];
                return JensenShannonCore(p, q, a, b, len);
            }

            if (len <= 512)
            {
                var (aArr, bArr) = _jsdScratch.Value;
                return JensenShannonCore(p, q,
                    aArr.AsSpan(0, len),
                    bArr.AsSpan(0, len), len);
            }

            var pool = ArrayPool<double>.Shared;
            double[] rentedA = pool.Rent(len);
            double[] rentedB = pool.Rent(len);
            try
            {
                return JensenShannonCore(p, q,
                    rentedA.AsSpan(0, len),
                    rentedB.AsSpan(0, len), len);
            }
            finally
            {
                pool.Return(rentedA);
                pool.Return(rentedB);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double JensenShannonCore(
            double[] p, double[] q,
            Span<double> mBuf, Span<double> termP, int len)
        {
            var pSpan = p.AsSpan(0, len);
            var qSpan = q.AsSpan(0, len);

            // m = (p+q)/2, clamped above 1e-15 to make division safe.
            TensorPrimitives.Add<double>(pSpan, qSpan, mBuf);
            TensorPrimitives.Multiply<double>(mBuf, 0.5, mBuf);
            TensorPrimitives.Max<double>(mBuf, 1e-15, mBuf);

            // termP = p * log(max(p/m, 1e-15))
            // Clamping the ratio is mathematically safe: where p ≈ 0 the
            // outer *p multiply collapses the term to 0, regardless of the
            // (clamped) log value. This replaces the original `if p[i] > 1e-15`
            // guard, which would otherwise kill SIMD vectorization.
            TensorPrimitives.Divide<double>(pSpan, mBuf, termP);
            TensorPrimitives.Max<double>(termP, 1e-15, termP);
            TensorPrimitives.Log<double>(termP, termP);
            TensorPrimitives.Multiply<double>(pSpan, termP, termP);
            double sumP = TensorPrimitives.Sum<double>(termP);

            // mBuf is repurposed for the q-side term. Order matters — m must
            // remain valid through both Divide calls; q-side runs second so
            // we can overwrite m in place.
            TensorPrimitives.Divide<double>(qSpan, mBuf, mBuf);
            TensorPrimitives.Max<double>(mBuf, 1e-15, mBuf);
            TensorPrimitives.Log<double>(mBuf, mBuf);
            TensorPrimitives.Multiply<double>(qSpan, mBuf, mBuf);
            double sumQ = TensorPrimitives.Sum<double>(mBuf);

            // JSD = 0.5 * (KL(p||m) + KL(q||m)); return sqrt for proper metric.
            return Math.Sqrt(0.5 * (sumP + sumQ));
        }
    }
}
