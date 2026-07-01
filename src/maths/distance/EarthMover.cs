using System;
using System.Buffers;

namespace Maths.Distance
{
    public static class EarthMover
    {
        /// <summary>Mismatched lengths: <c>Array.Copy(b, sb, a.Length)</c> throws if <c>b</c> is shorter; no explicit guard needed.</summary>
        public static double Distance1D(double[] a, double[] b)
        {
            // 1D EMD = integral of |CDF_a(x) - CDF_b(x)| dx
            // Arrays must be sorted ascending for the CDF formulation to be correct.
            //
            // Sort scratch comes from ArrayPool — graph construction calls this
            // O(N²) times. Two heap allocations per call (the previous Clone()
            // pattern) dominated GC pressure for medium-large N.
            int len = a.Length;
            if (len == 0) return 0.0;

            var pool = ArrayPool<double>.Shared;
            double[] sa = pool.Rent(len);
            double[] sb = pool.Rent(len);
            try
            {
                Array.Copy(a, sa, len);
                Array.Copy(b, sb, len);
                Array.Sort(sa, 0, len);
                Array.Sort(sb, 0, len);

                double cumA = 0, cumB = 0, sum = 0;
                for (int i = 0; i < len; i++)
                {
                    cumA += sa[i];
                    cumB += sb[i];
                    sum += Math.Abs(cumA - cumB);
                }
                return sum / len;
            }
            finally
            {
                pool.Return(sa);
                pool.Return(sb);
            }
        }
    }
}
