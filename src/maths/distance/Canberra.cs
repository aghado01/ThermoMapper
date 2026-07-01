using System;

namespace Maths.Distance
{
    public static class Canberra
    {
        /// <summary>Mismatched lengths throw at <c>b[i]</c> via array indexing; no explicit guard needed.</summary>
        public static double Distance(double[] a, double[] b)
        {
            double sum = 0;
            for (int i = 0; i < a.Length; i++)
            {
                double denom = Math.Abs(a[i]) + Math.Abs(b[i]);
                if (denom < 1e-15) continue;
                sum += Math.Abs(a[i] - b[i]) / denom;
            }
            return sum;
        }
    }
}
