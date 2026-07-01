using System;

namespace Maths.Distance
{
    public static class Jaccard
    {
        /// <summary>Jaccard distance on thresholded real-valued feature vectors.</summary>
        public static double Distance(double[] a, double[] b, double threshold = 0.5)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));
            if (a.Length != b.Length)
                throw new ArgumentException("Jaccard distance requires vectors of equal length.", nameof(b));
            if (threshold < 0.0 || threshold > 1.0)
                throw new ArgumentOutOfRangeException(nameof(threshold), threshold, "Threshold must be in [0, 1].");

            int intersection = 0, union = 0;
            for (int i = 0; i < a.Length; i++)
            {
                bool ai = a[i] > threshold;
                bool bi = b[i] > threshold;
                if (ai || bi) union++;
                if (ai && bi) intersection++;
            }
            if (union == 0) return 0.0; // both empty → identical
            return 1.0 - (double)intersection / union;
        }
    }
}
