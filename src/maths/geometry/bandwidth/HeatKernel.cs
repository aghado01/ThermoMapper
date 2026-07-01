using System;

namespace Maths.Geometry.Bandwidth
{
    public static class HeatKernel
    {
        public static double BandwidthFromBeta(double beta)
        {
            if (!double.IsFinite(beta) || beta <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(beta), beta, "Beta must be finite and positive.");

            return 1.0 / Math.Sqrt(2.0 * beta);
        }
    }
}
