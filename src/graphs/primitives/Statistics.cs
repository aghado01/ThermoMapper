using System;

namespace Graphs.Primitives
{
    /// <summary>
    /// Shared statistical primitives used across bandwidth estimation and
    /// graph diagnostics. Hoists what used to be 4 independent
    /// <c>MedianOfSorted</c> implementations (in <c>BandwidthEstimation</c>,
    /// <c>EdgeWeights</c>, <c>MstBridge</c>, <c>NeighborhoodScale</c>) and 2
    /// inline skewness implementations (<c>Hubness</c>, <c>MstBridge</c>)
    /// into one source of truth so the constants and edge cases don't
    /// drift apart over time.
    /// </summary>
    public static class Statistics
    {
        /// <summary>
        /// Median of a pre-sorted span. Even-length spans return the
        /// arithmetic mean of the two middle values; odd-length spans
        /// return the single middle. Empty span → 0.
        /// </summary>
        public static double MedianOfSorted(ReadOnlySpan<double> sorted)
        {
            int n = sorted.Length;
            if (n == 0) return 0.0;
            int mid = n / 2;
            return (n & 1) == 0
                ? (sorted[mid - 1] + sorted[mid]) / 2.0
                : sorted[mid];
        }

        /// <summary>
        /// Population skewness (γ₁) of a value sample:
        /// <c>m₃ / m₂^(3/2)</c> where <c>m_k = E[(x − x̄)^k]</c>.
        /// Returns 0 when the sample has fewer than 2 values or when
        /// variance is below the numerical floor.
        /// </summary>
        public static double Skewness(ReadOnlySpan<double> values)
        {
            int n = values.Length;
            if (n < 2) return 0.0;

            double mean = 0.0;
            for (int i = 0; i < n; i++) mean += values[i];
            mean /= n;

            double m2 = 0.0, m3 = 0.0;
            for (int i = 0; i < n; i++)
            {
                double d = values[i] - mean;
                double d2 = d * d;
                m2 += d2;
                m3 += d2 * d;
            }
            m2 /= n;
            m3 /= n;

            const double VarianceFloor = 1e-12;
            if (m2 < VarianceFloor) return 0.0;
            return m3 / Math.Pow(m2, 1.5);
        }

        /// <summary>
        /// Skewness over an integer sample (in-degree histograms, count
        /// vectors). Promotes to <see cref="double"/> internally; same
        /// definition + variance floor as the double overload.
        /// </summary>
        public static double Skewness(ReadOnlySpan<int> values)
        {
            int n = values.Length;
            if (n < 2) return 0.0;

            double mean = 0.0;
            for (int i = 0; i < n; i++) mean += values[i];
            mean /= n;

            double m2 = 0.0, m3 = 0.0;
            for (int i = 0; i < n; i++)
            {
                double d = values[i] - mean;
                double d2 = d * d;
                m2 += d2;
                m3 += d2 * d;
            }
            m2 /= n;
            m3 /= n;

            const double VarianceFloor = 1e-12;
            if (m2 < VarianceFloor) return 0.0;
            return m3 / Math.Pow(m2, 1.5);
        }
    }
}
