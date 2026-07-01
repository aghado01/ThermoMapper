using System;
using System.Collections.Generic;
using Graphs.Primitives;

namespace Graphs.Diagnostics
{
    public readonly record struct NeighborhoodScaleReport(
        int K,
        double Median1NN,
        double MedianKthNN,
        double ScaleRatio);

    public static class NeighborhoodScale
    {
        /// <summary>
        /// Median 1NN distance vs. median k-th NN distance, both read
        /// from the directed (pre-symmetrization) selection. The
        /// previously-required <c>mutual</c> selection argument was
        /// vestigial — the directed selection's
        /// <c>NearestNeighborDistances</c> already carries the same 1NN
        /// information without depending on whether MutualKnn was used
        /// upstream.
        /// </summary>
        public static NeighborhoodScaleReport Compute(NeighborSelection directed, int k)
        {
            if (k <= 0)
                throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

            double median1NN = ComputeMedian(directed.NearestNeighborDistances);

            var kthDistances = new List<double>(directed.AllNeighbors.Length);
            for (int i = 0; i < directed.AllNeighbors.Length; i++)
            {
                Neighbor[] row = directed.AllNeighbors[i];
                if (row.Length == 0)
                    continue;

                int index = Math.Min(k, row.Length) - 1;
                kthDistances.Add(row[index].Distance);
            }

            double medianKth  = ComputeMedian(kthDistances);
            double scaleRatio = median1NN > 0.0 ? medianKth / median1NN : 0.0;

            return new NeighborhoodScaleReport(
                K: k,
                Median1NN: median1NN,
                MedianKthNN: medianKth,
                ScaleRatio: scaleRatio);
        }

        private static double ComputeMedian(IReadOnlyList<double> values)
        {
            if (values.Count == 0)
                return 0.0;

            var sorted = new double[values.Count];
            for (int i = 0; i < values.Count; i++)
                sorted[i] = values[i];

            Array.Sort(sorted);
            return Statistics.MedianOfSorted(sorted);
        }
    }
}
