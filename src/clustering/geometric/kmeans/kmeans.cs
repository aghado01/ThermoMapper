using System;
using System.Collections.Generic;
using System.Linq;

namespace Clustering.Geometric.KMeans
{
    public static class KMeansPlusPlus
    {
        /// <summary>
        /// Full KMeans++ with smart initialization + Lloyd iterations.
        /// Returns labels and final centroids.
        /// </summary>
        public static KMeansResult Cluster(
        double[][] data,
        int k,
        int maxIterations = 100,
        int seed = 42,
        double tol = 1e-6)
        {
            if (data == null || data.Length == 0)
                throw new ArgumentException("Data cannot be empty.");
            if (k < 1 || k > data.Length)
                throw new ArgumentException("k must be between 1 and n.");

            var rng = new Random(seed);
            int n = data.Length;
            int dim = data[0].Length;

            // Step 1: KMeans++ initialization
            var centroids = InitializeCentroids(data, k, rng);

            int[] labels = new int[n];
            double prevDistortion = double.MaxValue;

            for (int iter = 0; iter < maxIterations; iter++)
            {
                // Assignment step
                bool changed = AssignLabels(data, centroids, labels);

                // Update step
                UpdateCentroids(data, labels, centroids, k, dim);

                // Convergence check
                double distortion = ComputeDistortion(data, centroids, labels);
                if (Math.Abs(prevDistortion - distortion) < tol || !changed)
                    break;

                prevDistortion = distortion;
            }

            return new KMeansResult(labels, centroids, k);
        }

        private static double[][] InitializeCentroids(double[][] data, int k, Random rng)
        {
            int n = data.Length;
            var centroids = new double[k][];
            var distances = new double[n];

            // First centroid: random point
            int idx = rng.Next(n);
            centroids[0] = (double[])data[idx].Clone();

            for (int c = 1; c < k; c++)
            {
                // Compute squared distance to nearest centroid
                double sumDist = 0;
                for (int i = 0; i < n; i++)
                {
                    distances[i] = SquaredDistanceToNearest(data[i], centroids, c);
                    sumDist += distances[i];
                }

                // Roulette wheel selection
                double r = rng.NextDouble() * sumDist;
                double cum = 0;
                for (int i = 0; i < n; i++)
                {
                    cum += distances[i];
                    if (cum >= r)
                    {
                        centroids[c] = (double[])data[i].Clone();
                        break;
                    }
                }
            }

            return centroids;
        }

        private static double SquaredDistanceToNearest(double[] point, double[][] centroids, int currentCount)
        {
            double minSq = double.MaxValue;
            for (int c = 0; c < currentCount; c++)
            {
                double d = SquaredEuclidean(point, centroids[c]);
                if (d < minSq) minSq = d;
            }
            return minSq;
        }

        private static bool AssignLabels(double[][] data, double[][] centroids, int[] labels)
        {
            bool changed = false;
            for (int i = 0; i < data.Length; i++)
            {
                int best = 0;
                double bestDist = SquaredEuclidean(data[i], centroids[0]);

                for (int c = 1; c < centroids.Length; c++)
                {
                    double d = SquaredEuclidean(data[i], centroids[c]);
                    if (d < bestDist)
                    {
                        bestDist = d;
                        best = c;
                    }
                }

                if (labels[i] != best)
                {
                    labels[i] = best;
                    changed = true;
                }
            }
            return changed;
        }

        private static void UpdateCentroids(double[][] data, int[] labels, double[][] centroids, int k, int dim)
        {
            var sums = new double[k][];
            var counts = new int[k];

            for (int c = 0; c < k; c++)
                sums[c] = new double[dim];

            for (int i = 0; i < data.Length; i++)
            {
                int c = labels[i];
                counts[c]++;
                for (int d = 0; d < dim; d++)
                    sums[c][d] += data[i][d];
            }

            for (int c = 0; c < k; c++)
            {
                if (counts[c] == 0) continue; // empty cluster — keep old centroid

                for (int d = 0; d < dim; d++)
                    centroids[c][d] = sums[c][d] / counts[c];
            }
        }

        private static double ComputeDistortion(double[][] data, double[][] centroids, int[] labels)
        {
            double sum = 0;
            for (int i = 0; i < data.Length; i++)
                sum += SquaredEuclidean(data[i], centroids[labels[i]]);
            return sum;
        }

        [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
        private static double SquaredEuclidean(double[] a, double[] b)
        {
            double sum = 0;
            for (int i = 0; i < a.Length; i++)
            {
                double diff = a[i] - b[i];
                sum += diff * diff;
            }
            return sum;
        }
    }

    public record KMeansResult(
        int[] Labels,
        double[][] Centroids,
        int K
    );

}
