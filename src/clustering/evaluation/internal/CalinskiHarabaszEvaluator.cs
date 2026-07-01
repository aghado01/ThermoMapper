using System;

namespace Clustering.Evaluation.Internal;

/// <summary>
/// Calinski-Harabasz index: ratio of between-cluster scatter to
/// within-cluster scatter, F-statistic-flavored:
/// <c>(BSS / (K−1)) / (WSS / (N−K))</c>. Higher is better.
/// </summary>
/// <remarks>
/// <para><b>Edge cases.</b> Returns <see cref="double.NaN"/> when the
/// index is undefined (fewer than 2 clusters, or N ≤ K). Returns
/// <see cref="double.PositiveInfinity"/> when within-cluster scatter
/// is zero (a degenerate but mathematically clean limit).</para>
///
/// <para><b>Distance metric.</b> Hardcoded to Euclidean — same caveat
/// as <see cref="SilhouetteEvaluator"/>.</para>
/// </remarks>
public sealed class CalinskiHarabaszEvaluator : IInternalClusterEvaluator
{
    public double Evaluate(double[][] data, int[] labels)
    {
        EvaluationHelpers.ValidateInputs(data, labels);
        // Drop unassigned points (label == Unassigned) before scoring — they
        // must not enter centroids, scatter, or the overall mean.
        (data, labels) = EvaluationHelpers.AssignedSubset(data, labels);

        int n = data.Length;
        int[] mappedLabels = EvaluationHelpers.MapLabelsToDense(labels, out _);
        int clusterCount = 0;
        foreach (int label in mappedLabels)
            clusterCount = Math.Max(clusterCount, label + 1);

        if (clusterCount < 2)
            return double.NaN;
        if (n <= clusterCount)
            return double.NaN;

        int[] counts = EvaluationHelpers.CountClusters(mappedLabels, clusterCount);
        double[][] centroids = EvaluationHelpers.ComputeCentroids(data, mappedLabels, clusterCount);
        double[] overallMean = ComputeOverallMean(data);

        double betweenSum = 0.0;
        double withinSum = 0.0;

        for (int c = 0; c < clusterCount; c++)
        {
            if (counts[c] == 0)
                continue;

            double[] clusterMean = centroids[c];
            double distance = SquaredEuclideanDistance(clusterMean, overallMean);
            betweenSum += counts[c] * distance;
        }

        for (int i = 0; i < n; i++)
        {
            int cluster = mappedLabels[i];
            double distance = SquaredEuclideanDistance(data[i], centroids[cluster]);
            withinSum += distance;
        }

        if (withinSum <= 0.0)
            return double.PositiveInfinity;

        double numerator = betweenSum / (clusterCount - 1);
        double denominator = withinSum / (n - clusterCount);
        return denominator <= 0.0 ? double.PositiveInfinity : numerator / denominator;
    }

    private static double[] ComputeOverallMean(double[][] data)
    {
        int n = data.Length;
        int dim = data[0].Length;
        var mean = new double[dim];

        for (int i = 0; i < n; i++)
        {
            for (int d = 0; d < dim; d++)
                mean[d] += data[i][d];
        }

        double inv = 1.0 / n;
        for (int d = 0; d < dim; d++)
            mean[d] *= inv;

        return mean;
    }

    private static double SquaredEuclideanDistance(double[] x, double[] y)
    {
        double sum = 0.0;
        for (int d = 0; d < x.Length; d++)
        {
            double diff = x[d] - y[d];
            sum += diff * diff;
        }
        return sum;
    }
}
