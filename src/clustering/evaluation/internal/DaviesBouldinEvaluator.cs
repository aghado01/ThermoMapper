using System;

namespace Clustering.Evaluation.Internal;

/// <summary>
/// Davies-Bouldin index: average over clusters of the worst-case
/// <c>(scatter_i + scatter_j) / separation_{ij}</c> ratio. Lower is
/// better; 0 indicates perfectly compact, well-separated clusters.
/// </summary>
/// <remarks>
/// <para><b>Edge cases.</b> Returns <see cref="double.NaN"/> when there
/// are fewer than 2 clusters (DB index is undefined). Returns
/// <see cref="double.PositiveInfinity"/> when any two cluster centroids
/// coincide (separation 0).</para>
///
/// <para><b>Distance metric.</b> Hardcoded to Euclidean — same caveat
/// as <see cref="SilhouetteEvaluator"/>.</para>
/// </remarks>
public sealed class DaviesBouldinEvaluator : IInternalClusterEvaluator
{
    public double Evaluate(double[][] data, int[] labels)
    {
        EvaluationHelpers.ValidateInputs(data, labels);
        // Drop unassigned points (label == Unassigned) before scoring — they
        // must not enter centroids or per-cluster scatter.
        (data, labels) = EvaluationHelpers.AssignedSubset(data, labels);

        int n = data.Length;
        int[] mappedLabels = EvaluationHelpers.MapLabelsToDense(labels, out _);
        int clusterCount = 0;
        foreach (int label in mappedLabels)
            clusterCount = Math.Max(clusterCount, label + 1);

        if (clusterCount < 2)
            return double.NaN;

        int[] counts = EvaluationHelpers.CountClusters(mappedLabels, clusterCount);
        double[][] centroids = EvaluationHelpers.ComputeCentroids(data, mappedLabels, clusterCount);
        double[] scatters = new double[clusterCount];

        for (int i = 0; i < n; i++)
        {
            int cluster = mappedLabels[i];
            double distance = SquaredEuclideanDistance(data[i], centroids[cluster]);
            scatters[cluster] += Math.Sqrt(distance);
        }

        for (int c = 0; c < clusterCount; c++)
        {
            if (counts[c] > 0)
                scatters[c] /= counts[c];
        }

        double sum = 0.0;
        for (int c = 0; c < clusterCount; c++)
        {
            double worstRatio = 0.0;
            for (int d = 0; d < clusterCount; d++)
            {
                if (c == d)
                    continue;

                double separation = Math.Sqrt(SquaredEuclideanDistance(centroids[c], centroids[d]));
                if (separation <= 0.0)
                {
                    worstRatio = double.PositiveInfinity;
                    break;
                }

                double ratio = (scatters[c] + scatters[d]) / separation;
                if (ratio > worstRatio)
                    worstRatio = ratio;
            }
            sum += worstRatio;
        }

        return sum / clusterCount;
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
