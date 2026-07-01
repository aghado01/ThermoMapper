using System;
using System.Collections.Generic;

namespace Clustering.Evaluation.Internal;

/// <summary>
/// Silhouette coefficient: <c>s_i = (b_i − a_i) / max(a_i, b_i)</c>
/// averaged over points, where <c>a_i</c> is the mean intra-cluster
/// distance and <c>b_i</c> is the mean distance to the nearest other
/// cluster. Range <c>[−1, 1]</c>; higher is better.
/// </summary>
/// <remarks>
/// <para><b>Edge cases.</b> Singletons score 0 (no intra-cluster pair
/// to average over). Partitions with fewer than 2 clusters return 0.0
/// (silhouette is undefined).</para>
///
/// <para><b>Distance metric.</b> Currently hardcoded to Euclidean.
/// A future revision should accept an injectable
/// <c>IDistanceMetric</c> so non-Euclidean data (Mahalanobis, product
/// manifolds, etc.) can be scored honestly.</para>
/// </remarks>
public sealed class SilhouetteEvaluator : IInternalClusterEvaluator
{
    public double Evaluate(double[][] data, int[] labels)
    {
        EvaluationHelpers.ValidateInputs(data, labels);
        // Drop unassigned points (label == Unassigned) before scoring — they
        // belong to no cluster and must not enter any intra/inter distance.
        (data, labels) = EvaluationHelpers.AssignedSubset(data, labels);

        int n = data.Length;
        int[] mappedLabels = EvaluationHelpers.MapLabelsToDense(labels, out _);
        int clusterCount = 0;
        foreach (int label in mappedLabels)
            clusterCount = Math.Max(clusterCount, label + 1);

        if (clusterCount < 2)
            return 0.0;

        var clusterMembers = new List<int>[clusterCount];
        for (int c = 0; c < clusterCount; c++)
            clusterMembers[c] = new List<int>();

        for (int i = 0; i < n; i++)
            clusterMembers[mappedLabels[i]].Add(i);

        double total = 0.0;
        for (int i = 0; i < n; i++)
        {
            int cluster = mappedLabels[i];
            var ownCluster = clusterMembers[cluster];

            double a = 0.0;
            if (ownCluster.Count > 1)
            {
                for (int j = 0; j < ownCluster.Count; j++)
                {
                    int other = ownCluster[j];
                    if (other == i) continue;
                    a += EuclideanDistance(data[i], data[other]);
                }
                a /= ownCluster.Count - 1;
            }

            double b = double.PositiveInfinity;
            for (int c = 0; c < clusterCount; c++)
            {
                if (c == cluster)
                    continue;

                var otherCluster = clusterMembers[c];
                if (otherCluster.Count == 0)
                    continue;

                double avgDistance = 0.0;
                for (int j = 0; j < otherCluster.Count; j++)
                    avgDistance += EuclideanDistance(data[i], data[otherCluster[j]]);

                avgDistance /= otherCluster.Count;
                if (avgDistance < b)
                    b = avgDistance;
            }

            double score;
            if (ownCluster.Count == 1)
            {
                score = 0.0;
            }
            else
            {
                double max = Math.Max(a, b);
                score = max <= 0.0 ? 0.0 : (b - a) / max;
            }

            total += score;
        }

        return total / n;
    }

    private static double EuclideanDistance(double[] x, double[] y)
    {
        double sum = 0.0;
        for (int d = 0; d < x.Length; d++)
        {
            double diff = x[d] - y[d];
            sum += diff * diff;
        }
        return Math.Sqrt(sum);
    }
}
