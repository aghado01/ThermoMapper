using System;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Evaluation.Internal;

/// <summary>
/// Newman–Girvan modularity on the edge-weighted graph:
/// <c>Q = Σ_c [ within_c / m − (D_c / 2m)² ]</c>, where <c>within_c</c> is the
/// weight inside cluster <c>c</c>, <c>D_c</c> its weighted-degree sum, and
/// <c>m</c> the total edge weight. Range <c>[-0.5, 1]</c>; <b>higher is better</b>
/// (the partition captures more within-cluster weight than a degree-preserving
/// null model expects).
/// </summary>
public sealed class BondModularity : IGraphPartitionEvaluator
{
    public string Name => "BondModularity";

    public double Evaluate(CsrGraph graph, double[] edgeWeight, int[] labels, int clusterCount)
    {
        ArgumentNullException.ThrowIfNull(edgeWeight);
        ArgumentNullException.ThrowIfNull(labels);

        int n = graph.NodeCount;
        if (n == 0 || clusterCount <= 0 || labels.Length != n)
            return 0.0;
        if (edgeWeight.Length != graph.Targets.Length)
            throw new InvalidOperationException(
                $"edgeWeight length ({edgeWeight.Length}) does not match CSR slot count ({graph.Targets.Length}).");

        var degree         = new double[n];
        var withinWeight   = new double[clusterCount];
        double totalWeight = 0.0;

        for (int i = 0; i < n; i++)
        {
            int rowEnd = graph.RowPointers[i + 1];
            for (int e = graph.RowPointers[i]; e < rowEnd; e++)
            {
                int j = graph.Targets[e];
                if (j <= i) continue;
                // Assigned-edges-only: an edge with an unassigned endpoint
                // contributes to neither degrees, totals, nor within-sums.
                if (labels[i] == Assignment.Unassigned || labels[j] == Assignment.Unassigned) continue;

                double w = edgeWeight[e];
                degree[i] += w;
                degree[j] += w;
                totalWeight += w;
                if (labels[i] == labels[j])
                    withinWeight[labels[i]] += w;
            }
        }

        if (totalWeight <= 0.0)
            return 0.0;

        var degreeSum = new double[clusterCount];
        for (int i = 0; i < n; i++)
        {
            if (labels[i] == Assignment.Unassigned) continue;
            degreeSum[labels[i]] += degree[i];
        }

        double twoM    = 2.0 * totalWeight;
        double fourMSq = twoM * twoM;

        double q = 0.0;
        for (int c = 0; c < clusterCount; c++)
            q += withinWeight[c] / totalWeight - (degreeSum[c] * degreeSum[c]) / fourMSq;
        return q;
    }
}
