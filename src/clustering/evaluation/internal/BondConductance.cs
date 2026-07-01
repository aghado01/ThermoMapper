using System;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Evaluation.Internal;

/// <summary>
/// Mean per-cluster conductance on the edge-weighted graph:
/// <c>φ(c) = cut(c) / min(vol(c), 2m − vol(c))</c>, averaged over clusters with
/// a defined denominator. Range <c>[0, 1]</c>; <b>lower is better</b> (clusters
/// keep their weight internal, leaking little across the boundary). The
/// per-cluster, boundary-tightness complement to <see cref="BondModularity"/>
/// (null-model density) and <see cref="BondCoverage"/> (raw internal fraction).
/// </summary>
public sealed class BondConductance : IGraphPartitionEvaluator
{
    public string Name => "BondConductance";

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

        var degree = new double[n];
        var cutPerCluster = new double[clusterCount];
        double totalWeight = 0.0;

        for (int i = 0; i < n; i++)
        {
            int rowEnd = graph.RowPointers[i + 1];
            for (int e = graph.RowPointers[i]; e < rowEnd; e++)
            {
                int j = graph.Targets[e];
                if (j <= i) continue;
                // Assigned-edges-only: an edge with an unassigned endpoint
                // contributes to neither degrees, volumes, nor cuts.
                if (labels[i] == Assignment.Unassigned || labels[j] == Assignment.Unassigned) continue;

                double w = edgeWeight[e];
                degree[i] += w;
                degree[j] += w;
                totalWeight += w;

                if (labels[i] != labels[j])
                {
                    // Cross-boundary edge: contributes to BOTH clusters' cuts.
                    cutPerCluster[labels[i]] += w;
                    cutPerCluster[labels[j]] += w;
                }
            }
        }

        if (totalWeight <= 0.0)
            return 0.0;

        var volumePerCluster = new double[clusterCount];
        for (int i = 0; i < n; i++)
        {
            if (labels[i] == Assignment.Unassigned) continue;
            volumePerCluster[labels[i]] += degree[i];
        }

        double twoM = 2.0 * totalWeight;
        double conductanceSum = 0.0;
        int contributingClusters = 0;

        for (int c = 0; c < clusterCount; c++)
        {
            double vol = volumePerCluster[c];
            double denom = Math.Min(vol, twoM - vol);
            if (denom <= 0.0) continue; // undefined for this cluster — skip
            conductanceSum += cutPerCluster[c] / denom;
            contributingClusters++;
        }

        return contributingClusters > 0 ? conductanceSum / contributingClusters : 0.0;
    }
}
