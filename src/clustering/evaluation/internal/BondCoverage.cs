using System;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Evaluation.Internal;

/// <summary>
/// Fraction of total edge weight lying within clusters:
/// <c>Coverage = (Σ_{c_i=c_j} w_ij) / (Σ w_ij)</c>. Range <c>[0, 1]</c>;
/// <b>higher is better</b>. Raw internal density — the unadjusted counterpart to
/// <see cref="BondModularity"/> (rewards lumping into one cluster, which
/// modularity penalizes; reading both separates "trivial big cluster" from
/// "genuinely community-shaped").
/// </summary>
public sealed class BondCoverage : IGraphPartitionEvaluator
{
    public string Name => "BondCoverage";

    public double Evaluate(CsrGraph graph, double[] edgeWeight, int[] labels, int clusterCount)
    {
        ArgumentNullException.ThrowIfNull(edgeWeight);
        ArgumentNullException.ThrowIfNull(labels);

        int n = graph.NodeCount;
        if (n == 0 || labels.Length != n)
            return 0.0;
        if (edgeWeight.Length != graph.Targets.Length)
            throw new InvalidOperationException(
                $"edgeWeight length ({edgeWeight.Length}) does not match CSR slot count ({graph.Targets.Length}).");

        double totalWeight = 0.0;
        double withinWeight = 0.0;

        for (int i = 0; i < n; i++)
        {
            int rowEnd = graph.RowPointers[i + 1];
            for (int e = graph.RowPointers[i]; e < rowEnd; e++)
            {
                int j = graph.Targets[e];
                if (j <= i) continue;
                // Assigned-edges-only: an edge with an unassigned endpoint
                // counts toward neither the within-weight nor the total.
                if (labels[i] == Assignment.Unassigned || labels[j] == Assignment.Unassigned) continue;

                double w = edgeWeight[e];
                totalWeight += w;
                if (labels[i] == labels[j])
                    withinWeight += w;
            }
        }

        return totalWeight > 0.0 ? withinWeight / totalWeight : 0.0;
    }
}
