using System;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Evaluation.External;

/// <summary>
/// Link overlap: the fraction of a graph's edges on which two partitions agree —
/// the edge's endpoints together in both, or apart in both.
/// <code>Q = Σ_e w_e·[same_A(e) = same_B(e)] / Σ_e w_e</code>
/// The spin-glass link overlap of two Potts states (Hed and Domany 2007), with
/// cluster labels as the states: it measures agreement on boundaries
/// (domain-wall area) where all-pairs indices measure agreement on volume. On the
/// complete graph it equals <see cref="RandIndex"/>.
/// </summary>
/// <remarks>
/// <para><b>Range.</b> <c>[0, 1]</c>; <c>1.0</c> = the partitions cut exactly the
/// same edges. <b>Higher is better.</b> Not chance-corrected.</para>
///
/// <para><b>Weighting.</b> Uniform by default (every edge counts once, the
/// spin-glass form). With <c>weighted</c> each edge counts by
/// <see cref="CsrGraph.Weights"/>, read as couplings (larger = stronger), so pass
/// a coupling-weighted graph.</para>
///
/// <para><b>Symmetry.</b> Symmetric in its two label arguments.</para>
///
/// <para><b>Edge cases.</b> Edges with a predicted-unassigned endpoint are not
/// scored. Returns 0.0 when no edge is scored.</para>
/// </remarks>
public sealed class LinkOverlap : IGraphExternalClusterEvaluator
{
    private readonly bool _weighted;

    public LinkOverlap(bool weighted = false) => _weighted = weighted;

    public string Name => _weighted ? "LinkOverlapWeighted" : "LinkOverlap";

    public double Evaluate(CsrGraph graph, int[] predictedLabels, int[] referenceLabels)
    {
        ArgumentNullException.ThrowIfNull(predictedLabels);
        ArgumentNullException.ThrowIfNull(referenceLabels);
        if (predictedLabels.Length != graph.NodeCount || referenceLabels.Length != graph.NodeCount)
            throw new ArgumentException(
                $"Label lengths ({predictedLabels.Length}, {referenceLabels.Length}) must match the node count ({graph.NodeCount}).");

        double agreeing = 0.0, total = 0.0;
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
        {
            int i = edge.Source, j = edge.Target;
            if (predictedLabels[i] == Assignment.Unassigned || predictedLabels[j] == Assignment.Unassigned)
                continue;

            double w = _weighted ? graph.Weights[edge.Slot] : 1.0;
            bool sameA = predictedLabels[i] == predictedLabels[j];
            bool sameB = referenceLabels[i] == referenceLabels[j];
            total += w;
            if (sameA == sameB) agreeing += w;
        }

        return total > 0.0 ? agreeing / total : 0.0;
    }
}
