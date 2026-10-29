using System;
using Graphs.Primitives;

namespace Graphs.Observables;

/// <summary>
/// Link overlap of two node labelings on a graph: the fraction of edges on which they agree — the
/// edge's endpoints together in both, or apart in both — reducing the per-edge agreement 1-form to
/// a degree-0 value.
/// <code>Q = Σ_e w_e·[same_A(e) = same_B(e)] / Σ_e w_e</code>
/// The spin-glass link overlap of two Potts states (Hed and Domany 2007), with cluster labels as the
/// states: it compares where two partitions draw their boundaries (domain-wall area), where
/// all-pairs indices compare volume. On the complete graph it equals the Rand index.
/// </summary>
/// <remarks>
/// <para><b>Uses.</b> Any two partitions of the same nodes: two methods, two settings, two replicas
/// of a sampler at one temperature (the replica overlap), or a result against a reference.</para>
///
/// <para><b>Weighting.</b> Uniform by default (each edge counts once, the spin-glass form). With
/// <c>weighted</c> each edge counts by <see cref="CsrGraph.Weights"/>, read as couplings (larger =
/// stronger), so pass a coupling-weighted graph.</para>
///
/// <para><b>Unplaced nodes.</b> A negative label marks a node a labeling did not place (clustering's
/// <c>Assignment.Unassigned</c> is −1). An edge touching such a node in either labeling is not
/// counted, which keeps the overlap symmetric. Returns 0 when no edge counts.</para>
/// </remarks>
public static class LinkOverlap
{
    /// <summary>Link overlap of labelings <paramref name="a"/> and <paramref name="b"/>, each of
    /// length <see cref="CsrGraph.NodeCount"/>, in <c>[0, 1]</c>.</summary>
    public static double Compute(CsrGraph graph, int[] a, int[] b, bool weighted = false)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        if (a.Length != graph.NodeCount || b.Length != graph.NodeCount)
            throw new ArgumentException(
                $"Label lengths ({a.Length}, {b.Length}) must match the node count ({graph.NodeCount}).");

        double agreeing = 0.0, total = 0.0;
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
        {
            int i = edge.Source, j = edge.Target;
            if (a[i] < 0 || a[j] < 0 || b[i] < 0 || b[j] < 0) continue;

            double w = weighted ? graph.Weights[edge.Slot] : 1.0;
            total += w;
            if ((a[i] == a[j]) == (b[i] == b[j])) agreeing += w;
        }

        return total > 0.0 ? agreeing / total : 0.0;
    }
}
