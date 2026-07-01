using System;
using Graphs.Primitives;

namespace Graphs.Observables;

/// <summary>
/// Per-node marginals of the per-edge affinity currency — degree-0 node
/// fields summing the edge field over each node's incident edges.
/// Producer-agnostic by the currency seam: the same definitions read SW's
/// sampled <see cref="Affinities"/> and PKWang's closed-form ones.
/// </summary>
public static class AffinityNodeMarginals
{
    /// <summary>
    /// Bond mass per node: Σ_{j∈N(i)} G_ij — the attachment field (ascend:
    /// high where the neighborhood's bonds survive).
    /// </summary>
    public static double[] BondMass(CsrGraph graph, Affinities affinities)
    {
        Validate(graph, affinities);
        var marginal = new double[graph.NodeCount];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
        {
            double g = affinities.G[edge.Slot];
            marginal[edge.Source] += g;
            marginal[edge.Target] += g;
        }
        return marginal;
    }

    /// <summary>
    /// Local energy per node: Σ_{j∈N(i)} J_ij·(1 − G_ij) — the frustration
    /// field (descend: low inside coherent cores, high on contested
    /// boundaries). The equilibrium recapture of PKWang's <c>Hcum</c> idea:
    /// for the solver it is exact; for SW it carries the sampled correlations.
    /// </summary>
    public static double[] LocalEnergy(CsrGraph graph, Affinities affinities)
    {
        Validate(graph, affinities);
        var marginal = new double[graph.NodeCount];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
        {
            double j = graph.Weights[edge.Slot];
            double frustration = j * (1.0 - affinities.G[edge.Slot]);
            marginal[edge.Source] += frustration;
            marginal[edge.Target] += frustration;
        }
        return marginal;
    }

    private static void Validate(CsrGraph graph, Affinities affinities)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(affinities);
        if (affinities.G.Length != graph.Targets.Length)
            throw new InvalidOperationException(
                $"Affinities.G length ({affinities.G.Length}) does not match CSR slot count ({graph.Targets.Length}).");
    }
}
