using System;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Core.Solver;

/// <summary>
/// Wang 2020 mean field: one global ascending sort of the undirected couplings,
/// cumulative-summed into a single energy ladder, scattered back to each edge's
/// CSR slot. The cut <c>Hcum &gt; T·ln2</c> reduces to thermal single-linkage
/// (Lemma B; see the spc-samplers note).
/// </summary>
/// <remarks>
/// Parity variant: ranks each undirected edge once (the <c>j &gt; i</c> half of
/// the symmetric CSR), not the paper's directed <c>N·K</c> sort which
/// double-counts every edge and rescales <c>T</c> by ~2. The clustering family
/// is identical either way — only the temperature axis is reparametrized.
/// </remarks>
internal readonly struct MeanField : IField
{
    public static bool DirectedSymmetrize => false;

    public static double[] BuildHcum(CsrGraph graph)
    {
        int[] rowPtr = graph.RowPointers;
        int[] targets = graph.Targets;
        double[] weights = graph.Weights;
        int n = graph.NodeCount;

        // Pass 1: count unique undirected edges (the j > i half).
        int edgeCount = 0;
        for (int i = 0; i < n; i++)
        {
            int rowEnd = rowPtr[i + 1];
            for (int e = rowPtr[i]; e < rowEnd; e++)
                if (targets[e] > i) edgeCount++;
        }

        var hcum = new double[targets.Length];
        if (edgeCount == 0) return hcum;

        // Pass 2: gather (CSR slot, coupling) for each unique edge.
        var slot = new int[edgeCount];
        var coupling = new double[edgeCount];
        int k = 0;
        for (int i = 0; i < n; i++)
        {
            int rowEnd = rowPtr[i + 1];
            for (int e = rowPtr[i]; e < rowEnd; e++)
            {
                if (targets[e] <= i) continue;
                slot[k] = e;
                coupling[k] = weights[e];
                k++;
            }
        }

        // Ascending sort by coupling, permuting the slot map in tandem, then a
        // single cumulative sum. Hcum is rank-inclusive: the rank-r edge's value
        // is the sum of couplings 0..r, so the strongest edge carries the total.
        Array.Sort(coupling, slot);
        double cum = 0.0;
        for (int r = 0; r < edgeCount; r++)
        {
            cum += coupling[r];
            hcum[slot[r]] = cum;
        }
        return hcum;
    }
}
