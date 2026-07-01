using System;

namespace Graphs.Primitives;

/// <summary>
/// Reconciles a <em>directed</em> per-edge field over a symmetric CSR graph — one
/// value per directed slot — into a single undirected value written to the
/// <c>j &gt; i</c> slot, per a <see cref="SymmetrizationRule"/>. The mirror map
/// (<see cref="CsrGraph.BuildReverseSlotMap"/>) pairs each CSR slot with its
/// reverse direction. The value-level twin of the set-level kNN symmetrization
/// in <c>Graphs.Neighbors</c>; lifted out of PKWang so any directed-field method
/// shares it.
/// </summary>
public static class EdgeFieldSymmetrization
{
    /// <summary>Collapse two directed values into one per the rule.</summary>
    public static double Combine(double a, double b, SymmetrizationRule rule) => rule switch
    {
        SymmetrizationRule.Mutual    => Math.Min(a, b),
        SymmetrizationRule.Inclusive => Math.Max(a, b),
        SymmetrizationRule.Mean      => 0.5 * (a + b),
        _ => throw new ArgumentOutOfRangeException(nameof(rule), rule, "Unsupported symmetrization rule."),
    };

    /// <summary>
    /// In-place: for each undirected edge (the <c>j &gt; i</c> CSR half), write the
    /// reconciled value into its slot. The mirror half is left as the raw directed
    /// value — undefined by contract. <paramref name="mirror"/>[e] is the reverse
    /// slot of edge slot <c>e</c> (see <see cref="CsrGraph.BuildReverseSlotMap"/>).
    /// </summary>
    public static void Symmetrize(CsrGraph graph, double[] g, int[] mirror, SymmetrizationRule rule)
    {
        for (int i = 0; i < graph.NodeCount; i++)
        {
            int rowEnd = graph.RowPointers[i + 1];
            for (int e = graph.RowPointers[i]; e < rowEnd; e++)
            {
                if (graph.Targets[e] <= i) continue;
                g[e] = Combine(g[e], g[mirror[e]], rule);
            }
        }
    }
}
