using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Core.Solver;

/// <summary>
/// Prepared, temperature-independent state for a PKWang run: the cumulative-
/// energy ladder built once from the graph, the field's symmetrize flag, and —
/// for directed fields — the chosen <see cref="SymmetrizationRule"/> plus the
/// mirror-slot map that pairs each CSR slot with its reverse direction. Reused
/// across every temperature in a sweep; <c>Solve</c> applies only the
/// closed-form kernel and never rebuilds.
/// </summary>
public sealed class PKWangContext
{
    internal CsrGraph Graph { get; }
    internal double[] Hcum { get; }
    internal bool DirectedSymmetrize { get; }
    internal SymmetrizationRule Rule { get; }

    /// <summary>For directed fields: <c>Mirror[e]</c> is the CSR slot of the
    /// reverse direction of edge slot <c>e</c>. Null for symmetric fields.</summary>
    internal int[]? Mirror { get; }

    internal PKWangContext(
        CsrGraph graph,
        double[] hcum,
        bool directedSymmetrize,
        SymmetrizationRule rule,
        int[]? mirror)
    {
        Graph = graph;
        Hcum = hcum;
        DirectedSymmetrize = directedSymmetrize;
        Rule = rule;
        Mirror = mirror;
    }
}
