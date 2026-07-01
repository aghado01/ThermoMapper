using System;
using Graphs.Primitives;

namespace Graphs.Observables;

/// <summary>
/// Per-node sum of incident affinities:
/// <c>D_i = Σ_{j ∈ N(i)} G(i,j)</c>. The expected weighted degree of node
/// <c>i</c> in the affinity-weighted subgraph at the equilibrium temperature —
/// for Swendsen–Wang the random subgraph of formed bonds, for PKWang the
/// closed-form survival-weighted subgraph.
/// </summary>
/// <remarks>
/// <para><b>Interpretation.</b> Low <c>D_i</c> = node sits on fragile bridges
/// where few incident edges carry affinity at this T; high <c>D_i</c> = node
/// sits in a structurally dense core where most incident edges bond. Sometimes
/// described in TDA contexts as the "expected structural degree" or a
/// "bottleneck score," but the math is plain weighted degree on the
/// affinity-weighted graph.</para>
///
/// <para><b>Subject / Op decomposition.</b>
/// Subject: <c>Affinity</c> (the per-edge currency <see cref="Affinities.G"/>).
/// Op: <c>Degree</c> (sum over node's incident edges).</para>
///
/// <para><b>CSR-walk convention.</b> Only the upper-triangular CSR slots
/// (<c>j &gt; i</c>) carry currency; this implementation walks that pattern and
/// accumulates each edge's affinity into <i>both</i> endpoints' degree (each
/// undirected edge contributes to both nodes' expected degree).</para>
/// </remarks>
public sealed class AffinityDegree : IGraphSignal<Affinities>
{
    public double[] Compute(Affinities affinities, CsrGraph graph)
    {
        ArgumentNullException.ThrowIfNull(affinities);
        ArgumentNullException.ThrowIfNull(graph);

        double[] g = affinities.G;
        if (g.Length != graph.Targets.Length)
            throw new ArgumentException(
                $"Affinities.G length ({g.Length}) does not match CSR slot count " +
                $"({graph.Targets.Length}).", nameof(affinities));

        int n = graph.NodeCount;
        var result = new double[n];

        for (int i = 0; i < n; i++)
        {
            int rowEnd = graph.RowPointers[i + 1];
            for (int e = graph.RowPointers[i]; e < rowEnd; e++)
            {
                int j = graph.Targets[e];
                if (j <= i) continue;   // walk only the upper triangle (where the currency is defined)
                double a = g[e];
                result[i] += a;
                result[j] += a;
            }
        }
        return result;
    }
}
