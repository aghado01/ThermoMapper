using System;
using Graphs.Primitives;

namespace Graphs.Observables;

/// <summary>
/// Per-node sum of binary entropies of incident affinities:
/// <c>V_i = Σ_{j ∈ N(i)} H_2(G(i,j))</c>
/// where <c>H_2(p) = −p log₂ p − (1−p) log₂(1−p)</c>.
/// </summary>
/// <remarks>
/// <para><b>Interpretation.</b> Binary entropy is zero when the affinity
/// saturates at 0 or 1 (the edge is permanently broken or permanently fused —
/// empty space or deep core respectively) and peaks at <c>p = 0.5</c> (the edge
/// is "flickering" — half-bonded at equilibrium). Summed over a node's incident
/// edges, <c>V_i</c> is high for nodes whose neighborhoods sit on volatile
/// boundaries — the structural fracture lines of the dataset.</para>
///
/// <para><b>Subject / Op decomposition.</b>
/// Subject: <c>Affinity</c> (per-edge <see cref="Affinities.G"/>).
/// Op: <c>BinaryEntropySum</c> — "pointwise binary entropy on each edge, summed
/// over the node's incident edges to produce a per-node scalar."</para>
///
/// <para><b>Units.</b> Result is in <i>bits</i>. The natural-log binary entropy
/// <c>−(p ln p + (1−p) ln(1−p))</c> is computed internally for numerical
/// stability, then divided by <c>ln 2</c>.</para>
///
/// <para><b>CSR-walk convention.</b> Same upper-triangular pattern as
/// <see cref="AffinityDegree"/>; each edge's binary entropy is accumulated into
/// both endpoints.</para>
/// </remarks>
public sealed class AffinityBinaryEntropySum : IGraphSignal<Affinities>
{
    private static readonly double Ln2 = Math.Log(2.0);

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
                if (j <= i) continue;
                double hBits = BinaryEntropyNats(g[e]) / Ln2;
                result[i] += hBits;
                result[j] += hBits;
            }
        }
        return result;
    }

    /// <summary>
    /// <c>H_2(p)</c> in nats: <c>−p ln p − (1−p) ln(1−p)</c>. Zero at the
    /// endpoints (handled explicitly to avoid <c>log(0)</c>); <c>ln 2</c> at
    /// <c>p = 0.5</c>.
    /// </summary>
    private static double BinaryEntropyNats(double p)
    {
        if (p <= 0.0 || p >= 1.0) return 0.0;
        return -(p * Math.Log(p) + (1.0 - p) * Math.Log(1.0 - p));
    }
}
