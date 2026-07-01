using System;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Strategies;

/// <summary>
/// Threshold-and-connected-components partition keyed on per-edge
/// equilibrium bond frequency <c>⟨ω(e)⟩</c>. The FK-direct counterpart
/// to <see cref="ThresholdSpinAgreement"/>.
/// </summary>
/// <remarks>
/// <para><b>Algorithm.</b> For each undirected edge <c>(i,j)</c>,
/// declare <c>(i,j)</c> an active bond iff <c>Affinities.G[e] &gt; Theta</c>.
/// Compute connected components on the active-bond subgraph via union-find.</para>
///
/// <para><b>Relationship to <see cref="ThresholdSpinAgreement"/>.</b>
/// In the joint Fortuin-Kasteleyn / Potts representation,
/// <c>⟨ω(e)⟩ = p_e · ⟨1_{s_i=s_j}⟩</c> where
/// <c>p_e = 1 − exp(−J_e/T)</c>. Spin agreement is what the canonical
/// Blatt cut uses; bond frequency is what FK theory and percolation
/// diagnostics consume directly. For homogeneous J the two ratios are
/// proportional, but for heterogeneous J they diverge — bond frequency
/// is the cleaner FK-aligned choice on weighted proximity graphs where
/// the coupling varies across edges.</para>
///
/// <para><b>Requires.</b> The <see cref="Affinities"/> currency
/// (<c>BondFormedCount / DrawCount</c>), always present at the chosen-T
/// equilibrium pass.</para>
/// </remarks>
public sealed class ThresholdBondFrequency : IPartitionStrategy
{
    /// <summary>
    /// Bond-activity threshold; edges with <c>G[e] &gt; Theta</c> are
    /// unioned. Defaults to 0.5.
    /// </summary>
    public double Theta { get; init; } = 0.5;

    /// <summary>
    /// When <see langword="true"/>, also unions each node with its single
    /// highest-affinity neighbor regardless of <see cref="Theta"/> (Domany1999
    /// step 2 — peripheral capture). Default off (strict BWD1995 parity).
    /// </summary>
    public bool PeripheralCapture { get; init; } = false;

    public Assignment Apply(CsrGraph graph, Affinities affinities, Alignments? alignments, CoMembership? coMembership = null)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(affinities);
        if (affinities.G.Length != graph.Targets.Length)
            throw new InvalidOperationException(
                $"Affinities.G length ({affinities.G.Length}) does not match " +
                $"CSR slot count ({graph.Targets.Length}).");
        if (Theta < 0.0 || Theta > 1.0)
            throw new InvalidOperationException($"Theta ({Theta}) must lie in [0, 1].");

        return AffinityThreshold.Connect(graph, affinities.G, Theta, PeripheralCapture);
    }
}
