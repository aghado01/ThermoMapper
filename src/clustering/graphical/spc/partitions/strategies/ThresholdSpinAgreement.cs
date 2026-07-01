using System;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Strategies;

/// <summary>
/// Threshold-and-connected-components partition keyed on per-edge
/// equilibrium spin agreement <c>⟨1_{s_i=s_j}⟩</c>. The canonical Blatt
/// 1996 friends-of-friends cut.
/// </summary>
/// <remarks>
/// <para><b>Algorithm.</b> For each undirected edge <c>(i,j)</c>,
/// declare <c>(i,j)</c> a friend iff <c>Alignments.G[e] &gt; Theta</c>.
/// Compute connected components on the friend subgraph via union-find;
/// each component is one output cluster.</para>
///
/// <para><b>Why it works.</b> At equilibrium in the superparamagnetic
/// phase, pairs that belong to the same physical cluster spend most
/// cycles in the same Potts color (<c>⟨1_{s_i=s_j}⟩ ≈ 1</c>); pairs in
/// different clusters align only by chance
/// (<c>⟨1_{s_i=s_j}⟩ ≈ 1/q</c>). A threshold between those two modes
/// sharply separates them; Blatt finds the cut is insensitive to θ
/// across roughly <c>0.2 ≤ θ ≤ 0.9</c>.</para>
///
/// <para><b>Peripheral capture.</b> Enable <see cref="PeripheralCapture"/> to
/// also union each node with its single highest-spin-agreement neighbor
/// regardless of θ (Domany1999 step 2). Default off (strict BWD1995 parity);
/// switch on for real-data runs where cluster density decreases toward the
/// perimeter (the 1999 Iris 25→2 unclassified improvement).</para>
///
/// <para><b>Requires.</b> The <see cref="Alignments"/> currency
/// (<c>SpinAgreementCount / DrawCount</c>); configure the sweep with
/// <see cref="Clustering.Graphical.SPC.Runtime.Core.Sampler.AccumulationSpec.Currencies"/>
/// to ensure it is collected.</para>
/// </remarks>
public sealed class ThresholdSpinAgreement : IPartitionStrategy
{
    /// <summary>
    /// Friends-of-friends threshold; edges with
    /// <c>G[e] &gt; Theta</c> get unioned. Defaults to 0.5;
    /// canonically insensitive across <c>[0.2, 0.9]</c>.
    /// </summary>
    public double Theta { get; init; } = 0.5;

    /// <summary>
    /// When <see langword="true"/>, also unions each node with its single
    /// highest-spin-agreement neighbor regardless of <see cref="Theta"/>
    /// (Domany1999 step 2 — peripheral capture). Default off.
    /// </summary>
    public bool PeripheralCapture { get; init; } = false;

    public Assignment Apply(CsrGraph graph, Affinities affinities, Alignments? alignments, CoMembership? coMembership = null)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(affinities);
        if (alignments is null)
            throw new InvalidOperationException(
                $"{nameof(ThresholdSpinAgreement)} requires the Alignments currency — " +
                "configure the sweep with AccumulationSpec.Currencies.");
        if (alignments.G.Length != graph.Targets.Length)
            throw new InvalidOperationException(
                $"Alignments.G length ({alignments.G.Length}) does not match " +
                $"CSR slot count ({graph.Targets.Length}).");
        if (Theta < 0.0 || Theta > 1.0)
            throw new InvalidOperationException($"Theta ({Theta}) must lie in [0, 1].");

        return AffinityThreshold.Connect(graph, alignments.G, Theta, PeripheralCapture);
    }
}
