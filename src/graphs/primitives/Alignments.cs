namespace Graphs.Primitives;

/// <summary>
/// Per-edge <b>alignment</b> currency — the bundle of pairwise spin alignments
/// <c>G_ij = ⟨δ(sᵢ,sⱼ)⟩ ∈ [0,1]</c> at one temperature: the equilibrium probability that
/// endpoints <c>i</c> and <c>j</c> occupy the <i>same Potts state</i> (are ferromagnetically
/// aligned), indexed by CSR slot (parallel to <see cref="CsrGraph.Targets"/>). Domany's
/// pair-correlation function for neighbouring spins — distinct from <see cref="Affinities"/>
/// (bond-survival strength): alignment counts state agreement, affinity counts bond formation.
/// </summary>
/// <remarks>
/// Alignment is the permutation-invariant order signal: the <c>q!</c> label symmetry of the Potts
/// model leaves a point→state map unidentified, but the pairwise <c>δ(sᵢ,sⱼ)</c> is invariant under
/// relabelling (it is equally McLachlan's fuzzy classification matrix / Bishop's responsibilities made
/// label-switch-safe). Swendsen–Wang mints it from <c>SpinAgreementCount / draws</c>; it is
/// <b>SW-native</b> — a forward solver (PKWang) draws no spins and has no analog, so unlike
/// <see cref="Affinities"/> this currency is not universal across samplers. Only the <c>j &gt; i</c>
/// CSR slots carry meaning; the mirror half stays zero.
/// </remarks>
public sealed record Alignments
{
    /// <summary>Temperature these alignments were evaluated at.</summary>
    public required double Temperature { get; init; }

    /// <summary>Per-CSR-slot alignment G; meaningful at the <c>j &gt; i</c> slots.</summary>
    public required double[] G { get; init; }

    /// <summary>Replica index for metadata/seed provenance (defaults 0).</summary>
    public int ReplicaIndex { get; init; }
}
