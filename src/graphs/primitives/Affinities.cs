namespace Graphs.Primitives;

/// <summary>
/// Universal cross-sampler currency — the per-edge bundle of <b>affinities</b>
/// <c>G_ij ∈ [0,1]</c> at one temperature: the bond-activity / edge-survival
/// strength between two endpoints, indexed by CSR slot (parallel to
/// <see cref="CsrGraph.Targets"/>) and thresholded to cut clusters.
/// "Affinity" is the categorical currency; this plural type is the concrete bundle a
/// sampler emits (paralleling the accumulators it implements).
/// </summary>
/// <remarks>
/// Every sampler emits this tier, computed however it likes: Swendsen–Wang derives
/// it from FK bond/cycle counts (a same-cluster bond frequency); PKWang writes it
/// directly from the closed-form survival kernel <c>G = 1 − exp(−Hcum/T)</c> — no
/// Monte-Carlo draws. It is <i>not</i> the spin–spin correlation
/// <c>⟨δ(sᵢ,sⱼ)⟩</c> (that is the separate <see cref="Alignments"/> channel);
/// "affinity" claims only the per-edge bonding strength the partition step consumes,
/// not how it was produced. Only the <c>j &gt; i</c> CSR slots carry meaning; the
/// mirror half stays zero, matching the upper-triangular walk every partition consumer uses.
/// </remarks>
public sealed record Affinities
{
    /// <summary>Temperature these affinities were evaluated at.</summary>
    public required double Temperature { get; init; }

    /// <summary>Per-CSR-slot affinity G; meaningful at the <c>j &gt; i</c> slots.</summary>
    public required double[] G { get; init; }

    /// <summary>Replica index for metadata/seed provenance (PKWang is replica-free; defaults 0).</summary>
    public int ReplicaIndex { get; init; }
}
