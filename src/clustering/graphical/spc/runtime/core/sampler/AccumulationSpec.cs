namespace Clustering.Graphical.SPC.Runtime.Core.Sampler;

/// <summary>
/// Declarative specification of <b>what sufficient-statistics a Swendsen–Wang run accumulates</b> —
/// the single config surface that replaces the old <c>PottsObservableTier</c> enum and the separate
/// per-node track-bools. Each gate names the currency or landscape it produces, not a procedural role:
/// "what data is collected" is a property of the run's config, never of whether it is a "probe" or a
/// "final" pass.
/// </summary>
/// <remarks>
/// <para>Scalar moments (FK χ, specific heat, magnetization) and the cluster-size histogram are free
/// byproducts of the union-find pass — always collected, not gated here. The four optional dimensions:</para>
/// <list type="bullet">
///   <item><see cref="Affinities"/> / <see cref="Alignments"/> — the per-edge currency precursors;
///   hot-loop scatter-writes, so they map to an <see cref="ISwConfig"/> specialization
///   (compile-time monomorphization, independently gated).</item>
///   <item><see cref="ClusterSizeLandscape"/> / <see cref="OrderLandscape"/> — the per-node landscapes
///   (un-reduced χ and order parameter); an <c>O(N)</c> post-pass, so runtime-gated.</item>
/// </list>
/// <para>The boundary mapping (spec → specialization + runtime flags) lives at the cold run-setup seam
/// (<c>SwRunner</c>); the hot loop never sees this type.</para>
/// </remarks>
public readonly record struct AccumulationSpec
{
    /// <summary>Accumulate per-edge bond-survival counts → the <c>Affinities</c> currency.</summary>
    public bool Affinities { get; init; }

    /// <summary>Accumulate per-edge spin-agreement counts → the <c>Alignments</c> currency.</summary>
    public bool Alignments { get; init; }

    /// <summary>Accumulate the per-node mean-cluster-size landscape (un-reduced χ).</summary>
    public bool ClusterSizeLandscape { get; init; }

    /// <summary>Accumulate the per-node giant-participation landscape (un-reduced order parameter).</summary>
    public bool OrderLandscape { get; init; }

    /// <summary>
    /// Accumulate per-edge co-membership counts → the <c>CoMembership</c> currency
    /// (<c>⟨n_ij⟩</c> — fraction of draws where i and j are in the same bond cluster,
    /// regardless of whether the direct bond froze; lower variance than <c>Affinities</c>).
    /// O(E) post-pass per draw; does NOT earn a JIT-specialization gate (no hot-loop scatter).
    /// </summary>
    public bool CoMembership { get; init; }

    /// <summary>Moments + cluster-size histogram only (nothing optional) — the lightweight default.</summary>
    public static AccumulationSpec None => default;

    /// <summary>Both per-edge currencies + co-membership: all three estimators for chosen-T comparison.</summary>
    public static AccumulationSpec Currencies => new() { Affinities = true, Alignments = true, CoMembership = true };
}
