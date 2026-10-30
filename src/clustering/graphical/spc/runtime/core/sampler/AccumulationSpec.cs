namespace Clustering.Graphical.SPC.Runtime.Core.Sampler;

/// <summary>
/// Declarative specification of <b>what sufficient-statistics a Swendsen–Wang run accumulates</b> —
/// the single config surface that replaces the old <c>PottsObservableTier</c> enum and the separate
/// per-node track-bools. Each gate names the currency or landscape it produces, not a procedural role:
/// "what data is collected" is a property of the run's config, never of whether it is a "probe" or a
/// "final" pass.
/// </summary>
/// <remarks>
/// <para>Scalar moments (FK χ, specific heat, magnetization), the per-draw record (E, m, Σ|c|², c_max)
/// and the cluster-size histogram are free byproducts of the union-find pass — always collected, not
/// gated here. The optional dimensions:</para>
/// <list type="bullet">
///   <item><see cref="Affinities"/> / <see cref="Alignments"/> — the per-edge currency precursors;
///   hot-loop scatter-writes, so they map to an <see cref="ISwConfig"/> specialization
///   (compile-time monomorphization, independently gated).</item>
///   <item><see cref="ClusterSizeLandscape"/> / <see cref="OrderLandscape"/> — the per-node landscapes
///   (un-reduced χ and order parameter); an <c>O(N)</c> post-pass, so runtime-gated.</item>
///   <item><see cref="CoMembership"/> — an <c>O(E)</c> post-pass, runtime-gated.</item>
///   <item><see cref="BurnInTrace"/> — a diagnostic copy of the burn-in draws' records, runtime-gated.</item>
///   <item><see cref="ExactClusterSizeHistogram"/> — counts at each exact cluster size; an
///   <c>O(clusters)</c> scatter beside the log-binned one, runtime-gated.</item>
///   <item><see cref="PartitionEvery"/> / <see cref="SpinsEvery"/> — every k-th draw's FK partition and spins;
///   an <c>O(N)</c> pass on kept draws only, runtime-gated.</item>
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

    /// <summary>
    /// Keep the per-draw record (E, m, Σ|c|², c_max) of the burn-in draws, which the reset after
    /// burn-in would otherwise discard, so the burn-in length can be checked rather than assumed.
    /// Diagnostic; off by default.
    /// </summary>
    public bool BurnInTrace { get; init; }

    /// <summary>
    /// Accumulate 64-bit counts of the FK clusters at each exact size 1..N beside the log-binned
    /// histogram, whose bins hold counts only: every cluster functional Σ_k g(c_k) is linear in these
    /// counts, so its mean is exact (the <c>Graphs.Observables.ClusterSizeSums</c> reductions). Off by default.
    /// </summary>
    public bool ExactClusterSizeHistogram { get; init; }

    /// <summary>
    /// Keep the FK partition of every k-th measured draw (k = this value; 0 = off), canonically labelled by each
    /// cluster's smallest node index: the partition draws that VI-loss estimates, credible balls and replica link
    /// overlaps need. N 32-bit labels per kept draw, so thin to about one draw per autocorrelation time. Off by default.
    /// </summary>
    public int PartitionEvery { get; init; }

    /// <summary>
    /// Keep the spins of every k-th measured draw (k = this value; 0 = off), one byte per node, so q ≤ 256: the
    /// Markov state, for Potts overlaps on spins rather than on FK clusters. Off by default.
    /// </summary>
    public int SpinsEvery { get; init; }

    /// <summary>Moments, the per-draw record and the cluster-size histogram only (nothing optional) — the lightweight default.</summary>
    public static AccumulationSpec None => default;

    /// <summary>Both per-edge currencies + co-membership: all three estimators for chosen-T comparison.</summary>
    public static AccumulationSpec Currencies => new() { Affinities = true, Alignments = true, CoMembership = true };
}
