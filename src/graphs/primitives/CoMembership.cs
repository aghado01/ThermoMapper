namespace Graphs.Primitives;

/// <summary>
/// Per-edge co-membership frequency at one temperature: <c>⟨n_ij⟩ = P(same cluster)</c>
/// — the fraction of SW draws in which nodes <c>i</c> and <c>j</c> fell in the same bond
/// cluster, regardless of whether the direct (i,j) bond froze (Niedermayer 1988 / eq. 4 of
/// Blatt–Wiedemann–Domany 1996). This is the improved estimator: lower variance than the raw
/// bond-frequency <see cref="Affinities"/> because it captures transitive co-clustering via
/// multi-hop paths, not just the single direct bond.
/// </summary>
/// <remarks>
/// Indexed by CSR slot parallel to <see cref="CsrGraph.Targets"/>; meaningful at
/// <c>j &gt; i</c> slots only (mirror slots carry zero, matching the upper-triangle
/// convention). Because the post-pass walks <see cref="CsrGraph.UndirectedEdges"/> (j &gt; i),
/// lower-index neighbors are captured correctly via the higher-index endpoint's row.
/// </remarks>
public sealed record CoMembership
{
    /// <summary>Temperature these co-membership rates were evaluated at.</summary>
    public required double Temperature { get; init; }

    /// <summary>Potts q used in the SW run that produced these counts.</summary>
    public int Q { get; init; }

    /// <summary>Per-CSR-slot co-membership fraction G; meaningful at the <c>j &gt; i</c> slots.</summary>
    public required double[] G { get; init; }

    /// <summary>Replica index for metadata/seed provenance (defaults 0).</summary>
    public int ReplicaIndex { get; init; }
}
