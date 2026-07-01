using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Core.Sampler;

/// <summary>
/// Configuration for a single Potts sampler execution.
/// </summary>
public sealed record SwRunSpec
{
    /// <summary>
    /// Weighted CSR graph for the sampler.
    /// </summary>
    public required CsrGraph Graph { get; init; }

    /// <summary>
    /// Temperature in the Potts model energy units.
    /// </summary>
    public required double Temperature { get; init; }

    /// <summary>
    /// Number of Potts colors (q) used by the sampler.
    /// </summary>
    public required int Q { get; init; }

    /// <summary>
    /// Declares what sufficient-statistics this run accumulates — the per-edge
    /// currencies (<c>Affinities</c>/<c>Alignments</c>) and the per-node landscapes.
    /// Scalar moments and the cluster-size histogram are always collected. Defaults
    /// to <see cref="AccumulationSpec.None"/>.
    /// </summary>
    public AccumulationSpec Accumulation { get; init; }

    /// <summary>
    /// Optional deterministic seed for reproducible runs.
    /// </summary>
    public int? Seed { get; init; }

    /// <summary>
    /// MC budget: burn-in cycles discarded before accumulation, then measurement cycles collected.
    /// </summary>
    public required RunBudget Budget { get; init; }

    /// <summary>
    /// Replica index for metadata/seed derivation.
    /// </summary>
    public int ReplicaIndex { get; init; }
}
