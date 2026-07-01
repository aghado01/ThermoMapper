using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Execution;

/// <summary>
/// One atomic unit of work for the SPC executor.
/// </summary>
public sealed record SpcTaskSpec
{
    public required double          Temperature  { get; init; }
    public required int             ReplicaIndex { get; init; }
    public required int             Q            { get; init; }

    /// <summary>
    /// Declares what sufficient-statistics this task's run accumulates — the per-edge
    /// currencies (<c>Affinities</c>/<c>Alignments</c>) and the per-node landscapes.
    /// Scalar moments + the cluster-size histogram are always collected.
    /// </summary>
    public AccumulationSpec Accumulation { get; init; }

    /// <summary>MC budget bundling burn-in and measurement cycles.</summary>
    public RunBudget Budget { get; init; }

    /// <summary>
    /// Path the persistent sink writes the checkpoint to. Required by
    /// disk-backed sinks (e.g.
    /// <see cref="Sinks.SpcxDiskFrameSink"/>); may be left null when the
    /// executor is wired to a non-persistent sink
    /// (<see cref="Sinks.NullFrameSink"/>,
    /// <see cref="Sinks.InMemoryFrameSink"/>) — those sinks ignore the
    /// path entirely.
    /// </summary>
    public string? CheckpointPath { get; init; }

    public int? BaseSeed { get; init; }
}
