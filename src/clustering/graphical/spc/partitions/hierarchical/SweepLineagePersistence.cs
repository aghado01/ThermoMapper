using System;
using System.Collections.Generic;
using Clustering.Dendrograms;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Hierarchical;

/// <summary>
/// Swendsen–Wang adapter for the shared <see cref="LineagePersistence"/> core:
/// pools a rich sweep's frames into per-T co-membership columns
/// (<see cref="SweepEdgeCurves.CoMembershipDelta"/>), cuts the dense stack, and
/// resolves lineages over it.
/// </summary>
public static class SweepLineagePersistence
{
    /// <summary>
    /// Resolve a rich sweep's frames by lineage persistence. Requires
    /// <c>AccumulationSpec.CoMembership</c> (the per-edge co-membership counts —
    /// the partition columns AND the decidedness field). No landscape needed.
    /// </summary>
    /// <param name="theta">The δ̄ cut of the dense stack.</param>
    /// <param name="config">What to read and how. Null reads the longest-lived lineage on
    /// each branch of the full grid, with <see cref="LineagePersistenceConfig"/>'s defaults.</param>
    public static LineagePersistenceResult Resolve(
        CsrGraph graph,
        IReadOnlyList<Accumulator> frames,
        double theta = 0.5,
        LineagePersistenceConfig? config = null)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(frames);

        int q = frames[0].Q;
        var (temperatures, delta) = SweepEdgeCurves.CoMembershipDelta(frames);
        PartitionHierarchy stack = DenseTStack.Build(graph, temperatures, delta, theta);

        return LineagePersistence.Resolve(graph, stack, delta, q, config);
    }
}
