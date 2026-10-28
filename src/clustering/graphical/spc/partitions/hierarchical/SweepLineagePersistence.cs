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
    /// <param name="minClusterSize">Clusters below this size are ignored when
    /// building lineages (the overclustering-tail filter; default 3).</param>
    /// <param name="splitShare">Second-successor share above which a step is a
    /// split, not erosion. Null ⇒ gap-based on the observed distribution
    /// (fallback 0.25).</param>
    /// <param name="temperatureWindow">Restrict to [Lo, Hi] — the SP-plateau
    /// bound (regime border: T_fs/T_ps from <see cref="SpcProfileAnalysis.SpPlateau"/>).
    /// Null ⇒ the full grid.</param>
    /// <param name="selectTopK">Force exactly this many lineages instead of the
    /// gap rule. Null ⇒ gap-based selection.</param>
    public static LineagePersistenceResult Resolve(
        CsrGraph graph,
        IReadOnlyList<Accumulator> frames,
        double theta = 0.5,
        int minClusterSize = 3,
        double? splitShare = null,
        (double Lo, double Hi)? temperatureWindow = null,
        int? selectTopK = null)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(frames);

        int q = frames[0].Q;
        var (temperatures, delta) = SweepEdgeCurves.CoMembershipDelta(frames);
        PartitionHierarchy stack = DenseTStack.Build(graph, temperatures, delta, theta);

        return LineagePersistence.Resolve(graph, stack, delta, q, minClusterSize, splitShare, temperatureWindow, selectTopK);
    }
}
