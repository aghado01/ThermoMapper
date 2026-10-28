using System;
using System.Collections.Generic;
using Clustering.Dendrograms;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Hierarchical;

/// <summary>
/// Swendsen–Wang adapter for the shared <see cref="DenseTStack"/> cut: pools a
/// rich sweep's frames into per-T co-membership columns
/// (<see cref="SweepEdgeCurves.CoMembershipDelta"/>) and cuts the dense stack.
/// </summary>
public static class SweepDenseTStack
{
    /// <summary>
    /// Build the dense per-T partition stack from rich sweep frames. Requires
    /// <c>AccumulationSpec.CoMembership</c> (the per-edge co-membership counts);
    /// frames are pooled across replicas per temperature inside
    /// <see cref="SweepEdgeCurves.CoMembershipDelta"/>.
    /// </summary>
    /// <param name="graph">The CSR graph the frames were sampled on; its edge
    /// slots index the δ̄ columns.</param>
    /// <param name="frames">A rich sweep's accumulators (one or more per T).</param>
    /// <param name="theta">Bond threshold on δ̄ (BWD step g; default 0.5).</param>
    /// <param name="peripheralCapture">Domany step 2 — union each node with its
    /// max-δ̄ neighbor regardless of θ, applied per level. Default off.</param>
    public static PartitionHierarchy Build(
        CsrGraph graph,
        IReadOnlyList<Accumulator> frames,
        double theta = 0.5,
        bool peripheralCapture = false)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(frames);

        var (temperatures, deltaByGridPoint) = SweepEdgeCurves.CoMembershipDelta(frames);
        return DenseTStack.Build(graph, temperatures, deltaByGridPoint, theta, peripheralCapture);
    }
}
