using System;
using System.Collections.Generic;
using Clustering.Dendrograms;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Hierarchical;

/// <summary>
/// Swendsen–Wang adapter for the shared <see cref="HierarchyEom"/> core: pools a
/// rich sweep's frames into the dense partition stack
/// (<see cref="SweepDenseTStack"/>) and the per-node landscape
/// (<see cref="SweepLandscapes"/>), then resolves.
/// </summary>
public static class SweepHierarchyEom
{
    /// <summary>
    /// Resolve a rich sweep's frames. Requires <c>AccumulationSpec.CoMembership</c>
    /// (the partition columns) and the landscape sink's per-node accumulation
    /// (<c>ClusterSizeLandscape</c> for the default sink).
    /// </summary>
    /// <param name="restoreNesting">Restore the FKG-monotone premise with the
    /// cheap refinement pre-pass before bridging (default true). Off = strict
    /// gate: bridge only if the raw stack nests, else flag + abstain.</param>
    /// <param name="minClusterSize">Selection eligibility floor (taming the
    /// persistence-selected micro-clump tail); displaced members fall to the
    /// periphery completion.</param>
    public static HierarchyEomResult Resolve(
        CsrGraph graph,
        IReadOnlyList<Accumulator> frames,
        double theta = 0.5,
        SwLandscapeSink sink = SwLandscapeSink.MeanClusterSize,
        string graphId = "unspecified",
        int minClusterSize = 1,
        bool restoreNesting = true,
        ThermalPeripheryCompletion completion = ThermalPeripheryCompletion.None)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(frames);

        PartitionHierarchy raw = SweepDenseTStack.Build(graph, frames, theta);
        Landscape landscape = SweepLandscapes.FromFrames(frames, sink, graphId);
        return HierarchyEom.Resolve(graph, raw, landscape, minClusterSize, restoreNesting, completion);
    }
}
