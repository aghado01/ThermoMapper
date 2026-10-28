using System;
using System.Collections.Generic;
using Clustering.Dendrograms;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Thermal;

/// <summary>
/// Swendsen–Wang adapter for the shared <see cref="ThermalEom"/> core: pools a
/// rich sweep's frames into the two currencies the walk consumes — the thermal
/// dendrogram from stacked per-T co-membership columns
/// (<see cref="SweepEdgeCurves"/> → <see cref="ThermalDendrogram"/>) and the
/// per-node landscape (<see cref="SweepLandscapes"/>).
/// </summary>
public static class SweepThermalEom
{
    /// <summary>
    /// Resolve a rich sweep's frames. Requires <c>AccumulationSpec.CoMembership</c>
    /// (the structure columns) and the landscape sink's per-node accumulation
    /// (<c>ClusterSizeLandscape</c> for the default sink).
    /// </summary>
    /// <param name="minClusterSize">Selection eligibility floor (taming the
    /// persistence-selected micro-clump tail); displaced members fall to the
    /// periphery completion.</param>
    public static ThermalEomResult Resolve(
        CsrGraph graph,
        IReadOnlyList<Accumulator> frames,
        double theta = 0.5,
        SwLandscapeSink sink = SwLandscapeSink.MeanClusterSize,
        string graphId = "unspecified",
        int minClusterSize = 1,
        ThermalPeripheryCompletion completion = ThermalPeripheryCompletion.None)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(frames);

        var (temperatures, delta) = SweepEdgeCurves.CoMembershipDelta(frames);
        Dendrogram dendrogram = ThermalDendrogram.FromEdgeCurves(graph, temperatures, delta, theta);
        Landscape landscape = SweepLandscapes.FromFrames(frames, sink, graphId);

        return ThermalEom.Resolve(graph, dendrogram, landscape, minClusterSize, completion);
    }
}
