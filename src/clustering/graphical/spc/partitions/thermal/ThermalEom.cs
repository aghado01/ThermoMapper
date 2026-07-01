using System;
using System.Collections.Generic;
using Clustering.Dendrograms;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Thermal;

/// <summary>
/// Everything the thermodynamic-EOM resolution produced, kept whole for
/// inspection: the structure (thermal dendrogram), the height (thermal
/// landscape), the walk report, the selection, and the resolved currency.
/// </summary>
public sealed record ThermalEomResult(
    Dendrogram Dendrogram,
    Landscape Landscape,
    ClusterWalkReport Walk,
    bool[] Selected,
    Assignment Assignment);

/// <summary>
/// The live thermodynamic-EOM composition over a rich sweep: stacked per-T
/// eq-4 co-membership columns → thermal dendrogram (structure) + pooled
/// per-node landscape (height) → excess-of-mass walk → <see cref="Assignment"/>.
/// NOT an <c>IPartitionStrategy</c> — it consumes the whole sweep, not the
/// chosen-T slice; the flat threshold cut is its single-T degenerate sibling.
/// Unselected leaves come back <see cref="Assignment.Unassigned"/> — the
/// honest abstain a periphery policy may complete downstream.
/// </summary>
/// <summary>How (whether) to complete the abstains after selection.</summary>
public enum ThermalPeripheryCompletion
{
    /// <summary>Leave abstains as <see cref="Assignment.Unassigned"/>.</summary>
    None,

    /// <summary>Modal ascent on the COLDEST landscape column (the most-coupled
    /// slice): height-greedy, valley-respecting, gauge-free.</summary>
    Ascend,
}

public static class ThermalEom
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

        return Resolve(graph, dendrogram, landscape, minClusterSize, completion);
    }

    /// <summary>
    /// Producer-agnostic core: resolve a thermal merge tree + per-node landscape
    /// into an excess-of-mass assignment, optionally completing the periphery by
    /// modal ascent on the coldest (most-coupled) landscape column. The structure
    /// and height are currencies, not frames — SW builds them from sampled
    /// accumulators (the overload above); the PKWang solver builds them exactly
    /// from its closed form. Both share this one walk.
    /// </summary>
    public static ThermalEomResult Resolve(
        CsrGraph graph,
        Dendrogram dendrogram,
        Landscape landscape,
        int minClusterSize = 1,
        ThermalPeripheryCompletion completion = ThermalPeripheryCompletion.None)
    {
        ArgumentNullException.ThrowIfNull(dendrogram);
        ArgumentNullException.ThrowIfNull(landscape);

        ClusterWalkReport walk = LandscapeWalk.ClusterProfiles(dendrogram, landscape);
        bool[] selected = LandscapeWalk.SelectByExcessOfMass(
            dendrogram, walk.Mass, allowRoot: false, minClusterSize);
        Assignment assignment = LandscapeWalk.ToAssignment(dendrogram, selected);

        if (completion == ThermalPeripheryCompletion.Ascend)
            assignment = PeripheryPolicies.Ascend(assignment, graph, landscape.ValuesByGridPoint[0]);

        return new ThermalEomResult(dendrogram, landscape, walk, selected, assignment);
    }
}
