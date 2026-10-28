using System;
using System.Collections.Generic;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Dendrograms;

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

/// <summary>How (whether) to complete the abstains after selection.</summary>
public enum ThermalPeripheryCompletion
{
    /// <summary>Leave abstains as <see cref="Assignment.Unassigned"/>.</summary>
    None,

    /// <summary>Modal ascent on the COLDEST landscape column (the most-coupled
    /// slice): height-greedy, valley-respecting, gauge-free.</summary>
    Ascend,
}

/// <summary>
/// The thermodynamic-EOM composition: a thermal dendrogram (structure, e.g. from
/// stacked per-T edge columns via <see cref="ThermalDendrogram"/>) + a per-node
/// landscape (height) → excess-of-mass walk → <see cref="Assignment"/>.
/// Consumes the whole temperature axis, not a chosen-T slice; the flat threshold
/// cut is its single-T degenerate sibling. Unselected leaves come back
/// <see cref="Assignment.Unassigned"/> — the honest abstain a periphery policy
/// may complete downstream.
/// </summary>
public static class ThermalEom
{
    /// <summary>
    /// Producer-agnostic core: resolve a thermal merge tree + per-node landscape
    /// into an excess-of-mass assignment, optionally completing the periphery by
    /// modal ascent on the coldest (most-coupled) landscape column. The structure
    /// and height are currencies, not frames — SW builds them from sampled
    /// accumulators (SPC's <c>SweepThermalEom</c>); the cumulative-field solver builds them exactly
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
