using System;
using System.Collections.Generic;
using Clustering.Dendrograms;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Hierarchical;

/// <summary>
/// Everything the Blatt/Domany hierarchical resolution produced, kept whole for
/// inspection: the dense T-stack (structure read off the partitions), whether
/// the RAW sampled stack nested, the bridged dendrogram (null when a non-nested
/// stack was flagged rather than restored), the thermal landscape (height), the
/// walk report, the selection, and the resolved currency.
/// </summary>
/// <param name="Stack">The dense per-T partition stack actually walked (after
/// any nesting restoration).</param>
/// <param name="RawNestingHeld">Did the sampled stack nest BEFORE restoration —
/// the FKG-gate diagnostic. For ferromagnetic SPC the ideal sweep is FKG-
/// monotone ⇒ nested; a false here is (almost always) finite-sampling noise a
/// BARS monotonizer would remove.</param>
/// <param name="Restored">True when the cheap monotone nesting restoration was
/// applied to make the stack a refinement chain.</param>
public sealed record HierarchyEomResult(
    PartitionHierarchy Stack,
    bool RawNestingHeld,
    bool Restored,
    Dendrogram? Dendrogram,
    Landscape? Landscape,
    ClusterWalkReport? Walk,
    bool[]? Selected,
    Assignment Assignment);

/// <summary>
/// Track 1 — the classical BWD/Domany hierarchical T-stack resolver. Reads the
/// canonical dendrogram-across-T off the partition stack itself: dense per-T
/// cuts (<see cref="DenseTStack"/>) → bridge the nested case to a
/// <c>Clustering.Dendrograms.Dendrogram</c>
/// (<see cref="PartitionHierarchyDendrogram"/>) → excess-of-mass walk over the
/// thermal landscape (<see cref="LandscapeWalk"/>) → <see cref="Assignment"/>.
/// </summary>
/// <remarks>
/// <para><b>Distinct from <c>ThermalEom</c>.</b> Both end in the same EOM walk,
/// but the STRUCTURE differs at the source: thermal-EOM builds the merge tree
/// from per-edge co-membership curves G_e(T); this builds it from the discrete
/// per-T <i>partitions</i> — the Domany Fig-5 dendrogram. Kept a separate
/// resolver on purpose (the don't-warp discipline: <c>PartitionHierarchy</c> is
/// its own family member).</para>
///
/// <para><b>The FKG nesting gate.</b> For ferromagnetic SPC the ideal sweep is
/// FKG-monotone ⇒ nested ⇒ a dendrogram. A sampled stack that fails
/// <see cref="PartitionHierarchy.NestingHolds"/> is (almost always) a
/// finite-sampling artifact. The resolver therefore either (a) restores the
/// premise with a cheap monotone pre-pass — split any hot cluster straddling a
/// cold boundary, trusting the colder/more-stable partition (the
/// isotonic/PAVA-flavored fallback, NOT zigzag machinery) — and bridges; or
/// (b) with <paramref name="restoreNesting"/> off, bridges only when the raw
/// stack nests and otherwise FLAGS and abstains (genuine non-nesting is
/// telos-tier, not this baseline).</para>
/// </remarks>
public static class HierarchyEom
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
        Thermal.ThermalPeripheryCompletion completion = Thermal.ThermalPeripheryCompletion.None)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(frames);

        PartitionHierarchy raw = DenseTStack.Build(graph, frames, theta);
        Landscape landscape = SweepLandscapes.FromFrames(frames, sink, graphId);
        return Resolve(graph, raw, landscape, minClusterSize, restoreNesting, completion);
    }

    /// <summary>
    /// Producer-agnostic core: resolve a dense partition stack + per-node
    /// landscape. SW builds both from sampled frames (the overload above); a
    /// solver / BARS-monotonized columns build them upstream and feed this core.
    /// </summary>
    public static HierarchyEomResult Resolve(
        CsrGraph graph,
        PartitionHierarchy rawStack,
        Landscape landscape,
        int minClusterSize = 1,
        bool restoreNesting = true,
        Thermal.ThermalPeripheryCompletion completion = Thermal.ThermalPeripheryCompletion.None)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(rawStack);
        ArgumentNullException.ThrowIfNull(landscape);

        bool rawNestingHeld = rawStack.NestingHolds;

        PartitionHierarchy stack = rawStack;
        bool restored = false;
        if (!rawNestingHeld && restoreNesting)
        {
            stack = RestoreNesting(rawStack);
            restored = true;
        }

        if (!stack.NestingHolds)
        {
            // Strict gate, premise not restored: flag and abstain. Genuine
            // non-nesting is the lineage-persistence resolver's / telos-tier's territory.
            var abstain = new int[landscape.NodeCount];
            Array.Fill(abstain, Assignment.Unassigned);
            return new HierarchyEomResult(
                stack, rawNestingHeld, restored,
                Dendrogram: null, Landscape: landscape, Walk: null, Selected: null,
                Assignment: new Assignment { Labels = abstain, Count = 0 });
        }

        Dendrogram dendrogram = PartitionHierarchyDendrogram.ToDendrogram(stack);
        ClusterWalkReport walk = LandscapeWalk.ClusterProfiles(dendrogram, landscape);
        bool[] selected = LandscapeWalk.SelectByExcessOfMass(
            dendrogram, walk.Mass, allowRoot: false, minClusterSize);
        Assignment assignment = LandscapeWalk.ToAssignment(dendrogram, selected);

        if (completion == Thermal.ThermalPeripheryCompletion.Ascend)
            assignment = PeripheryPolicies.Ascend(assignment, graph, landscape.ValuesByGridPoint[0]);

        return new HierarchyEomResult(
            stack, rawNestingHeld, restored, dendrogram, landscape, walk, selected, assignment);
    }

    /// <summary>
    /// The cheap monotone premise-restoration: rebuild the stack as a strict
    /// refinement chain so it nests. Walk cold→hot; at each step force the hot
    /// partition to refine the cold one by splitting any hot cluster along the
    /// colder cluster boundaries (label = pair (coldLabel, hotLabel)). Trusts
    /// the colder/more-stable partition for the coarse structure and only ever
    /// SPLITS (never merges across a cold boundary), so it cannot fabricate
    /// coarse structure — the conservative isotonic-flavored denoiser the FKG
    /// gate licenses. Unassigned points stay unassigned.
    /// </summary>
    internal static PartitionHierarchy RestoreNesting(PartitionHierarchy stack)
    {
        IReadOnlyList<HierarchyLevel> levels = stack.Levels;
        if (levels.Count < 2) return stack;

        int n = levels[0].Partition.Labels.Length;
        var restored = new List<HierarchyLevel>(levels.Count) { levels[0] };
        int[] coldLabels = levels[0].Partition.Labels;

        for (int li = 1; li < levels.Count; li++)
        {
            int[] hotLabels = levels[li].Partition.Labels;
            var refinedLabels = new int[n];
            var keyToDense = new Dictionary<(int, int), int>();
            int next = 0;
            for (int p = 0; p < n; p++)
            {
                if (hotLabels[p] == Assignment.Unassigned)
                {
                    refinedLabels[p] = Assignment.Unassigned;
                    continue;
                }
                // A point with no cold cluster (cold-unassigned) keeps the hot
                // label alone — it imposes no coarse constraint.
                var key = (coldLabels[p], hotLabels[p]);
                if (!keyToDense.TryGetValue(key, out int dense))
                {
                    dense = next++;
                    keyToDense[key] = dense;
                }
                refinedLabels[p] = dense;
            }

            var refined = new Assignment { Labels = refinedLabels, Count = next };
            restored.Add(new HierarchyLevel(
                Temperature: levels[li].Temperature,
                Partition:   refined,
                Provenance:  (levels[li].Provenance ?? "dense T-stack") + " (nesting-restored)"));
            coldLabels = refinedLabels;   // the refined hot level is the next cold reference
        }

        return new PartitionHierarchy(restored, NestingHolds: PartitionNesting.Holds(restored),
            TopologyAxis: stack.TopologyAxis);
    }
}
