using System.Collections.Generic;
using Clustering.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Hierarchical;

/// <summary>
/// Shared nesting check for any cold→hot stack of partitions — the
/// don't-warp gate (<see cref="PartitionHierarchy.NestingHolds"/>) that
/// decides whether a T-stack collapses to a single-linkage
/// <c>Clustering.Dendrograms.Dendrogram</c> (the degenerate-nested
/// projection) or stays its own member (the lineage-persistence-resolved
/// <see cref="PartitionHierarchy"/>). Hoisted out of
/// <see cref="BlattPartitionStrategy"/> so the per-phase reading and the
/// dense per-T reading (<c>DenseTStack</c>) check nesting one way.
/// </summary>
public static class PartitionNesting
{
    /// <summary>
    /// True when every cluster at a hotter level is a subset of some cluster
    /// at the immediately colder level — the Blatt nesting invariant.
    /// <paramref name="levels"/> must be ordered cold→hot (ascending T).
    /// Empty or single-level stacks trivially nest. O(N) per consecutive pair.
    /// </summary>
    public static bool Holds(IReadOnlyList<HierarchyLevel> levels)
    {
        if (levels is null || levels.Count < 2) return true;
        for (int i = 0; i < levels.Count - 1; i++)
        {
            int[] coldLabels = levels[i].Partition.Labels;
            int[] hotLabels = levels[i + 1].Partition.Labels;
            if (coldLabels.Length != hotLabels.Length) return false;

            // Each hot cluster id must map to a single cold cluster id; a hot
            // cluster straddling two cold clusters breaks nesting. Unassigned
            // points carry no constraint (they belong to no cluster at that
            // level) and are skipped on whichever side abstains.
            var hotToCold = new Dictionary<int, int>();
            for (int k = 0; k < hotLabels.Length; k++)
            {
                int hot = hotLabels[k];
                int cold = coldLabels[k];
                if (hot == Assignment.Unassigned || cold == Assignment.Unassigned) continue;
                if (!hotToCold.TryAdd(hot, cold) && hotToCold[hot] != cold)
                    return false;
            }
        }
        return true;
    }
}

/// <summary>
/// One stable super-paramagnetic phase from the Blatt picture: the
/// representative temperature, the partition read off at that
/// temperature, and a short provenance string documenting how the level
/// was produced (e.g. <c>"phase-midpoint between χ_m peaks at T=0.05
/// and T=0.18"</c>).
/// </summary>
public sealed record HierarchyLevel(
    double     Temperature,
    Assignment Partition,
    string?    Provenance = null);

/// <summary>
/// Ordered sequence of <see cref="HierarchyLevel"/>s spanning the
/// temperature axis from cold (ordered phase, few large clusters) to
/// hot (paramagnetic, many small or singleton clusters). The
/// hierarchical output of an <see cref="IHierarchicalPartitionStrategy"/>.
/// </summary>
/// <remarks>
/// <para><b>Ordering.</b> <see cref="Levels"/> is sorted ascending in
/// <see cref="HierarchyLevel.Temperature"/>. The cold-end level is
/// <c>Levels[0]</c>; the hot-end level is <c>Levels[^1]</c>. Consumers
/// that want to walk merges (cold→hot) should iterate ascending;
/// consumers that want to walk splits (hot→cold) iterate descending.</para>
///
/// <para><b>Nesting.</b> The Blatt picture predicts strict nesting:
/// every cluster at a hotter level is a subset of some cluster at a
/// colder level. <see cref="NestingHolds"/> records whether that
/// invariant actually held for this run; an undersampled or noisy
/// sweep may produce non-nesting levels that are still individually
/// useful as per-phase partitions, just not stackable as a single
/// merge tree.</para>
///
/// <para><b>Not a <c>Clustering.Dendrograms.Dendrogram</c>.</b> Both
/// shapes are hierarchical trees of clusters with a scalar cost axis,
/// but the semantics differ at the foundation. A
/// <c>Clustering.Dendrograms.Dendrogram</c> records pairwise
/// single-linkage merges between subtrees of data points; its cost
/// axis is a dissimilarity (mutual-reachability distance, ΔH, …) that
/// is intrinsic to the geometry of the data. A
/// <see cref="PartitionHierarchy"/> records phase-resolved partitions
/// along an SPC temperature sweep; its cost axis is the thermodynamic
/// control parameter and each level is the equilibrium clustering at
/// that T, not a pairwise merge event. Don't convert one into the
/// other — the shared "tree of clusters with a y-axis" is a napkin-
/// sketch coincidence, not a structural identity.</para>
///
/// <para><b>Topology axis.</b> <see cref="TopologyAxis"/> defaults to
/// <c>"temperature"</c> — the y-axis label for plot consumers. Custom
/// detectors that fold T into a transformed coordinate (e.g. log T,
/// inverse temperature β) should override.</para>
/// </remarks>
public sealed record PartitionHierarchy(
    IReadOnlyList<HierarchyLevel> Levels,
    bool                          NestingHolds,
    string                        TopologyAxis = "temperature")
{
    /// <summary>Number of phases / hierarchy levels.</summary>
    public int Count => Levels.Count;

    /// <summary>True when no phases were detected (typically a single-
    /// component graph or a sweep that didn't bracket any
    /// transition).</summary>
    public bool IsEmpty => Levels.Count == 0;
}
