using System;
using System.Collections.Generic;
using Clustering.Dendrograms;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Hierarchical;

/// <summary>
/// The <b>nested-degenerate bridge</b>: projects a strictly-nested SPC T-stack
/// (<see cref="PartitionHierarchy"/>) onto an N-leaf single-linkage-style
/// <see cref="Dendrogram"/> so the shared excess-of-mass walk
/// (<see cref="LandscapeWalk"/>) can resolve it — the canonical Domany Fig-5
/// dendrogram-across-T, read off the partition stack itself rather than the
/// per-edge co-membership curves (<c>ThermalDendrogram.FromEdgeCurves</c>).
/// </summary>
/// <remarks>
/// <para><b>Don't-warp gate.</b> This is a DECLARED lossy projection valid only
/// when <see cref="PartitionHierarchy.NestingHolds"/>. A non-nested stack
/// carries contest/frustration a dendrogram cannot represent (members regroup
/// as T varies); collapsing it would fabricate a tree. Such a stack is the
/// lineage-persistence resolver's territory (overlap-linked lineages), not this bridge —
/// so a non-nested input THROWS rather than warping. The genuinely-non-nested
/// handling (zigzag/formigram) is telos-tier, not this classical baseline.</para>
///
/// <para><b>Construction.</b> Process levels hot→cold (the cooling direction):
/// at each temperature, points sharing a cluster label are unioned, emitting a
/// binary merge node at that level's temperature (height = the decoupling T at
/// which the sub-clusters coalesce; an N-ary join becomes a tie cascade). Merge
/// heights therefore DESCEND in build order — the thermal orientation
/// <see cref="LandscapeWalk"/> detects; <c>CostAxis = "temperature"</c>, so the
/// axis-alignment law holds against a thermal landscape. Heights are
/// grid-quantized to the observed temperatures.</para>
///
/// <para><b>Forests are first-class.</b> A point never grouped with another at
/// any level (a permanent singleton / always-unassigned) stays an isolated leaf
/// — the result is then a forest (fewer than N−1 internal nodes). The walk
/// handles forests natively and such leaves resolve to
/// <see cref="Assignment.Unassigned"/>. A stack whose coldest level is a single
/// component yields a spanning tree (N−1 merges), on which
/// <see cref="Dendrogram.CutToK"/> reproduces the intermediate stages.</para>
/// </remarks>
public static class PartitionHierarchyDendrogram
{
    /// <summary>
    /// Bridge a nested stack to a thermal merge tree. Throws when the stack is
    /// not strictly nested (<see cref="PartitionHierarchy.NestingHolds"/> false)
    /// — the contest a dendrogram cannot hold.
    /// </summary>
    public static Dendrogram ToDendrogram(PartitionHierarchy hierarchy)
    {
        ArgumentNullException.ThrowIfNull(hierarchy);
        if (hierarchy.IsEmpty)
            throw new ArgumentException("Cannot bridge an empty hierarchy to a dendrogram.", nameof(hierarchy));
        if (!hierarchy.NestingHolds)
            throw new InvalidOperationException(
                "PartitionHierarchyDendrogram.ToDendrogram: the T-stack is NOT strictly nested. The " +
                "dendrogram bridge is the declared nested-degenerate projection only; a non-nested stack " +
                "carries contest a single-linkage tree cannot represent — resolve it with the " +
                "lineage-persistence resolver (or, telos-tier, zigzag/formigram), not this bridge.");

        IReadOnlyList<HierarchyLevel> levels = hierarchy.Levels;
        int n = levels[0].Partition.Labels.Length;
        if (n < 2)
            throw new ArgumentException("A dendrogram needs at least 2 leaves.", nameof(hierarchy));
        for (int li = 1; li < levels.Count; li++)
            if (levels[li].Partition.Labels.Length != n)
                throw new ArgumentException(
                    $"Level {li} has {levels[li].Partition.Labels.Length} points; expected {n} (all levels " +
                    "must partition the same point set).",
                    nameof(hierarchy));

        var uf = new UnionFind(n);
        var compNodeId = new int[n];
        for (int i = 0; i < n; i++) compNodeId[i] = i;   // each point is leaf id i
        int nextId = n;
        var merges = new List<DendrogramNode>(n - 1);

        // Hot → cold: the hottest level a pair shares a label is the temperature
        // at which they coalesce as the system cools (the decoupling T).
        var anchorByLabel = new Dictionary<int, int>();
        for (int li = levels.Count - 1; li >= 0; li--)
        {
            double temperature = levels[li].Temperature;
            int[] labels = levels[li].Partition.Labels;
            anchorByLabel.Clear();

            for (int p = 0; p < n; p++)
            {
                int lab = labels[p];
                if (lab == Assignment.Unassigned) continue;
                if (!anchorByLabel.TryGetValue(lab, out int anchor))
                {
                    anchorByLabel[lab] = p;
                    continue;
                }

                int ra = uf.Find(anchor);
                int rb = uf.Find(p);
                if (ra == rb) continue;   // already in the same subtree

                int idA = compNodeId[ra];
                int idB = compNodeId[rb];
                int size = uf.Size(ra) + uf.Size(rb);
                uf.Union(ra, rb);
                compNodeId[uf.Find(ra)] = nextId;
                merges.Add(new DendrogramNode(idA, idB, temperature, size));
                nextId++;
            }
        }

        return new Dendrogram(merges.ToArray(), n, CostAxis: "temperature");
    }
}
