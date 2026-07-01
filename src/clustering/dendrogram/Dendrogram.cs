using System;
using System.Collections.Generic;
using Graphs.Primitives;

namespace Clustering.Dendrograms;

/// <summary>
/// Immutable binary merge tree over <see cref="LeafCount"/> leaves.
/// Carries the <see cref="DendrogramNode"/> sequence plus a free-text
/// <see cref="CostAxis"/> label so downstream renderers and analyses can
/// label the y-axis correctly (e.g. <c>"mutual_reachability_distance"</c>
/// for HDBSCAN, <c>"entropy_delta"</c> for GMM agglomerative,
/// <c>"single_linkage_distance"</c> for vanilla SLINK).
/// </summary>
/// <remarks>
/// <para><b>Build invariant.</b> <see cref="Merges"/> is expected to be
/// in build order — equivalently, monotone non-decreasing in
/// <see cref="DendrogramNode.Distance"/>. The <see cref="CutAt"/>
/// helper relies on this; <see cref="CutToK"/> does not.</para>
///
/// <para><b>Internal-node ids.</b> Each merge at index <c>i</c> in
/// <see cref="Merges"/> is implicitly assigned id <c>LeafCount + i</c>,
/// matching the build-order convention HDBSCAN's condensation pass and
/// scipy/MATLAB linkage expect.</para>
/// </remarks>
public sealed record Dendrogram(
    DendrogramNode[] Merges,
    int              LeafCount,
    string           CostAxis = "distance")
{
    /// <summary>Number of internal merge nodes. Equals <c>LeafCount - 1</c>
    /// for a fully-merged tree; smaller when the tree is truncated.</summary>
    public int InternalNodeCount => Merges.Length;

    /// <summary>Total node ids in <c>[0, TotalNodeCount)</c>: leaves +
    /// internal nodes.</summary>
    public int TotalNodeCount => LeafCount + Merges.Length;

    /// <summary>
    /// Cuts to <paramref name="k"/> clusters by applying the first
    /// <c>LeafCount - k</c> merges. Returns a dense
    /// <c>int[LeafCount]</c> of cluster labels in <c>[0, k)</c>.
    /// </summary>
    /// <remarks>
    /// Requires <c>1 ≤ k ≤ LeafCount</c> and a fully-merged tree
    /// (<c>InternalNodeCount ≥ LeafCount - k</c>). Does not depend on
    /// the build-order monotonicity invariant — only on the merge
    /// sequence itself.
    /// </remarks>
    public int[] CutToK(int k)
    {
        if (k < 1 || k > LeafCount)
            throw new ArgumentOutOfRangeException(nameof(k),
                $"k ({k}) must be in [1, LeafCount={LeafCount}].");

        int mergesToApply = LeafCount - k;
        if (mergesToApply > Merges.Length)
            throw new InvalidOperationException(
                $"Tree is truncated: needs {mergesToApply} merges to reach " +
                $"k={k}, but only {Merges.Length} are available.");

        var uf = new UnionFind(LeafCount);
        for (int i = 0; i < mergesToApply; i++)
            UnionLeaves(uf, Merges[i]);

        return Densify(uf, LeafCount);
    }

    /// <summary>
    /// Cuts at cost level <paramref name="level"/>: applies all merges
    /// with <see cref="DendrogramNode.Distance"/> ≤ <paramref name="level"/>.
    /// Returns a dense <c>int[LeafCount]</c> of cluster labels.
    /// </summary>
    /// <remarks>
    /// Assumes <see cref="Merges"/> is monotone non-decreasing in
    /// <see cref="DendrogramNode.Distance"/> (the build invariant) and
    /// short-circuits the scan at the first merge exceeding the level.
    /// For consumers needing λ-style "merges with λ ≥ level" cuts,
    /// convert <c>level → 1/level</c> at the call site.
    /// </remarks>
    public int[] CutAt(double level)
    {
        var uf = new UnionFind(LeafCount);
        for (int i = 0; i < Merges.Length; i++)
        {
            if (Merges[i].Distance > level) break;
            UnionLeaves(uf, Merges[i]);
        }
        return Densify(uf, LeafCount);
    }

    /// <summary>
    /// Resolves a node id (leaf or internal merge id) down to its
    /// representative leaf for union-find purposes. For a leaf id
    /// <c>x &lt; LeafCount</c>, returns <c>x</c>; for an internal id
    /// <c>LeafCount + i</c>, walks <c>Merges[i].LeftChild</c> down to
    /// a leaf.
    /// </summary>
    private int FirstLeafBelow(int id)
    {
        while (id >= LeafCount)
            id = Merges[id - LeafCount].LeftChild;
        return id;
    }

    private void UnionLeaves(UnionFind uf, DendrogramNode m)
    {
        uf.Union(FirstLeafBelow(m.LeftChild), FirstLeafBelow(m.RightChild));
    }

    /// <summary>
    /// Densifies a <see cref="UnionFind"/>'s sparse root ids into dense
    /// labels in <c>[0, clusterCount)</c>. Mirrors the SPC partition
    /// helper but kept local to avoid the cross-namespace dependency.
    /// </summary>
    private static int[] Densify(UnionFind uf, int n)
    {
        var labels   = new int[n];
        var labelMap = new Dictionary<int, int>();
        int next = 0;
        for (int i = 0; i < n; i++)
        {
            int root = uf.Find(i);
            if (!labelMap.TryGetValue(root, out int dense))
            {
                dense = next++;
                labelMap[root] = dense;
            }
            labels[i] = dense;
        }
        return labels;
    }
}
