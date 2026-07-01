using System;
using Graphs.Primitives;
using Graphs.Primitives.Mst;

namespace Clustering.Dendrograms;

/// <summary>
/// Shared helpers for building and traversing single-linkage-style
/// dendrograms.
/// </summary>
public static class DendrogramBuilder
{
    /// <summary>
    /// Builds a single-linkage dendrogram from a sorted MST edge list.
    /// The caller is responsible for ensuring <paramref name="sortedEdges"/>
    /// is ordered by non-decreasing <see cref="MstEdge.Weight"/>.
    /// </summary>
    /// <param name="sortedEdges">The MST edges sorted by weight.</param>
    /// <param name="leafCount">Number of leaves in the tree.</param>
    /// <param name="uf">Caller-owned UnionFind scratch whose capacity is at
    /// least <c>2 * leafCount - 1</c>.</param>
    public static DendrogramNode[] BuildSingleLinkageDendrogram(
        ReadOnlySpan<MstEdge> sortedEdges,
        int                   leafCount,
        UnionFind             uf)
    {
        if (leafCount < 2)
            throw new ArgumentOutOfRangeException(nameof(leafCount), "leafCount must be >= 2.");
        if (uf is null)
            throw new ArgumentNullException(nameof(uf));
        if (sortedEdges.Length != leafCount - 1)
            throw new ArgumentException("sortedEdges length must be leafCount - 1.", nameof(sortedEdges));

        uf.Reset();
        var tree = new DendrogramNode[leafCount - 1];
        int nextId = leafCount;

        for (int i = 0; i < sortedEdges.Length; i++)
        {
            var e = sortedEdges[i];
            int ra = uf.Find(e.U);
            int rb = uf.Find(e.V);

            int sizeA = uf.Size(ra);
            int sizeB = uf.Size(rb);

            uf.Union(ra, rb);
            uf.Reroot(ra, nextId, sizeA + sizeB);

            tree[i] = new DendrogramNode(ra, rb, e.Weight, sizeA + sizeB);
            nextId++;
        }

        return tree;
    }

    /// <summary>
    /// Visits every leaf in the subtree rooted at <paramref name="subtreeRootId"/>.
    /// </summary>
    /// <param name="tree">Dendrogram merge nodes.</param>
    /// <param name="subtreeRootId">Leaf or internal node id.</param>
    /// <param name="leafCount">Number of leaves in the tree.</param>
    /// <param name="stack">Caller-owned scratch stack sized to at least leafCount.
    /// Used to avoid recursion.</param>
    /// <param name="onLeaf">Action invoked for each visited leaf id.</param>
    public static void VisitLeaves(
        DendrogramNode[] tree,
        int subtreeRootId,
        int leafCount,
        Span<int> stack,
        Action<int> onLeaf)
    {
        if (tree is null)
            throw new ArgumentNullException(nameof(tree));
        if (onLeaf is null)
            throw new ArgumentNullException(nameof(onLeaf));
        if (subtreeRootId < 0 || subtreeRootId >= leafCount + tree.Length)
            throw new ArgumentOutOfRangeException(nameof(subtreeRootId));
        if (leafCount < 1)
            throw new ArgumentOutOfRangeException(nameof(leafCount), "leafCount must be positive.");

        if (subtreeRootId < leafCount)
        {
            onLeaf(subtreeRootId);
            return;
        }

        int top = 0;
        stack[top++] = subtreeRootId - leafCount;

        while (top > 0)
        {
            int cur = stack[--top];
            int l = tree[cur].LeftChild;
            int r = tree[cur].RightChild;

            if (l < leafCount)
                onLeaf(l);
            else
                stack[top++] = l - leafCount;

            if (r < leafCount)
                onLeaf(r);
            else
                stack[top++] = r - leafCount;
        }
    }
}
