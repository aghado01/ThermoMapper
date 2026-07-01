using System;
using System.Collections.Generic;
using Clustering.Primitives;

namespace Clustering.Dendrograms;

/// <summary>
/// A condensed cluster tree: the HDBSCAN condensation of a binary single-linkage
/// <see cref="Dendrogram"/> into a tree of "real" clusters (those whose both
/// sides clear <c>minClusterSize</c> at a split). Producer-agnostic — any
/// dendrogram on any cost axis condenses here, given the cost→λ persistence map.
/// </summary>
/// <remarks>
/// <para>This is the shared form of the machinery currently inlined in
/// <c>HdbscanRunner.ExtractClusters</c>; HDBSCAN's stability
/// <c>Σ_p (λ_fallout(p) − λ_birth(C))</c> is exactly the per-leaf walk mass with
/// L ≡ 1 over each point's condensed membership interval (landscape-carrier note
/// §2). Landed ALONGSIDE the runner — the rewire/rip-out is fresh-look cleanup
/// (post-unification ledger #7).</para>
/// <para>Cluster ids: 0 = root, then two per real split, in top-down order.
/// <see cref="Death"/> is 0 for a cluster that never splits again (a leaf
/// cluster persisting to λ→∞), matching the runner's membership-window
/// convention.</para>
/// </remarks>
public sealed record CondensedTree(
    int      LeafCount,
    int      ClusterCount,
    int[]    Parent,
    double[] Birth,
    double[] Death,
    double[] Stability,
    int[]    LeafCluster,
    double[] LeafFalloutLambda)
{
    /// <summary>
    /// Excess-of-mass selection over the condensed tree: bottom-up, select C iff
    /// its own stability ≥ Σ subtree-stability of its children; top-down,
    /// deselect any cluster under a selected ancestor. Root (id 0) is eligible
    /// only when <paramref name="allowSingleCluster"/>.
    /// </summary>
    public bool[] SelectByExcessOfMass(bool allowSingleCluster)
    {
        int loopStart = allowSingleCluster ? 0 : 1;
        var childSum    = new double[ClusterCount];
        var subtreeStab = new double[ClusterCount];
        var selected    = new bool[ClusterCount];

        for (int i = ClusterCount - 1; i >= loopStart; i--)
        {
            double own = Stability[i];
            if (own >= childSum[i]) { subtreeStab[i] = own; selected[i] = true; }
            else                     subtreeStab[i] = childSum[i];
            int p = Parent[i];
            if (p >= 0) childSum[p] += subtreeStab[i];
        }

        var hasSelectedAncestor = new bool[ClusterCount];
        for (int i = 1; i < ClusterCount; i++)
        {
            int p = Parent[i];
            hasSelectedAncestor[i] = hasSelectedAncestor[p] || selected[p];
            if (hasSelectedAncestor[i]) selected[i] = false;
        }
        return selected;
    }

    /// <summary>
    /// Leaf selection: select every condensed cluster that has no condensed child
    /// (a <i>leaf</i> of the condensed tree — the finest stable clusters). The
    /// "Leaf" sibling of <see cref="SelectByExcessOfMass"/> over the same tree —
    /// sklearn's <c>cluster_selection_method='leaf'</c>. Recovers structure EOM
    /// collapses (EOM prefers few large persistent clusters; leaf takes the
    /// terminal ones). Leaf clusters are mutually non-ancestral by construction,
    /// so no top-down suppression is needed. Root (id 0) is a leaf only when the
    /// tree never split, and is then selected only if
    /// <paramref name="allowSingleCluster"/>.
    /// </summary>
    public bool[] SelectByLeaf(bool allowSingleCluster)
    {
        var isParent = new bool[ClusterCount];
        for (int i = 1; i < ClusterCount; i++)   // i=0 is the root (no parent)
            if (Parent[i] >= 0) isParent[Parent[i]] = true;

        var selected = new bool[ClusterCount];
        int start = allowSingleCluster ? 0 : 1;
        for (int i = start; i < ClusterCount; i++)
            if (!isParent[i]) selected[i] = true;
        return selected;
    }

    /// <summary>
    /// Epsilon post-filter (sklearn's <c>cluster_selection_epsilon</c>) over a
    /// base selection (EOM or leaf). For each selected cluster born at a
    /// mutual-reachability distance below <paramref name="epsilon"/> (i.e. split
    /// off too finely), walk up to the nearest ancestor born <i>above</i>
    /// <paramref name="epsilon"/> and select that instead — a DBSCAN-like distance
    /// floor that merges over-fine splits (e.g. tames leaf's fragmentation). A
    /// cluster whose parent is the root resolves to the root when
    /// <paramref name="allowSingleCluster"/>, else stays put. Selections may end up
    /// ancestor/descendant-overlapping; <see cref="ResolveLabeled"/> /
    /// <see cref="ToAssignment"/> resolve each leaf to its <i>deepest</i> selected
    /// ancestor, which is the intended epsilon semantics. <paramref name="epsilon"/>
    /// ≤ 0 is a no-op (returns a copy of the base selection).
    /// </summary>
    public bool[] SelectByEpsilon(bool[] baseSelected, double epsilon, bool allowSingleCluster)
    {
        ArgumentNullException.ThrowIfNull(baseSelected);
        if (baseSelected.Length != ClusterCount)
            throw new ArgumentException($"baseSelected length ({baseSelected.Length}) != cluster count ({ClusterCount}).", nameof(baseSelected));
        if (epsilon <= 0.0) return (bool[])baseSelected.Clone();

        var children = new List<int>[ClusterCount];
        for (int c = 0; c < ClusterCount; c++) children[c] = new List<int>();
        for (int c = 1; c < ClusterCount; c++)
            if (Parent[c] >= 0) children[Parent[c]].Add(c);

        var result    = new bool[ClusterCount];
        var processed = new bool[ClusterCount];

        for (int leaf = 0; leaf < ClusterCount; leaf++)
        {
            if (!baseSelected[leaf]) continue;
            if (BirthDistance(leaf) >= epsilon) { result[leaf] = true; continue; }
            if (processed[leaf]) continue;

            int merged = TraverseUpwards(leaf, epsilon, allowSingleCluster);
            result[merged] = true;
            MarkDescendants(merged, children, processed);   // siblings under `merged` fold in
        }
        return result;
    }

    /// <summary>Mutual-reachability distance at which cluster <paramref name="c"/>
    /// was born (1/λ_birth); +∞ for the root (born at λ=0).</summary>
    private double BirthDistance(int c) => Birth[c] > 0.0 ? 1.0 / Birth[c] : double.PositiveInfinity;

    /// <summary>Climb from <paramref name="node"/> to the nearest ancestor born
    /// above <paramref name="epsilon"/>; stop at the root (returns it when
    /// <paramref name="allowSingleCluster"/>, else the node closest to it).</summary>
    private int TraverseUpwards(int node, double epsilon, bool allowSingleCluster)
    {
        int parent = Parent[node];
        if (parent < 0) return node;                                  // node is the root
        if (parent == 0) return allowSingleCluster ? 0 : node;        // parent is the root
        return BirthDistance(parent) > epsilon
            ? parent
            : TraverseUpwards(parent, epsilon, allowSingleCluster);
    }

    private static void MarkDescendants(int a, List<int>[] children, bool[] processed)
    {
        var stack = new Stack<int>(children[a]);
        while (stack.Count > 0)
        {
            int c = stack.Pop();
            processed[c] = true;
            foreach (int ch in children[c]) stack.Push(ch);
        }
    }

    /// <summary>
    /// Resolves a selection to dense labels <b>and</b> HDBSCAN membership
    /// probabilities in one leaf walk. Labels match <see cref="ToAssignment"/>
    /// (each leaf → its nearest selected ancestor; unselected → −1). The
    /// probability for leaf x assigned to selected cluster L is the condensed-tree
    /// soft score λ_fallout(x)/λ_death(L) when L is x's <i>deepest</i> condensed
    /// cluster (x may have left L early via a small-side fallout), 1.0 when L is a
    /// strict ancestor of x's deepest cluster (x stayed until L died), and 0.0 for
    /// noise. This is the λ-ratio membership the in-silo HDBSCAN extraction
    /// produced, lifted onto the shared tree so the runner reads it here rather
    /// than recomputing condensation.
    /// </summary>
    public (int[] Labels, double[] MembershipProbabilities, int Count) ResolveLabeled(bool[] selected)
    {
        ArgumentNullException.ThrowIfNull(selected);
        if (selected.Length != ClusterCount)
            throw new ArgumentException($"selected length ({selected.Length}) != cluster count ({ClusterCount}).", nameof(selected));

        var cidLabel = new int[ClusterCount];
        Array.Fill(cidLabel, -1);
        int count = 0;
        for (int i = 0; i < ClusterCount; i++)
            if (selected[i]) cidLabel[i] = count++;

        var labels = new int[LeafCount];
        var prob   = new double[LeafCount];
        for (int x = 0; x < LeafCount; x++)
        {
            int deepest = LeafCluster[x];
            int cid     = deepest;
            while (cid >= 0 && cidLabel[cid] < 0) cid = Parent[cid];

            if (cid < 0)
            {
                labels[x] = Assignment.Unassigned;
                prob[x]   = 0.0;
                continue;
            }

            labels[x] = cidLabel[cid];
            double lambdaMax = Death[cid];
            if (cid == deepest)
            {
                double lambdaIn = LeafFalloutLambda[x];
                prob[x] = lambdaMax > 0.0 && !double.IsInfinity(lambdaMax)
                    ? Math.Min(1.0, lambdaIn / lambdaMax)
                    : 1.0;
            }
            else
            {
                prob[x] = 1.0;
            }
        }
        return (labels, prob, count);
    }

    /// <summary>
    /// Densely labels each selected cluster and assigns each leaf to its nearest
    /// selected ancestor; unselected leaves are <see cref="Assignment.Unassigned"/>.
    /// </summary>
    public Assignment ToAssignment(bool[] selected)
    {
        ArgumentNullException.ThrowIfNull(selected);
        if (selected.Length != ClusterCount)
            throw new ArgumentException($"selected length ({selected.Length}) != cluster count ({ClusterCount}).", nameof(selected));

        var cidLabel = new int[ClusterCount];
        Array.Fill(cidLabel, -1);
        int count = 0;
        for (int i = 0; i < ClusterCount; i++)
            if (selected[i]) cidLabel[i] = count++;

        var labels = new int[LeafCount];
        for (int x = 0; x < LeafCount; x++)
        {
            int cid = LeafCluster[x];
            while (cid >= 0 && cidLabel[cid] < 0) cid = Parent[cid];
            labels[x] = cid < 0 ? Assignment.Unassigned : cidLabel[cid];
        }
        return new Assignment { Labels = labels, Count = count };
    }
}

/// <summary>
/// Condenses a binary single-linkage <see cref="Dendrogram"/> into a
/// <see cref="CondensedTree"/> — the producer-agnostic lift of HDBSCAN's
/// condensation pass.
/// </summary>
public static class Condensation
{
    /// <summary>λ = 1/cost — the HDBSCAN persistence map (mutual-reachability
    /// distance → density-like λ). The cost-axis interpretation a producer
    /// stamps; passed explicitly so the condensation stays cost-axis-neutral.</summary>
    public static double InverseLambda(double cost) => cost > 0.0 ? 1.0 / cost : double.PositiveInfinity;

    /// <summary>
    /// Condense <paramref name="dendrogram"/> (build-ordered binary merges) using
    /// <paramref name="costToLambda"/> for the persistence axis. Defaults to the
    /// HDBSCAN map <see cref="InverseLambda"/>.
    /// </summary>
    public static CondensedTree Condense(
        Dendrogram dendrogram,
        int minClusterSize,
        Func<double, double>? costToLambda = null)
    {
        ArgumentNullException.ThrowIfNull(dendrogram);
        if (minClusterSize < 2)
            throw new ArgumentOutOfRangeException(nameof(minClusterSize), "Must be >= 2.");
        costToLambda ??= InverseLambda;

        DendrogramNode[] tree = dendrogram.Merges;
        int n = dendrogram.LeafCount;
        int numMerges = tree.Length;

        int      maxCondensed = 2 * n;
        var      cParent      = new int[maxCondensed];
        var      cBirth       = new double[maxCondensed];
        var      cDeath       = new double[maxCondensed];
        var      cStab        = new double[maxCondensed];
        cParent[0] = -1;
        int numCondensed = 1; // root = id 0

        var clusterAtMerge = new int[numMerges];
        if (numMerges > 0)
        {
            Array.Fill(clusterAtMerge, -1);
            clusterAtMerge[numMerges - 1] = 0;
        }

        var clusterAtLeaf     = new int[n];
        var leafFalloutLambda = new double[n];
        var dfsStack          = new int[n];

        for (int i = numMerges - 1; i >= 0; i--)
        {
            int parentCid = clusterAtMerge[i];
            int leftId    = tree[i].LeftChild;
            int rightId   = tree[i].RightChild;

            if (parentCid < 0)
            {
                if (leftId  >= n) clusterAtMerge[leftId  - n] = -1;
                if (rightId >= n) clusterAtMerge[rightId - n] = -1;
                continue;
            }

            int leftSize  = leftId  < n ? 1 : tree[leftId  - n].Size;
            int rightSize = rightId < n ? 1 : tree[rightId - n].Size;

            double cost        = tree[i].Distance;
            double splitLambda = cost > 0.0 ? costToLambda(cost) : double.PositiveInfinity;
            double pBirth      = cBirth[parentCid];

            bool leftBig  = leftSize  >= minClusterSize;
            bool rightBig = rightSize >= minClusterSize;

            if (leftBig && rightBig)
            {
                cStab[parentCid] += (leftSize + rightSize) * (splitLambda - pBirth);
                cDeath[parentCid] = splitLambda;

                int newLeftCid  = numCondensed++;
                int newRightCid = numCondensed++;
                cParent[newLeftCid]  = parentCid; cBirth[newLeftCid]  = splitLambda;
                cParent[newRightCid] = parentCid; cBirth[newRightCid] = splitLambda;

                DendrogramBuilder.VisitLeaves(tree, leftId,  n, dfsStack, leaf => clusterAtLeaf[leaf] = newLeftCid);
                DendrogramBuilder.VisitLeaves(tree, rightId, n, dfsStack, leaf => clusterAtLeaf[leaf] = newRightCid);

                if (leftId  >= n) clusterAtMerge[leftId  - n] = newLeftCid;
                if (rightId >= n) clusterAtMerge[rightId - n] = newRightCid;
            }
            else if (leftBig)
            {
                cStab[parentCid] += rightSize * (splitLambda - pBirth);
                DendrogramBuilder.VisitLeaves(tree, rightId, n, dfsStack, leaf => leafFalloutLambda[leaf] = splitLambda);
                if (leftId  >= n) clusterAtMerge[leftId  - n] = parentCid;
                if (rightId >= n) clusterAtMerge[rightId - n] = -1;
            }
            else if (rightBig)
            {
                cStab[parentCid] += leftSize * (splitLambda - pBirth);
                DendrogramBuilder.VisitLeaves(tree, leftId, n, dfsStack, leaf => leafFalloutLambda[leaf] = splitLambda);
                if (rightId >= n) clusterAtMerge[rightId - n] = parentCid;
                if (leftId  >= n) clusterAtMerge[leftId  - n] = -1;
            }
            else
            {
                cStab[parentCid] += (leftSize + rightSize) * (splitLambda - pBirth);
                cDeath[parentCid] = splitLambda;
                DendrogramBuilder.VisitLeaves(tree, leftId,  n, dfsStack, leaf => leafFalloutLambda[leaf] = splitLambda);
                DendrogramBuilder.VisitLeaves(tree, rightId, n, dfsStack, leaf => leafFalloutLambda[leaf] = splitLambda);
                if (leftId  >= n) clusterAtMerge[leftId  - n] = -1;
                if (rightId >= n) clusterAtMerge[rightId - n] = -1;
            }
        }

        return new CondensedTree(
            LeafCount:         n,
            ClusterCount:      numCondensed,
            Parent:            cParent[..numCondensed],
            Birth:             cBirth[..numCondensed],
            Death:             cDeath[..numCondensed],
            Stability:         cStab[..numCondensed],
            LeafCluster:       clusterAtLeaf,
            LeafFalloutLambda: leafFalloutLambda);
    }
}
