using System;
using Clustering.Dendrograms;
using Clustering.Statistical.GMM;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// Unification cross-check: GMM's agglomerative merge sequence, emitted as a
/// shared <see cref="Dendrogram"/> (BuildDendrogram), reproduces
/// EntropyMergeStrategy.Merge's component→cluster partition at every cut level.
/// Validation independence: the oracle is the bespoke flatten path (Baudry
/// entropy merging, a separate code path), not the dendrogram's own assumptions.
/// Landed alongside the flatten path; the rewire is fresh-look cleanup.
/// </summary>
public sealed class GmmDendrogramTests
{
    // EntropyMergeStrategy reads only components.Length; minimal dummies suffice.
    private static GaussianComponent[] Dummies(int k)
    {
        var arr = new GaussianComponent[k];
        for (int i = 0; i < k; i++) arr[i] = new GaussianComponent(1);
        return arr;
    }

    // A 6-point × 5-component responsibility matrix with structure: components
    // {0,1} share points, {2,3} share points, 4 is distinct — so the greedy
    // entropy merges are non-trivial (not a degenerate chain).
    private static double[,] Responsibilities() => new double[,]
    {
        { 0.55, 0.40, 0.02, 0.02, 0.01 },
        { 0.45, 0.50, 0.02, 0.02, 0.01 },
        { 0.02, 0.03, 0.50, 0.43, 0.02 },
        { 0.02, 0.02, 0.46, 0.48, 0.02 },
        { 0.10, 0.05, 0.10, 0.05, 0.70 },
        { 0.30, 0.30, 0.15, 0.15, 0.10 },
    };

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void BuildDendrogram_CutToK_ReproducesMergePartition(int k)
    {
        var components = Dummies(5);
        var resp = Responsibilities();

        int[] oracle = new EntropyMergeStrategy(targetClusters: k).Merge(components, resp);

        Dendrogram tree = EntropyMergeStrategy.BuildDendrogram(components, resp);
        int[] viaTree = tree.CutToK(k);

        Assert.Equal("entropy_reduction", tree.CostAxis);
        Assert.Equal(5, tree.LeafCount);
        AssertSamePartition(oracle, viaTree);
    }

    [Fact]
    public void BuildDendrogram_HeightsAreMonotoneNonDecreasing()
    {
        // Cumulative entropy reduction only grows as components merge.
        Dendrogram tree = EntropyMergeStrategy.BuildDendrogram(Dummies(5), Responsibilities());
        for (int i = 1; i < tree.Merges.Length; i++)
            Assert.True(tree.Merges[i].Distance >= tree.Merges[i - 1].Distance - 1e-12,
                $"merge {i} height {tree.Merges[i].Distance} < {tree.Merges[i - 1].Distance}");
        Assert.Equal(5, tree.Merges[^1].Size); // root spans all components
    }

    private static void AssertSamePartition(int[] expected, int[] actual)
    {
        Assert.Equal(expected.Length, actual.Length);
        var fwd = new System.Collections.Generic.Dictionary<int, int>();
        var rev = new System.Collections.Generic.Dictionary<int, int>();
        for (int i = 0; i < expected.Length; i++)
        {
            if (fwd.TryGetValue(expected[i], out int m)) Assert.Equal(m, actual[i]);
            else fwd[expected[i]] = actual[i];
            if (rev.TryGetValue(actual[i], out int b)) Assert.Equal(b, expected[i]);
            else rev[actual[i]] = expected[i];
        }
    }
}
