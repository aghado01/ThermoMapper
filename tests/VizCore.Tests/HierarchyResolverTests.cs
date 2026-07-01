using System;
using System.Linq;
using Clustering.Dendrograms;
using Clustering.Graphical.SPC.Partitions.Hierarchical;
using Clustering.Primitives;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// Track 1 — the Blatt/Domany hierarchical T-stack resolver: the dense per-T
/// cut substrate (<see cref="DenseTStack"/>), the nested-degenerate dendrogram
/// bridge (<see cref="PartitionHierarchyDendrogram"/>), and the EOM composition
/// (<see cref="HierarchyEom"/>). Oracles are hand-computed (validation
/// independence — the oracle is the math, never the resolver's own output).
/// </summary>
public sealed class HierarchyResolverTests
{
    // Path graph 0—1—2: two undirected edges, slots filled per (source,target).
    private static CsrGraph PathGraph() => CsrGraph.FromEdges(
        new[] { new Edge(0, 1, 2.0), new Edge(1, 2, 4.0) }, 3);

    // A per-edge δ̄ column for the path graph at one T: value01 on (0,1),
    // value12 on (1,2). Indexed by CSR slot, as AffinityThreshold.Connect reads.
    private static double[] Column(CsrGraph graph, double value01, double value12)
    {
        var col = new double[graph.Targets.Length];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
            col[edge.Slot] = (edge.Source, edge.Target) is (0, 1) or (1, 0) ? value01 : value12;
        return col;
    }

    /// <summary>
    /// A clean nested stack: cold one-cluster → split tail → all singletons.
    ///   T=1: both bonds hot  → {0,1,2}
    ///   T=2: (0,1) hot only  → {0,1} {2}
    ///   T=3: both cold       → {0} {1} {2}
    /// </summary>
    private static PartitionHierarchy NestedStack(CsrGraph graph)
    {
        var temps   = new[] { 1.0, 2.0, 3.0 };
        var columns = new[]
        {
            Column(graph, 0.9, 0.9),  // T=1
            Column(graph, 0.9, 0.1),  // T=2
            Column(graph, 0.1, 0.1),  // T=3
        };
        return DenseTStack.Build(graph, temps, columns, theta: 0.5);
    }

    [Fact]
    public void DenseTStack_CutsEveryTemperature_AndNests()
    {
        var graph = PathGraph();
        PartitionHierarchy stack = NestedStack(graph);

        Assert.Equal(3, stack.Count);
        Assert.True(stack.NestingHolds);

        // Cluster counts cool → hot: 1, 2, 3.
        Assert.Equal(1, stack.Levels[0].Partition.Count);
        Assert.Equal(2, stack.Levels[1].Partition.Count);
        Assert.Equal(3, stack.Levels[2].Partition.Count);

        // Cold level: all three together.
        int[] cold = stack.Levels[0].Partition.Labels;
        Assert.Equal(cold[0], cold[1]);
        Assert.Equal(cold[1], cold[2]);

        // Mid level: {0,1} together, 2 apart.
        int[] mid = stack.Levels[1].Partition.Labels;
        Assert.Equal(mid[0], mid[1]);
        Assert.NotEqual(mid[0], mid[2]);
    }

    [Fact]
    public void Bridge_NestedStack_ProducesSpanningThermalTree()
    {
        var graph = PathGraph();
        PartitionHierarchy stack = NestedStack(graph);

        Dendrogram tree = PartitionHierarchyDendrogram.ToDendrogram(stack);

        Assert.Equal("temperature", tree.CostAxis);
        Assert.Equal(3, tree.LeafCount);
        Assert.Equal(2, tree.InternalNodeCount);   // spanning (cold = one cluster)

        // Hot→cold build order ⇒ DESCENDING heights: {0,1} couples at T=2, the
        // tail joins at T=1.
        Assert.Equal(2.0, tree.Merges[0].Distance, precision: 12);
        Assert.Equal(1.0, tree.Merges[1].Distance, precision: 12);

        // CutToK reproduces the intermediate stages a single-T cut cannot show
        // jointly: k=2 = the {0,1}|{2} stage, k=1 = the merged cold root.
        int[] k2 = tree.CutToK(2);
        Assert.Equal(k2[0], k2[1]);
        Assert.NotEqual(k2[0], k2[2]);

        int[] k1 = tree.CutToK(1);
        Assert.Equal(k1[0], k1[1]);
        Assert.Equal(k1[1], k1[2]);
    }

    [Fact]
    public void Bridge_NonNestedStack_Throws_DontWarpGuard()
    {
        var graph = PathGraph();
        // Cold {0,1}|{2}; hot {0}|{1,2}: the hot cluster {1,2} straddles two cold
        // clusters — non-nested contest a single-linkage tree cannot represent.
        var temps   = new[] { 1.0, 2.0 };
        var columns = new[]
        {
            Column(graph, 0.9, 0.1),  // T=1: {0,1} {2}
            Column(graph, 0.1, 0.9),  // T=2: {0} {1,2}
        };
        PartitionHierarchy stack = DenseTStack.Build(graph, temps, columns, theta: 0.5);

        Assert.False(stack.NestingHolds);
        var ex = Assert.Throws<InvalidOperationException>(
            () => PartitionHierarchyDendrogram.ToDendrogram(stack));
        Assert.Contains("not strictly nested", ex.Message, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void HierarchyEom_NestedStack_SelectsPersistentClusterAndAbstainsTail()
    {
        var graph = PathGraph();
        PartitionHierarchy stack = NestedStack(graph);

        // L ≡ 1 landscape on the same axis/grid: masses are |C| × lifetime.
        //   {0,1} born T=2, dies T=1 → alive (1,2], mass 2·w(T=2)=2.
        //   root {0,1,2} born T=1 → ineligible (root); EOM selects {0,1}, 2 abstains.
        var grid = new[] { 1.0, 2.0, 3.0 };
        var columns = new[]
        {
            new[] { 1.0, 1.0, 1.0 },
            new[] { 1.0, 1.0, 1.0 },
            new[] { 1.0, 1.0, 1.0 },
        };
        var landscape = Landscape.Create("temperature", grid, columns,
            new LandscapeProvenance("UnitLandscape", "test"));

        HierarchyEomResult result = HierarchyEom.Resolve(graph, stack, landscape);

        Assert.True(result.RawNestingHeld);
        Assert.False(result.Restored);
        Assert.NotNull(result.Dendrogram);
        Assert.Equal(1, result.Assignment.Count);
        Assert.Equal(result.Assignment.Labels[0], result.Assignment.Labels[1]);
        Assert.Equal(Assignment.Unassigned, result.Assignment.Labels[2]);
    }

    [Fact]
    public void HierarchyEom_NonNestedStack_RestoresPremiseThenResolves()
    {
        var graph = PathGraph();
        var temps   = new[] { 1.0, 2.0 };
        var columns = new[]
        {
            Column(graph, 0.9, 0.1),  // {0,1} {2}
            Column(graph, 0.1, 0.9),  // {0} {1,2}  — straddles, non-nested
        };
        PartitionHierarchy raw = DenseTStack.Build(graph, temps, columns, theta: 0.5);
        Assert.False(raw.NestingHolds);

        var grid = new[] { 1.0, 2.0 };
        var land = new[] { new[] { 1.0, 1.0, 1.0 }, new[] { 1.0, 1.0, 1.0 } };
        var landscape = Landscape.Create("temperature", grid, land,
            new LandscapeProvenance("UnitLandscape", "test"));

        HierarchyEomResult result = HierarchyEom.Resolve(graph, raw, landscape, restoreNesting: true);

        // Premise restored ⇒ the stack now nests and the bridge succeeds.
        Assert.False(result.RawNestingHeld);
        Assert.True(result.Restored);
        Assert.True(result.Stack.NestingHolds);
        Assert.NotNull(result.Dendrogram);
    }

    [Fact]
    public void HierarchyEom_NonNestedStack_StrictGate_FlagsAndAbstains()
    {
        var graph = PathGraph();
        var temps   = new[] { 1.0, 2.0 };
        var columns = new[]
        {
            Column(graph, 0.9, 0.1),
            Column(graph, 0.1, 0.9),
        };
        PartitionHierarchy raw = DenseTStack.Build(graph, temps, columns, theta: 0.5);

        var grid = new[] { 1.0, 2.0 };
        var land = new[] { new[] { 1.0, 1.0, 1.0 }, new[] { 1.0, 1.0, 1.0 } };
        var landscape = Landscape.Create("temperature", grid, land);

        HierarchyEomResult result = HierarchyEom.Resolve(graph, raw, landscape, restoreNesting: false);

        // Strict gate, premise not restored: no dendrogram, honest all-abstain.
        Assert.False(result.RawNestingHeld);
        Assert.False(result.Restored);
        Assert.Null(result.Dendrogram);
        Assert.All(result.Assignment.Labels, l => Assert.Equal(Assignment.Unassigned, l));
    }
}
