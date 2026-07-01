using System;
using Clustering.Graphical.SPC.Partitions.Thermal;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Primitives;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// R1: the live thermodynamic-EOM composition over synthetic rich-sweep
/// frames — pooled eq-4 columns → thermal tree → pooled landscape → walk →
/// Assignment. Oracles hand-computed (the same path fixture as
/// ThermalDendrogramTests, expressed through real Accumulator frames).
/// </summary>
public sealed class ThermalEomTests
{
    private static CsrGraph PathGraph() => CsrGraph.FromEdges(
        new[] { new Edge(0, 1, 2.0), new Edge(1, 2, 4.0) }, 3);

    private static readonly double[] Temps = { 0.5, 1.0, 2.0, 3.0 };

    // Per-edge co-membership counts over draws=20, q=20: rates → eq-4 δ̄.
    //   edge (0,1): 19,19,16,12 → δ̄ = .9525,.9525,.81,.62 (≥ .5 through T=3) → T_e = 3
    //   edge (1,2): 19,12, 4, 2 → δ̄ = .9525,.62, .24,.145 (≥ .5 through T=1) → T_e = 1
    private static Accumulator MakeFrame(CsrGraph graph, double temperature, int c01, int c12)
    {
        var counts = new int[graph.Targets.Length];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
            counts[edge.Slot] = (edge.Source, edge.Target) is (0, 1) or (1, 0) ? c01 : c12;

        return new Accumulator
        {
            Temperature = temperature,
            Q = 20,
            DrawCount = 20,
            Spins = new int[3],
            ClusterSizeHistogram = new int[3],
            RngState0 = 1, RngState1 = 2, RngState2 = 3, RngState3 = 4,
            RunningSumSqClusterSizes = 0.0,
            RunningSumSqClusterSizesExcl = 0.0,
            RunningSumEnergy = 0.0,
            RunningSumEnergySq = 0.0,
            RunningSumMag = 0.0,
            RunningSumMagSq = 0.0,
            CoMembershipCount = counts,
            SumClusterSizePerNode = new[] { 20.0, 20.0, 20.0 }, // L ≡ 1 after the draw divide
        };
    }

    private static Accumulator[] Frames(CsrGraph graph) => new[]
    {
        MakeFrame(graph, 0.5, 19, 19),
        MakeFrame(graph, 1.0, 19, 12),
        MakeFrame(graph, 2.0, 16,  4),
        MakeFrame(graph, 3.0, 12,  2),
    };

    [Fact]
    public void CoMembershipDelta_PoolsCountsAndAppliesEqFour()
    {
        var graph = PathGraph();
        // Two replicas at one T: pooled rate (10+0)/(20+20) = 0.25 → δ̄ = (19·0.25+1)/20.
        var frames = new[]
        {
            MakeFrame(graph, 1.0, 10, 10),
            MakeFrame(graph, 1.0, 0, 0),
        };

        var (temps, delta) = SweepEdgeCurves.CoMembershipDelta(frames);

        Assert.Equal(new[] { 1.0 }, temps);
        double expected = (19.0 * 0.25 + 1.0) / 20.0;
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
            Assert.Equal(expected, delta[0][edge.Slot], precision: 12);
    }

    [Fact]
    public void CoMembershipDelta_MissingCounts_ThrowsWithAccumulationHint()
    {
        var graph = PathGraph();
        Accumulator frame = MakeFrame(graph, 1.0, 10, 10) with { CoMembershipCount = null };

        var ex = Assert.Throws<InvalidOperationException>(
            () => SweepEdgeCurves.CoMembershipDelta(new[] { frame }));
        Assert.Contains("comembership", ex.Message, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void Resolve_EndToEnd_ThermalTreeWalkAndHonestAbstain()
    {
        // Same physics as ThermalDendrogramTests, through real frames:
        //   tree heights 3 (pair {0,1}) and 1 (tail joins); L ≡ 1 →
        //   masses 4 (pair) and 4.5 (root); EOM selects the pair; node 2
        //   abstains.
        var graph = PathGraph();

        ThermalEomResult result = ThermalEom.Resolve(graph, Frames(graph), theta: 0.5, graphId: "test");

        Assert.Equal("temperature", result.Dendrogram.CostAxis);
        Assert.Equal(3.0, result.Dendrogram.Merges[0].Distance, precision: 12);
        Assert.Equal(1.0, result.Dendrogram.Merges[1].Distance, precision: 12);

        Assert.Equal("MeanClusterSize", result.Landscape.Provenance!.Sink);
        Assert.Equal(4.0, result.Walk.Mass[0], precision: 12);
        Assert.Equal(4.5, result.Walk.Mass[1], precision: 12);

        Assert.Equal(new[] { true, false }, result.Selected);
        Assert.Equal(1, result.Assignment.Count);
        Assert.Equal(result.Assignment.Labels[0], result.Assignment.Labels[1]);
        Assert.Equal(Assignment.Unassigned, result.Assignment.Labels[2]);
    }

    [Fact]
    public void Resolve_WithAscendCompletion_AdoptsTheAbstainedTail()
    {
        // Same fixture; completion = Ascend on the coldest landscape column
        // (flat here, so the index tie-break carries node 2 to its assigned
        // neighbor): the abstain is completed, count unchanged.
        var graph = PathGraph();

        ThermalEomResult result = ThermalEom.Resolve(
            graph, Frames(graph), theta: 0.5,
            completion: ThermalPeripheryCompletion.Ascend);

        Assert.Equal(1, result.Assignment.Count);
        Assert.Equal(result.Assignment.Labels[0], result.Assignment.Labels[2]);
        Assert.DoesNotContain(Assignment.Unassigned, result.Assignment.Labels);
    }
}
