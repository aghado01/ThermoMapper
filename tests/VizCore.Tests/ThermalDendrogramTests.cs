using System;
using Clustering.Dendrograms;
using Clustering.Primitives;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// R1: the thermal dendrogram producer (decoupling-temperature heights) and
/// the thermodynamic-EOM pipeline end to end on hand-computed curves. The
/// orientation contract (descending heights; windows (death, birth]; root
/// open at the cold end) is pinned here.
/// </summary>
public sealed class ThermalDendrogramTests
{
    // Path 0 —— 1 —— 2; the (0,1) bond survives to T=3, the (1,2) bond to T=1.
    private static CsrGraph PathGraph() => CsrGraph.FromEdges(
        new[] { new Edge(0, 1, 2.0), new Edge(1, 2, 4.0) }, 3);

    private static readonly double[] Temps = { 0.5, 1.0, 2.0, 3.0 };

    private static double[] Column(CsrGraph graph, double g01, double g12)
    {
        var column = new double[graph.Targets.Length];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
            column[edge.Slot] = (edge.Source, edge.Target) is (0, 1) or (1, 0) ? g01 : g12;
        return column;
    }

    private static double[][] Curves(CsrGraph graph) => new[]
    {
        Column(graph, 0.9, 0.9), // T = 0.5
        Column(graph, 0.9, 0.6), // T = 1.0
        Column(graph, 0.8, 0.2), // T = 2.0
        Column(graph, 0.6, 0.1), // T = 3.0
    };

    [Fact]
    public void FromEdgeCurves_HeightsAreDecouplingTemperatures()
    {
        var graph = PathGraph();

        Dendrogram tree = ThermalDendrogram.FromEdgeCurves(graph, Temps, Curves(graph), theta: 0.5);

        Assert.Equal("temperature", tree.CostAxis);
        Assert.Equal(2, tree.InternalNodeCount);
        Assert.Equal(3.0, tree.Merges[0].Distance, precision: 12); // (0,1) couples at T=3
        Assert.Equal(1.0, tree.Merges[1].Distance, precision: 12); // tail joins at T=1
        Assert.Equal(2, tree.Merges[0].Size);
        Assert.Equal(3, tree.Merges[1].Size);
    }

    [Fact]
    public void FromEdgeCurves_NeverCoupledEdge_YieldsForestWithThermalOutlier()
    {
        // The (1,2) bond never reaches θ within the grid: node 2 is a thermal
        // outlier in the observed window — the result is a forest, and the
        // outlier belongs to no merge; with nothing selected it resolves
        // Unassigned (the honest abstain, not an error).
        var graph = PathGraph();
        var curves = new[]
        {
            Column(graph, 0.9, 0.2),
            Column(graph, 0.9, 0.2),
            Column(graph, 0.8, 0.1),
            Column(graph, 0.6, 0.1),
        };

        Dendrogram forest = ThermalDendrogram.FromEdgeCurves(graph, Temps, curves, theta: 0.5);

        Assert.Equal(1, forest.InternalNodeCount);
        Assert.Equal(3.0, forest.Merges[0].Distance, precision: 12);
        Assert.Equal(2, forest.Merges[0].Size);

        Assignment none = LandscapeWalk.ToAssignment(forest, new[] { false });
        Assert.All(none.Labels, label => Assert.Equal(Assignment.Unassigned, label));
    }

    [Fact]
    public void ThermodynamicEom_EndToEnd_SelectsTheHotPairAndAbstainsTheTail()
    {
        // The first instrument: thermal dendrogram (structure) walked with a
        // thermal landscape (height). L ≡ 1 on grid [0.5,1,2,3]
        // (widths [0.5,1,1,1]); descending windows:
        //   id3 (birth 3, death 1): (1, 3]  → cells at 2, 3   → M = 2·(1+1) = 4
        //   root (birth 1, cold-open): (−∞, 1] → cells at 0.5, 1 → M = 3·(0.5+1) = 4.5
        // EOM (root ineligible): the hot pair {0,1} is selected; node 2 only
        // ever belongs to the root → Unassigned, the honest abstain.
        var graph = PathGraph();
        Dendrogram tree = ThermalDendrogram.FromEdgeCurves(graph, Temps, Curves(graph), theta: 0.5);
        Landscape landscape = Landscape.Create(
            "temperature",
            Temps,
            new[] { new[] { 1.0, 1.0, 1.0 }, new[] { 1.0, 1.0, 1.0 }, new[] { 1.0, 1.0, 1.0 }, new[] { 1.0, 1.0, 1.0 } });

        ClusterWalkReport report = LandscapeWalk.ClusterProfiles(tree, landscape);

        Assert.Equal(4.0, report.Mass[0], precision: 12);
        Assert.Equal(4.5, report.Mass[1], precision: 12);
        Assert.Equal(3.0, report.Birth[0], precision: 12);
        Assert.Equal(1.0, report.Death[0], precision: 12);
        Assert.True(double.IsNegativeInfinity(report.Death[1]));

        bool[] selected = LandscapeWalk.SelectByExcessOfMass(tree, report.Mass);
        Assert.Equal(new[] { true, false }, selected);

        Assignment assignment = LandscapeWalk.ToAssignment(tree, selected);
        Assert.Equal(1, assignment.Count);
        Assert.Equal(assignment.Labels[0], assignment.Labels[1]);
        Assert.Equal(Assignment.Unassigned, assignment.Labels[2]);
    }
}
