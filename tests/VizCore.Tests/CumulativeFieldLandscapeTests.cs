using System;
using Clustering.Graphical.SPC.Runtime.Core.Solver;
using Clustering.Primitives;
using Graphs;
using Graphs.Observables;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// R1: the solver-side landscape mint (PKWang closed form → per-node
/// marginals → Landscape) and the producer-agnostic affinity node marginals.
/// Marginal oracles are hand-computed; the OverGrid test pins wiring,
/// provenance, and the closed form's monotonicity in T.
/// </summary>
public sealed class PKWangLandscapeTests
{
    // Path 0 —(J=2)— 1 —(J=4)— 2.
    private static CsrGraph PathGraph() => CsrGraph.FromEdges(
        new[] { new Edge(0, 1, 2.0), new Edge(1, 2, 4.0) }, 3);

    private static Affinities HandAffinities(CsrGraph graph)
    {
        // G(0,1) = 0.25, G(1,2) = 0.5, assigned via the canonical slot walk.
        var g = new double[graph.Targets.Length];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
            g[edge.Slot] = (edge.Source, edge.Target) switch
            {
                (0, 1) or (1, 0) => 0.25,
                (1, 2) or (2, 1) => 0.5,
                _ => throw new InvalidOperationException("unexpected edge"),
            };
        return new Affinities { Temperature = 1.0, G = g, ReplicaIndex = 0 };
    }

    [Fact]
    public void BondMass_SumsAffinityOverIncidentEdges()
    {
        var graph = PathGraph();
        double[] marginal = AffinityNodeMarginals.BondMass(graph, HandAffinities(graph));

        Assert.Equal(0.25, marginal[0], precision: 12);
        Assert.Equal(0.75, marginal[1], precision: 12);
        Assert.Equal(0.50, marginal[2], precision: 12);
    }

    [Fact]
    public void LocalEnergy_SumsCouplingWeightedFrustration()
    {
        // J·(1−G): edge (0,1): 2·0.75 = 1.5; edge (1,2): 4·0.5 = 2.0.
        var graph = PathGraph();
        double[] marginal = AffinityNodeMarginals.LocalEnergy(graph, HandAffinities(graph));

        Assert.Equal(1.5, marginal[0], precision: 12);
        Assert.Equal(3.5, marginal[1], precision: 12);
        Assert.Equal(2.0, marginal[2], precision: 12);
    }

    [Fact]
    public void OverGrid_MintsThermalLandscapeFromTheClosedForm()
    {
        var graph = PathGraph();
        PKWangContext context = PKWang.Prepare(graph, EdgeWeightKind.Coupling, Field.Mean);
        var grid = new[] { 0.5, 1.0, 2.0 };

        Landscape landscape = PKWangLandscapes.OverGrid(context, grid, PKWangLandscapeSink.BondMass, graphId: "test");

        Assert.Equal("temperature", landscape.Axis);
        Assert.Equal(grid, landscape.Grid);
        Assert.Equal(3, landscape.NodeCount);
        Assert.Equal("BondMass", landscape.Provenance!.Sink);
        Assert.Equal("pkwang:closed-form", landscape.Provenance.GaugeNote);

        // G = 1 − exp(−Hcum/T) falls strictly as T rises wherever Hcum > 0,
        // so every node's bond mass is strictly decreasing along the grid.
        for (int node = 0; node < 3; node++)
        {
            Assert.True(landscape.ValuesByGridPoint[0][node] > landscape.ValuesByGridPoint[1][node]);
            Assert.True(landscape.ValuesByGridPoint[1][node] > landscape.ValuesByGridPoint[2][node]);
        }

        // Column consistency: the mint's column IS the marginal of Solve at that T.
        double[] direct = AffinityNodeMarginals.BondMass(graph, PKWang.Solve(context, 1.0));
        Assert.Equal(direct, landscape.ValuesByGridPoint[1]);
    }

    [Fact]
    public void OverGrid_LocalEnergyRisesWithTemperature()
    {
        var graph = PathGraph();
        PKWangContext context = PKWang.Prepare(graph, EdgeWeightKind.Coupling, Field.Mean);

        Landscape landscape = PKWangLandscapes.OverGrid(
            context, new[] { 0.5, 2.0 }, PKWangLandscapeSink.LocalEnergy);

        for (int node = 0; node < 3; node++)
            Assert.True(landscape.ValuesByGridPoint[1][node] > landscape.ValuesByGridPoint[0][node]);
    }
}
