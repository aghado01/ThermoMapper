#nullable enable
using System;
using System.Linq;
using Graphs.Primitives;
using Graphs.Spectral;
using Maths.LinAlg;
using Maths.Topology;
using Xunit;

using TDA.Ph;
namespace TDA.Mapper.Tests;

/// <summary>
/// <see cref="RipsFiltration.RipsFromGraph"/> → <see cref="PersistentHomology.Compute"/>
/// → <see cref="DiagramMetrics.Wasserstein"/> pipeline tests.
/// </summary>
public sealed class RipsFiltrationTests
{
    [Fact]
    public void SingleEdge_MatchesHandBuiltH0Barcode()
    {
        var g = CsrGraph.FromEdges(new[] { new Edge(0, 1, 1.0) }, nodeCount: 2);
        var filtration = RipsFiltration.RipsFromGraph(g, FiltrationWeights.RawDistance, maxDimension: 1);

        Barcode barcode = PersistentHomology.Compute(filtration);

        var h0 = barcode.Bars.Where(b => b.Dimension == 0).ToList();
        Assert.Equal(2, h0.Count);
        Assert.Single(h0, b => b.IsInfinite);
        Assert.Single(h0, b => !b.IsInfinite && b.Death == 1.0);
    }

    [Fact]
    public void TriangleSkeleton_NoFillers_InfiniteH1()
    {
        var g = TriangleGraph(edgeWeight: 1.0);
        var filtration = RipsFiltration.RipsFromGraph(g, FiltrationWeights.RawDistance, maxDimension: 1);

        Barcode barcode = PersistentHomology.Compute(filtration);

        var h1 = barcode.Bars.Where(b => b.Dimension == 1).ToList();
        var bar = Assert.Single(h1);
        Assert.True(bar.IsInfinite);
        Assert.Equal(1.0, bar.Birth);
    }

    [Fact]
    public void TriangleGraph_WithFillers_FiniteH1()
    {
        var g = TriangleGraph(edgeWeight: 1.0);
        var filtration = RipsFiltration.RipsFromGraph(g, FiltrationWeights.RawDistance, maxDimension: 2);

        Barcode barcode = PersistentHomology.Compute(filtration);

        var h1 = barcode.Bars.Where(b => b.Dimension == 1).ToList();
        var bar = Assert.Single(h1);
        // All edges and triangle share filtration value 1 — loop killed at same scale.
        Assert.False(bar.IsInfinite);
        Assert.Equal(1.0, bar.Birth);
        Assert.Equal(1.0, bar.Death);
    }

    [Fact]
    public void IdenticalGraphs_PipelineWassersteinIsZero()
    {
        var g = TriangleGraph(edgeWeight: 1.0);
        var f = RipsFiltration.RipsFromGraph(g, FiltrationWeights.RawDistance);

        Barcode a = PersistentHomology.Compute(f);
        Barcode b = PersistentHomology.Compute(f);

        double w1 = DiagramMetrics.Wasserstein(a, b, dimension: 1, p: 1.0);
        Assert.Equal(0.0, w1);
    }

    [Fact]
    public void DifferentEdgeWeights_NonzeroWassersteinOnH0()
    {
        var gSmall = CsrGraph.FromEdges(new[] { new Edge(0, 1, 1.0) }, nodeCount: 2);
        var gLarge = CsrGraph.FromEdges(new[] { new Edge(0, 1, 2.0) }, nodeCount: 2);

        Barcode small = PersistentHomology.Compute(
            RipsFiltration.RipsFromGraph(gSmall, FiltrationWeights.RawDistance, maxDimension: 1));
        Barcode large = PersistentHomology.Compute(
            RipsFiltration.RipsFromGraph(gLarge, FiltrationWeights.RawDistance, maxDimension: 1));

        double w1 = DiagramMetrics.Wasserstein(small, large, dimension: 0, p: 1.0);
        Assert.Equal(1.0, w1, precision: 12);
    }

    [Fact]
    public void EffectiveResistanceWeights_BuildsFiltrationAndRunsPh()
    {
        var g = TriangleGraph(edgeWeight: 1.0);
        var weights = EffectiveResistanceWeights.FromGraph(g, kMax: 8, solverKind: SolverKind.Dense);

        var filtration = RipsFiltration.RipsFromGraph(g, weights, maxDimension: 2);
        Barcode barcode = PersistentHomology.Compute(filtration, maxDimension: 1);

        Assert.NotEmpty(barcode.Bars);
        Assert.All(filtration.Simplices.Where(s => s.Dimension == 1), s => Assert.True(s.FiltrationValue > 0.0));
    }

    static CsrGraph TriangleGraph(double edgeWeight)
    {
        return CsrGraph.FromEdges(new[]
        {
            new Edge(0, 1, edgeWeight),
            new Edge(0, 2, edgeWeight),
            new Edge(1, 2, edgeWeight),
        }, nodeCount: 3);
    }
}
