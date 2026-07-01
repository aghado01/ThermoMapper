#nullable enable
using System.Linq;
using Xunit;

using Maths.Topology;
namespace TDA.Ph.Tests;

public sealed class BarCycleEdgesTests
{
    [Fact]
    public void Circle_InfiniteH1_ExtractsThreeEdges()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(0.0, 2),
            new Simplex(1.0, 0, 1),
            new Simplex(1.0, 0, 2),
            new Simplex(1.0, 1, 2),
        }, "t");

        Barcode barcode = PersistentInvolutedHomology.Compute(filtration, representatives: true);
        var h1 = barcode.Bars.Single(b => b.Dimension == 1);

        var edges = BarCycleEdges.GetEdgePairs(h1, filtration);
        Assert.Equal(3, edges.Count);
        Assert.Equal(
            new[] { new UndirectedEdge(0, 1), new UndirectedEdge(0, 2), new UndirectedEdge(1, 2) },
            edges.OrderBy(e => e.Lo).ThenBy(e => e.Hi).ToArray());
        Assert.All(edges, e => Assert.True(e.Lo < e.Hi));
    }

    [Fact]
    public void FilledTriangle_FiniteH1_ExtractsThreeEdgesNoTriangle()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(0.0, 2),
            new Simplex(1.0, 0, 1),
            new Simplex(1.0, 0, 2),
            new Simplex(1.0, 1, 2),
            new Simplex(2.0, 0, 1, 2),
        }, "t");

        Barcode barcode = PersistentInvolutedHomology.Compute(filtration, representatives: true);
        var h1 = barcode.Bars.Single(b => b.Dimension == 1);

        var edges = BarCycleEdges.GetEdgePairs(h1, filtration);
        Assert.Equal(3, edges.Count);
        Assert.Equal(
            new[] { new UndirectedEdge(0, 1), new UndirectedEdge(0, 2), new UndirectedEdge(1, 2) },
            edges.OrderBy(e => e.Lo).ThenBy(e => e.Hi).ToArray());
    }

    [Fact]
    public void TwoLoopWedge_TwoH1Loops()
    {
        var filtration = SimplicialFiltrationFixtures.TwoLoopWedge();
        Barcode barcode = PersistentInvolutedHomology.Compute(filtration, representatives: true);

        var loops = BarCycleEdges.H1Loops(barcode, filtration).ToList();
        Assert.Equal(2, loops.Count);
        Assert.All(loops, l => Assert.NotEmpty(l.Edges));
        Assert.Contains(loops, l => l.Edges.Contains(new UndirectedEdge(0, 2)));
        Assert.Contains(loops, l => l.Edges.Contains(new UndirectedEdge(0, 4)));
        Assert.All(loops, l => Assert.Equal(3, l.Edges.Count));
    }

    [Fact]
    public void NoCycle_ReturnsEmpty()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(1.0, 0, 1),
        }, "t");

        Barcode barcode = PersistentInvolutedHomology.Compute(filtration, representatives: false);
        var bar = barcode.Bars.First(b => b.Dimension == 0 && !b.IsInfinite);

        Assert.Empty(BarCycleEdges.GetEdgePairs(bar, filtration));
    }
}
