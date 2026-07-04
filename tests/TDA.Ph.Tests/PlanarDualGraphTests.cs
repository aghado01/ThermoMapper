#nullable enable
using System.Collections.Generic;
using System.Linq;
using Xunit;

using Maths.Topology;
namespace TDA.Ph.Tests;

/// <summary>
/// Z5d (A2, p=2) — the planar face-tracer (<see cref="PlanarDualGraph"/>) builds the dual graph from
/// coordinates, and the full geometry→dual→duality pipeline must reproduce the (p−1)=1 barcode the direct
/// engines compute on F. No hand-supplied dual: this proves the geometry, not just the reduction core.
/// </summary>
public sealed class PlanarDualGraphTests
{
    static IEnumerable<(double, double, int, int, int)> Sig(Barcode bc, int dim) =>
        bc.Bars.Where(b => b.Dimension == dim)
              .Select(b => (b.Birth, b.Death, b.Dimension, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item1).ThenBy(x => x.Item2).ThenBy(x => x.Item4).ThenBy(x => x.Item5);

    // p=2, FILLED: triangle 0-1-2 (edges 3,4,5; 2-simplex 6) — loop forms → fills → tears down (∅→∅).
    static (ZigzagFiltration f, DualGraphSpec dual) TriangleFill()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 }); f.Add(5, new[] { 0, 2 });
        f.Add(6, new[] { 3, 4, 5 });
        f.Delete(6); f.Delete(5); f.Delete(4); f.Delete(3);
        f.Delete(2); f.Delete(1); f.Delete(0);

        var dual = PlanarDualGraph.Build(
            vertexCoords: new Dictionary<int, (double, double)> { { 0, (0, 0) }, { 1, (1, 0) }, { 2, (0, 1) } },
            edges: new Dictionary<int, (int, int)> { { 3, (0, 1) }, { 4, (1, 2) }, { 5, (0, 2) } },
            triangleVertices: new Dictionary<int, IReadOnlyList<int>> { { 6, new[] { 0, 1, 2 } } },
            firstDualVertexId: 100);
        return (f, dual);
    }

    // p=2, PURE GRAPH (no 2-simplex): square 0-1-2-3 (edges 4,5,6,7) — loop forms then opens (∅→∅).
    static (ZigzagFiltration f, DualGraphSpec dual) SquareLoop()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]); f.Add(3, new int[0]);
        f.Add(4, new[] { 0, 1 }); f.Add(5, new[] { 1, 2 }); f.Add(6, new[] { 2, 3 }); f.Add(7, new[] { 0, 3 });
        f.Delete(7); f.Delete(6); f.Delete(5); f.Delete(4);
        f.Delete(3); f.Delete(2); f.Delete(1); f.Delete(0);

        var dual = PlanarDualGraph.Build(
            vertexCoords: new Dictionary<int, (double, double)>
                { { 0, (0, 0) }, { 1, (1, 0) }, { 2, (1, 1) }, { 3, (0, 1) } },
            edges: new Dictionary<int, (int, int)>
                { { 4, (0, 1) }, { 5, (1, 2) }, { 6, (2, 3) }, { 7, (0, 3) } },
            triangleVertices: new Dictionary<int, IReadOnlyList<int>>(),
            firstDualVertexId: 100);
        return (f, dual);
    }

    [Fact]
    public void TriangleFill_GeometryToDualityMatchesGeneralOracles()
    {
        var (f, dual) = TriangleFill();
        var duality = Sig(CodimensionOneZigzag.Compute(f, dual), 1).ToList();
        Assert.Equal(Sig(ZigzagBarcodeNaive.Compute(f, 1), 1), duality);
        Assert.Equal(Sig(FastZigzag.Compute(f, 1), 1), duality);
        // Two H1 features (form→fill, then unfill→open) — same as the hand-supplied-dual skeleton.
        Assert.Equal(2, duality.Count);
    }

    [Fact]
    public void TriangleFill_DualStructure()
    {
        var (_, dual) = TriangleFill();
        // One outside void, one 2-simplex dual vertex, three parallel dual edges between them.
        Assert.Single(dual.VoidDualVertexIds);
        Assert.Single(dual.PSimplexDualVertex);
        int vid = dual.VoidDualVertexIds[0], tid = dual.PSimplexDualVertex[6];
        Assert.Equal(3, dual.PMinus1SimplexDualEdge.Count);
        Assert.All(dual.PMinus1SimplexDualEdge.Values,
            e => Assert.True((e.A == vid && e.B == tid) || (e.A == tid && e.B == vid)));
    }

    [Fact]
    public void SquareLoop_GeometryToDualityMatchesAllEngines()
    {
        var (f, dual) = SquareLoop();
        var duality = Sig(CodimensionOneZigzag.Compute(f, dual), 1).ToList();
        Assert.Equal(Sig(ZigzagBarcodeNaive.Compute(f, 1), 1), duality);
        Assert.Equal(Sig(FastZigzag.Compute(f, 1), 1), duality);
        // Pure graph (no 2-simplex), so GraphZigzagH1 IS a valid independent route here (regime 1).
        Assert.Equal(Sig(GraphZigzagH1.Compute(f), 1), duality);
        // The single loop: born when the square closes (step 7, Add), dies when an edge leaves (step 8, Delete).
        Assert.Equal((7.0, 8.0, 1, (int)IntervalEnd.Closed, (int)IntervalEnd.Closed), Assert.Single(duality));
    }

    [Fact]
    public void SquareLoop_DualStructure()
    {
        var (_, dual) = SquareLoop();
        // Two voids (inside / outside), no 2-simplex duals, four parallel dual edges between the two voids.
        Assert.Equal(2, dual.VoidDualVertexIds.Count);
        Assert.Empty(dual.PSimplexDualVertex);
        Assert.Equal(4, dual.PMinus1SimplexDualEdge.Count);
        int a = dual.VoidDualVertexIds[0], b = dual.VoidDualVertexIds[1];
        Assert.All(dual.PMinus1SimplexDualEdge.Values,
            e => Assert.True((e.A == a && e.B == b) || (e.A == b && e.B == a)));
    }
}
