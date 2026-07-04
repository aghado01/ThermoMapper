#nullable enable
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// Z5d (A2, general p) — the dimension-general void-boundary builder (<see cref="CodimensionOneDualGraph"/>)
/// builds the dual graph from ℝᵖ coordinates, and the geometry→dual→duality pipeline must reproduce the
/// (p−1)-th barcode the direct engines compute on F. The p=3 hollow tetrahedron is the milestone (first
/// nontrivial transverse plane); p=2 is exercised through the same general code as the degenerate case.
/// </summary>
public sealed class CodimensionOneDualGraphTests
{
    static IEnumerable<(double, double, int, int, int)> Sig(Barcode bc, int dim) =>
        bc.Bars.Where(b => b.Dimension == dim)
              .Select(b => (b.Birth, b.Death, b.Dimension, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item1).ThenBy(x => x.Item2).ThenBy(x => x.Item4).ThenBy(x => x.Item5);

    // p=3: hollow tetrahedron 0-1-2-3 (triangles 10-13, 3-simplex 14). Same filtration as the hand-dual
    // fixture in CodimensionOneZigzagTests, but the dual is built FROM GEOMETRY.
    static (ZigzagFiltration f, DualGraphSpec dual) HollowTetra()
    {
        var f = new ZigzagFiltration();
        for (int v = 0; v < 4; v++) f.Add(v, new int[0]);
        f.Add(4, new[] { 0, 1 }); f.Add(5, new[] { 0, 2 }); f.Add(6, new[] { 0, 3 });
        f.Add(7, new[] { 1, 2 }); f.Add(8, new[] { 1, 3 }); f.Add(9, new[] { 2, 3 });
        f.Add(10, new[] { 4, 5, 7 }); f.Add(11, new[] { 4, 6, 8 });
        f.Add(12, new[] { 5, 6, 9 }); f.Add(13, new[] { 7, 8, 9 });   // 2-sphere closes (step 13)
        f.Add(14, new[] { 10, 11, 12, 13 });                          // 3-simplex fills cavity (step 14)
        f.Delete(14);
        f.Delete(13); f.Delete(12); f.Delete(11); f.Delete(10);
        f.Delete(9); f.Delete(8); f.Delete(7); f.Delete(6); f.Delete(5); f.Delete(4);
        for (int v = 3; v >= 0; v--) f.Delete(v);

        var dual = CodimensionOneDualGraph.Build(
            p: 3,
            vertexCoords: new Dictionary<int, double[]>
            {
                { 0, new[] { 0.0, 0.0, 0.0 } }, { 1, new[] { 1.0, 0.0, 0.0 } },
                { 2, new[] { 0.0, 1.0, 0.0 } }, { 3, new[] { 0.0, 0.0, 1.0 } },
            },
            pMinus1Simplices: new Dictionary<int, IReadOnlyList<int>>
            {
                { 10, new[] { 0, 1, 2 } }, { 11, new[] { 0, 1, 3 } },
                { 12, new[] { 0, 2, 3 } }, { 13, new[] { 1, 2, 3 } },
            },
            pSimplices: new Dictionary<int, IReadOnlyList<int>> { { 14, new[] { 0, 1, 2, 3 } } },
            firstDualVertexId: 100);
        return (f, dual);
    }

    [Fact]
    public void HollowTetra_GeometryToDualityMatchesGeneralOracles()
    {
        var (f, dual) = HollowTetra();
        var duality = Sig(CodimensionOneZigzag.Compute(f, dual), 2).ToList();   // H₂ = p−1 for p=3
        Assert.Equal(Sig(ZigzagBarcodeNaive.Compute(f, 2), 2), duality);
        Assert.Equal(Sig(FastZigzag.Compute(f, 2), 2), duality);
        Assert.Equal(2, duality.Count);
    }

    [Fact]
    public void HollowTetra_DualStructure()
    {
        var (_, dual) = HollowTetra();
        // One outer void, one 3-simplex dual vertex, four parallel dual edges between them.
        Assert.Single(dual.VoidDualVertexIds);
        Assert.Single(dual.PSimplexDualVertex);
        int vid = dual.VoidDualVertexIds[0], tid = dual.PSimplexDualVertex[14];
        Assert.Equal(4, dual.PMinus1SimplexDualEdge.Count);
        Assert.All(dual.PMinus1SimplexDualEdge.Values,
            e => Assert.True((e.A == vid && e.B == tid) || (e.A == tid && e.B == vid)));
    }

    // p=2 through the GENERAL builder: triangle (filled) and square (pure graph) — the degenerate d=1 case.
    static (ZigzagFiltration f, DualGraphSpec dual) TriangleFill()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 }); f.Add(5, new[] { 0, 2 });
        f.Add(6, new[] { 3, 4, 5 });
        f.Delete(6); f.Delete(5); f.Delete(4); f.Delete(3);
        f.Delete(2); f.Delete(1); f.Delete(0);

        var dual = CodimensionOneDualGraph.Build(
            p: 2,
            vertexCoords: new Dictionary<int, double[]>
                { { 0, new[] { 0.0, 0.0 } }, { 1, new[] { 1.0, 0.0 } }, { 2, new[] { 0.0, 1.0 } } },
            pMinus1Simplices: new Dictionary<int, IReadOnlyList<int>>
                { { 3, new[] { 0, 1 } }, { 4, new[] { 1, 2 } }, { 5, new[] { 0, 2 } } },
            pSimplices: new Dictionary<int, IReadOnlyList<int>> { { 6, new[] { 0, 1, 2 } } },
            firstDualVertexId: 100);
        return (f, dual);
    }

    static (ZigzagFiltration f, DualGraphSpec dual) SquareLoop()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]); f.Add(3, new int[0]);
        f.Add(4, new[] { 0, 1 }); f.Add(5, new[] { 1, 2 }); f.Add(6, new[] { 2, 3 }); f.Add(7, new[] { 0, 3 });
        f.Delete(7); f.Delete(6); f.Delete(5); f.Delete(4);
        f.Delete(3); f.Delete(2); f.Delete(1); f.Delete(0);

        var dual = CodimensionOneDualGraph.Build(
            p: 2,
            vertexCoords: new Dictionary<int, double[]>
                { { 0, new[] { 0.0, 0.0 } }, { 1, new[] { 1.0, 0.0 } }, { 2, new[] { 1.0, 1.0 } }, { 3, new[] { 0.0, 1.0 } } },
            pMinus1Simplices: new Dictionary<int, IReadOnlyList<int>>
                { { 4, new[] { 0, 1 } }, { 5, new[] { 1, 2 } }, { 6, new[] { 2, 3 } }, { 7, new[] { 0, 3 } } },
            pSimplices: new Dictionary<int, IReadOnlyList<int>>(),
            firstDualVertexId: 100);
        return (f, dual);
    }

    [Fact]
    public void TriangleFill_GeneralBuilderMatchesOracles()
    {
        var (f, dual) = TriangleFill();
        var duality = Sig(CodimensionOneZigzag.Compute(f, dual), 1).ToList();
        Assert.Equal(Sig(ZigzagBarcodeNaive.Compute(f, 1), 1), duality);
        Assert.Equal(Sig(FastZigzag.Compute(f, 1), 1), duality);
        // One outer void + the 2-simplex dual, three parallel edges (matches PlanarDualGraph).
        Assert.Single(dual.VoidDualVertexIds);
        Assert.Single(dual.PSimplexDualVertex);
        Assert.Equal(3, dual.PMinus1SimplexDualEdge.Count);
    }

    [Fact]
    public void SquareLoop_GeneralBuilderMatchesAllEngines()
    {
        var (f, dual) = SquareLoop();
        var duality = Sig(CodimensionOneZigzag.Compute(f, dual), 1).ToList();
        Assert.Equal(Sig(ZigzagBarcodeNaive.Compute(f, 1), 1), duality);
        Assert.Equal(Sig(FastZigzag.Compute(f, 1), 1), duality);
        Assert.Equal(Sig(GraphZigzagH1.Compute(f), 1), duality);   // pure graph: independent H1 route
        // Two voids (inside / outside), no 2-simplex duals, four parallel edges.
        Assert.Equal(2, dual.VoidDualVertexIds.Count);
        Assert.Empty(dual.PSimplexDualVertex);
        Assert.Equal(4, dual.PMinus1SimplexDualEdge.Count);
    }

    [Fact]
    public void TwoTetraSharedFace_InteriorSimplexDualsTwoCells()
    {
        // Two tetrahedra 0-1-2-3 (above) and 0-1-2-4 (below) sharing triangle 0-1-2 (a bipyramid).
        // The shared triangle is INTERIOR (two p-cofaces) → its dual edge joins the two 3-simplex duals;
        // the other six triangles are boundary → joined to the single outer void.
        var coords = new Dictionary<int, double[]>
        {
            { 0, new[] { 0.0, 0.0, 0.0 } }, { 1, new[] { 1.0, 0.0, 0.0 } }, { 2, new[] { 0.0, 1.0, 0.0 } },
            { 3, new[] { 0.0, 0.0, 1.0 } }, { 4, new[] { 0.0, 0.0, -1.0 } },
        };
        var tris = new Dictionary<int, IReadOnlyList<int>>
        {
            { 20, new[] { 0, 1, 2 } },                                   // shared (interior)
            { 21, new[] { 0, 1, 3 } }, { 22, new[] { 0, 2, 3 } }, { 23, new[] { 1, 2, 3 } },
            { 24, new[] { 0, 1, 4 } }, { 25, new[] { 0, 2, 4 } }, { 26, new[] { 1, 2, 4 } },
        };
        var tetras = new Dictionary<int, IReadOnlyList<int>>
            { { 30, new[] { 0, 1, 2, 3 } }, { 31, new[] { 0, 1, 2, 4 } } };

        var dual = CodimensionOneDualGraph.Build(3, coords, tris, tetras, firstDualVertexId: 100);

        Assert.Equal(2, dual.PSimplexDualVertex.Count);
        Assert.Single(dual.VoidDualVertexIds);
        int dA = dual.PSimplexDualVertex[30], dB = dual.PSimplexDualVertex[31];
        int vid = dual.VoidDualVertexIds[0];

        // shared triangle 20 → edge between the two tetra duals (interior, no void)
        var shared = dual.PMinus1SimplexDualEdge[20];
        Assert.True((shared.A == dA && shared.B == dB) || (shared.A == dB && shared.B == dA));

        // the six outer triangles → each joins a tetra dual to the outer void
        foreach (int t in new[] { 21, 22, 23, 24, 25, 26 })
        {
            var e = dual.PMinus1SimplexDualEdge[t];
            Assert.True(e.A == vid || e.B == vid);
            Assert.True(e.A == dA || e.A == dB || e.B == dA || e.B == dB);
        }
    }
}
