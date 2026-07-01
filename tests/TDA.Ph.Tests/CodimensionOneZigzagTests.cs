#nullable enable
using System.Collections.Generic;
using System.Linq;
using Xunit;

using Maths.Topology;
namespace TDA.Ph.Tests;

/// <summary>
/// Z5d (A1 skeleton) — codimension-one zigzag via Alexander duality (<see cref="CodimensionOneZigzag"/>,
/// Dey–Hou §5). The duality route reduces H_{p−1}(F) to H̃₀ of a dual filtration and runs the H0 engine;
/// its output must equal the (p−1)-dimensional barcode that every direct engine computes on F itself.
/// p=2 fixture: a triangle whose loop forms then fills, then tears back down to ∅ — exercising both the
/// split-by-add (loop close) and split-by-delete (loop re-open) sides of the end-flip.
/// </summary>
public sealed class CodimensionOneZigzagTests
{
    static IEnumerable<(double, double, int, int, int)> Sig(Barcode bc, int dim) =>
        bc.Bars.Where(b => b.Dimension == dim)
              .Select(b => (b.Birth, b.Death, b.Dimension, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item1).ThenBy(x => x.Item2).ThenBy(x => x.Item4).ThenBy(x => x.Item5);

    // p=2: triangle 0-1-2, edges 3,4,5, 2-simplex 6. Loop closes at step 5, fills at step 6, then
    // everything is deleted (∅→∅, the zigzag convention). Dual graph: one outside void (id 100) and the
    // 2-simplex's dual vertex (id 101), joined by three dual edges (one per triangle edge).
    static (ZigzagFiltration f, DualGraphSpec dual) TriangleFillUpDown()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 }); f.Add(5, new[] { 0, 2 }); // loop closes (step 5)
        f.Add(6, new[] { 3, 4, 5 });                                                   // 2-simplex fills (step 6)
        f.Delete(6);                                                                   // step 7: loop reappears
        f.Delete(5); f.Delete(4); f.Delete(3);                                         // steps 8-10
        f.Delete(2); f.Delete(1); f.Delete(0);                                         // steps 11-13

        var dual = new DualGraphSpec(
            p: 2,
            voidDualVertexIds: new[] { 100 },
            pSimplexDualVertex: new Dictionary<int, int> { { 6, 101 } },
            pMinus1SimplexDualEdge: new Dictionary<int, (int, int)>
            {
                { 3, (100, 101) }, { 4, (100, 101) }, { 5, (100, 101) },
            });
        return (f, dual);
    }

    [Fact]
    public void TriangleFillUpDown_DualityMatchesGeneralOracles()
    {
        var (f, dual) = TriangleFillUpDown();
        var duality = Sig(CodimensionOneZigzag.Compute(f, dual), 1).ToList();

        // The general (all-dimensional) oracles see the 2-simplex fill, so they are the correct ground
        // truth for the duality. NOT GraphZigzagH1 — it is a graph-only engine that cannot process the
        // 2-simplex, so it reports the uncut [5,8] loop; the duality (correctly) cuts it at the fill.
        // A GraphZigzagH1 cross-check needs a pure-graph (no 2-simplex) fixture — deferred.
        Assert.Equal(Sig(ZigzagBarcodeNaive.Compute(f, 1), 1), duality);
        Assert.Equal(Sig(FastZigzag.Compute(f, 1), 1), duality);
    }

    [Fact]
    public void TriangleFillUpDown_AnchorBars()
    {
        var (f, dual) = TriangleFillUpDown();
        var bars = CodimensionOneZigzag.Compute(f, dual).Bars.OrderBy(b => b.Birth).ToList();

        // Anti-circularity anchor: two H1 features. The fill is born when the loop closes (step 5, Add ⇒
        // Closed) and dies when the 2-simplex fills it (step 6, Add ⇒ Open); the teardown re-opens a loop
        // born at the 2-simplex delete (step 7, Delete ⇒ Open) that dies when an edge is removed (step 8,
        // Delete ⇒ Closed). The dual computes both with the ends reversed; the engine flips them back.
        Assert.Equal(2, bars.Count);
        Assert.Equal((5.0, 6.0, 1, IntervalEnd.Closed, IntervalEnd.Open),
            (bars[0].Birth, bars[0].Death, bars[0].Dimension, bars[0].BirthEnd, bars[0].DeathEnd));
        Assert.Equal((7.0, 8.0, 1, IntervalEnd.Open, IntervalEnd.Closed),
            (bars[1].Birth, bars[1].Death, bars[1].Dimension, bars[1].BirthEnd, bars[1].DeathEnd));
    }

    // p=3: hollow tetrahedron on 0-1-2-3 (6 edges 4-9, 4 triangles 10-13, 3-simplex 14). The 2-sphere of four
    // triangles closes → H₂ born; the tetrahedron fills the cavity → dies; teardown re-opens it (∅→∅). Dual:
    // outside void (100) + the 3-simplex's dual vertex (101), joined by four parallel edges (the 4 triangles).
    static (ZigzagFiltration f, DualGraphSpec dual) HollowTetraFillUpDown()
    {
        var f = new ZigzagFiltration();
        for (int v = 0; v < 4; v++) f.Add(v, new int[0]);
        f.Add(4, new[] { 0, 1 }); f.Add(5, new[] { 0, 2 }); f.Add(6, new[] { 0, 3 });
        f.Add(7, new[] { 1, 2 }); f.Add(8, new[] { 1, 3 }); f.Add(9, new[] { 2, 3 });
        f.Add(10, new[] { 4, 5, 7 }); f.Add(11, new[] { 4, 6, 8 });
        f.Add(12, new[] { 5, 6, 9 }); f.Add(13, new[] { 7, 8, 9 });   // 2-sphere closes (step 13)
        f.Add(14, new[] { 10, 11, 12, 13 });                          // 3-simplex fills cavity (step 14)
        f.Delete(14);                                                 // step 15: cavity reopens
        f.Delete(13); f.Delete(12); f.Delete(11); f.Delete(10);
        f.Delete(9); f.Delete(8); f.Delete(7); f.Delete(6); f.Delete(5); f.Delete(4);
        for (int v = 3; v >= 0; v--) f.Delete(v);

        var dual = new DualGraphSpec(
            p: 3,
            voidDualVertexIds: new[] { 100 },
            pSimplexDualVertex: new Dictionary<int, int> { { 14, 101 } },
            pMinus1SimplexDualEdge: new Dictionary<int, (int, int)>
            {
                { 10, (100, 101) }, { 11, (100, 101) }, { 12, (100, 101) }, { 13, (100, 101) },
            });
        return (f, dual);
    }

    [Fact]
    public void HollowTetraFillUpDown_DualityMatchesGeneralOracles()
    {
        var (f, dual) = HollowTetraFillUpDown();
        var duality = Sig(CodimensionOneZigzag.Compute(f, dual), 2).ToList();   // H₂ = p−1 for p=3

        // Proves the duality core is dimension-general (not p=2-specific) against the general oracles at dim 2.
        Assert.Equal(Sig(ZigzagBarcodeNaive.Compute(f, 2), 2), duality);
        Assert.Equal(Sig(FastZigzag.Compute(f, 2), 2), duality);

        // Two H₂ features (sphere closes→fills, then unfill→break) — the p=3 analogue of the triangle's two H₁ bars.
        Assert.Equal(2, duality.Count);
        Assert.Equal((13.0, 14.0, 2, (int)IntervalEnd.Closed, (int)IntervalEnd.Open), duality[0]);
        Assert.Equal((15.0, 16.0, 2, (int)IntervalEnd.Open, (int)IntervalEnd.Closed), duality[1]);
    }
}
