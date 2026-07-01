#nullable enable
using System.Collections.Generic;
using System.Linq;
using Xunit;

using TDA.Primitives;
namespace TDA.Ph.Tests;

/// <summary>
/// Z5d (item C/F) — the arbitrary-complex front door <see cref="EmbeddedComplexZigzag"/>: decompose K into
/// (p−1)-connected pieces, build each piece's dual from geometry, reduce, and union the barcodes (Prop 28).
/// Squared away on <b>both</b> H0 reduction paths (Reference oracle and the multigraph-safe Fast/HDT engine):
/// every fixture asserts Reference == Fast == the general oracles (<see cref="ZigzagBarcodeNaive"/>,
/// <see cref="FastZigzag"/>) at dimension p−1.
/// </summary>
public sealed class EmbeddedComplexZigzagTests
{
    static IEnumerable<(double, double, int, int, int)> Sig(Barcode bc, int dim) =>
        bc.Bars.Where(b => b.Dimension == dim)
              .Select(b => (b.Birth, b.Death, b.Dimension, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item1).ThenBy(x => x.Item2).ThenBy(x => x.Item4).ThenBy(x => x.Item5);

    static void AssertBothPathsMatchOracles(
        int p,
        Dictionary<int, double[]> coords,
        Dictionary<int, IReadOnlyList<int>> tris,
        Dictionary<int, IReadOnlyList<int>> tetras,
        ZigzagFiltration f,
        int expectedCount)
    {
        var naive = Sig(ZigzagBarcodeNaive.Compute(f, p - 1), p - 1).ToList();
        var fast = Sig(FastZigzag.Compute(f, p - 1), p - 1).ToList();
        Assert.Equal(naive, fast);   // the two general oracles agree (ground truth)

        var refPath = Sig(EmbeddedComplexZigzag.Compute(p, coords, tris, tetras, f, GraphZigzagAlgorithm.Reference), p - 1).ToList();
        var fastPath = Sig(EmbeddedComplexZigzag.Compute(p, coords, tris, tetras, f, GraphZigzagAlgorithm.Fast), p - 1).ToList();

        Assert.Equal(naive, refPath);     // slow path == oracles
        Assert.Equal(refPath, fastPath);  // fast path == slow path (the multigraph-safe HDT engine)
        Assert.Equal(expectedCount, refPath.Count);
    }

    // Two disjoint hollow tetrahedra (filled 2-spheres): A on vertices 0-3 near the origin, B on 20-23
    // translated +10 in x. Each closes its 2-sphere, is filled by a 3-simplex, then everything tears down
    // (∅→∅). The two triangle sets share no edge → decomposition must split K into two (p−1)-connected pieces.
    static (ZigzagFiltration f,
            Dictionary<int, double[]> coords,
            Dictionary<int, IReadOnlyList<int>> tris,
            Dictionary<int, IReadOnlyList<int>> tetras) TwoDisjointSpheres()
    {
        var f = new ZigzagFiltration();

        for (int v = 0; v < 4; v++) f.Add(v, new int[0]);
        f.Add(4, new[] { 0, 1 }); f.Add(5, new[] { 0, 2 }); f.Add(6, new[] { 0, 3 });
        f.Add(7, new[] { 1, 2 }); f.Add(8, new[] { 1, 3 }); f.Add(9, new[] { 2, 3 });
        f.Add(10, new[] { 4, 5, 7 }); f.Add(11, new[] { 4, 6, 8 });
        f.Add(12, new[] { 5, 6, 9 }); f.Add(13, new[] { 7, 8, 9 });
        f.Add(14, new[] { 10, 11, 12, 13 });                          // A filled

        for (int v = 20; v < 24; v++) f.Add(v, new int[0]);
        f.Add(24, new[] { 20, 21 }); f.Add(25, new[] { 20, 22 }); f.Add(26, new[] { 20, 23 });
        f.Add(27, new[] { 21, 22 }); f.Add(28, new[] { 21, 23 }); f.Add(29, new[] { 22, 23 });
        f.Add(30, new[] { 24, 25, 27 }); f.Add(31, new[] { 24, 26, 28 });
        f.Add(32, new[] { 25, 26, 29 }); f.Add(33, new[] { 27, 28, 29 });
        f.Add(34, new[] { 30, 31, 32, 33 });                          // B filled

        f.Delete(34);
        f.Delete(33); f.Delete(32); f.Delete(31); f.Delete(30);
        f.Delete(29); f.Delete(28); f.Delete(27); f.Delete(26); f.Delete(25); f.Delete(24);
        for (int v = 23; v >= 20; v--) f.Delete(v);
        f.Delete(14);
        f.Delete(13); f.Delete(12); f.Delete(11); f.Delete(10);
        f.Delete(9); f.Delete(8); f.Delete(7); f.Delete(6); f.Delete(5); f.Delete(4);
        for (int v = 3; v >= 0; v--) f.Delete(v);

        var coords = new Dictionary<int, double[]>
        {
            { 0, new[] { 0.0, 0.0, 0.0 } }, { 1, new[] { 1.0, 0.0, 0.0 } },
            { 2, new[] { 0.0, 1.0, 0.0 } }, { 3, new[] { 0.0, 0.0, 1.0 } },
            { 20, new[] { 10.0, 0.0, 0.0 } }, { 21, new[] { 11.0, 0.0, 0.0 } },
            { 22, new[] { 10.0, 1.0, 0.0 } }, { 23, new[] { 10.0, 0.0, 1.0 } },
        };
        var tris = new Dictionary<int, IReadOnlyList<int>>
        {
            { 10, new[] { 0, 1, 2 } }, { 11, new[] { 0, 1, 3 } }, { 12, new[] { 0, 2, 3 } }, { 13, new[] { 1, 2, 3 } },
            { 30, new[] { 20, 21, 22 } }, { 31, new[] { 20, 21, 23 } }, { 32, new[] { 20, 22, 23 } }, { 33, new[] { 21, 22, 23 } },
        };
        var tetras = new Dictionary<int, IReadOnlyList<int>>
        {
            { 14, new[] { 0, 1, 2, 3 } }, { 34, new[] { 20, 21, 22, 23 } },
        };
        return (f, coords, tris, tetras);
    }

    // Two CONCENTRIC unfilled hollow tetrahedra (nested 2-spheres): outer on 0-3 (large), inner on 20-23
    // (small), both centered at the origin, neither filled. The complement of full K has three regions —
    // outside-outer, the between-shell, inside-inner — and the between-void's boundary is DISCONNECTED (outer's
    // inner face ∪ inner's outer face). That is the regime where (p−1)-connectedness is load-bearing for
    // correctness (Thm 4.1: one union-find class per void): Build-on-whole would mis-split the between-void, so
    // the per-piece decomposition is *required*, not merely consistent. Each piece is an unfilled hollow tetra
    // → a 2-void multigraph dual (4 parallel edges), which also exercises the multigraph-safe Fast path.
    static (ZigzagFiltration f,
            Dictionary<int, double[]> coords,
            Dictionary<int, IReadOnlyList<int>> tris,
            Dictionary<int, IReadOnlyList<int>> tetras) NestedSpheres()
    {
        var f = new ZigzagFiltration();

        // outer 2-sphere (no 3-simplex)
        for (int v = 0; v < 4; v++) f.Add(v, new int[0]);
        f.Add(4, new[] { 0, 1 }); f.Add(5, new[] { 0, 2 }); f.Add(6, new[] { 0, 3 });
        f.Add(7, new[] { 1, 2 }); f.Add(8, new[] { 1, 3 }); f.Add(9, new[] { 2, 3 });
        f.Add(10, new[] { 4, 5, 7 }); f.Add(11, new[] { 4, 6, 8 });
        f.Add(12, new[] { 5, 6, 9 }); f.Add(13, new[] { 7, 8, 9 });   // outer sphere closes

        // inner 2-sphere (no 3-simplex)
        for (int v = 20; v < 24; v++) f.Add(v, new int[0]);
        f.Add(24, new[] { 20, 21 }); f.Add(25, new[] { 20, 22 }); f.Add(26, new[] { 20, 23 });
        f.Add(27, new[] { 21, 22 }); f.Add(28, new[] { 21, 23 }); f.Add(29, new[] { 22, 23 });
        f.Add(30, new[] { 24, 25, 27 }); f.Add(31, new[] { 24, 26, 28 });
        f.Add(32, new[] { 25, 26, 29 }); f.Add(33, new[] { 27, 28, 29 });   // inner sphere closes

        // tear down (∅→∅)
        f.Delete(33); f.Delete(32); f.Delete(31); f.Delete(30);
        f.Delete(29); f.Delete(28); f.Delete(27); f.Delete(26); f.Delete(25); f.Delete(24);
        for (int v = 23; v >= 20; v--) f.Delete(v);
        f.Delete(13); f.Delete(12); f.Delete(11); f.Delete(10);
        f.Delete(9); f.Delete(8); f.Delete(7); f.Delete(6); f.Delete(5); f.Delete(4);
        for (int v = 3; v >= 0; v--) f.Delete(v);

        var coords = new Dictionary<int, double[]>
        {
            // outer tetra (radius ~3), centered at origin
            { 0, new[] { 3.0, 3.0, 3.0 } }, { 1, new[] { 3.0, -3.0, -3.0 } },
            { 2, new[] { -3.0, 3.0, -3.0 } }, { 3, new[] { -3.0, -3.0, 3.0 } },
            // inner tetra (radius ~1), same centroid → strictly inside the outer solid
            { 20, new[] { 1.0, 1.0, 1.0 } }, { 21, new[] { 1.0, -1.0, -1.0 } },
            { 22, new[] { -1.0, 1.0, -1.0 } }, { 23, new[] { -1.0, -1.0, 1.0 } },
        };
        var tris = new Dictionary<int, IReadOnlyList<int>>
        {
            { 10, new[] { 0, 1, 2 } }, { 11, new[] { 0, 1, 3 } }, { 12, new[] { 0, 2, 3 } }, { 13, new[] { 1, 2, 3 } },
            { 30, new[] { 20, 21, 22 } }, { 31, new[] { 20, 21, 23 } }, { 32, new[] { 20, 22, 23 } }, { 33, new[] { 21, 22, 23 } },
        };
        var tetras = new Dictionary<int, IReadOnlyList<int>>();   // unfilled
        return (f, coords, tris, tetras);
    }

    [Fact]
    public void TwoDisjointSpheres_DecomposeUnionMatchesOraclesBothPaths()
    {
        var (f, coords, tris, tetras) = TwoDisjointSpheres();
        // Four H₂ features: each sphere closes→fills (2 bars) and unfills→breaks on teardown (2 bars).
        AssertBothPathsMatchOracles(3, coords, tris, tetras, f, 4);
    }

    [Fact]
    public void NestedSpheres_DecompositionIsLoadBearingBothPaths()
    {
        var (f, coords, tris, tetras) = NestedSpheres();
        // Two H₂ features (the two independent 2-cycles), each born when its sphere closes and dying on teardown.
        // Requires decomposition (the between-void has a disconnected boundary); also exercises the 2-void
        // multigraph dual on both paths.
        AssertBothPathsMatchOracles(3, coords, tris, tetras, f, 2);
    }

    [Fact]
    public void SinglePiece_FrontDoorAgreesWithPerPieceOnBothPaths()
    {
        // One (p−1)-connected piece must route through decomposition unchanged: the front door equals the
        // direct per-piece builder (CodimensionOneDualGraph + CodimensionOneZigzag) on both engines.
        var (f, coords, tris, tetras) = TwoDisjointSpheres();
        var aTris = tris.Where(kv => kv.Key < 20).ToDictionary(kv => kv.Key, kv => kv.Value);
        var aTetras = tetras.Where(kv => kv.Key < 20).ToDictionary(kv => kv.Key, kv => kv.Value);

        foreach (var alg in new[] { GraphZigzagAlgorithm.Reference, GraphZigzagAlgorithm.Fast })
        {
            var frontDoor = Sig(EmbeddedComplexZigzag.Compute(3, coords, aTris, aTetras, f, alg), 2).ToList();
            var perPiece = Sig(CodimensionOneZigzag.Compute(
                f, CodimensionOneDualGraph.Build(3, coords, aTris, aTetras, firstDualVertexId: 100), alg), 2).ToList();
            Assert.Equal(perPiece, frontDoor);
            Assert.Equal(2, frontDoor.Count);
        }
    }
}
