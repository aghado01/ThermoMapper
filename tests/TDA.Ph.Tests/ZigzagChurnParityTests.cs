#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// Z5d hardening — random-churn parity for the codimension-one duality. The monotone build-up-then-teardown
/// fixtures exercise only one shape of zigzag; these sweeps drive each fixed embedded complex through many
/// random <i>valid</i> zigzag filtrations (interleaved adds/deletes — faces before a coface, never a coface
/// before its face) and assert the duality on <b>both</b> H0 paths reproduces the (p−1)-th barcode the general
/// oracles compute on F. Stresses the reversal-flip, the reduced-homology drop, the no-op collapse, the union
/// of (p−1)-connected pieces, and the multigraph Fast engine under non-monotone churn. A negative test pins the
/// endpoint-presence guard against a malformed filtration.
/// <para><b>Input class.</b> The generator exercises full cell <i>re-entry</i> (delete-then-re-add). The churn
/// sweep originally surfaced that <see cref="FastZigzag"/> (Z2) conflated a re-added cell with its prior
/// incarnation and diverged from the rank-definitional <see cref="ZigzagBarcodeNaive"/> (Z1); that cone bug is
/// fixed (per-incarnation accounting), so all four engines now agree on re-entry filtrations too.</para>
/// </summary>
public sealed class ZigzagChurnParityTests
{
    // A fixed embedded complex: every cell by vertex set (the full closure), vertex coordinates, embedding dim p.
    sealed class Complex
    {
        public int P;
        public Dictionary<int, double[]> Coords = new();
        public (int Id, int[] Verts)[] Cells = Array.Empty<(int, int[])>();

        public readonly Dictionary<int, int[]> FaceCells = new();             // cell id -> (dim-1)-face cell ids
        public readonly Dictionary<int, IReadOnlyList<int>> PMinus1 = new();  // dim == p-1
        public readonly Dictionary<int, IReadOnlyList<int>> PSimplices = new(); // dim == p

        public Complex Derive()
        {
            var byKey = new Dictionary<string, int>();
            foreach (var (id, verts) in Cells) byKey[Key(verts)] = id;
            foreach (var (id, verts) in Cells)
            {
                var v = verts.OrderBy(x => x).ToArray();
                var faces = new List<int>();
                if (v.Length >= 2)
                    for (int drop = 0; drop < v.Length; drop++)
                        faces.Add(byKey[Key(Omit(v, drop))]);   // closure: every face must be listed
                FaceCells[id] = faces.ToArray();
                int dim = v.Length - 1;
                if (dim == P - 1) PMinus1[id] = v;
                if (dim == P) PSimplices[id] = v;
            }
            return this;
        }

        static string Key(int[] verts) => string.Join(",", verts.OrderBy(x => x));
        static int[] Omit(int[] v, int drop) { var r = new int[v.Length - 1]; for (int i = 0, j = 0; i < v.Length; i++) if (i != drop) r[j++] = v[i]; return r; }
    }

    // A random valid ∅→∅ zigzag over the complex: each churn step adds a cell whose faces are all present (and
    // it absent) or deletes a present cell with no present coface; then tears everything down. Deterministic in seed.
    static ZigzagFiltration RandomZigzag(Complex k, int seed, int churnSteps, bool allowReentry = false)
    {
        var rng = new Random(seed);
        var present = new HashSet<int>();
        var deletedEver = new HashSet<int>();   // each cell gets at most one lifetime (simplex-wise zigzag, no re-entry)
        var f = new ZigzagFiltration();

        var cofaces = k.Cells.ToDictionary(c => c.Id, _ => new List<int>());
        foreach (var (id, _) in k.Cells) foreach (int fc in k.FaceCells[id]) cofaces[fc].Add(id);

        bool Addable(int id) => !present.Contains(id) && (allowReentry || !deletedEver.Contains(id)) && k.FaceCells[id].All(present.Contains);
        bool Deletable(int id) => present.Contains(id) && cofaces[id].All(c => !present.Contains(c));
        void DoAdd(int id) { f.Add(id, k.FaceCells[id]); present.Add(id); }
        void DoDel(int id) { f.Delete(id); present.Remove(id); deletedEver.Add(id); }

        for (int s = 0; s < churnSteps; s++)
        {
            var adds = k.Cells.Where(c => Addable(c.Id)).Select(c => c.Id).ToList();
            var dels = present.Where(Deletable).ToList();
            if (adds.Count == 0 && dels.Count == 0) break;
            bool doAdd = dels.Count == 0 || (adds.Count > 0 && rng.NextDouble() < 0.6);
            if (doAdd) DoAdd(adds[rng.Next(adds.Count)]);
            else DoDel(dels[rng.Next(dels.Count)]);
        }
        // teardown to ∅: repeatedly delete a maximal present cell (a present cell with no present coface exists
        // whenever present is non-empty).
        while (present.Count > 0) DoDel(present.First(Deletable));
        return f;
    }

    static IEnumerable<(double, double, int, int, int)> Sig(Barcode bc, int dim) =>
        bc.Bars.Where(b => b.Dimension == dim)
              .Select(b => (b.Birth, b.Death, b.Dimension, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item1).ThenBy(x => x.Item2).ThenBy(x => x.Item4).ThenBy(x => x.Item5);

    static void AssertChurnParity(Complex k, int seeds, int churnSteps)
    {
        k.Derive();
        int d = k.P - 1;
        for (int seed = 0; seed < seeds; seed++)
        {
            var f = RandomZigzag(k, seed, churnSteps, allowReentry: true);
            var naive = Sig(ZigzagBarcodeNaive.Compute(f, d), d).ToList();
            var fast = Sig(FastZigzag.Compute(f, d), d).ToList();
            var refPath = Sig(EmbeddedComplexZigzag.Compute(k.P, k.Coords, k.PMinus1, k.PSimplices, f, GraphZigzagAlgorithm.Reference), d).ToList();
            var fastPath = Sig(EmbeddedComplexZigzag.Compute(k.P, k.Coords, k.PMinus1, k.PSimplices, f, GraphZigzagAlgorithm.Fast), d).ToList();

            string ctx = $"p={k.P} seed={seed} steps={f.Count}";
            Assert.True(naive.SequenceEqual(fast), $"oracles disagree ({ctx})");
            Assert.True(naive.SequenceEqual(refPath), $"Reference duality != oracle ({ctx})");
            Assert.True(naive.SequenceEqual(fastPath), $"Fast duality != oracle ({ctx})");
        }
    }

    // ── Fixed complexes (cells by vertex set; faces derived) ──────────────────────────────────────────

    static Complex HollowTetra() => new()
    {
        P = 3,
        Coords = new() { { 0, new[] { 0.0, 0, 0 } }, { 1, new[] { 1.0, 0, 0 } }, { 2, new[] { 0.0, 1, 0 } }, { 3, new[] { 0.0, 0, 1 } } },
        Cells = new (int, int[])[]
        {
            (0, new[]{0}), (1, new[]{1}), (2, new[]{2}), (3, new[]{3}),
            (4, new[]{0,1}), (5, new[]{0,2}), (6, new[]{0,3}), (7, new[]{1,2}), (8, new[]{1,3}), (9, new[]{2,3}),
            (10, new[]{0,1,2}), (11, new[]{0,1,3}), (12, new[]{0,2,3}), (13, new[]{1,2,3}),
            (14, new[]{0,1,2,3}),
        },
    };

    static Complex TriangleFill() => new()
    {
        P = 2,
        Coords = new() { { 0, new[] { 0.0, 0 } }, { 1, new[] { 1.0, 0 } }, { 2, new[] { 0.0, 1 } } },
        Cells = new (int, int[])[]
        {
            (0, new[]{0}), (1, new[]{1}), (2, new[]{2}),
            (3, new[]{0,1}), (4, new[]{1,2}), (5, new[]{0,2}),
            (6, new[]{0,1,2}),
        },
    };

    // Bipyramid: two tetrahedra sharing triangle {0,1,2}. The shared triangle is an INTERIOR (p−1)-simplex
    // (two p-cofaces → a dual edge between the two p-simplex duals), so churn exercises interior-edge re-entry.
    static Complex Bipyramid() => new()
    {
        P = 3,
        Coords = new()
        {
            { 0, new[] { 0.0, 0, 0 } }, { 1, new[] { 1.0, 0, 0 } }, { 2, new[] { 0.0, 1, 0 } },
            { 3, new[] { 0.0, 0, 1 } }, { 4, new[] { 0.0, 0, -1 } },
        },
        Cells = new (int, int[])[]
        {
            (0, new[]{0}), (1, new[]{1}), (2, new[]{2}), (3, new[]{3}), (4, new[]{4}),
            (5, new[]{0,1}), (6, new[]{0,2}), (7, new[]{0,3}), (8, new[]{0,4}), (9, new[]{1,2}),
            (10, new[]{1,3}), (11, new[]{1,4}), (12, new[]{2,3}), (13, new[]{2,4}),
            (20, new[]{0,1,2}), (21, new[]{0,1,3}), (22, new[]{0,2,3}), (23, new[]{1,2,3}),
            (24, new[]{0,1,4}), (25, new[]{0,2,4}), (26, new[]{1,2,4}),
            (30, new[]{0,1,2,3}), (31, new[]{0,1,2,4}),
        },
    };

    // Two disjoint filled triangles (p=2) — share no vertex, so decomposition yields two (p−1)-connected pieces;
    // churn exercises the per-piece union under non-monotone filtrations.
    static Complex TwoTriangles() => new()
    {
        P = 2,
        Coords = new()
        {
            { 0, new[] { 0.0, 0 } }, { 1, new[] { 1.0, 0 } }, { 2, new[] { 0.0, 1 } },
            { 7, new[] { 10.0, 0 } }, { 8, new[] { 11.0, 0 } }, { 9, new[] { 10.0, 1 } },
        },
        Cells = new (int, int[])[]
        {
            (0, new[]{0}), (1, new[]{1}), (2, new[]{2}),
            (3, new[]{0,1}), (4, new[]{1,2}), (5, new[]{0,2}), (6, new[]{0,1,2}),
            (7, new[]{7}), (8, new[]{8}), (9, new[]{9}),
            (10, new[]{7,8}), (11, new[]{8,9}), (12, new[]{7,9}), (13, new[]{7,8,9}),
        },
    };

    [Fact] public void HollowTetra_ChurnParityBothPaths() => AssertChurnParity(HollowTetra(), seeds: 25, churnSteps: 20);
    [Fact] public void TriangleFill_ChurnParityBothPaths() => AssertChurnParity(TriangleFill(), seeds: 25, churnSteps: 16);
    [Fact] public void Bipyramid_ChurnParityBothPaths() => AssertChurnParity(Bipyramid(), seeds: 25, churnSteps: 24);
    [Fact] public void TwoTriangles_ChurnUnionParityBothPaths() => AssertChurnParity(TwoTriangles(), seeds: 25, churnSteps: 18);

    [Fact]
    public void EndpointGuard_RejectsMalformedFiltration()
    {
        // Build a hollow tetra, fill it (tetra 14), then DELETE a triangle (10) while 14 is still present —
        // illegal (10 is a face of present 14). The duality's E-Add of triangle 10's dual edge must throw,
        // because tetra 14's dual vertex is absent from G while we re-add an edge incident to it.
        var f = new ZigzagFiltration();
        for (int v = 0; v < 4; v++) f.Add(v, new int[0]);
        f.Add(4, new[] { 0, 1 }); f.Add(5, new[] { 0, 2 }); f.Add(6, new[] { 0, 3 });
        f.Add(7, new[] { 1, 2 }); f.Add(8, new[] { 1, 3 }); f.Add(9, new[] { 2, 3 });
        f.Add(10, new[] { 4, 5, 7 }); f.Add(11, new[] { 4, 6, 8 }); f.Add(12, new[] { 5, 6, 9 }); f.Add(13, new[] { 7, 8, 9 });
        f.Add(14, new[] { 10, 11, 12, 13 });
        f.Delete(10);   // illegal: triangle 10 is a face of the present tetra 14

        var dual = new DualGraphSpec(
            p: 3,
            voidDualVertexIds: new[] { 100 },
            pSimplexDualVertex: new Dictionary<int, int> { { 14, 101 } },
            pMinus1SimplexDualEdge: new Dictionary<int, (int, int)>
            {
                { 10, (100, 101) }, { 11, (100, 101) }, { 12, (100, 101) }, { 13, (100, 101) },
            });

        Assert.Throws<InvalidOperationException>(() => CodimensionOneZigzag.Compute(f, dual));
    }

    // ── Minimal re-entry adjudication (the churn finding) ─────────────────────────────────────────────
    // B = +v0 +v1 +e(0,1) −e +e −e −v1 −v0 : the edge (cell 2) is deleted and re-added. The H0 component
    // containing v0 persists steps 0..7; the isolated-{v1} class is born and dies three separate times — it
    // exists only while v0,v1 are disconnected (after +v1, after each −e) and is killed at each merge. So the
    // true H0 barcode (engine convention birth=b−1, death=d; Closed=0, Open=1) is four bars:
    //   (0,7) C/C, (1,2) C/O, (3,4) O/O, (5,6) O/C.  Z1 (ZigzagBarcodeNaive) matches; Z2 (FastZigzag)
    //   previously merged (1,2)+(3,4)→(1,4) and dropped a bar — now fixed via per-incarnation accounting.
    static ZigzagFiltration MinimalH0Reentry()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]);
        f.Add(2, new[] { 0, 1 }); f.Delete(2); f.Add(2, new[] { 0, 1 }); f.Delete(2);
        f.Delete(1); f.Delete(0);
        return f;
    }

    static readonly (double, double, int, int, int)[] TrueH0Reentry =
    {
        (0, 7, 0, 0, 0), (1, 2, 0, 0, 1), (3, 4, 0, 1, 1), (5, 6, 0, 1, 0),
    };

    [Fact]
    public void MinimalReentry_NaiveMatchesHandDerivedTruth()
    {
        // Z1 (rank-definitional) is correct on cell re-entry: it reproduces the hand-derived true barcode.
        Assert.Equal(TrueH0Reentry, Sig(ZigzagBarcodeNaive.Compute(MinimalH0Reentry(), 0), 0).ToArray());
    }

    [Fact]
    public void MinimalReentry_FastZigzagShouldMatchTruth()
    {
        Assert.Equal(TrueH0Reentry, Sig(FastZigzag.Compute(MinimalH0Reentry(), 0), 0).ToArray());
    }
}
