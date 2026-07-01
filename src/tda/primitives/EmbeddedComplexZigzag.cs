#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

using Graphs.Primitives;
using Maths.Topology;
using TDA.Ph;
namespace TDA.Primitives;

/// <summary>
/// Z5d item C/F — the full Dey–Hou §5 Algorithm 3 front door for an <b>arbitrary</b> ℝᵖ-embedded simplicial
/// complex (no (p−1)-connectedness assumption). Lifts the per-piece pipeline
/// (<see cref="CodimensionOneDualGraph"/> geometry → <see cref="DualGraphSpec"/> → the pure
/// <see cref="CodimensionOneZigzag"/> duality engine) to all of K by decomposing into (p−1)-connected pieces
/// and unioning their barcodes — the (p−1)-th zigzag barcode of <c>f</c>:
/// <c>Pers(H_{p−1}(F)) = ⋃ℓ Pers(H_{p−1}(Xℓ))</c> (Prop 28).
///
/// <para><b>Pipeline.</b>
/// <list type="number">
/// <item><b>Decompose</b> (Def 27): union-find the (p−1)-simplices, merging any two that share a (p−2)-face.
/// Each class is one (p−1)-connected piece — the granularity at which <see cref="CodimensionOneDualGraph"/>'s
/// local void reconstruction is provably one-class-per-void (Thm 4.1).</item>
/// <item><b>Assign p-simplices</b> to the piece of their (p−1)-faces (all p+1 faces are pairwise (p−2)-adjacent,
/// hence one piece).</item>
/// <item><b>Build + reduce</b> per piece: the restriction <c>Xℓ : K_i ∩ Cℓ</c> is left <i>implicit</i> — the full
/// <c>f</c> is fed to <see cref="CodimensionOneZigzag.Compute"/> with the piece's
/// <see cref="DualGraphSpec"/>; cells outside the piece are absent from the spec and so collapse to no-op arrows
/// (the same no-op collapse the engine already applies to dims ∉ {p−1, p}).</item>
/// <item><b>Union</b> (item F): concatenate the per-piece bars — they already share F's step axis and dimension p−1,
/// and the pieces share no (p−1)-simplices, so the homologies are independent (no double-count).</item>
/// </list></para>
///
/// <para>The geometry-coupled companion of the pure engine: this is where the embedding (coordinates + simplices
/// by dimension) is consumed; nothing in <c>TDA.Ph</c> changes. Assumes a valid embedded simplicial complex in
/// general position (the per-piece builder's degeneracy guards still apply).</para>
/// </summary>
public static class EmbeddedComplexZigzag
{
    /// <param name="p">Embedding dimension (≥ 2); computes the H_{p−1} barcode.</param>
    /// <param name="vertexCoords">Vertex cell-id → its ℝᵖ coordinates (length p).</param>
    /// <param name="pMinus1Simplices">(p−1)-simplex cell-id → its p vertex ids (the dual edges).</param>
    /// <param name="pSimplices">p-simplex cell-id → its p+1 vertex ids (the filled cells).</param>
    /// <param name="f">The zigzag filtration of all of K (∅ → … → ∅), in F's step axis.</param>
    /// <param name="algorithm">H0 reduction engine for each piece: <see cref="GraphZigzagAlgorithm.Reference"/>
    /// (default, the proven oracle) or <see cref="GraphZigzagAlgorithm.Fast"/> (near-linear); both multigraph-safe.</param>
    public static Barcode Compute(
        int p,
        IReadOnlyDictionary<int, double[]> vertexCoords,
        IReadOnlyDictionary<int, IReadOnlyList<int>> pMinus1Simplices,
        IReadOnlyDictionary<int, IReadOnlyList<int>> pSimplices,
        ZigzagFiltration f,
        GraphZigzagAlgorithm algorithm = GraphZigzagAlgorithm.Reference)
    {
        if (p < 2) throw new ArgumentOutOfRangeException(nameof(p), "Codimension-one duality requires p ≥ 2.");
        ArgumentNullException.ThrowIfNull(vertexCoords);
        ArgumentNullException.ThrowIfNull(pMinus1Simplices);
        ArgumentNullException.ThrowIfNull(pSimplices);
        ArgumentNullException.ThrowIfNull(f);

        // No (p−1)-simplices → no (p−1)-cycles → empty barcode (and UnionFind(0) would be degenerate).
        if (pMinus1Simplices.Count == 0) return new Barcode(new List<Bar>(), "Zigzag Step");

        // Dense index + canonical vertex sets of the (p−1)-simplices.
        var t1Keys = pMinus1Simplices.Keys.OrderBy(x => x).ToList();
        int M = t1Keys.Count;
        var t1Idx = new Dictionary<int, int>();
        for (int i = 0; i < M; i++) t1Idx[t1Keys[i]] = i;

        var t1Verts = new Dictionary<int, int[]>();
        var t1ByKey = new Dictionary<string, int>();
        foreach (int c in t1Keys)
        {
            var v = pMinus1Simplices[c].OrderBy(x => x).ToArray();
            if (v.Length != p) throw new ArgumentException($"(p−1)-simplex {c} must have p={p} vertices, has {v.Length}.");
            t1Verts[c] = v;
            t1ByKey[Key(v)] = c;
        }

        // (1) Decompose: union (p−1)-simplices sharing a (p−2)-face.
        var uf = new UnionFind(M);
        var firstByFace = new Dictionary<string, int>();   // (p−2)-face key -> first (p−1)-simplex dense idx
        foreach (int c in t1Keys)
        {
            var v = t1Verts[c];
            for (int drop = 0; drop < v.Length; drop++)
            {
                string fk = Key(Omit(v, drop));
                if (firstByFace.TryGetValue(fk, out int first)) uf.Union(first, t1Idx[c]);
                else firstByFace[fk] = t1Idx[c];
            }
        }

        var pieceT1 = new Dictionary<int, List<int>>();    // root -> (p−1)-simplex cell-ids
        foreach (int c in t1Keys)
        {
            int r = uf.Find(t1Idx[c]);
            if (!pieceT1.TryGetValue(r, out var lst)) { lst = new List<int>(); pieceT1[r] = lst; }
            lst.Add(c);
        }

        // (2) Assign each p-simplex to the piece of its (p−1)-faces (all in one piece by (p−2)-adjacency).
        var pieceP = new Dictionary<int, List<int>>();
        foreach (var kv in pSimplices)
        {
            var pv = kv.Value.OrderBy(x => x).ToArray();
            if (pv.Length != p + 1) throw new ArgumentException($"p-simplex {kv.Key} must have p+1={p + 1} vertices, has {pv.Length}.");
            int root = -1;
            for (int drop = 0; drop < pv.Length; drop++)
                if (t1ByKey.TryGetValue(Key(Omit(pv, drop)), out int fc)) { root = uf.Find(t1Idx[fc]); break; }
            if (root < 0)
                throw new ArgumentException($"p-simplex {kv.Key} has no (p−1)-face among the supplied (p−1)-simplices.");
            if (!pieceP.TryGetValue(root, out var lst)) { lst = new List<int>(); pieceP[root] = lst; }
            lst.Add(kv.Key);
        }

        // (3+4) Per piece: build the dual, reduce via the pure engine, union the bars.
        var allBars = new List<Bar>();
        int dualBase = 0;
        foreach (int root in pieceT1.Keys.OrderBy(x => x))
        {
            var t1Sub = pieceT1[root].ToDictionary(c => c, c => pMinus1Simplices[c]);
            var pSub = (pieceP.TryGetValue(root, out var ps) ? ps : new List<int>())
                       .ToDictionary(c => c, c => pSimplices[c]);

            // firstDualVertexId need only keep ids distinct within one spec; a generous per-piece stride is
            // belt-and-suspenders, since each piece's Compute remaps dual ids to fresh E-cell ids anyway.
            var dual = CodimensionOneDualGraph.Build(p, vertexCoords, t1Sub, pSub, firstDualVertexId: dualBase);
            dualBase += 2 * t1Sub.Count + pSub.Count + 1;

            allBars.AddRange(CodimensionOneZigzag.Compute(f, dual, algorithm).Bars);
        }

        return new Barcode(allBars, "Zigzag Step");
    }

    static string Key(int[] sortedVerts) => string.Join(",", sortedVerts);

    static int[] Omit(int[] v, int drop)
    {
        var r = new int[v.Length - 1];
        for (int i = 0, j = 0; i < v.Length; i++) if (i != drop) r[j++] = v[i];
        return r;
    }
}
