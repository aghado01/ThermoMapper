#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

namespace TDA.Ph;

/// <summary>
/// p=2 planar specialization of Dey–Hou–Mandal §4.1 void-boundary reconstruction (`1907.04889`): builds the
/// codimension-one dual graph of a planar embedded complex K ⊂ ℝ² by tracing the faces of K's straight-line
/// 1-skeleton (DCEL-style angular sort + face walk). Produces the combinatorial <see cref="DualGraphSpec"/>
/// the pure <see cref="CodimensionOneZigzag"/> engine consumes — the **A2 seam**: geometry in, dual out,
/// duality unchanged downstream (A1's caller-supplied dual becomes A2's computed dual).
///
/// <para>Voids = faces of ℝ²∖|K|: every bounded face of the edge arrangement that is NOT a filled triangle is
/// a distinct void (adjacent faces are separated by K-edges), plus the one unbounded outer void. A 2-simplex
/// present in K fills its triangular face → that face is the triangle's dual vertex, not a void. Each K-edge's
/// two half-edges border its two faces; the dual edge joins those two cells. For d=1, §4.1's "transverse-plane
/// pairing around a (d−1)-simplex" collapses to angular sort of half-edges around each vertex.</para>
///
/// <para>Assumes K's 1-skeleton is connected and in general position (no degenerate equal angles). The general
/// p≥3 reconstruction (natural orientation + transverse-plane pairing) and the (p−1)-connected decomposition
/// are the deferred §4.1 work; this is the planar case that closes the geometry→dual→duality loop end-to-end.</para>
/// </summary>
public static class PlanarDualGraph
{
    /// <param name="firstDualVertexId">Allocation base for fresh dual-vertex ids (voids + 2-simplex duals);
    /// only needs to keep those ids mutually distinct.</param>
    public static DualGraphSpec Build(
        IReadOnlyDictionary<int, (double X, double Y)> vertexCoords,
        IReadOnlyDictionary<int, (int U, int V)> edges,
        IReadOnlyDictionary<int, IReadOnlyList<int>> triangleVertices,
        int firstDualVertexId = 0)
    {
        ArgumentNullException.ThrowIfNull(vertexCoords);
        ArgumentNullException.ThrowIfNull(edges);
        ArgumentNullException.ThrowIfNull(triangleVertices);

        // --- half-edges: each K-edge (u,v) -> directed u->v and v->u, both tagged with the edge cell-id. ---
        var heFrom = new List<int>();
        var heTo = new List<int>();
        var heIndex = new Dictionary<(int, int), int>();
        void AddHalf(int a, int b) { heIndex[(a, b)] = heFrom.Count; heFrom.Add(a); heTo.Add(b); }
        foreach (var kv in edges) { var (u, v) = kv.Value; AddHalf(u, v); AddHalf(v, u); }

        // --- at each vertex, CCW angular order of outgoing half-edges (by atan2 of the edge direction). ---
        var ccwAt = new Dictionary<int, List<int>>();
        for (int he = 0; he < heFrom.Count; he++)
        {
            if (!ccwAt.TryGetValue(heFrom[he], out var lst)) { lst = new List<int>(); ccwAt[heFrom[he]] = lst; }
            lst.Add(he);
        }
        double Angle(int he)
        {
            var (sx, sy) = vertexCoords[heFrom[he]];
            var (tx, ty) = vertexCoords[heTo[he]];
            return Math.Atan2(ty - sy, tx - sx);
        }
        var posInCcw = new int[heFrom.Count];
        foreach (var lst in ccwAt.Values)
        {
            lst.Sort((a, b) => Angle(a).CompareTo(Angle(b)));
            for (int i = 0; i < lst.Count; i++) posInCcw[lst[i]] = i;
        }

        // --- next(u->v) = the CCW-predecessor of v->u at v (= CW-next): traces faces with the interior on the left. ---
        int Next(int he)
        {
            int back = heIndex[(heTo[he], heFrom[he])];     // v->u
            var lst = ccwAt[heTo[he]];
            return lst[(posInCcw[back] - 1 + lst.Count) % lst.Count];
        }

        // --- trace faces (each half-edge belongs to exactly one face). ---
        var faceOf = new int[heFrom.Count];
        Array.Fill(faceOf, -1);
        var faces = new List<List<int>>();
        for (int he = 0; he < heFrom.Count; he++)
        {
            if (faceOf[he] != -1) continue;
            int fid = faces.Count;
            var cycle = new List<int>();
            for (int cur = he; faceOf[cur] == -1; cur = Next(cur)) { faceOf[cur] = fid; cycle.Add(cur); }
            faces.Add(cycle);
        }

        // --- classify faces: a bounded (CCW, +area) triangular face whose vertex set is a 2-simplex in K is
        //     filled (= that triangle's dual vertex); every other face is a distinct void (incl. the outer face). ---
        var triByVerts = new Dictionary<string, int>();
        foreach (var kv in triangleVertices) triByVerts[VertKey(kv.Value)] = kv.Key;

        int nextDual = firstDualVertexId;
        var triDual = new Dictionary<int, int>();      // 2-simplex cell-id -> dual vertex id
        var voidIds = new List<int>();
        var faceCell = new int[faces.Count];           // face id -> dual vertex id
        for (int f = 0; f < faces.Count; f++)
        {
            var cyc = faces[f];
            double area2 = 0;
            foreach (int h in cyc)
            {
                var (x1, y1) = vertexCoords[heFrom[h]];
                var (x2, y2) = vertexCoords[heTo[h]];
                area2 += x1 * y2 - x2 * y1;
            }
            int matched = -1;
            if (area2 > 0 && cyc.Count == 3 &&
                triByVerts.TryGetValue(VertKey(cyc.Select(h => heFrom[h]).ToList()), out int tcell))
                matched = tcell;

            if (matched >= 0)
            {
                if (!triDual.TryGetValue(matched, out int dv)) { dv = nextDual++; triDual[matched] = dv; }
                faceCell[f] = dv;
            }
            else { int dv = nextDual++; voidIds.Add(dv); faceCell[f] = dv; }
        }

        // --- each K-edge's dual joins the cells of the two faces its half-edges bound. ---
        var edgeDual = new Dictionary<int, (int A, int B)>();
        foreach (var kv in edges)
        {
            var (u, v) = kv.Value;
            edgeDual[kv.Key] = (faceCell[faceOf[heIndex[(u, v)]]], faceCell[faceOf[heIndex[(v, u)]]]);
        }

        return new DualGraphSpec(2, voidIds, triDual, edgeDual);
    }

    static string VertKey(IReadOnlyList<int> verts)
    {
        var s = verts.ToList();
        s.Sort();
        return string.Join(",", s);
    }
}
