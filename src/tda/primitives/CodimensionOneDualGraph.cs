#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

using Graphs.Primitives;
using TDA.Ph;
namespace TDA.Primitives;

/// <summary>
/// General p≥3 codimension-one dual-graph builder — Dey–Hou–Mandal §4.1 <c>VoidBoundary</c>
/// (<c>1907.04889</c>) for an ℝᵖ-embedded simplicial complex. Reconstructs the voids of ℝᵖ∖|K| and
/// assembles the <see cref="DualGraphSpec"/> the pure <see cref="CodimensionOneZigzag"/> engine consumes
/// (the A2 seam: geometry in, dual out, duality unchanged downstream). The dimension-general companion of
/// <see cref="PlanarDualGraph"/> (its d=1 specialization — the transverse-2-plane pairing collapses to the
/// planar angular sort around a vertex).
///
/// <para><b>The algorithm (§4.1), per (p−2)-simplex σ.</b> The (p−1)-simplices τ ⊃ σ form a fan in the
/// transverse 2-plane Δ ⟂ σ; sort them circularly. Each incident p-simplex σ_p ⊃ σ fills the wedge between
/// its two σ-containing (p−1)-faces. The unmarked gaps are <i>void wedges</i>: the two sides facing each void
/// wedge bound the same local void and are paired. A union-find over <b>half-(p−1)-simplices</b> (a (p−1)-simplex
/// + a global co-orientation sign) groups paired sides into voids — one class per void, since (p−1)-connectedness
/// (the caller's responsibility, Algorithm 3 item C) forbids voids with disconnected boundaries (Thm 4.1).</para>
///
/// <para><b>Co-orientation is global, not local (the load-bearing subtlety).</b> A half-simplex's sign is the
/// sign of the natural-orientation determinant <c>det(v₁−v₀,…,v_{p−1}−v₀, P−v₀)</c> for a point <c>P</c> in the
/// void wedge (Def 4.1) — computed in real ℝᵖ coordinates, so it names the same physical side of τ no matter
/// which of τ's (p−2)-faces is being processed (each face's transverse-plane basis has arbitrary handedness;
/// the determinant does not). This is what keeps the union-find consistent across faces.</para>
///
/// <para><b>Self-loop drop.</b> A naked (p−1)-simplex whose two sides bound the <i>same</i> void (a
/// non-separating flap — Prop 4.1's "two sides bound different voids" applies only to a (p−1)-simplex <i>in a
/// cycle</i>, not to every naked one) duals to a self-loop — topologically inert in H0. We drop it at assembly:
/// barcode-neutral, and it keeps the pure engine free of inert edges (≡ the optional §4.1 Prune, done dual-side).</para>
///
/// <para>Assumes K is a single (p−1)-connected piece in general position (distinct transverse angles, non-degenerate
/// orientation determinants). The (p−1)-connected decomposition + per-component union (item C/F) is the
/// orchestration layer above this per-piece builder. Degeneracies throw rather than silently mis-sort.</para>
/// </summary>
public static class CodimensionOneDualGraph
{
    const double Eps = 1e-9;

    /// <param name="p">Embedding dimension (≥ 2); the dual is for H_{p−1}.</param>
    /// <param name="vertexCoords">Vertex cell-id → its ℝᵖ coordinates (length p).</param>
    /// <param name="pMinus1Simplices">(p−1)-simplex cell-id → its p vertex ids (the dual edges).</param>
    /// <param name="pSimplices">p-simplex cell-id → its p+1 vertex ids (the filled cells).</param>
    /// <param name="firstDualVertexId">Allocation base for fresh dual-vertex ids (p-simplex duals + voids);
    /// only needs to keep those ids mutually distinct.</param>
    public static DualGraphSpec Build(
        int p,
        IReadOnlyDictionary<int, double[]> vertexCoords,
        IReadOnlyDictionary<int, IReadOnlyList<int>> pMinus1Simplices,
        IReadOnlyDictionary<int, IReadOnlyList<int>> pSimplices,
        int firstDualVertexId = 0)
    {
        if (p < 2) throw new ArgumentOutOfRangeException(nameof(p), "Codimension-one duality requires p ≥ 2.");
        ArgumentNullException.ThrowIfNull(vertexCoords);
        ArgumentNullException.ThrowIfNull(pMinus1Simplices);
        ArgumentNullException.ThrowIfNull(pSimplices);

        // --- dense indexing + canonical vertex sets ---
        var t1Keys = pMinus1Simplices.Keys.OrderBy(x => x).ToList();      // (p−1)-simplex cell-ids
        var t1Idx = new Dictionary<int, int>();
        for (int i = 0; i < t1Keys.Count; i++) t1Idx[t1Keys[i]] = i;
        int M = t1Keys.Count;

        var t1Verts = new Dictionary<int, int[]>();                       // cell-id -> sorted vertex ids
        var t1ByKey = new Dictionary<string, int>();                      // vertex-set key -> (p−1) cell-id
        foreach (int c in t1Keys)
        {
            var v = pMinus1Simplices[c].OrderBy(x => x).ToArray();
            if (v.Length != p) throw new ArgumentException($"(p−1)-simplex {c} must have p={p} vertices, has {v.Length}.");
            t1Verts[c] = v;
            t1ByKey[Key(v)] = c;
        }

        var pKeys = pSimplices.Keys.OrderBy(x => x).ToList();
        var pVerts = new Dictionary<int, int[]>();
        foreach (int c in pKeys)
        {
            var v = pSimplices[c].OrderBy(x => x).ToArray();
            if (v.Length != p + 1) throw new ArgumentException($"p-simplex {c} must have p+1={p + 1} vertices, has {v.Length}.");
            pVerts[c] = v;
        }

        // --- p-cofaces of each (p−1)-simplex: a p-simplex contributes to each of its p+1 (p−1)-faces ---
        var cofaces = new Dictionary<int, List<int>>();
        foreach (int c in t1Keys) cofaces[c] = new List<int>();
        foreach (int pc in pKeys)
        {
            var pv = pVerts[pc];
            for (int drop = 0; drop < pv.Length; drop++)
            {
                var face = Omit(pv, drop);
                if (t1ByKey.TryGetValue(Key(face), out int fc)) cofaces[fc].Add(pc);
            }
        }
        foreach (int c in t1Keys)
            if (cofaces[c].Count > 2)
                throw new ArgumentException($"(p−1)-simplex {c} has {cofaces[c].Count} p-cofaces (> 2): not a valid codim-one embedding.");

        // --- group (p−1)-simplices by shared (p−2)-face (drop one vertex of each τ) ---
        // key -> (faceVerts, list of (t1 cell-id, extra vertex of τ not in σ))
        var groups = new Dictionary<string, (int[] Face, List<(int Cell, int Extra)> Incident)>();
        foreach (int c in t1Keys)
        {
            var v = t1Verts[c];
            for (int drop = 0; drop < v.Length; drop++)
            {
                var face = Omit(v, drop);                                 // (p−2)-face, p−1 vertices
                string k = Key(face);
                if (!groups.TryGetValue(k, out var g)) { g = (face, new List<(int, int)>()); groups[k] = g; }
                g.Incident.Add((c, v[drop]));
            }
        }

        // --- void-boundary reconstruction: union-find over half-(p−1)-simplices (2 sides each) ---
        var uf = new UnionFind(2 * M);
        var appeared = new HashSet<int>();                               // half-ids that face some void wedge
        int Half(int cell, int sign) => 2 * t1Idx[cell] + (sign > 0 ? 0 : 1);

        foreach (var g in groups.Values)
        {
            int[] sigma = g.Face;
            var (basisA, e1, e2) = TransversePlane(sigma, vertexCoords, p);
            double[] cSigma = Centroid(sigma, vertexCoords, p);
            double[] w0 = vertexCoords[sigma[0]];

            // circular order of incident (p−1)-simplices in Δ
            int n = g.Incident.Count;
            var ray = new (int Cell, double Cx, double Cy, double Ang)[n];
            for (int i = 0; i < n; i++)
            {
                var (cell, extra) = g.Incident[i];
                var (cx, cy) = ProjectDir(Sub(vertexCoords[extra], w0, p), basisA, e1, e2, p);
                ray[i] = (cell, cx, cy, Norm2(Math.Atan2(cy, cx)));
            }
            Array.Sort(ray, (a, b) => a.Ang.CompareTo(b.Ang));
            var posOf = new Dictionary<int, int>();
            for (int i = 0; i < n; i++) posOf[ray[i].Cell] = i;

            // mark filled gaps: each incident p-simplex fills the wedge between its two σ-faces
            var filled = new bool[n];
            foreach (int pc in CofacesContaining(sigma, g.Incident, cofaces))
            {
                var extras = pVerts[pc].Where(x => Array.BinarySearch(sigma, x) < 0).ToArray();   // 2 vertices ∉ σ
                if (extras.Length != 2) throw new InvalidOperationException($"p-simplex {pc} does not contain (p−2)-face {Key(sigma)} cleanly.");
                int fa = posOf[t1ByKey[Key(Insert(sigma, extras[0]))]];   // ray of σ∪{a}
                int fb = posOf[t1ByKey[Key(Insert(sigma, extras[1]))]];   // ray of σ∪{b}
                // bisector of the two face-rays points into the p-simplex's (filled) cone
                var (sx, sy) = Bisector(ray[fa].Cx, ray[fa].Cy, ray[fb].Cx, ray[fb].Cy);
                filled[FindGap(ray, n, Norm2(Math.Atan2(sy, sx)))] = true;
            }

            // void wedges = unfilled gaps; pair the two sides facing each
            for (int gi = 0; gi < n; gi++)
            {
                if (filled[gi]) continue;
                int li = gi, ri = (gi + 1) % n;
                int tauL = ray[li].Cell, tauR = ray[ri].Cell;
                double thL = ray[li].Ang;
                double width = (gi < n - 1) ? ray[ri].Ang - ray[li].Ang : (2 * Math.PI - ray[n - 1].Ang + ray[0].Ang);
                double thR = thL + width;
                double delta = Math.Min(width, Math.PI) * 0.5;

                double[] pL = Add(cSigma, Dir(thL + delta, e1, e2, p), p);
                double[] pR = Add(cSigma, Dir(thR - delta, e1, e2, p), p);
                int sL = OrientationSign(t1Verts[tauL], pL, vertexCoords, p);
                int sR = OrientationSign(t1Verts[tauR], pR, vertexCoords, p);

                int hL = Half(tauL, sL), hR = Half(tauR, sR);
                appeared.Add(hL); appeared.Add(hR);
                uf.Union(hL, hR);
            }
        }

        // --- assembly: allocate dual-vertex ids (p-simplices first, then void classes), then edges ---
        int nextDual = firstDualVertexId;
        var pDual = new Dictionary<int, int>();
        foreach (int pc in pKeys) pDual[pc] = nextDual++;

        var rootVoid = new Dictionary<int, int>();
        var voidIds = new List<int>();
        int VoidOf(int half)
        {
            int r = uf.Find(half);
            if (!rootVoid.TryGetValue(r, out int id)) { id = nextDual++; rootVoid[r] = id; voidIds.Add(id); }
            return id;
        }
        foreach (int h in appeared.OrderBy(x => x)) VoidOf(h);          // deterministic void-id allocation

        var edgeDual = new Dictionary<int, (int A, int B)>();
        foreach (int c in t1Keys)
        {
            var cof = cofaces[c];
            if (cof.Count == 2)
            {
                edgeDual[c] = (pDual[cof[0]], pDual[cof[1]]);            // interior: two p-cofaces
            }
            else if (cof.Count == 1)
            {
                bool plus = appeared.Contains(Half(c, +1)), minus = appeared.Contains(Half(c, -1));
                if (plus == minus)
                    throw new InvalidOperationException($"Boundary (p−1)-simplex {c} must face exactly one void; faced {(plus ? 2 : 0)}.");
                edgeDual[c] = (pDual[cof[0]], VoidOf(Half(c, plus ? +1 : -1)));
            }
            else
            {
                // naked: both sides void-facing; drop if the same void (inert self-loop)
                int vP = VoidOf(Half(c, +1)), vM = VoidOf(Half(c, -1));
                if (vP == vM) continue;                                 // self-loop drop
                edgeDual[c] = (vP, vM);
            }
        }

        return new DualGraphSpec(p, voidIds, pDual, edgeDual);
    }

    // ---- p-simplices that contain a given (p−2)-face σ (dedup of incident (p−1)-simplices' cofaces) ----
    static IEnumerable<int> CofacesContaining(int[] sigma, List<(int Cell, int Extra)> incident,
        Dictionary<int, List<int>> cofaces)
    {
        var seen = new HashSet<int>();
        foreach (var (cell, _) in incident)
            foreach (int pc in cofaces[cell])
                if (seen.Add(pc)) yield return pc;
    }

    // ---- geometry primitives (small p; self-contained double[] ops) ----

    // Transverse 2-plane Δ ⟂ σ: Gram–Schmidt σ's edge vectors → B_A (dim p−2), complete to ℝᵖ → leftover (e1,e2).
    static (List<double[]> BasisA, double[] E1, double[] E2) TransversePlane(
        int[] sigma, IReadOnlyDictionary<int, double[]> coords, int p)
    {
        var basisA = new List<double[]>();
        double[] w0 = coords[sigma[0]];
        for (int i = 1; i < sigma.Length; i++)
        {
            double[] u = Sub(coords[sigma[i]], w0, p);
            foreach (var b in basisA) AxpyInto(u, -Dot(u, b, p), b, p);
            double nrm = Math.Sqrt(Dot(u, u, p));
            if (nrm > Eps) { Scale(u, 1.0 / nrm, p); basisA.Add(u); }
        }
        if (basisA.Count != p - 2)
            throw new InvalidOperationException($"(p−2)-face {Key(sigma)} is degenerate: transverse span has rank {basisA.Count}, expected {p - 2}.");

        var extra = new List<double[]>();
        for (int k = 0; k < p && extra.Count < 2; k++)
        {
            double[] u = new double[p]; u[k] = 1.0;
            foreach (var b in basisA) AxpyInto(u, -Dot(u, b, p), b, p);
            foreach (var b in extra) AxpyInto(u, -Dot(u, b, p), b, p);
            double nrm = Math.Sqrt(Dot(u, u, p));
            if (nrm > Eps) { Scale(u, 1.0 / nrm, p); extra.Add(u); }
        }
        if (extra.Count != 2)
            throw new InvalidOperationException($"Could not complete transverse plane for (p−2)-face {Key(sigma)}.");
        return (basisA, extra[0], extra[1]);
    }

    // Project a direction into Δ and return its (e1,e2) coordinates.
    static (double Cx, double Cy) ProjectDir(double[] v, List<double[]> basisA, double[] e1, double[] e2, int p)
    {
        double[] d = (double[])v.Clone();
        foreach (var b in basisA) AxpyInto(d, -Dot(d, b, p), b, p);
        return (Dot(d, e1, p), Dot(d, e2, p));
    }

    // Sign of det(v₁−v₀,…,v_{p−1}−v₀, P−v₀): which global side of τ's hyperplane P is on (Def 4.1).
    static int OrientationSign(int[] tau, double[] pPoint, IReadOnlyDictionary<int, double[]> coords, int p)
    {
        double[] v0 = coords[tau[0]];
        var m = new double[p][];
        for (int i = 1; i < p; i++) m[i - 1] = Sub(coords[tau[i]], v0, p);
        m[p - 1] = Sub(pPoint, v0, p);
        double det = Determinant(m, p);
        if (Math.Abs(det) < Eps)
            throw new InvalidOperationException($"Degenerate orientation for (p−1)-simplex [{string.Join(",", tau)}] (det≈0): not in general position.");
        return det > 0 ? +1 : -1;
    }

    static double Determinant(double[][] a, int n)
    {
        var m = new double[n][];
        for (int i = 0; i < n; i++) m[i] = (double[])a[i].Clone();
        double det = 1.0;
        for (int col = 0; col < n; col++)
        {
            int piv = col;
            for (int r = col + 1; r < n; r++) if (Math.Abs(m[r][col]) > Math.Abs(m[piv][col])) piv = r;
            if (Math.Abs(m[piv][col]) < 1e-300) return 0.0;
            if (piv != col) { (m[piv], m[col]) = (m[col], m[piv]); det = -det; }
            det *= m[col][col];
            for (int r = col + 1; r < n; r++)
            {
                double f = m[r][col] / m[col][col];
                for (int k = col; k < n; k++) m[r][k] -= f * m[col][k];
            }
        }
        return det;
    }

    // Gap (consecutive-ray interval) of the circular order containing angle θ ∈ [0,2π).
    static int FindGap((int Cell, double Cx, double Cy, double Ang)[] ray, int n, double theta)
    {
        for (int i = 0; i < n - 1; i++)
            if (ray[i].Ang <= theta && theta < ray[i + 1].Ang) return i;
        return n - 1;   // wrap gap: θ ≥ last or θ < first
    }

    static double[] Dir(double theta, double[] e1, double[] e2, int p)
    {
        double c = Math.Cos(theta), s = Math.Sin(theta);
        var r = new double[p];
        for (int k = 0; k < p; k++) r[k] = c * e1[k] + s * e2[k];
        return r;
    }

    static double[] Centroid(int[] verts, IReadOnlyDictionary<int, double[]> coords, int p)
    {
        var r = new double[p];
        foreach (int v in verts) { var c = coords[v]; for (int k = 0; k < p; k++) r[k] += c[k]; }
        for (int k = 0; k < p; k++) r[k] /= verts.Length;
        return r;
    }

    // sum of the two unit vectors (ax,ay),(bx,by) — the bisector direction of the cone they span
    static (double Sx, double Sy) Bisector(double ax, double ay, double bx, double by)
    {
        double na = Math.Sqrt(ax * ax + ay * ay), nb = Math.Sqrt(bx * bx + by * by);
        return (ax / na + bx / nb, ay / na + by / nb);
    }

    static double Norm2(double a) { a %= 2 * Math.PI; return a < 0 ? a + 2 * Math.PI : a; }

    static double[] Sub(double[] a, double[] b, int p) { var r = new double[p]; for (int k = 0; k < p; k++) r[k] = a[k] - b[k]; return r; }
    static double[] Add(double[] a, double[] b, int p) { var r = new double[p]; for (int k = 0; k < p; k++) r[k] = a[k] + b[k]; return r; }
    static double Dot(double[] a, double[] b, int p) { double s = 0; for (int k = 0; k < p; k++) s += a[k] * b[k]; return s; }
    static void Scale(double[] a, double f, int p) { for (int k = 0; k < p; k++) a[k] *= f; }
    static void AxpyInto(double[] y, double f, double[] x, int p) { for (int k = 0; k < p; k++) y[k] += f * x[k]; }

    static int[] Omit(int[] v, int drop)
    {
        var r = new int[v.Length - 1];
        for (int i = 0, j = 0; i < v.Length; i++) if (i != drop) r[j++] = v[i];
        return r;
    }

    static int[] Insert(int[] sorted, int x)
    {
        var r = new int[sorted.Length + 1];
        int i = 0; while (i < sorted.Length && sorted[i] < x) { r[i] = sorted[i]; i++; }
        r[i] = x; for (int j = i; j < sorted.Length; j++) r[j + 1] = sorted[j];
        return r;
    }

    static string Key(int[] sortedVerts) => string.Join(",", sortedVerts);
}
