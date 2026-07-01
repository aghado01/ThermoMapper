#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// §4 of 1809.10945 — zigzag persistence via strong-collapse core-assembly. Each complex in the
/// zigzag is collapsed to its core; the induced maps become retractions (f^c_j = r_{j+1}∘f_j∘i_j),
/// and by Theorem 2 the barcode is unchanged. The smaller core sequence is evaluated by the general
/// simplicial-zigzag oracle <see cref="ZigzagMapBarcode"/>.
/// <para>Requires a SIMPLICIAL zigzag: a cell's vertices are recovered from its transitive dim-0
/// boundary (vertex labels = dim-0 cell ids).</para>
/// </summary>
public static class StrongCollapseZigzag
{
    public static Barcode Compute(ZigzagFiltration f, int maxDimension = int.MaxValue)
    {
        ArgumentNullException.ThrowIfNull(f);
        var (complexes, dimOf, bndOf, n, steps) = BuildCoreSequence(f);
        return ZigzagMapBarcode.Compute(complexes, dimOf, bndOf, n, steps, maxDimension);
    }

    static (List<HashSet<int>> Complexes, int[] DimOf, int[][] BndOf, int N, List<ZigzagMapBarcode.Step> Steps)
        BuildCoreSequence(ZigzagFiltration f)
    {
        int numSteps = f.Count;

        // --- 1. Materialize the zigzag complexes; recover each cell's vertex set ---
        int maxCellId = -1;
        foreach (var s in f) if (s.GlobalCellId > maxCellId) maxCellId = s.GlobalCellId;
        int U = maxCellId + 1;
        var zdim = new int[U];
        var zbnd = new int[U][];
        for (int i = 0; i < U; i++) zbnd[i] = Array.Empty<int>();

        var zComplexes = new List<HashSet<int>> { new HashSet<int>() };
        var cur = new HashSet<int>();
        for (int i = 0; i < numSteps; i++)
        {
            var s = f[i];
            if (s.Direction == ZigzagDirection.Add)
            {
                cur.Add(s.GlobalCellId);
                var b = s.BoundaryAtAdd ?? Array.Empty<int>();
                zdim[s.GlobalCellId] = b.Length > 0 ? zdim[b[0]] + 1 : 0;
                zbnd[s.GlobalCellId] = b;
            }
            else cur.Remove(s.GlobalCellId);
            zComplexes.Add(new HashSet<int>(cur));
        }

        var vmemo = new Dictionary<int, int[]>();
        int[] Verts(int cell)
        {
            if (vmemo.TryGetValue(cell, out var v)) return v;
            int[] r;
            if (zdim[cell] == 0) r = new[] { cell };
            else
            {
                var set = new SortedSet<int>();
                foreach (int face in zbnd[cell]) set.UnionWith(Verts(face));
                r = set.ToArray();
            }
            vmemo[cell] = r;
            return r;
        }

        // --- 2. Collapse each complex: core maximal simplices (vertex tuples) + retraction r_i ---
        int m = numSteps; // complexes indexed 0..m
        var coreMaximal = new List<int[]>[m + 1];
        var retraction = new IReadOnlyDictionary<int, int>[m + 1];
        for (int i = 0; i <= m; i++)
        {
            var complex = zComplexes[i];
            var nonMaximal = new HashSet<int>();
            foreach (int c in complex)
                foreach (int face in zbnd[c])
                    nonMaximal.Add(face);

            var maxVerts = new List<int[]>();
            foreach (int c in complex)
                if (!nonMaximal.Contains(c))
                    maxVerts.Add(Verts(c));

            var (core, retr) = StrongCollapse.CoreWithRetraction(maxVerts);
            coreMaximal[i] = core.ToList();
            retraction[i] = retr;
        }

        // --- 3. Global-id all core simplices (by vertex tuple) -> cells with dim + boundary ---
        var idOf = new Dictionary<string, int>();
        var vertsList = new List<int[]>();
        int Id(int[] verts)
        {
            string key = string.Join(",", verts);
            if (idOf.TryGetValue(key, out int id)) return id;
            id = vertsList.Count;
            idOf[key] = id;
            vertsList.Add(verts);
            return id;
        }

        var coreComplexes = new List<HashSet<int>>();
        for (int i = 0; i <= m; i++)
        {
            var cells = new HashSet<int>();
            foreach (var mx in coreMaximal[i])
                foreach (var face in Faces(mx))
                    cells.Add(Id(face));
            coreComplexes.Add(cells);
        }

        int N = vertsList.Count;
        var dimOf = new int[N];
        var bndOf = new int[N][];
        for (int c = 0; c < N; c++)
        {
            var verts = vertsList[c];
            dimOf[c] = verts.Length - 1;
            if (verts.Length <= 1) { bndOf[c] = Array.Empty<int>(); continue; }
            var b = new int[verts.Length];
            for (int k = 0; k < verts.Length; k++) b[k] = Id(Omit(verts, k));
            bndOf[c] = b;
        }

        // --- 4. Steps: cell maps induced by the retractions (f^c_j = r on vertices) ---
        var steps = new List<ZigzagMapBarcode.Step>(numSteps);
        for (int i = 0; i < numSteps; i++)
        {
            bool fwd = f[i].Direction == ZigzagDirection.Add;
            int srcIdx = fwd ? i : i + 1;                 // source complex of the step's arrow
            var r = fwd ? retraction[i + 1] : retraction[i];

            var cellMap = new int[N];
            for (int c = 0; c < N; c++) cellMap[c] = -1;
            foreach (int c in coreComplexes[srcIdx])
            {
                var verts = vertsList[c];
                var image = new SortedSet<int>();
                foreach (int v in verts) image.Add(r[v]);
                if (image.Count < verts.Length) continue;  // r collapses vertices -> maps to 0
                cellMap[c] = Id(image.ToArray());
            }
            steps.Add(new ZigzagMapBarcode.Step(fwd, cellMap));
        }

        return (coreComplexes, dimOf, bndOf, N, steps);
    }

    static IEnumerable<int[]> Faces(int[] simplex)
    {
        int n = simplex.Length;
        for (int mask = 1; mask < (1 << n); mask++)
        {
            var verts = new List<int>();
            for (int b = 0; b < n; b++) if ((mask & (1 << b)) != 0) verts.Add(simplex[b]);
            yield return verts.ToArray(); // simplex is sorted, so the subset is sorted
        }
    }

    static int[] Omit(int[] verts, int k)
    {
        var r = new int[verts.Length - 1];
        for (int i = 0, w = 0; i < verts.Length; i++) if (i != k) r[w++] = verts[i];
        return r;
    }
}
