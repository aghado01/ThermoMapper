#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Strong collapse (Boissonnat–Pritam, arXiv:1809.10945): iteratively delete dominated vertices to
/// reach the <em>core</em> — a minimal complex with the same strong homotopy type, hence the same
/// homology. By Remark 1, a vertex <c>v</c> is dominated iff every maximal simplex containing <c>v</c>
/// also contains some other vertex <c>v'</c> (i.e. the intersection of all maximal simplices
/// containing <c>v</c> has a vertex besides <c>v</c>).
/// <para>This is the slow, obviously-correct ground truth; the nerve-based fast algorithm (§3
/// Algorithm 1) and the §4 sequence/persistence version are perf follow-ons validated against it.</para>
/// </summary>
public static class StrongCollapse
{
    /// <summary>Collapse a complex (given by its maximal simplices) to its core's maximal simplices.</summary>
    public static int[][] Core(IReadOnlyList<int[]> maximalSimplices) =>
        CoreWithRetraction(maximalSimplices).Core;

    /// <summary>
    /// Collapse to the core, also returning the retraction map r: every original vertex → its image
    /// in the core (core vertices map to themselves; a deleted vertex maps to a dominator, chased to
    /// the core). This is the r_j needed for §4 core-assembly of a sequence.
    /// </summary>
    public static (int[][] Core, IReadOnlyDictionary<int, int> Retraction) CoreWithRetraction(
        IReadOnlyList<int[]> maximalSimplices)
    {
        ArgumentNullException.ThrowIfNull(maximalSimplices);
        var maximals = Maximalize(maximalSimplices.Select(s => new HashSet<int>(s)));

        var allVertices = new HashSet<int>();
        foreach (var m in maximals) allVertices.UnionWith(m);

        // dom[v] = a vertex that dominated v at the moment v was deleted (resolved to the core below).
        var dom = new Dictionary<int, int>();

        bool changed = true;
        while (changed)
        {
            changed = false;
            var vertices = new HashSet<int>();
            foreach (var m in maximals) vertices.UnionWith(m);

            foreach (int v in vertices)
            {
                // Intersection of every maximal simplex containing v; v is dominated iff something
                // other than v survives that intersection (Remark 1).
                HashSet<int>? inter = null;
                foreach (var m in maximals)
                    if (m.Contains(v))
                    {
                        if (inter is null) inter = new HashSet<int>(m);
                        else inter.IntersectWith(m);
                    }
                if (inter is null) continue;
                inter.Remove(v);
                if (inter.Count == 0) continue; // v not dominated

                dom[v] = inter.Min(); // deterministic dominator choice

                // Delete v: each maximal simplex containing v drops to its v-omitting face.
                var next = new List<HashSet<int>>();
                foreach (var m in maximals)
                {
                    if (!m.Contains(v)) { next.Add(m); continue; }
                    var face = new HashSet<int>(m);
                    face.Remove(v);
                    if (face.Count > 0) next.Add(face);
                }
                maximals = Maximalize(next);
                changed = true;
                break; // re-scan from scratch after a deletion
            }
        }

        var core = maximals
            .Select(m => { var a = m.ToArray(); Array.Sort(a); return a; })
            .OrderBy(a => a.Length).ThenBy(a => a[0])
            .ToArray();

        // Resolve each vertex to its core image by chasing dominators.
        var retraction = new Dictionary<int, int>(allVertices.Count);
        foreach (int v in allVertices)
        {
            int x = v;
            while (dom.TryGetValue(x, out int next)) x = next;
            retraction[v] = x;
        }

        return (core, retraction);
    }

    /// <summary>Keep only inclusion-maximal sets (dedups equal sets).</summary>
    static List<HashSet<int>> Maximalize(IEnumerable<HashSet<int>> sets)
    {
        var sorted = sets.OrderByDescending(s => s.Count).ToList();
        var result = new List<HashSet<int>>();
        foreach (var s in sorted)
            if (!result.Any(r => s.IsSubsetOf(r)))
                result.Add(s);
        return result;
    }
}
