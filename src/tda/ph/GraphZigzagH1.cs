#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

namespace TDA.Ph;

/// <summary>
/// Slow-correct <b>H1</b> graph-zigzag persistence — Z5b of the zigzag engine (Dey–Hou,
/// <c>2103.07353</c>, <b>Algorithm 2, §4</b>). On graphs the only events that move H1 are edges: an
/// edge added between already-connected endpoints <b>births</b> a 1-cycle (a positive index), and an
/// edge deleted whose endpoints stay connected <b>kills</b> one (a negative index). A death pairs with
/// the <b>earliest</b> compatible unpaired birth — the smallest <c>j ∈ U</c> for which a 1-cycle
/// through both edges lives in every intermediate graph (Pairing Principle).
/// <para>This faithful version realizes the Pairing Principle directly via §4.1's reduction (Prop 19–21):
/// weight each present edge by its latest add-step; the pairing condition <c>j &gt; max{w*, w(σ)}</c>
/// uses <c>w*</c> = the max edge-weight on the <i>minimax</i> (bottleneck) path joining the deleted
/// edge's endpoints in <c>G_{i+1}</c> — computed here by a from-scratch Kruskal union-find, not the
/// near-linear dynamic-MSF data structure (that is Z5c). Connectivity is inline union-find (pure — no
/// <c>CsrGraph</c>). Level-intervals <c>[b,d]</c> map to <see cref="Bar"/> by <see cref="ZigzagBarcodeNaive"/>'s
/// convention, dimension 1. Oracle: <see cref="FastZigzag"/> (Z2) / <see cref="ZigzagBarcodeNaive"/> (Z1).</para>
/// <para>All H1 births on graphs are edge additions (closed birth); deaths are edge deletions (closed)
/// or survive to the end (open at <c>m</c>). Vertex add/delete never touch H1, so they are skipped.</para>
/// </summary>
public static class GraphZigzagH1
{
    /// <param name="representatives">When true, every H1 bar carries a representative 1-cycle in
    /// <see cref="Bar.Cycle"/> as edge cell ids — Dey–Hou Prop 17: each unpaired birth holds a cycle
    /// (the closing edge plus its endpoints' MSF path), and a death symmetric-differences the cycles
    /// that contain the dying edge. The cycle lives across the whole bar and contains both the birth
    /// and death edges.</param>
    public static Barcode Compute(ZigzagFiltration f, bool representatives = false)
    {
        int m = f.Count;
        var bars = new List<Bar>();
        if (m == 0) return new Barcode(bars, "Zigzag Step");

        var isForward = new bool[m];
        for (int s = 0; s < m; s++) isForward[s] = f[s].Direction == ZigzagDirection.Add;

        // Present edges by global id: endpoints + weight (= the step of its latest addition).
        var present = new Dictionary<int, (int U, int V, int Weight)>();
        // Unpaired positive (birth) indices, kept sorted so the smallest qualifying one is cheap.
        var unpaired = new SortedSet<int>();
        var intervals = new List<(int B, int D, int[]? Cyc)>();
        // birth index -> representative 1-cycle (edge cell ids), maintained only when requested.
        var cyc = representatives ? new Dictionary<int, HashSet<int>>() : null;

        for (int s = 0; s < m; s++)
        {
            var step = f[s];
            int id = step.GlobalCellId;

            if (step.Direction == ZigzagDirection.Add)
            {
                int[] bnd = step.BoundaryAtAdd!;
                if (bnd.Length != 2) continue;          // vertex (or higher cell) — no H1 change
                int u = bnd[0], v = bnd[1];
                // Birth iff the endpoints are already connected in G_s (the new edge closes a cycle).
                bool birth = Bottleneck(u, v, present).Connected;
                present[id] = (u, v, s);
                if (birth)
                {
                    unpaired.Add(s + 1);
                    if (cyc != null) { var z = new HashSet<int> { id }; z.SymmetricExceptWith(MsfPath(u, v, present)); cyc[s + 1] = z; }
                }
            }
            else
            {
                if (!present.TryGetValue(id, out var e)) continue;   // vertex deletion — no H1 change
                present.Remove(id);                                   // now G_{s+1}
                var (conn, wStar) = Bottleneck(e.U, e.V, present);
                if (!conn) continue;                                  // bridge — the edge carried no cycle

                int jStar = -1;
                if (cyc != null)
                {
                    // Prop 17: pair with the smallest unpaired birth whose cycle contains the dying edge.
                    foreach (int j in unpaired) if (cyc[j].Contains(id)) { jStar = j; break; }
                    var rep = cyc[jStar];
                    intervals.Add((jStar, s, rep.ToArray()));
                    foreach (int j in unpaired) if (j != jStar && cyc[j].Contains(id)) cyc[j].SymmetricExceptWith(rep);
                    cyc.Remove(jStar);
                }
                else
                {
                    // Bottleneck pairing (§4.1): smallest birth above max{w*, w(σ)}.
                    int threshold = Math.Max(wStar, e.Weight);
                    foreach (int j in unpaired) if (j > threshold) { jStar = j; break; }
                    intervals.Add((jStar, s, null));
                }
                unpaired.Remove(jStar);
            }
        }

        foreach (int j in unpaired) intervals.Add((j, m, cyc?[j].ToArray()));

        foreach (var (b, d, c) in intervals)
        {
            IntervalEnd bEnd = (b > 0 && isForward[b - 1]) ? IntervalEnd.Closed : IntervalEnd.Open;
            IntervalEnd dEnd = (d < m && !isForward[d]) ? IntervalEnd.Closed : IntervalEnd.Open;
            bars.Add(new Bar(b - 1, d, 1, null, null, c, bEnd, dEnd));
        }
        return new Barcode(bars, "Zigzag Step");
    }

    // Edge cell ids on the u–v path of the current MSF (Kruskal by weight, then BFS the tree path).
    static List<int> MsfPath(int u, int v, Dictionary<int, (int U, int V, int Weight)> present)
    {
        var parent = new Dictionary<int, int>();
        int Find(int x) { if (!parent.TryGetValue(x, out int p)) { parent[x] = x; return x; } while (p != x) { parent[x] = parent[p]; x = p; p = parent[x]; } return x; }
        var adj = new Dictionary<int, List<(int To, int Cell)>>();
        void Edge(int a, int b, int cell) { (adj.TryGetValue(a, out var la) ? la : adj[a] = new()).Add((b, cell)); }
        foreach (var kv in present.OrderBy(p => p.Value.Weight))
        {
            var (a, b, _) = kv.Value;
            if (Find(a) != Find(b)) { parent[Find(a)] = Find(b); Edge(a, b, kv.Key); Edge(b, a, kv.Key); }
        }
        var prev = new Dictionary<int, (int From, int Cell)> { [u] = (-1, -1) };
        var q = new Queue<int>(); q.Enqueue(u);
        while (q.Count > 0)
        {
            int x = q.Dequeue();
            if (x == v) break;
            if (!adj.TryGetValue(x, out var lst)) continue;
            foreach (var (y, cell) in lst) if (!prev.ContainsKey(y)) { prev[y] = (x, cell); q.Enqueue(y); }
        }
        var path = new List<int>();
        if (prev.ContainsKey(v)) for (int c = v; c != u; c = prev[c].From) path.Add(prev[c].Cell);
        return path;
    }

    // Minimax (bottleneck) query over the present weighted graph: are u, v connected, and if so what
    // is the largest edge-weight on the lightest-bottleneck u–v path? Kruskal over edges by ascending
    // weight — the weight that first joins u and v is exactly that bottleneck (the max-weight edge on
    // the u–v path of the minimum spanning forest, Prop 21). Rebuilt per query (slow but correct).
    static (bool Connected, int WStar) Bottleneck(int u, int v, Dictionary<int, (int U, int V, int Weight)> present)
    {
        if (u == v) return (true, int.MinValue);
        var parent = new Dictionary<int, int>();

        int Find(int x)
        {
            if (!parent.TryGetValue(x, out int p)) { parent[x] = x; return x; }
            while (p != x) { parent[x] = parent[p]; x = p; p = parent[x]; }
            return x;
        }

        foreach (var edge in present.Values.OrderBy(x => x.Weight))
        {
            int ru = Find(edge.U), rv = Find(edge.V);
            if (ru != rv) parent[ru] = rv;
            if (Find(u) == Find(v)) return (true, edge.Weight);
        }
        return (false, 0);
    }
}
