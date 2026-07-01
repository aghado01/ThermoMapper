#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using TDA.Ph.Dynamic;

namespace TDA.Ph;

/// <summary>
/// Near-linear <b>H1</b> graph-zigzag persistence — Z5c of the zigzag engine (Dey–Hou,
/// <c>2103.07353</c>, <b>§4.1</b>). Same Algorithm-2 pairing as the slow <see cref="GraphZigzagH1"/>
/// (Z5b), but the per-deletion Kruskal rebuild is replaced by a <see cref="DynamicMsf"/>: each edge is
/// weighted by its add-step, an edge added between already-connected endpoints births a cycle, and a
/// cycle-killing deletion pairs with the earliest unpaired birth <c>j &gt; max{w*, w(σ)}</c> where
/// <c>w*</c> is the <see cref="DynamicMsf.PathMax">MSF path-max</see> between the deleted edge's
/// endpoints in <c>G_{i+1}</c> — the §4.1 minimum-spanning-forest bottleneck, now answered in
/// O(log n) instead of O(E log E).
/// <para>Slow-correct oracle: <see cref="GraphZigzagH1"/> (Z5b). Vertex events never touch H1 and are
/// skipped. Level-intervals map to <see cref="Bar"/> by <see cref="ZigzagBarcodeNaive"/>'s convention,
/// dimension 1. Pure; the dynamic structure is integer-id only.</para>
/// </summary>
public static class GraphZigzagH1Fast
{
    /// <param name="representatives">When true, each H1 bar carries a representative 1-cycle (edge cell
    /// ids) in <see cref="Bar.Cycle"/> — Prop 17, identical to <see cref="GraphZigzagH1"/>, but the MSF
    /// path is read from the dynamic MSF rather than rebuilt by Kruskal.</param>
    public static Barcode Compute(ZigzagFiltration f, bool representatives = false)
    {
        int m = f.Count;
        var bars = new List<Bar>();
        if (m == 0) return new Barcode(bars, "Zigzag Step");

        var isForward = new bool[m];
        for (int s = 0; s < m; s++) isForward[s] = f[s].Direction == ZigzagDirection.Add;

        var vid = new Dictionary<int, int>();
        foreach (var step in f)
            if (step.Direction == ZigzagDirection.Add && step.BoundaryAtAdd!.Length == 0 && !vid.ContainsKey(step.GlobalCellId))
                vid[step.GlobalCellId] = vid.Count;
        int n = Math.Max(1, vid.Count);

        var msf = new DynamicMsf(n);
        var present = new Dictionary<int, (int U, int V, int W)>();   // edge cell id -> dense ends + weight
        var unpaired = new SortedSet<int>();
        var intervals = new List<(int B, int D, int[]? Cyc)>();
        // Representatives (Prop 17, as in GraphZigzagH1): per-birth cycle of edge cell ids; the MSF path
        // edges come from the dynamic MSF (weights) mapped back to cells via weightToCell.
        var cyc = representatives ? new Dictionary<int, HashSet<int>>() : null;
        var weightToCell = representatives ? new Dictionary<int, int>() : null;

        for (int s = 0; s < m; s++)
        {
            var step = f[s];
            int id = step.GlobalCellId;

            if (step.Direction == ZigzagDirection.Add)
            {
                int[] bnd = step.BoundaryAtAdd!;
                if (bnd.Length != 2) continue;            // vertex — no H1 change
                int u = vid[bnd[0]], v = vid[bnd[1]];
                if (weightToCell != null) weightToCell[s] = id;
                bool cycle = msf.Connected(u, v);         // closing a cycle => birth
                msf.Insert(u, v, s);                      // weight = add-step (strictly increasing)
                present[id] = (u, v, s);
                if (cycle)
                {
                    unpaired.Add(s + 1);
                    if (cyc != null)
                    {
                        var z = new HashSet<int> { id };  // the closing edge + its endpoints' MSF path
                        foreach (int w in msf.PathEdgeWeights(u, v)) z.Add(weightToCell![w]);
                        cyc[s + 1] = z;
                    }
                }
            }
            else
            {
                if (!present.TryGetValue(id, out var e)) continue;   // vertex deletion — no H1 change
                present.Remove(id);
                msf.Delete(e.U, e.V);                                 // now MSF of G_{s+1}
                if (!msf.Connected(e.U, e.V)) continue;               // bridge — carried no cycle

                int jStar = -1;
                if (cyc != null)
                {
                    foreach (int j in unpaired) if (cyc[j].Contains(id)) { jStar = j; break; }   // Prop 17 pairing
                    var rep = cyc[jStar];
                    intervals.Add((jStar, s, rep.ToArray()));
                    foreach (int j in unpaired) if (j != jStar && cyc[j].Contains(id)) cyc[j].SymmetricExceptWith(rep);
                    cyc.Remove(jStar);
                }
                else
                {
                    int threshold = Math.Max(msf.PathMax(e.U, e.V), e.W);
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
}
