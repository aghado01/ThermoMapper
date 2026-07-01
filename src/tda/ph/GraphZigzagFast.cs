#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using TDA.Ph.Dynamic;

namespace TDA.Ph;

/// <summary>
/// Near-linear graph-zigzag persistence (the <see cref="GraphZigzagAlgorithm.Fast"/> path of
/// <see cref="GraphZigzag.Compute"/>) — Z5c of the zigzag engine (Dey–Hou, <c>2103.07353</c>, <b>§3.1</b>).
/// Computes H0 here and folds in H1 (<see cref="GraphZigzagH1Fast"/>) when <c>maxDimension ≥ 1</c>, so it
/// covers the same dimensions as the reference path. Same Algorithm-1 barcode forest as the slow <see cref="GraphZigzag"/>
/// (Z5a), but the dense per-level copy-up is gone: connectivity events are classified by an
/// <see cref="DynamicConnectivity"/> (HDT, O(log² n) amortized) and the forest is kept <i>sparse</i>
/// (one node per entrance / split / merge — O(m) nodes total) with leaves advanced to the current level
/// only when touched. Merge gluing zips the two leaf paths by level: the current-level leaves coalesce,
/// the split nodes (at distinct past steps) interleave — the same identification the dense union-find
/// performed, without materializing the fillers. A φ-map (component token → leaf) tracks which leaf
/// belongs to which live component.
/// <para>Slow-correct oracle: <see cref="GraphZigzag"/> (Z5a) — every fast path ships with its ground
/// truth. Level-intervals map to <see cref="Bar"/> by <see cref="ZigzagBarcodeNaive"/>'s convention.
/// Pure; the dynamic structures are integer-id only (no <c>CsrGraph</c>).</para>
/// </summary>
public static class GraphZigzagFast
{
    sealed class FNode
    {
        public int Level;
        public FNode? Par;
        public bool Splitting;
        public bool Paired;     // a split consumed by a merge/departure, or an NCA resolved
        public FNode? Glue;     // union-find rep link (coalesced nodes)
        public FNode(int level, FNode? par) { Level = level; Par = par; }
    }

    static FNode Find(FNode n) { while (n.Glue != null) n = n.Glue; return n; }
    static FNode? ParOf(FNode n) { FNode? p = Find(n).Par; return p == null ? null : Find(p); }
    static FNode RootOf(FNode n) { n = Find(n); FNode? p; while ((p = ParOf(n)) != null) n = p; return n; }

    static FNode Hca(FNode a, FNode b)
    {
        var anc = new HashSet<FNode>();
        for (FNode? x = Find(a); x != null; x = ParOf(x)) anc.Add(x);
        for (FNode? y = Find(b); y != null; y = ParOf(y)) if (anc.Contains(y)) return y;
        return RootOf(a);
    }

    static void Coalesce(FNode keep, FNode absorb)
    {
        FNode k = Find(keep), a = Find(absorb);
        if (ReferenceEquals(k, a)) return;
        a.Glue = k;
        if (a.Splitting) k.Splitting = true;
        if (a.Paired) k.Paired = true;
    }

    static void SetPar(FNode node, FNode? par) => Find(node).Par = par == null ? null : Find(par);

    /// <param name="representatives">Same H0 birth-chain reps as <see cref="GraphZigzag"/>: entrance →
    /// <c>{v}</c>, split → the kernel of the deleted edge's two endpoints (one per piece). Emitted only
    /// when present across the whole bar; null otherwise (sound omission).</param>
    /// <param name="strictRepresentatives">Throw instead of leaving an H0 rep null (see
    /// <see cref="GraphZigzag.Compute"/>).</param>
    public static Barcode Compute(ZigzagFiltration f, int maxDimension = 0, bool representatives = false, bool strictRepresentatives = false)
    {
        int m = f.Count;
        var bars = new List<Bar>();
        if (m == 0) return new Barcode(bars, "Zigzag Step");

        var isForward = new bool[m];
        for (int s = 0; s < m; s++) isForward[s] = f[s].Direction == ZigzagDirection.Add;

        // Dense vertex universe: map every vertex cell id to 0..n-1 (HDT needs a fixed vertex set).
        var vid = new Dictionary<int, int>();
        foreach (var step in f)
            if (step.Direction == ZigzagDirection.Add && step.BoundaryAtAdd!.Length == 0 && !vid.ContainsKey(step.GlobalCellId))
                vid[step.GlobalCellId] = vid.Count;
        int n = Math.Max(1, vid.Count);

        var hdt = new DynamicConnectivity(n);
        var edgeEnds = new Dictionary<int, (int U, int W)>();   // edge cell id -> dense endpoints
        var edgeMult = new Dictionary<(int, int), int>();       // canonical vertex pair -> live parallel count (multigraph-safe)
        var edgeOrig = representatives ? new Dictionary<int, (int U, int W)>() : null;  // -> original endpoints
        var isVertexCell = new Dictionary<int, bool>();
        var compLeaf = new Dictionary<object, FNode>();          // component token -> leaf
        var intervals = new List<(int B, int D)>();
        // H0 birth 0-chains (original cell ids) + per-vertex present-interval, as in GraphZigzag.
        var birthChain = representatives ? new Dictionary<int, int[]>() : null;
        var vAdd = representatives ? new Dictionary<int, int>() : null;
        var vDel = representatives ? new Dictionary<int, int>() : null;

        for (int s = 0; s < m; s++)
        {
            var step = f[s];
            int id = step.GlobalCellId;

            if (step.Direction == ZigzagDirection.Add)
            {
                int[] bnd = step.BoundaryAtAdd!;
                if (bnd.Length == 0)
                {
                    // Entrance: isolated vertex already exists in HDT; give its component a root leaf.
                    isVertexCell[id] = true;
                    int v = vid[id];
                    compLeaf[hdt.ComponentToken(v)] = new FNode(s + 1, null);
                    if (birthChain != null) { birthChain[s + 1] = new[] { id }; vAdd![id] = s; }
                }
                else
                {
                    isVertexCell[id] = false;
                    int u = vid[bnd[0]], w = vid[bnd[1]];
                    edgeEnds[id] = (u, w);
                    if (edgeOrig != null) edgeOrig[id] = (bnd[0], bnd[1]);
                    var pair = u < w ? (u, w) : (w, u);
                    int prior = edgeMult.GetValueOrDefault(pair);
                    edgeMult[pair] = prior + 1;
                    if (prior == 0)
                    {
                        // First edge across this vertex pair — the connectivity-relevant insert.
                        bool merge = !hdt.Connected(u, w);
                        FNode? leafU = merge ? compLeaf[hdt.ComponentToken(u)] : null;
                        FNode? leafW = merge ? compLeaf[hdt.ComponentToken(w)] : null;
                        object oldU = hdt.ComponentToken(u), oldW = hdt.ComponentToken(w);
                        hdt.Insert(u, w);
                        if (merge)
                        {
                            FNode newLeaf = Merge(leafU!, leafW!, s, intervals);
                            compLeaf.Remove(oldU); compLeaf.Remove(oldW);
                            compLeaf[hdt.ComponentToken(u)] = newLeaf;
                        }
                        // else: no H0 change (token unchanged for a non-tree edge).
                    }
                    // else: parallel edge — endpoints already connected, no HDT op, no H0 change.
                }
            }
            else if (isVertexCell[id])
            {
                // Departure of an isolated vertex.
                vDel?.Add(id, s);
                int v = vid[id];
                object tok = hdt.ComponentToken(v);
                FNode leaf = compLeaf[tok];
                Departure(leaf, s, intervals);
                compLeaf.Remove(tok);
            }
            else
            {
                // Delete edge.
                var (u, w) = edgeEnds[id];
                var pair = u < w ? (u, w) : (w, u);
                int prior = edgeMult[pair];
                edgeMult[pair] = prior - 1;
                if (prior == 1)
                {
                    // Last edge across this vertex pair removed — the connectivity-relevant delete.
                    object oldTok = hdt.ComponentToken(u);
                    FNode leaf = compLeaf[oldTok];
                    hdt.Delete(u, w);
                    bool split = !hdt.Connected(u, w);
                    if (!split)
                    {
                        // Still connected: same component, possibly a new token after replacement.
                        compLeaf.Remove(oldTok);
                        compLeaf[hdt.ComponentToken(u)] = leaf;
                    }
                    else
                    {
                        var (childU, childW) = Split(leaf, s);
                        compLeaf.Remove(oldTok);
                        compLeaf[hdt.ComponentToken(u)] = childU;
                        compLeaf[hdt.ComponentToken(w)] = childW;
                        if (birthChain != null) { var (ou, ow) = edgeOrig![id]; birthChain[s + 1] = new[] { ou, ow }; }  // kernel = deleted edge's endpoints
                    }
                }
                // else: a parallel edge remains — endpoints stay connected, no HDT op, no H0 change.
            }
        }

        // End: each root -> [level, m]; each active splitting node -> [level+1, m].
        var seen = new HashSet<FNode>();
        foreach (var leaf in compLeaf.Values)
            for (FNode? x = Find(leaf); x != null; x = ParOf(x)) seen.Add(x);
        foreach (var node in seen)
        {
            if (ParOf(node) == null) intervals.Add((node.Level, m));
            if (node.Splitting && !node.Paired) intervals.Add((node.Level + 1, m));
        }

        foreach (var (b, d) in intervals)
        {
            IntervalEnd bEnd = (b > 0 && isForward[b - 1]) ? IntervalEnd.Closed : IntervalEnd.Open;
            IntervalEnd dEnd = (d < m && !isForward[d]) ? IntervalEnd.Closed : IntervalEnd.Open;
            int[]? cyc = null;
            if (birthChain != null)
            {
                int[] chain = birthChain[b];
                bool ok = chain.All(v => vAdd![v] < b && (!vDel!.TryGetValue(v, out int del) || del >= d));
                if (!ok && strictRepresentatives)
                    throw new NotSupportedException(
                        $"H0 representative for the bar born at level {b} cannot be formed as a single " +
                        "persistent 0-chain (a generating vertex departs mid-interval); use FastZigzag " +
                        "for filtrations with vertex deletions, or drop strictRepresentatives.");
                cyc = ok ? chain : null;
            }
            bars.Add(new Bar(b - 1, d, 0, null, null, cyc, bEnd, dEnd));
        }
        // Full near-linear engine: fold in the H1 cycle barcode (Z5c) when asked, mirroring the
        // reference GraphZigzag entry so both paths cover the same dimensions identically.
        if (maxDimension >= 1) bars.AddRange(GraphZigzagH1Fast.Compute(f, representatives).Bars);
        return new Barcode(bars, "Zigzag Step");
    }

    // Advance a leaf to the current level by hanging a fresh node above it (only if it lags).
    static FNode Advance(FNode leaf, int s)
    {
        FNode l = Find(leaf);
        if (l.Level == s) return l;
        return new FNode(s, l);
    }

    static List<FNode> CollectPath(FNode top, int lowLevel)
    {
        var path = new List<FNode>();
        for (FNode? x = Find(top); x != null && x.Level >= lowLevel; x = ParOf(x)) path.Add(x);
        return path;
    }

    static FNode Merge(FNode leafA, FNode leafB, int s, List<(int, int)> intervals)
    {
        FNode rA = RootOf(leafA), rB = RootOf(leafB);
        bool same = ReferenceEquals(rA, rB);
        int jLow;
        if (!same)
        {
            int j = Math.Max(rA.Level, rB.Level);
            intervals.Add((j, s));
            jLow = j;
        }
        else
        {
            FNode v = Hca(leafA, leafB);
            int j = v.Level;
            intervals.Add((j + 1, s));
            jLow = j + 1;
            Find(v).Paired = true;
        }

        FNode advA = Advance(leafA, s);
        FNode advB = Advance(leafB, s);
        var pathA = CollectPath(advA, jLow);
        var pathB = CollectPath(advB, jLow);

        // Continuation below jLow: the shared ancestor (same-tree NCA) or the elder's deeper path.
        FNode? contA = ParOf(pathA[^1]);
        FNode? contB = ParOf(pathB[^1]);
        FNode? cont = contA == null ? contB
                    : contB == null ? contA
                    : (contA.Level <= contB.Level ? contA : contB);

        // Zip the two level-descending paths; coalesce equal levels (the current-level leaves, and any
        // diff-tree level-j ancestor), interleave the rest.
        var merged = new List<FNode>();
        int ia = 0, ib = 0;
        while (ia < pathA.Count && ib < pathB.Count)
        {
            FNode a = Find(pathA[ia]), b = Find(pathB[ib]);
            if (a.Level == b.Level) { Coalesce(a, b); merged.Add(a); ia++; ib++; }
            else if (a.Level > b.Level) { merged.Add(a); ia++; }
            else { merged.Add(b); ib++; }
        }
        while (ia < pathA.Count) merged.Add(Find(pathA[ia++]));
        while (ib < pathB.Count) merged.Add(Find(pathB[ib++]));

        for (int i = 0; i + 1 < merged.Count; i++) SetPar(merged[i], merged[i + 1]);
        SetPar(merged[^1], cont);

        return new FNode(s + 1, Find(merged[0]));   // new leaf hangs off the top (level s) node
    }

    static (FNode, FNode) Split(FNode leaf, int s)
    {
        FNode adv = Advance(leaf, s);
        adv.Splitting = true;
        return (new FNode(s + 1, adv), new FNode(s + 1, adv));
    }

    static void Departure(FNode leaf, int s, List<(int, int)> intervals)
    {
        FNode? chosen = null;
        for (FNode? x = Find(leaf); x != null; x = ParOf(x)) if (x.Splitting && !x.Paired) { chosen = x; break; }
        if (chosen != null) { intervals.Add((chosen.Level + 1, s)); chosen.Paired = true; }
        else intervals.Add((RootOf(leaf).Level, s));
    }
}
