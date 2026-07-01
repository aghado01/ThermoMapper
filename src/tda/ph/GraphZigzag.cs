#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

namespace TDA.Ph;

/// <summary>
/// Slow-correct <b>H0</b> graph-zigzag persistence — Z5a of the zigzag engine (Dey–Hou,
/// <c>2103.07353</c>, <b>Algorithm 1, §3</b>). Faithful <i>barcode forest</i>: leaves track current
/// connected components, splits add a splitting node + two children, and merges <b>glue the two paths
/// to their level-j ancestors</b> (union-find — siblings of an earlier split share the lower path, so
/// gluing must identify nodes, not discard them). Dense per-level nodes keep the level-j ancestors
/// present; the near-linear data structures are Z5c.
/// <para>Cells are dim-0 (vertices, empty boundary) and dim-1 (edges, two-vertex boundary). Connectivity
/// is an inline component map + BFS (pure — no <c>CsrGraph</c>). Level-intervals <c>[b,d]</c> map to
/// <see cref="Bar"/> by <see cref="ZigzagBarcodeNaive"/>'s convention. Oracle: <see cref="FastZigzag"/>
/// (Z2) / <see cref="ZigzagBarcodeNaive"/> (Z1).</para>
/// <para><see cref="Compute"/> is the single graph-zigzag entry point: it always returns H0, and folds
/// in the H1 cycle barcode from <see cref="GraphZigzagH1"/> (Z5b) when <c>maxDimension &gt;= 1</c>.
/// Graphs are 1-complexes, so those are the only two dimensions.</para>
/// <para><b>STATUS — parity-complete (H0).</b> Agrees with both oracles across the hand zoo and random
/// edge-/vertex-churn sweeps. Two subtleties the forest transcription must respect: a merge must
/// <i>capture</i> both glue paths before unioning (see <see cref="GluePaths"/>) so every shared level is
/// identified, not just the leaf; and a departure pairs with the <i>highest-level</i> (first-encountered)
/// splitting ancestor — Dey–Hou's "highest", where higher means higher level — leaving lower splits on
/// the path unpaired. Next: Z5b (H1, §4 + Prop 19) and Z5c (near-linear data structures).</para>
/// </summary>

/// <summary>Which implementation the graph-zigzag entry runs: <see cref="Reference"/> is the simple,
/// obviously-correct oracle (forest + Kruskal); <see cref="Fast"/> is the near-linear path (HDT + dynamic
/// MSF). They are independent so they cross-validate each other and the Z1/Z2 oracles.</summary>
public enum GraphZigzagAlgorithm { Reference, Fast }

/// <summary>
/// Single graph-zigzag entry point. Returns one <see cref="Barcode"/> covering the requested dimensions
/// (H0 always; H1 folded when <paramref name="maxDimension"/> ≥ 1), with the implementation chosen by
/// <paramref name="algorithm"/>. The Reference path lives here; the Fast path is composed from the
/// near-linear components. Mirrors ripserer's one-entry / algorithm-dispatch / uniform-output shape.
/// </summary>
public static class GraphZigzag
{
    sealed class Node
    {
        public int Level;
        public Node? Par;       // toward root (strictly lower level); null = root
        public bool Splitting;
        public Node? Glue;      // union-find representative link; null = self is representative
        public Node(int level, Node? par) { Level = level; Par = par; }
    }

    static Node Find(Node n) { while (n.Glue != null) n = n.Glue; return n; }
    static Node? ParOf(Node n) { Node? p = Find(n).Par; return p == null ? null : Find(p); }
    static Node RootOf(Node n) { n = Find(n); Node? p; while ((p = ParOf(n)) != null) n = p; return n; }
    static Node AtLevel(Node leaf, int level) { Node n = Find(leaf); while (n.Level > level) n = ParOf(n)!; return n; }

    static Node Hca(Node a, Node b)
    {
        var anc = new HashSet<Node>();
        for (Node? x = Find(a); x != null; x = ParOf(x)) anc.Add(x);
        for (Node? y = Find(b); y != null; y = ParOf(y)) if (anc.Contains(y)) return y;
        return RootOf(a);
    }

    static void Union(Node keep, Node absorb)
    {
        Node k = Find(keep), a = Find(absorb);
        if (ReferenceEquals(k, a)) return;
        a.Glue = k;
        if (a.Splitting) k.Splitting = true;
    }

    // Glue two dense paths to their common level by unioning them node-for-node over [low..high].
    // The two paths must be captured BEFORE any union: unioning rewires Find, so re-evaluating
    // AtLevel(absorbLeaf, L) after the top union would re-walk the *keep* path and silently glue
    // only the top level — orphaning the absorbed path's lower nodes (and any siblings hanging off
    // them). Capturing first keeps every level identified, so HCA/RootOf stay correct downstream.
    static void GluePaths(Node keepLeaf, Node absorbLeaf, int low, int high)
    {
        int n = high - low + 1;
        var keep = new Node[n];
        var absorb = new Node[n];
        for (int L = high; L >= low; L--) { keep[high - L] = AtLevel(keepLeaf, L); absorb[high - L] = AtLevel(absorbLeaf, L); }
        for (int i = 0; i < n; i++) Union(keep[i], absorb[i]);
    }

    /// <param name="representatives">When true, every bar carries a representative cycle in
    /// <see cref="Bar.Cycle"/> (vertex cell ids for H0, edge cell ids for the folded-in H1). The H0 rep
    /// of a bar is its <i>birth 0-chain</i>: an entrance births the singleton <c>{v}</c>; a split births
    /// the kernel <c>{minA, minB}</c> (the dual of the H1 cycle — a vertex from each side, symmetric so no
    /// "which piece is new" guess). The forest's same-tree/diff-tree pairing then assigns each bar its
    /// birth chain (<c>birthChain[B]</c>). The chain is emitted only when all its vertices live across the
    /// whole bar; if a vertex departs mid-interval the rep is left null (that case needs the full evolving
    /// 0-chain reduction — see <see cref="FastZigzag"/>). For a fixed vertex set (e.g. nerve-over-T) every
    /// H0 bar gets a sound rep.</param>
    /// <param name="strictRepresentatives">Opt-in: when true (with <paramref name="representatives"/>),
    /// throw <see cref="NotSupportedException"/> instead of leaving an H0 rep null — for callers that need
    /// every bar represented and would rather fail loud than receive a silent partial. Default false keeps
    /// the sound lenient behaviour (omit the unrepresentable rep). H1 reps are always complete, so this
    /// only ever fires on H0 bars over filtrations with vertex departures.</param>
    /// <param name="algorithm">Fast (near-linear, default) or Reference (the simple oracle). Both produce
    /// the identical barcode and equally-sound reps; Reference exists to cross-validate Fast and is what
    /// the oracle-style tests pin.</param>
    public static Barcode Compute(ZigzagFiltration f, int maxDimension = 0, bool representatives = false,
                                  bool strictRepresentatives = false,
                                  GraphZigzagAlgorithm algorithm = GraphZigzagAlgorithm.Fast)
    {
        if (algorithm == GraphZigzagAlgorithm.Fast)
            return GraphZigzagFast.Compute(f, maxDimension, representatives, strictRepresentatives);

        int m = f.Count;
        var bars = new List<Bar>();
        if (m == 0) return new Barcode(bars, "Zigzag Step");

        var isForward = new bool[m];
        for (int s = 0; s < m; s++) isForward[s] = f[s].Direction == ZigzagDirection.Add;

        var vertexComp = new Dictionary<int, int>();
        var adjacency = new Dictionary<int, Dictionary<int, int>>();  // neighbor -> edge multiplicity (multigraph-safe)
        var edgeEndpoints = new Dictionary<int, (int U, int W)>();
        var isVertexCell = new Dictionary<int, bool>();
        var compNode = new Dictionary<int, Node>();
        int nextComp = 0;
        var intervals = new List<(int B, int D)>();
        // H0 birth 0-chains, write-once per birth level (entrance singleton / split kernel), plus the
        // present-interval of each vertex cell so a rep is only emitted when present across its whole bar.
        var birthChain = representatives ? new Dictionary<int, int[]>() : null;
        var vAdd = representatives ? new Dictionary<int, int>() : null;
        var vDel = representatives ? new Dictionary<int, int>() : null;

        // Continuation: every component not directly touched by the step copies its leaf up to level+1,
        // keeping all paths dense (a node at every level, so merge gluing always finds level-j ancestors).
        void CopyUpExcept(HashSet<int> involved, int level)
        {
            foreach (int c in compNode.Keys.ToList())
                if (!involved.Contains(c)) compNode[c] = new Node(level, compNode[c]);
        }
        var none = new HashSet<int>();

        for (int s = 0; s < m; s++)
        {
            var step = f[s];
            int id = step.GlobalCellId;

            if (step.Direction == ZigzagDirection.Add)
            {
                int[] bnd = step.BoundaryAtAdd!;
                if (bnd.Length == 0)
                {
                    // Entrance.
                    CopyUpExcept(none, s + 1);
                    isVertexCell[id] = true;
                    int c = nextComp++;
                    compNode[c] = new Node(s + 1, null);
                    vertexComp[id] = c;
                    adjacency[id] = new Dictionary<int, int>();
                    if (birthChain != null) { birthChain[s + 1] = new[] { id }; vAdd![id] = s; }   // entrance: singleton class
                }
                else
                {
                    isVertexCell[id] = false;
                    int u = bnd[0], w = bnd[1];
                    edgeEndpoints[id] = (u, w);
                    adjacency[u][w] = adjacency[u].GetValueOrDefault(w) + 1;
                    adjacency[w][u] = adjacency[w].GetValueOrDefault(u) + 1;
                    int cu = vertexComp[u], cw = vertexComp[w];
                    if (cu == cw)
                    {
                        CopyUpExcept(none, s + 1);     // no H0 change
                    }
                    else
                    {
                        // Merge — one class dies; glue the two level-s leaves' paths.
                        CopyUpExcept(new HashSet<int> { cu, cw }, s + 1);
                        Node u1 = compNode[cu], u2 = compNode[cw];
                        Node r1 = RootOf(u1), r2 = RootOf(u2);
                        Node mergedBase;
                        if (!ReferenceEquals(r1, r2))
                        {
                            // Different trees: the higher (younger) root's class dies; glue down to it.
                            int j = Math.Max(r1.Level, r2.Level);
                            intervals.Add((j, s));
                            Node elder = r1.Level < r2.Level ? u1 : u2;
                            Node younger = r1.Level < r2.Level ? u2 : u1;
                            GluePaths(elder, younger, j, s);
                            mergedBase = elder;
                        }
                        else
                        {
                            // Same tree: the split-born class at the highest common ancestor dies.
                            Node v = Hca(u1, u2);
                            int j = v.Level;
                            intervals.Add((j + 1, s));
                            GluePaths(u1, u2, j + 1, s);
                            Find(v).Splitting = false;     // the split is undone (one branch remains)
                            mergedBase = u1;
                        }
                        var mNode = new Node(s + 1, Find(mergedBase));
                        foreach (int vtx in vertexComp.Keys.ToList())
                            if (vertexComp[vtx] == cw) vertexComp[vtx] = cu;
                        compNode.Remove(cw);
                        compNode[cu] = mNode;
                    }
                }
            }
            else if (isVertexCell[id])
            {
                // Departure (the vertex is isolated -> a singleton component leaves).
                vDel?.Add(id, s);
                int c = vertexComp[id];
                CopyUpExcept(new HashSet<int> { c }, s + 1);
                Node u = compNode[c];
                // Pair with the *highest-level* (first-encountered, i.e. lowest-in-tree) unpaired
                // splitting ancestor — Dey–Hou's "highest splitting ancestor", where they define
                // higher = higher level. Walking u -> root, that is the FIRST splitting node hit.
                // Lower splitting ancestors on the path stay unpaired (they resolve on later events
                // or at the end); taking the last one (lowest level) misdates nested-split departures.
                Node? chosen = null;
                for (Node? x = Find(u); x != null; x = ParOf(x)) if (x.Splitting) { chosen = x; break; }
                if (chosen != null)
                {
                    intervals.Add((chosen.Level + 1, s));
                    chosen.Splitting = false;     // delete u's branch -> the split loses a branch
                }
                else
                {
                    intervals.Add((RootOf(u).Level, s));
                }
                compNode.Remove(c);
                vertexComp.Remove(id);
                adjacency.Remove(id);
            }
            else
            {
                // Delete edge.
                var (u, w) = edgeEndpoints[id];
                if (--adjacency[u][w] == 0) adjacency[u].Remove(w);
                if (--adjacency[w][u] == 0) adjacency[w].Remove(u);
                if (Connected(u, w, adjacency))
                {
                    CopyUpExcept(none, s + 1);     // still connected, no H0 change
                }
                else
                {
                    // Split: the level-s leaf becomes a splitting node with two level-(s+1) children.
                    int c = vertexComp[u];
                    CopyUpExcept(new HashSet<int> { c }, s + 1);
                    Node uLeaf = Find(compNode[c]);
                    uLeaf.Splitting = true;
                    var childU = new Node(s + 1, uLeaf);
                    var childW = new Node(s + 1, uLeaf);
                    var sideU = ComponentOf(u, adjacency);
                    var sideW = ComponentOf(w, adjacency);
                    int cU = nextComp++, cW = nextComp++;
                    foreach (int x in sideU) vertexComp[x] = cU;
                    foreach (int x in sideW) vertexComp[x] = cW;
                    compNode.Remove(c);
                    compNode[cU] = childU;
                    compNode[cW] = childW;
                    if (birthChain != null) birthChain[s + 1] = new[] { sideU.Min(), sideW.Min() };  // split: kernel [A]+[B]
                }
            }
        }

        // End: each root -> [level, m]; each splitting node -> [level+1, m].
        var seen = new HashSet<Node>();
        foreach (var leaf in compNode.Values)
            for (Node? x = Find(leaf); x != null; x = ParOf(x)) seen.Add(x);
        foreach (var node in seen)
        {
            if (ParOf(node) == null) intervals.Add((node.Level, m));
            if (node.Splitting) intervals.Add((node.Level + 1, m));
        }

        foreach (var (b, d) in intervals)
        {
            IntervalEnd bEnd = (b > 0 && isForward[b - 1]) ? IntervalEnd.Closed : IntervalEnd.Open;
            IntervalEnd dEnd = (d < m && !isForward[d]) ? IntervalEnd.Closed : IntervalEnd.Open;
            int[]? cyc = null;
            if (birthChain != null)
            {
                // Emit the birth chain only when every vertex of it is present across the whole bar
                // [b,d]; a vertex deleted mid-interval would make the chain a wrong representative (that
                // case needs the full evolving 0-chain reduction — left to FastZigzag for now).
                int[] chain = birthChain[b];
                bool ok = chain.All(v => vAdd![v] < b && (!vDel!.TryGetValue(v, out int del) || del >= d));
                if (!ok && strictRepresentatives)
                    throw new NotSupportedException(
                        $"H0 representative for the bar born at level {b} cannot be formed as a single " +
                        "persistent 0-chain because a generating vertex departs mid-interval. Use " +
                        "FastZigzag.Compute(.., representatives: true) for filtrations with vertex deletions, " +
                        "or call without strictRepresentatives to accept a null (omitted) rep.");
                cyc = ok ? chain : null;
            }
            bars.Add(new Bar(b - 1, d, 0, null, null, cyc, bEnd, dEnd));
        }
        // Graphs are 1-complexes, so H1 is the only higher dimension; fold in the cycle barcode
        // (Z5b) when asked. Keeps a single graph-zigzag entry point that mirrors the oracles.
        if (maxDimension >= 1) bars.AddRange(GraphZigzagH1.Compute(f, representatives).Bars);
        return new Barcode(bars, "Zigzag Step");
    }

    static bool Connected(int a, int b, Dictionary<int, Dictionary<int, int>> adj)
    {
        if (a == b) return true;
        var seen = new HashSet<int> { a };
        var stack = new Stack<int>();
        stack.Push(a);
        while (stack.Count > 0)
        {
            int x = stack.Pop();
            foreach (int y in adj[x].Keys) if (seen.Add(y)) { if (y == b) return true; stack.Push(y); }
        }
        return false;
    }

    static List<int> ComponentOf(int start, Dictionary<int, Dictionary<int, int>> adj)
    {
        var seen = new HashSet<int> { start };
        var stack = new Stack<int>();
        stack.Push(start);
        var result = new List<int> { start };
        while (stack.Count > 0)
        {
            int x = stack.Pop();
            foreach (int y in adj[x].Keys) if (seen.Add(y)) { result.Add(y); stack.Push(y); }
        }
        return result;
    }
}
