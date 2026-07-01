#nullable enable
using System;

namespace TDA.Ph.Dynamic;

/// <summary>
/// Splay-based link-cut tree over a dynamic forest, supporting O(log n) <see cref="Connected"/>,
/// link, cut, and — the reason it exists here — <see cref="PathMax"/>: the maximum edge weight on the
/// tree path between two vertices. Used by <see cref="DynamicMsf"/> (Z5c, H1) for the §4.1 bottleneck
/// query "max edge-weight of the u–v path in the MSF".
/// <para>Edges carry weights, so each represented edge is materialized as its own <i>edge node</i>
/// (vertices hold weight −∞); a path's max-weight aggregate is then the heaviest edge on it. Re-rooting
/// uses lazy subtree reversal. Pure, integer-id; internal to <c>TDA.Ph</c>.</para>
/// </summary>
internal sealed class LinkCutTree
{
    sealed class Node
    {
        public Node? L, R, P;
        public bool Rev;
        public int Val;     // edge weight, or int.MinValue for a vertex node
        public int Max;     // max Val in this splay subtree
        public Node(int val) { Val = val; Max = val; }
    }

    readonly Node[] _vertex;

    public LinkCutTree(int n)
    {
        _vertex = new Node[n];
        for (int i = 0; i < n; i++) _vertex[i] = new Node(int.MinValue);
    }

    static bool IsRoot(Node x) => x.P == null || (x.P.L != x && x.P.R != x);

    static int Mx(Node? x) => x?.Max ?? int.MinValue;

    static void Update(Node x) => x.Max = Math.Max(x.Val, Math.Max(Mx(x.L), Mx(x.R)));

    static void ApplyRev(Node x) { (x.L, x.R) = (x.R, x.L); x.Rev = !x.Rev; }

    static void Push(Node x)
    {
        if (!x.Rev) return;
        if (x.L != null) ApplyRev(x.L);
        if (x.R != null) ApplyRev(x.R);
        x.Rev = false;
    }

    static void Rotate(Node x)
    {
        Node p = x.P!, g = p.P!;
        if (!IsRoot(p)) { if (g.L == p) g.L = x; else g.R = x; }
        x.P = g;
        if (p.L == x) { p.L = x.R; if (x.R != null) x.R.P = p; x.R = p; }
        else { p.R = x.L; if (x.L != null) x.L.P = p; x.L = p; }
        p.P = x;
        Update(p); Update(x);
    }

    static void Splay(Node x)
    {
        while (!IsRoot(x))
        {
            Node p = x.P!, g = p.P!;
            if (!IsRoot(p)) Push(g);
            Push(p); Push(x);
            if (!IsRoot(p)) Rotate((g.L == p) == (p.L == x) ? p : x);
            Rotate(x);
        }
        Push(x);
        Update(x);
    }

    static Node Access(Node x)
    {
        Node? last = null;
        for (Node? y = x; y != null; y = y.P)
        {
            Splay(y);
            y.R = last;          // splice the previously-accessed path in as the preferred child
            Update(y);
            last = y;
        }
        Splay(x);
        return last!;
    }

    static void MakeRoot(Node x) { Access(x); ApplyRev(x); Push(x); }

    static Node FindRoot(Node x)
    {
        Access(x);
        while (x.L != null) { Push(x); x = x.L; }
        Splay(x);
        return x;
    }

    // ---- public API (vertices + represented edges) ----------------------------------------------
    static void Link(Node a, Node b) { MakeRoot(a); a.P = b; }

    static void Cut(Node a, Node b)
    {
        MakeRoot(a); Access(b);
        // a is b's left neighbour on the (a..b) path; detach.
        b.L!.P = null; b.L = null; Update(b);
    }

    public bool Connected(int u, int v)
    {
        if (u == v) return true;
        Node a = _vertex[u];
        MakeRoot(a);
        return ReferenceEquals(FindRoot(_vertex[v]), a);
    }

    /// <summary>Insert a represented edge (u,v,w) via a fresh weighted edge node (u and v must be in
    /// different trees).</summary>
    public void LinkEdge(int u, int v, int w)
    {
        var e = new Node(w);
        Link(e, _vertex[u]);
        Link(_vertex[v], e);
        _edge[(u, v)] = e;
        _edge[(v, u)] = e;
    }

    public void CutEdge(int u, int v)
    {
        Node e = _edge[(u, v)];
        Cut(_vertex[u], e);
        Cut(e, _vertex[v]);
        _edge.Remove((u, v));
        _edge.Remove((v, u));
    }

    /// <summary>Max edge weight on the u–v tree path (u,v must be connected); int.MinValue if u==v.</summary>
    public int PathMax(int u, int v)
    {
        if (u == v) return int.MinValue;
        MakeRoot(_vertex[u]);
        Access(_vertex[v]);
        return _vertex[v].Max;
    }

    /// <summary>Weights of the edges on the u–v tree path, in path order (empty if u==v).</summary>
    public System.Collections.Generic.List<int> PathEdgeWeights(int u, int v)
    {
        var acc = new System.Collections.Generic.List<int>();
        if (u == v) return acc;
        MakeRoot(_vertex[u]);
        Access(_vertex[v]);
        // In-order over v's splay tree (= the path u..v); collect edge nodes (Val != MinValue).
        var stack = new System.Collections.Generic.Stack<Node>();
        Node? cur = _vertex[v];
        while (cur != null || stack.Count > 0)
        {
            while (cur != null) { Push(cur); stack.Push(cur); cur = cur.L; }
            cur = stack.Pop();
            if (cur.Val != int.MinValue) acc.Add(cur.Val);
            cur = cur.R;
        }
        return acc;
    }

    readonly System.Collections.Generic.Dictionary<(int, int), Node> _edge = new();
}
