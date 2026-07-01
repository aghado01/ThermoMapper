#nullable enable
using System;
using System.Collections.Generic;

namespace TDA.Ph.Dynamic;

/// <summary>
/// Euler-Tour Tree (ETT) — one level of the Holm–de Lichtenberg–Thorup dynamic forest. Represents a
/// forest of <i>tree</i> edges over a fixed vertex set as a set of balanced sequences (one per tree),
/// each the Euler tour of its tree. Backed by an <b>implicit treap</b> (randomized BST keyed by
/// position; deterministic SplitMix64 priorities, no <c>Random</c>), giving O(log n) expected
/// <see cref="Link"/> / <see cref="Cut"/> / <see cref="Connected"/> / re-root.
/// <para>The tour interleaves <b>vertex occurrences</b> (one permanent node per vertex) with <b>arc
/// occurrences</b> (two per tree edge, <c>u→v</c> and <c>v→u</c>). Linking u,v re-roots both tours at
/// u and v and concatenates <c>tour(u) · arc(u,v) · tour(v) · arc(v,u)</c>; cutting splices out the two
/// arcs and the enclosed sub-tour. Subtree aggregates carried for HDT: vertex count (the n/2^i size
/// invariant), plus two finger-markable classes — vertices incident to a level non-tree edge, and
/// level tree-edge arcs — each with a subtree counter so a marked element in a tree is found in
/// O(log n). Pure (integer ids only); internal to <c>TDA.Ph</c>.</para>
/// </summary>
internal sealed class EulerTourTree
{
    sealed class Node
    {
        public long Priority;
        public Node? Left, Right, Parent;
        public int Size = 1;          // nodes in subtree (implicit-treap index)
        public int Vertices;          // vertex-occurrence nodes in subtree
        public int VMarkCount;        // marked vertex occurrences in subtree
        public int EMarkCount;        // marked arc occurrences in subtree

        public bool IsVertex;         // vertex occurrence vs arc occurrence
        public int A, B;              // vertex id (A, with IsVertex) or arc endpoints (A->B)
        public bool VMark;            // this vertex has an incident level non-tree edge
        public bool EMark;            // this arc is the marked representative of a level tree edge
    }

    static int Sz(Node? n) => n?.Size ?? 0;
    static int Vc(Node? n) => n?.Vertices ?? 0;
    static int Vm(Node? n) => n?.VMarkCount ?? 0;
    static int Em(Node? n) => n?.EMarkCount ?? 0;

    static void Update(Node n)
    {
        n.Size = 1 + Sz(n.Left) + Sz(n.Right);
        n.Vertices = (n.IsVertex ? 1 : 0) + Vc(n.Left) + Vc(n.Right);
        n.VMarkCount = (n.IsVertex && n.VMark ? 1 : 0) + Vm(n.Left) + Vm(n.Right);
        n.EMarkCount = (!n.IsVertex && n.EMark ? 1 : 0) + Em(n.Left) + Em(n.Right);
    }

    // Deterministic priority stream (SplitMix64) — distinct, reproducible, no Random.
    ulong _rngState;
    long NextPriority()
    {
        _rngState += 0x9E3779B97F4A7C15UL;
        ulong z = _rngState;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
        return (long)(z ^ (z >> 31));
    }

    readonly Node[] _vertex;                       // vertexNode[v]
    readonly Dictionary<(int, int), Node> _arc = new();  // arc[(u,v)]

    public EulerTourTree(int n)
    {
        _vertex = new Node[n];
        for (int v = 0; v < n; v++)
        {
            var node = new Node { Priority = NextPriority(), IsVertex = true, A = v };
            Update(node);
            _vertex[v] = node;
        }
    }

    static Node Root(Node n) { while (n.Parent != null) n = n.Parent; return n; }

    // Position (0-based index) of n within its treap, via left-subtree sizes walking to the root.
    static int Position(Node n)
    {
        int pos = Sz(n.Left);
        for (Node cur = n; cur.Parent != null; cur = cur.Parent)
            if (cur.Parent.Right == cur) pos += Sz(cur.Parent.Left) + 1;
        return pos;
    }

    // Split treap `t` into the first `k` nodes (L) and the rest (R). Detaches t from any parent.
    static (Node? L, Node? R) Split(Node? t, int k)
    {
        if (t == null) return (null, null);
        t.Parent = null;
        int leftSize = Sz(t.Left);
        if (leftSize >= k)
        {
            var (l, r) = Split(t.Left, k);
            t.Left = r; if (r != null) r.Parent = t;
            Update(t);
            if (l != null) l.Parent = null;
            return (l, t);
        }
        else
        {
            var (l, r) = Split(t.Right, k - leftSize - 1);
            t.Right = l; if (l != null) l.Parent = t;
            Update(t);
            if (r != null) r.Parent = null;
            return (t, r);
        }
    }

    static Node? Merge(Node? a, Node? b)
    {
        if (a == null) return b;
        if (b == null) return a;
        if (a.Priority > b.Priority)
        {
            a.Right = Merge(a.Right, b); a.Right!.Parent = a;
            Update(a); a.Parent = null; return a;
        }
        else
        {
            b.Left = Merge(a, b.Left); b.Left!.Parent = b;
            Update(b); b.Parent = null; return b;
        }
    }

    // Rotate a tree's tour so vertexNode[v] is first.
    void Reroot(int v)
    {
        Node node = _vertex[v];
        Node root = Root(node);
        int k = Position(node);
        var (l, r) = Split(root, k);
        Merge(r, l);
    }

    public bool Connected(int u, int v) => ReferenceEquals(Root(_vertex[u]), Root(_vertex[v]));

    /// <summary>Opaque per-component token (the tour's treap root); stable until v's tree changes.</summary>
    public object ComponentToken(int v) => Root(_vertex[v]);

    /// <summary>Vertices in v's tree (= the tree's size, for the n/2^i invariant).</summary>
    public int TreeSize(int v) => Root(_vertex[v]).Vertices;

    public void Link(int u, int v)
    {
        Reroot(u);
        Reroot(v);
        Node ru = Root(_vertex[u]);
        Node rv = Root(_vertex[v]);
        var arcUV = new Node { Priority = NextPriority(), IsVertex = false, A = u, B = v };
        var arcVU = new Node { Priority = NextPriority(), IsVertex = false, A = v, B = u };
        Update(arcUV); Update(arcVU);
        _arc[(u, v)] = arcUV;
        _arc[(v, u)] = arcVU;
        Merge(Merge(Merge(ru, arcUV), rv), arcVU);
    }

    public void Cut(int u, int v)
    {
        Node a1 = _arc[(u, v)], a2 = _arc[(v, u)];
        Node root = Root(a1);
        int p1 = Position(a1), p2 = Position(a2);
        if (p1 > p2) { (p1, p2) = (p2, p1); (a1, a2) = (a2, a1); }
        // tour = [A : p1][a1][B : p2-p1-1][a2][C]
        var (A, restA) = Split(root, p1);
        var (n1, restB) = Split(restA, 1);     // n1 == a1
        var (B, restC) = Split(restB, p2 - p1 - 1);
        var (n2, C) = Split(restC, 1);         // n2 == a2
        // B becomes its own tree; A and C rejoin. Discard the two arc nodes.
        if (n1 != null) { n1.Left = n1.Right = n1.Parent = null; Update(n1); }
        if (n2 != null) { n2.Left = n2.Right = n2.Parent = null; Update(n2); }
        if (B != null) B.Parent = null;
        Merge(A, C);
        _arc.Remove((u, v));
        _arc.Remove((v, u));
    }

    void Refresh(Node n) { for (Node? c = n; c != null; c = c.Parent) Update(c); }

    public void SetVertexMark(int v, bool mark)
    {
        Node n = _vertex[v];
        if (n.VMark == mark) return;
        n.VMark = mark;
        Refresh(n);
    }

    public void SetEdgeMark(int u, int v, bool mark)
    {
        Node n = _arc[(u, v)];
        if (n.EMark == mark) return;
        n.EMark = mark;
        Refresh(n);
    }

    /// <summary>Any vertex in v's tree carrying a set vertex-mark, or -1 if none.</summary>
    public int FindMarkedVertex(int v)
    {
        Node? n = Root(_vertex[v]);
        if (Vm(n) == 0) return -1;
        while (n != null)
        {
            if (Vm(n.Left) > 0) { n = n.Left; continue; }
            if (n.IsVertex && n.VMark) return n.A;
            n = n.Right;
        }
        return -1;
    }

    /// <summary>All vertices in v's tree carrying a set vertex-mark.</summary>
    public List<int> CollectMarkedVertices(int v)
    {
        var acc = new List<int>();
        CollectMarked(Root(_vertex[v]), acc);
        return acc;
    }

    static void CollectMarked(Node? n, List<int> acc)
    {
        if (n == null || n.VMarkCount == 0) return;
        CollectMarked(n.Left, acc);
        if (n.IsVertex && n.VMark) acc.Add(n.A);
        CollectMarked(n.Right, acc);
    }

    /// <summary>Any marked tree-edge arc (u,v) in v's tree, or null if none.</summary>
    public (int U, int V)? FindMarkedEdge(int v)
    {
        Node? n = Root(_vertex[v]);
        if (Em(n) == 0) return null;
        while (n != null)
        {
            if (Em(n.Left) > 0) { n = n.Left; continue; }
            if (!n.IsVertex && n.EMark) return (n.A, n.B);
            n = n.Right;
        }
        return null;
    }
}
