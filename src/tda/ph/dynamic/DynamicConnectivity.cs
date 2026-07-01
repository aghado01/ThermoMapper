#nullable enable
using System;
using System.Collections.Generic;

namespace TDA.Ph.Dynamic;

/// <summary>
/// Fully-dynamic connectivity (Holm–de Lichtenberg–Thorup, 2001) — edge insert/delete and
/// connectivity queries in O(log² n) amortized on a graph with a fixed vertex set. This is the
/// dynamic-connectivity structure §3.1 of Dey–Hou uses for H0 graph-zigzag (Z5c).
/// <para>Each edge carries a <i>level</i> in <c>[0, L)</c>. <c>F_i</c> = the spanning forest induced by
/// tree edges of level ≥ i, maintained as <see cref="EulerTourTree"/> <c>_forest[i]</c>; the invariants
/// are (1) <c>F_i</c> spans the level-≥i subgraph and (2) every tree in <c>F_i</c> has ≤ <c>n/2^i</c>
/// vertices. New edges enter at level 0 (tree edge if it joins two components, else a non-tree edge held
/// in per-level adjacency). Deleting a tree edge searches for a replacement from its level down to 0:
/// push the smaller side's level-i tree edges up to i+1 (legal by invariant 2), then scan that side's
/// level-i non-tree edges — the first one crossing to the other side reconnects; the rest are pushed to
/// i+1. Markable vertices (incident non-tree edge) and tree-edge arcs are found via the ETT's subtree
/// fingers. Pure, integer-id; internal to <c>TDA.Ph</c>.</para>
/// </summary>
internal sealed class DynamicConnectivity
{
    readonly int _n;
    readonly int _levels;
    readonly EulerTourTree[] _forest;
    readonly Dictionary<int, HashSet<int>>[] _adj;       // _adj[i][x] = level-i non-tree neighbors of x
    readonly Dictionary<(int, int), int> _edgeLevel = new();
    readonly HashSet<(int, int)> _treeEdges = new();     // canonical (min,max) tree edges

    public DynamicConnectivity(int n)
    {
        _n = Math.Max(1, n);
        int maxLevel = 1;
        while ((1 << maxLevel) < _n) maxLevel++;          // 2^maxLevel >= n
        _levels = maxLevel + 2;                           // headroom; edges never need beyond ~log2 n
        _forest = new EulerTourTree[_levels];
        _adj = new Dictionary<int, HashSet<int>>[_levels];
        for (int i = 0; i < _levels; i++)
        {
            _forest[i] = new EulerTourTree(_n);
            _adj[i] = new Dictionary<int, HashSet<int>>();
        }
    }

    static (int, int) Key(int u, int v) => u < v ? (u, v) : (v, u);

    public bool Connected(int u, int v) => _forest[0].Connected(u, v);

    /// <summary>Opaque token identifying u's current connected component (stable until it changes).</summary>
    public object ComponentToken(int u) => _forest[0].ComponentToken(u);

    HashSet<int> Adj(int level, int x)
    {
        if (!_adj[level].TryGetValue(x, out var s)) { s = new HashSet<int>(); _adj[level][x] = s; }
        return s;
    }

    void AddNonTree(int level, int u, int v)
    {
        Adj(level, u).Add(v);
        Adj(level, v).Add(u);
        _forest[level].SetVertexMark(u, true);
        _forest[level].SetVertexMark(v, true);
    }

    void RemoveNonTree(int level, int u, int v)
    {
        var su = Adj(level, u); su.Remove(v);
        var sv = Adj(level, v); sv.Remove(u);
        _forest[level].SetVertexMark(u, su.Count > 0);
        _forest[level].SetVertexMark(v, sv.Count > 0);
    }

    public void Insert(int u, int v)
    {
        var key = Key(u, v);
        _edgeLevel[key] = 0;
        if (!Connected(u, v))
        {
            _forest[0].Link(u, v);
            _treeEdges.Add(key);
            // tree edge lives at level 0: marked only in forest[0] (its top level).
            _forest[0].SetEdgeMark(key.Item1, key.Item2, true);
        }
        else
        {
            AddNonTree(0, u, v);
        }
    }

    public void Delete(int u, int v)
    {
        var key = Key(u, v);
        int lvl = _edgeLevel[key];
        _edgeLevel.Remove(key);
        if (!_treeEdges.Contains(key))
        {
            RemoveNonTree(lvl, u, v);
            return;
        }
        // Tree edge: remove from F_0..F_lvl, then look for a replacement.
        _treeEdges.Remove(key);
        _forest[lvl].SetEdgeMark(key.Item1, key.Item2, false);
        for (int i = 0; i <= lvl; i++) _forest[i].Cut(u, v);
        Replace(u, v, lvl);
    }

    bool Replace(int u, int v, int lvl)
    {
        for (int i = lvl; i >= 0; i--)
        {
            // Smaller side is `a`; the replacement must cross to `b`'s side.
            int a = _forest[i].TreeSize(u) <= _forest[i].TreeSize(v) ? u : v;
            int b = a == u ? v : u;

            // Step 1: push the smaller tree's level-i tree edges up to i+1 (legal: |T_a| <= n/2^{i+1}).
            while (_forest[i].FindMarkedEdge(a) is var fe && fe != null)
            {
                var (x, y) = fe.Value;
                var ek = Key(x, y);
                _forest[i].SetEdgeMark(ek.Item1, ek.Item2, false);
                _forest[i + 1].Link(x, y);
                _forest[i + 1].SetEdgeMark(ek.Item1, ek.Item2, true);
                _edgeLevel[ek] = i + 1;
            }

            // Step 2: scan the smaller tree's level-i non-tree edges.
            while (_forest[i].FindMarkedVertex(a) is int x && x >= 0)
            {
                var list = Adj(i, x);
                while (list.Count > 0)
                {
                    int y = First(list);
                    if (_forest[i].Connected(y, b))
                    {
                        // Replacement found — promote (x,y) to a level-i tree edge in F_0..F_i.
                        RemoveNonTree(i, x, y);
                        var ek = Key(x, y);
                        for (int j = 0; j <= i; j++) _forest[j].Link(x, y);
                        _treeEdges.Add(ek);
                        _edgeLevel[ek] = i;
                        _forest[i].SetEdgeMark(ek.Item1, ek.Item2, true);
                        return true;
                    }
                    else
                    {
                        // Internal to T_a — push up to level i+1.
                        RemoveNonTree(i, x, y);
                        AddNonTree(i + 1, x, y);
                        _edgeLevel[Key(x, y)] = i + 1;
                    }
                }
            }
        }
        return false;
    }

    static int First(HashSet<int> s) { foreach (int x in s) return x; return -1; }
}
