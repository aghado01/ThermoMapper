#nullable enable
using System;
using System.Collections.Generic;

namespace TDA.Ph.Dynamic;

/// <summary>
/// Fully-dynamic minimum spanning forest with a max-edge-weight path query — the structure §4.1 of
/// Dey–Hou uses for H1 graph-zigzag (Z5c). Distinct integer edge weights (the zigzag's add-steps).
/// <para>The MSF itself is held twice: an <see cref="EulerTourTree"/> for connectivity / tree size /
/// enumeration of the smaller side, and a <see cref="LinkCutTree"/> for <see cref="PathMax"/>. Non-tree
/// edges sit in per-vertex weighted adjacency; vertices incident to one are ETT-marked. <see cref="Insert"/>
/// adds a tree edge if it joins two components, else compares the new weight against the heaviest edge
/// on the cycle and swaps if lighter (kept general, though the zigzag inserts in increasing weight so a
/// connecting edge is always non-tree). <see cref="Delete"/> of a tree edge reconnects with the
/// <i>minimum</i>-weight crossing non-tree edge, found by scanning the smaller side's marked vertices.
/// Pure, integer-id; internal to <c>TDA.Ph</c>.</para>
/// </summary>
internal sealed class DynamicMsf
{
    readonly EulerTourTree _tree;
    readonly LinkCutTree _lct;
    readonly Dictionary<int, Dictionary<int, int>> _nonTree = new();   // x -> (y -> weight)
    readonly HashSet<(int, int)> _treeEdges = new();
    readonly Dictionary<(int, int), int> _weight = new();
    readonly Dictionary<int, (int U, int V)> _byWeight = new();         // tree-edge weight -> endpoints

    public DynamicMsf(int n)
    {
        _tree = new EulerTourTree(Math.Max(1, n));
        _lct = new LinkCutTree(Math.Max(1, n));
    }

    static (int, int) Key(int u, int v) => u < v ? (u, v) : (v, u);

    public bool Connected(int u, int v) => _tree.Connected(u, v);

    /// <summary>Max edge weight on the u–v MSF path (caller ensures they are connected).</summary>
    public int PathMax(int u, int v) => _lct.PathMax(u, v);

    /// <summary>Weights of the edges on the u–v MSF path (caller ensures they are connected). With the
    /// zigzag's distinct add-step weights, each maps back to a unique edge.</summary>
    public List<int> PathEdgeWeights(int u, int v) => _lct.PathEdgeWeights(u, v);

    Dictionary<int, int> Bag(int x)
    {
        if (!_nonTree.TryGetValue(x, out var d)) { d = new Dictionary<int, int>(); _nonTree[x] = d; }
        return d;
    }

    void Mark(int x) => _tree.SetVertexMark(x, _nonTree.TryGetValue(x, out var d) && d.Count > 0);

    void AddNonTree(int u, int v, int w) { Bag(u)[v] = w; Bag(v)[u] = w; Mark(u); Mark(v); }

    void RemoveNonTree(int u, int v)
    {
        if (_nonTree.TryGetValue(u, out var du)) du.Remove(v);
        if (_nonTree.TryGetValue(v, out var dv)) dv.Remove(u);
        Mark(u); Mark(v);
    }

    void AddTree(int u, int v, int w)
    {
        _tree.Link(u, v); _lct.LinkEdge(u, v, w);
        _treeEdges.Add(Key(u, v)); _byWeight[w] = (u, v);
    }

    void RemoveTree(int u, int v, int w)
    {
        _tree.Cut(u, v); _lct.CutEdge(u, v);
        _treeEdges.Remove(Key(u, v)); _byWeight.Remove(w);
    }

    public void Insert(int u, int v, int w)
    {
        _weight[Key(u, v)] = w;
        if (!_tree.Connected(u, v)) { AddTree(u, v, w); return; }
        int heaviest = _lct.PathMax(u, v);
        if (w < heaviest)
        {
            var (mu, mv) = _byWeight[heaviest];
            RemoveTree(mu, mv, heaviest);
            AddNonTree(mu, mv, heaviest);
            AddTree(u, v, w);
        }
        else AddNonTree(u, v, w);
    }

    public void Delete(int u, int v)
    {
        var k = Key(u, v);
        int w = _weight[k]; _weight.Remove(k);
        if (!_treeEdges.Contains(k)) { RemoveNonTree(u, v); return; }

        RemoveTree(u, v, w);
        // Reconnect with the minimum-weight non-tree edge crossing the cut (scan the smaller side).
        int a = _tree.TreeSize(u) <= _tree.TreeSize(v) ? u : v;
        int bestU = -1, bestV = -1, bestW = int.MaxValue;
        foreach (int x in _tree.CollectMarkedVertices(a))
        {
            if (!_nonTree.TryGetValue(x, out var dx)) continue;
            foreach (var kv in dx)
                if (!_tree.Connected(x, kv.Key) && kv.Value < bestW) { bestW = kv.Value; bestU = x; bestV = kv.Key; }
        }
        if (bestU >= 0)
        {
            RemoveNonTree(bestU, bestV);
            AddTree(bestU, bestV, bestW);
        }
    }
}
