using System;
using System.Collections.Generic;
using System.Linq;

namespace VizCore.Tests;

/// <summary>
/// A closure that won't cross a thin cut between two large sides, for the test in
/// <c>issues/spc/bayesian-clustering-design.md</c> ("A closure that won't cross a thin cut"). It reads only the edges
/// above the cut, and it lives with the tests until that question says it should be anything more.
/// The edges build a cut tree (Gomory &amp; Hu 1961, by Gusfield 1990's flows), which holds a minimum cut for every pair
/// of nodes. A cut weighs one per edge plus a fraction of its edges' values, so it has the fewest edges and, among those,
/// the weakest. <see cref="Partition"/> removes tree edges lightest first. Taking k = 0 gives closure, and m = 1 gives
/// the connected parts of the (k + 1)-edge-connected classes, which only coarsen as edges are added.
/// </summary>
internal sealed class ThinCutClosure
{
    // A cut of c edges weighs c units and less than c / 16 more, so counts order first for anything below 16 edges.
    private const int Unit = 1 << 20, ValueScale = 1 << 16;

    private readonly int _nodeCount;
    private readonly int[] _sources, _targets;
    private readonly int _thick;

    // Arc 2e runs from source to target and arc 2e + 1 back, each holding the edge's weight.
    private readonly int[] _arcStart, _arcs, _head, _capacity, _residual;
    private readonly int[] _level, _next, _path, _queue, _reached;
    private readonly List<int> _touched = [];
    private int _round;

    // The tree, rooted per component: each node's parent and the weight of the cut to it, and a preorder.
    private readonly int[] _parent, _weight, _order, _first, _size, _root, _lightest;

    /// <summary>The edges of <paramref name="column"/> above <paramref name="cut"/>, as closure reads them.</summary>
    internal static ThinCutClosure Above(int nodeCount, int[] sources, int[] targets, double[] column, double cut, int maxCut)
    {
        int[] kept = Enumerable.Range(0, column.Length).Where(e => column[e] > cut).ToArray();
        return new ThinCutClosure(nodeCount, kept.Select(e => sources[e]).ToArray(), kept.Select(e => targets[e]).ToArray(), kept.Select(e => column[e]).ToArray(), maxCut);
    }

    /// <param name="maxCut">The most edges a cut read here can hold. A pair no cut that small parts stays together.</param>
    internal ThinCutClosure(int nodeCount, int[] sources, int[] targets, double[] values, int maxCut)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(maxCut);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(maxCut, 14);

        _nodeCount = nodeCount;
        _sources = sources;
        _targets = targets;
        _thick = (maxCut + 1) * Unit;

        int edges = sources.Length;
        _arcStart = new int[nodeCount + 1];
        for (int e = 0; e < edges; e++)
        {
            _arcStart[sources[e] + 1]++;
            _arcStart[targets[e] + 1]++;
        }
        for (int v = 0; v < nodeCount; v++) _arcStart[v + 1] += _arcStart[v];

        _arcs = new int[2 * edges];
        _head = new int[2 * edges];
        _capacity = new int[2 * edges];
        int[] fill = (int[])_arcStart.Clone();
        for (int e = 0; e < edges; e++)
        {
            int weight = EdgeWeight(values[e]);
            _arcs[fill[sources[e]]++] = 2 * e;
            _arcs[fill[targets[e]]++] = 2 * e + 1;
            _head[2 * e] = targets[e];
            _head[2 * e + 1] = sources[e];
            _capacity[2 * e] = _capacity[2 * e + 1] = weight;
        }

        _level = new int[nodeCount];
        _next = new int[nodeCount];
        _path = new int[nodeCount];
        _queue = new int[nodeCount];
        _reached = new int[nodeCount];

        // A bridge is its own tree edge, and no flow between two nodes of one block uses it.
        bool[] bridge = Bridges();
        var links = new List<(int A, int B, int Weight)>();
        for (int e = 0; e < edges; e++)
        {
            if (!bridge[e]) continue;
            links.Add((sources[e], targets[e], Math.Min(_capacity[2 * e], _thick)));
            _capacity[2 * e] = _capacity[2 * e + 1] = 0;
        }
        _residual = (int[])_capacity.Clone();

        var above = new int[nodeCount];
        var cost = new int[nodeCount];
        foreach (int[] block in Blocks(bridge))
            CutTree(block, above, cost, links);

        (_parent, _weight, _order, _first, _size, _root) = Rooted(links);
        _lightest = Enumerable.Range(0, nodeCount).Where(v => _parent[v] != v && _weight[v] < _thick)
            .OrderBy(v => _weight[v]).ThenBy(v => v).ToArray();
    }

    /// <summary>What one edge adds to a cut's weight.</summary>
    internal static int EdgeWeight(double value) => Unit + (int)Math.Round(value * ValueScale);

    /// <summary>The weight at which a pair counts as joined by more than the largest cut read here.</summary>
    internal int Thick => _thick;

    /// <summary>The tree: each node's parent, itself at a root, and the weight of the cut between them, at most <see cref="Thick"/>.</summary>
    internal (int[] Parent, int[] Weight) Tree => (_parent, _weight);

    /// <summary>
    /// Clusters that no thin cut between two large sides runs through: going lightest first over the tree's cuts of at
    /// most <paramref name="k"/> edges, each is removed if both its sides within the current piece hold at least
    /// <paramref name="m"/> points. A cluster is then a piece's part connected above the cut.
    /// </summary>
    internal int[] Partition(int k, int m)
    {
        int[] size = (int[])_size.Clone();
        int[] top = (int[])_root.Clone();
        int limit = (k + 1) * Unit;

        foreach (int v in _lightest)
        {
            if (_weight[v] >= limit) break;
            int piece = top[v], side = size[v];
            if (side < m || size[piece] - side < m) continue;

            for (int a = _parent[v]; ; a = _parent[a])
            {
                size[a] -= side;
                if (a == piece) break;
            }
            for (int i = _first[v]; i < _first[v] + _size[v]; i++)
                if (top[_order[i]] == piece) top[_order[i]] = v;
        }

        int[] set = Enumerable.Range(0, _nodeCount).ToArray();
        for (int e = 0; e < _sources.Length; e++)
            if (top[_sources[e]] == top[_targets[e]])
                set[Find(set, _sources[e])] = Find(set, _targets[e]);

        var dense = new int[_nodeCount];
        Array.Fill(dense, -1);
        var labels = new int[_nodeCount];
        int count = 0;
        for (int v = 0; v < _nodeCount; v++)
        {
            int r = Find(set, v);
            if (dense[r] < 0) dense[r] = count++;
            labels[v] = dense[r];
        }
        return labels;
    }

    private static int Find(int[] set, int x)
    {
        while (set[x] != x)
        {
            set[x] = set[set[x]];
            x = set[x];
        }
        return x;
    }

    private bool[] Bridges()
    {
        var bridge = new bool[_sources.Length];
        var order = new int[_nodeCount];
        Array.Fill(order, -1);
        var low = new int[_nodeCount];
        var entry = new int[_nodeCount];
        int[] cursor = (int[])_arcStart.Clone();
        var stack = new int[_nodeCount];
        int clock = 0;

        for (int start = 0; start < _nodeCount; start++)
        {
            if (order[start] >= 0) continue;
            int depth = 0;
            stack[depth++] = start;
            order[start] = low[start] = clock++;
            entry[start] = -1;
            while (depth > 0)
            {
                int v = stack[depth - 1];
                if (cursor[v] < _arcStart[v + 1])
                {
                    int arc = _arcs[cursor[v]++];
                    if (arc == (entry[v] ^ 1)) continue;
                    int w = _head[arc];
                    if (order[w] < 0)
                    {
                        order[w] = low[w] = clock++;
                        entry[w] = arc;
                        stack[depth++] = w;
                    }
                    else low[v] = Math.Min(low[v], order[w]);
                }
                else
                {
                    depth--;
                    if (depth == 0) continue;
                    int up = stack[depth - 1];
                    low[up] = Math.Min(low[up], low[v]);
                    if (low[v] > order[up]) bridge[entry[v] / 2] = true;
                }
            }
        }
        return bridge;
    }

    private IEnumerable<int[]> Blocks(bool[] bridge)
    {
        int[] set = Enumerable.Range(0, _nodeCount).ToArray();
        for (int e = 0; e < _sources.Length; e++)
            if (!bridge[e])
                set[Find(set, _sources[e])] = Find(set, _targets[e]);
        return Enumerable.Range(0, _nodeCount).GroupBy(v => Find(set, v)).Where(g => g.Count() > 1).Select(g => g.ToArray());
    }

    // Gusfield's cut tree over one block. A pair too thick to read keeps its tree edge at the thick weight, and since no
    // cut read here can part that pair, the rest of the run goes on as if the two were contracted.
    private void CutTree(int[] block, int[] parent, int[] weight, List<(int A, int B, int Weight)> links)
    {
        int root = block[0];
        foreach (int v in block) parent[v] = root;

        for (int b = 1; b < block.Length; b++)
        {
            int s = block[b], t = parent[s];
            int flow = MaxFlow(s, t);
            if (flow >= _thick)
            {
                weight[s] = _thick;
                continue;
            }

            foreach (int i in block)
                if (i != s && _reached[i] == _round && parent[i] == t)
                    parent[i] = s;

            if (_reached[parent[t]] == _round)
            {
                parent[s] = parent[t];
                parent[t] = s;
                weight[s] = weight[t];
                weight[t] = flow;
            }
            else weight[s] = flow;
        }

        foreach (int v in block)
            if (v != root)
                links.Add((v, parent[v], weight[v]));
    }

    // Dinic's flow, stopping at the thick weight. Once it returns under that, _reached marks the side of s.
    private int MaxFlow(int s, int t)
    {
        int flow = 0;
        _touched.Clear();
        while (flow < _thick && Levels(s, t))
        {
            while (flow < _thick)
            {
                int pushed = Augment(s, t);
                if (pushed == 0) break;
                flow += pushed;
            }
        }
        foreach (int arc in _touched)
        {
            _residual[arc] = _capacity[arc];
            _residual[arc ^ 1] = _capacity[arc ^ 1];
        }
        return flow;
    }

    private bool Levels(int s, int t)
    {
        _round++;
        int read = 0, write = 0;
        _queue[write++] = s;
        _reached[s] = _round;
        _level[s] = 0;
        _next[s] = _arcStart[s];
        while (read < write)
        {
            int v = _queue[read++];
            for (int i = _arcStart[v]; i < _arcStart[v + 1]; i++)
            {
                int arc = _arcs[i], w = _head[arc];
                if (_residual[arc] <= 0 || _reached[w] == _round) continue;
                _reached[w] = _round;
                _level[w] = _level[v] + 1;
                _next[w] = _arcStart[w];
                if (w == t) return true;
                _queue[write++] = w;
            }
        }
        return false;
    }

    private int Augment(int s, int t)
    {
        int depth = 0, v = s;
        while (v != t)
        {
            bool moved = false;
            for (; _next[v] < _arcStart[v + 1]; _next[v]++)
            {
                int arc = _arcs[_next[v]], w = _head[arc];
                if (_residual[arc] > 0 && _reached[w] == _round && _level[w] == _level[v] + 1)
                {
                    _path[depth++] = arc;
                    v = w;
                    moved = true;
                    break;
                }
            }
            if (moved) continue;
            if (depth == 0) return 0;

            // A dead end: nothing comes back here this phase.
            _level[v] = -1;
            v = _head[_path[--depth] ^ 1];
            _next[v]++;
        }

        int pushed = int.MaxValue;
        for (int i = 0; i < depth; i++) pushed = Math.Min(pushed, _residual[_path[i]]);
        for (int i = 0; i < depth; i++)
        {
            _residual[_path[i]] -= pushed;
            _residual[_path[i] ^ 1] += pushed;
            _touched.Add(_path[i]);
        }
        return pushed;
    }

    private (int[] Parent, int[] Weight, int[] Order, int[] First, int[] Size, int[] Root) Rooted(List<(int A, int B, int Weight)> links)
    {
        var near = new List<(int To, int Weight)>[_nodeCount];
        foreach ((int a, int b, int w) in links)
        {
            (near[a] ??= []).Add((b, w));
            (near[b] ??= []).Add((a, w));
        }

        var parent = new int[_nodeCount];
        Array.Fill(parent, -1);
        var weight = new int[_nodeCount];
        var order = new int[_nodeCount];
        var first = new int[_nodeCount];
        var root = new int[_nodeCount];
        var stack = new int[_nodeCount];
        int count = 0;

        for (int start = 0; start < _nodeCount; start++)
        {
            if (parent[start] >= 0) continue;
            parent[start] = start;
            root[start] = start;
            int depth = 0;
            stack[depth++] = start;
            while (depth > 0)
            {
                int v = stack[--depth];
                first[v] = count;
                order[count++] = v;
                if (near[v] == null) continue;
                foreach ((int to, int w) in near[v])
                {
                    if (parent[to] >= 0) continue;
                    parent[to] = v;
                    weight[to] = w;
                    root[to] = start;
                    stack[depth++] = to;
                }
            }
        }

        var size = new int[_nodeCount];
        Array.Fill(size, 1);
        for (int i = _nodeCount - 1; i >= 0; i--)
            if (parent[order[i]] != order[i])
                size[parent[order[i]]] += size[order[i]];

        return (parent, weight, order, first, size, root);
    }
}
