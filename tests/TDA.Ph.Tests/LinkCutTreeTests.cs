#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using TDA.Ph.Dynamic;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// Z5c Phase C (a) — <see cref="LinkCutTree"/> connectivity and path-max-edge-weight must match a
/// brute-force forest under random weighted link/cut.
/// </summary>
public sealed class LinkCutTreeTests
{
    sealed class RefForest
    {
        readonly int _n;
        public readonly Dictionary<(int, int), int> W = new();   // canonical (min,max) -> weight
        readonly Dictionary<int, List<int>> _adj = new();
        public RefForest(int n) { _n = n; for (int i = 0; i < n; i++) _adj[i] = new List<int>(); }

        public void Add(int u, int v, int w) { W[(Math.Min(u, v), Math.Max(u, v))] = w; _adj[u].Add(v); _adj[v].Add(u); }
        public void Remove(int u, int v) { W.Remove((Math.Min(u, v), Math.Max(u, v))); _adj[u].Remove(v); _adj[v].Remove(u); }

        public bool Connected(int u, int v) => PathMax(u, v) != null || u == v;

        // max edge weight on the unique u-v tree path, or null if not connected.
        public int? PathMax(int u, int v)
        {
            if (u == v) return null;
            var prev = new Dictionary<int, (int from, int w)> { [u] = (-1, int.MinValue) };
            var q = new Queue<int>(); q.Enqueue(u);
            while (q.Count > 0)
            {
                int x = q.Dequeue();
                foreach (int y in _adj[x])
                {
                    if (prev.ContainsKey(y)) continue;
                    prev[y] = (x, W[(Math.Min(x, y), Math.Max(x, y))]);
                    if (y == v)
                    {
                        int mx = int.MinValue;
                        for (int c = v; c != u; c = prev[c].from) mx = Math.Max(mx, prev[c].w);
                        return mx;
                    }
                    q.Enqueue(y);
                }
            }
            return null;
        }
    }

    [Theory]
    [InlineData(1)]
    [InlineData(7)]
    [InlineData(42)]
    [InlineData(2024)]
    [InlineData(31337)]
    public void RandomLinkCut_PathMax_MatchesBruteForce(int seed)
    {
        var rng = new Random(seed);
        int n = rng.Next(2, 26);
        var lct = new LinkCutTree(n);
        var fref = new RefForest(n);
        int nextW = 1;

        for (int step = 0; step < 500; step++)
        {
            bool doLink = fref.W.Count == 0 || (fref.W.Count < n - 1 && rng.NextDouble() < 0.6);
            if (doLink)
            {
                int u = rng.Next(n), v = rng.Next(n);
                if (u == v || fref.Connected(u, v)) continue;
                int w = nextW++;
                lct.LinkEdge(u, v, w);
                fref.Add(u, v, w);
            }
            else
            {
                var e = fref.W.Keys.ElementAt(rng.Next(fref.W.Count));
                lct.CutEdge(e.Item1, e.Item2);
                fref.Remove(e.Item1, e.Item2);
            }

            for (int q = 0; q < 8; q++)
            {
                int a = rng.Next(n), b = rng.Next(n);
                bool conn = fref.Connected(a, b);
                Assert.Equal(conn, lct.Connected(a, b));
                if (conn && a != b) Assert.Equal(fref.PathMax(a, b), lct.PathMax(a, b));
            }
        }
    }
}
