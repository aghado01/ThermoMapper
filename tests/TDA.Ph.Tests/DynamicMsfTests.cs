#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using TDA.Ph.Dynamic;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// Z5c Phase C — <see cref="DynamicMsf"/> must match a Kruskal-rebuilt minimum spanning forest under
/// random weighted insert/delete (arbitrary distinct weights, so both the insertion swap and the
/// minimum-replacement-on-deletion paths are exercised): connectivity, and the max-edge-weight on the
/// MSF path.
/// </summary>
public sealed class DynamicMsfTests
{
    sealed class Brute
    {
        readonly int _n;
        public readonly Dictionary<(int, int), int> Edges = new();   // (min,max) -> weight

        public Brute(int n) => _n = n;

        Dictionary<int, List<(int to, int w)>> Msf()
        {
            var adj = new Dictionary<int, List<(int, int)>>();
            for (int i = 0; i < _n; i++) adj[i] = new List<(int, int)>();
            var p = Enumerable.Range(0, _n).ToArray();
            int Find(int x) { while (p[x] != x) x = p[x] = p[p[x]]; return x; }
            foreach (var (e, w) in Edges.OrderBy(kv => kv.Value))
            {
                int a = Find(e.Item1), b = Find(e.Item2);
                if (a == b) continue;
                p[a] = b;
                adj[e.Item1].Add((e.Item2, w));
                adj[e.Item2].Add((e.Item1, w));
            }
            return adj;
        }

        public int? PathMax(int u, int v)
        {
            if (u == v) return null;
            var adj = Msf();
            var prev = new Dictionary<int, (int from, int w)> { [u] = (-1, int.MinValue) };
            var q = new Queue<int>(); q.Enqueue(u);
            while (q.Count > 0)
            {
                int x = q.Dequeue();
                foreach (var (y, w) in adj[x])
                {
                    if (prev.ContainsKey(y)) continue;
                    prev[y] = (x, w);
                    if (y == v) { int mx = int.MinValue; for (int c = v; c != u; c = prev[c].from) mx = Math.Max(mx, prev[c].w); return mx; }
                    q.Enqueue(y);
                }
            }
            return null;
        }

        public bool Connected(int u, int v) => u == v || PathMax(u, v) != null;
    }

    [Theory]
    [InlineData(1)]
    [InlineData(7)]
    [InlineData(42)]
    [InlineData(99)]
    [InlineData(2024)]
    [InlineData(31337)]
    public void RandomWeightedInsertDelete_MatchesKruskal(int seed)
    {
        var rng = new Random(seed);
        int n = rng.Next(2, 22);
        var dmsf = new DynamicMsf(n);
        var brute = new Brute(n);
        var usedW = new HashSet<int>();

        for (int step = 0; step < 600; step++)
        {
            bool doInsert = brute.Edges.Count == 0 || (brute.Edges.Count < n * (n - 1) / 2 && rng.NextDouble() < 0.6);
            if (doInsert)
            {
                int u = rng.Next(n), v = rng.Next(n);
                if (u == v) continue;
                var key = (Math.Min(u, v), Math.Max(u, v));
                if (brute.Edges.ContainsKey(key)) continue;
                int w; do { w = rng.Next(1, 100000); } while (!usedW.Add(w));
                dmsf.Insert(u, v, w);
                brute.Edges[key] = w;
            }
            else
            {
                var key = brute.Edges.Keys.ElementAt(rng.Next(brute.Edges.Count));
                usedW.Remove(brute.Edges[key]);
                brute.Edges.Remove(key);
                dmsf.Delete(key.Item1, key.Item2);
            }

            for (int q = 0; q < 8; q++)
            {
                int a = rng.Next(n), b = rng.Next(n);
                bool conn = brute.Connected(a, b);
                Assert.Equal(conn, dmsf.Connected(a, b));
                if (conn && a != b) Assert.Equal(brute.PathMax(a, b), dmsf.PathMax(a, b));
            }
        }
    }
}
