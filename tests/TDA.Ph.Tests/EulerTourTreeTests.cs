#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using TDA.Ph.Dynamic;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// Z5c Phase A — <see cref="EulerTourTree"/> must behave exactly like a brute-force forest under random
/// link/cut: connectivity, tree size, and the two finger-marks (markable vertices / tree-edge arcs).
/// </summary>
public sealed class EulerTourTreeTests
{
    // Brute-force forest over the same tree-edge set.
    sealed class RefForest
    {
        readonly int _n;
        public readonly HashSet<(int, int)> Edges = new();   // canonical (min,max)
        public RefForest(int n) => _n = n;

        int[] Comp()
        {
            var p = Enumerable.Range(0, _n).ToArray();
            int Find(int x) { while (p[x] != x) x = p[x] = p[p[x]]; return x; }
            foreach (var (a, b) in Edges) p[Find(a)] = Find(b);
            var c = new int[_n];
            for (int i = 0; i < _n; i++) c[i] = Find(i);
            return c;
        }

        public bool Connected(int u, int v) { var c = Comp(); return c[u] == c[v]; }
        public int Size(int v) { var c = Comp(); return Enumerable.Range(0, _n).Count(i => c[i] == c[v]); }
        public IEnumerable<int> Component(int v) { var c = Comp(); return Enumerable.Range(0, _n).Where(i => c[i] == c[v]); }
    }

    [Theory]
    [InlineData(1)]
    [InlineData(7)]
    [InlineData(42)]
    [InlineData(2024)]
    [InlineData(31337)]
    public void RandomLinkCut_MatchesBruteForce(int seed)
    {
        var rng = new Random(seed);
        int n = rng.Next(2, 30);
        var ett = new EulerTourTree(n);
        var fref = new RefForest(n);

        for (int step = 0; step < 400; step++)
        {
            bool doLink = fref.Edges.Count == 0 || (fref.Edges.Count < n - 1 && rng.NextDouble() < 0.6);
            if (doLink)
            {
                // pick u,v in different components
                int u = rng.Next(n), v = rng.Next(n);
                if (fref.Connected(u, v)) continue;
                ett.Link(u, v);
                fref.Edges.Add((Math.Min(u, v), Math.Max(u, v)));
            }
            else
            {
                var e = fref.Edges.ElementAt(rng.Next(fref.Edges.Count));
                ett.Cut(e.Item1, e.Item2);
                fref.Edges.Remove(e);
            }

            // connectivity + size on random probes
            for (int q = 0; q < 6; q++)
            {
                int a = rng.Next(n), b = rng.Next(n);
                Assert.Equal(fref.Connected(a, b), ett.Connected(a, b));
                Assert.Equal(fref.Size(a), ett.TreeSize(a));
            }

            // mark round: set a fresh random subset, verify finger-finds, then clear.
            var markedV = new HashSet<int>();
            for (int v = 0; v < n; v++) if (rng.NextDouble() < 0.2) { ett.SetVertexMark(v, true); markedV.Add(v); }
            var markedE = new HashSet<(int, int)>();
            foreach (var e in fref.Edges) if (rng.NextDouble() < 0.3) { ett.SetEdgeMark(e.Item1, e.Item2, true); markedE.Add(e); }

            for (int v = 0; v < n; v++)
            {
                var comp = fref.Component(v).ToHashSet();
                bool anyV = comp.Any(markedV.Contains);
                int found = ett.FindMarkedVertex(v);
                if (!anyV) Assert.Equal(-1, found);
                else { Assert.NotEqual(-1, found); Assert.Contains(found, markedV); Assert.Contains(found, comp); }

                bool anyE = markedE.Any(e => comp.Contains(e.Item1));
                var fe = ett.FindMarkedEdge(v);
                if (!anyE) Assert.Null(fe);
                else
                {
                    Assert.NotNull(fe);
                    var canon = (Math.Min(fe!.Value.U, fe.Value.V), Math.Max(fe.Value.U, fe.Value.V));
                    Assert.Contains(canon, markedE);
                    Assert.Contains(fe.Value.U, comp);
                }
            }

            foreach (int v in markedV) ett.SetVertexMark(v, false);
            foreach (var e in markedE) ett.SetEdgeMark(e.Item1, e.Item2, false);
        }
    }
}
