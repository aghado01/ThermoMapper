#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using TDA.Ph.Dynamic;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// Z5c Phase B — <see cref="DynamicConnectivity"/> (HDT) must answer connectivity exactly like a DSU
/// rebuilt over the live edge set, under long random insert/delete sequences that force many
/// replacement searches (tree-edge deletions reconnected from non-tree edges).
/// </summary>
public sealed class DynamicConnectivityTests
{
    static int[] Components(int n, IEnumerable<(int, int)> edges)
    {
        var p = Enumerable.Range(0, n).ToArray();
        int Find(int x) { while (p[x] != x) x = p[x] = p[p[x]]; return x; }
        foreach (var (a, b) in edges) p[Find(a)] = Find(b);
        var c = new int[n];
        for (int i = 0; i < n; i++) c[i] = Find(i);
        return c;
    }

    [Theory]
    [InlineData(1)]
    [InlineData(7)]
    [InlineData(42)]
    [InlineData(99)]
    [InlineData(2024)]
    [InlineData(31337)]
    public void RandomInsertDelete_MatchesDsu(int seed)
    {
        var rng = new Random(seed);
        int n = rng.Next(2, 28);
        var dc = new DynamicConnectivity(n);
        var present = new HashSet<(int, int)>();

        for (int step = 0; step < 700; step++)
        {
            bool doInsert = present.Count == 0 || (present.Count < n * (n - 1) / 2 && rng.NextDouble() < 0.6);
            if (doInsert)
            {
                int u = rng.Next(n), v = rng.Next(n);
                if (u == v) continue;
                var key = (Math.Min(u, v), Math.Max(u, v));
                if (!present.Add(key)) continue;
                dc.Insert(u, v);
            }
            else
            {
                var key = present.ElementAt(rng.Next(present.Count));
                present.Remove(key);
                dc.Delete(key.Item1, key.Item2);
            }

            var c = Components(n, present);
            for (int q = 0; q < 8; q++)
            {
                int a = rng.Next(n), b = rng.Next(n);
                Assert.Equal(c[a] == c[b], dc.Connected(a, b));
            }
        }
    }
}
