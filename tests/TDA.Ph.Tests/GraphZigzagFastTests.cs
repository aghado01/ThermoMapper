#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// Z5c Phase E — the near-linear H0 engine <see cref="GraphZigzagFast"/> must agree with the slow
/// Z5a oracle <see cref="GraphZigzag"/> and with Z1/Z2 on dimension 0, across the hand zoo and the
/// random edge-/vertex-churn sweeps that stress merge gluing and nested-split departures.
/// </summary>
public sealed class GraphZigzagFastTests
{
    static IEnumerable<(double, double, int, int, int)> Sig(Barcode bc) =>
        bc.Bars.Where(b => b.Dimension == 0)
              .Select(b => (b.Birth, b.Death, b.Dimension, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item1).ThenBy(x => x.Item2).ThenBy(x => x.Item4).ThenBy(x => x.Item5);

    static void AssertParity(ZigzagFiltration f)
    {
        var fast = Sig(GraphZigzagFast.Compute(f, 0)).ToList();
        Assert.Equal(Sig(GraphZigzag.Compute(f, 0, algorithm: GraphZigzagAlgorithm.Reference)), fast);  // vs reference Z5a
        Assert.Equal(Sig(ZigzagBarcodeNaive.Compute(f, 0)), fast);     // vs Z1
        Assert.Equal(Sig(FastZigzag.Compute(f, 0)), fast);            // vs Z2
    }

    // All-dimension signature, for the unified-entry parity (H0 + folded H1).
    static IEnumerable<(double, double, int, int, int)> SigAll(Barcode bc) =>
        bc.Bars.Select(b => (b.Birth, b.Death, b.Dimension, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item3).ThenBy(x => x.Item1).ThenBy(x => x.Item2).ThenBy(x => x.Item4).ThenBy(x => x.Item5);

    static void AssertUnifiedEntry(ZigzagFiltration f)
    {
        // One entry, two interchangeable algorithms, both dimensions — agreeing with each other and Z1/Z2.
        var reference = SigAll(GraphZigzag.Compute(f, 1, algorithm: GraphZigzagAlgorithm.Reference)).ToList();
        var fast = SigAll(GraphZigzag.Compute(f, 1, algorithm: GraphZigzagAlgorithm.Fast)).ToList();
        Assert.Equal(reference, fast);
        Assert.Equal(SigAll(ZigzagBarcodeNaive.Compute(f, 1)), reference);
        Assert.Equal(SigAll(FastZigzag.Compute(f, 1)), reference);
    }

    [Theory]
    [InlineData(1)] [InlineData(7)] [InlineData(42)] [InlineData(99)] [InlineData(123)]
    public void UnifiedEntry_BothAlgorithms_AgreeAcrossDimensions(int seed)
    {
        AssertUnifiedEntry(BuildEdgeChurn(seed));
        AssertUnifiedEntry(BuildVertexChurn(seed));
    }

    [Fact]
    public void HandZoo()
    {
        // entrance/merge/split/departure mix
        var a = new ZigzagFiltration();
        a.Add(0, new int[0]); a.Add(1, new int[0]); a.Add(2, new int[0]);
        a.Add(3, new[] { 0, 1 }); a.Add(4, new[] { 1, 2 }); a.Delete(3);
        AssertParity(a);

        var b = new ZigzagFiltration();
        b.Add(0, new int[0]); b.Add(1, new int[0]); b.Add(2, new int[0]); b.Add(3, new int[0]);
        b.Add(4, new[] { 0, 1 }); b.Add(5, new[] { 1, 2 }); b.Add(6, new[] { 2, 3 });
        b.Delete(5); b.Delete(6); b.Delete(2);   // nested split then departure
        AssertParity(b);
    }

    static ZigzagFiltration BuildEdgeChurn(int seed)
    {
        var rng = new Random(seed);
        int V = rng.Next(3, 7);
        var f = new ZigzagFiltration();
        for (int v = 0; v < V; v++) f.Add(v, new int[0]);
        int nextId = V;
        var edges = new Dictionary<(int, int), int>();
        int steps = rng.Next(8, 22);
        for (int s = 0; s < steps; s++)
        {
            var addable = new List<(int, int)>();
            for (int u = 0; u < V; u++)
                for (int w = u + 1; w < V; w++)
                    if (!edges.ContainsKey((u, w))) addable.Add((u, w));
            bool doAdd = edges.Count == 0 || (addable.Count > 0 && rng.NextDouble() < 0.6);
            if (doAdd)
            {
                var e = addable[rng.Next(addable.Count)];
                int eid = nextId++;
                edges[e] = eid;
                f.Add(eid, new[] { e.Item1, e.Item2 });
            }
            else
            {
                var e = edges.Keys.ElementAt(rng.Next(edges.Count));
                f.Delete(edges[e]);
                edges.Remove(e);
            }
        }
        return f;
    }

    static ZigzagFiltration BuildVertexChurn(int seed)
    {
        var rng = new Random(seed);
        var f = new ZigzagFiltration();
        var present = new HashSet<int>();
        var edges = new Dictionary<(int, int), int>();
        var incident = new Dictionary<int, int>();
        int nextId = 0;
        int initV = rng.Next(3, 6);
        for (int i = 0; i < initV; i++) { int v = nextId++; f.Add(v, new int[0]); present.Add(v); incident[v] = 0; }
        int steps = rng.Next(12, 30);
        for (int s = 0; s < steps; s++)
        {
            var verts = present.ToList();
            var addable = new List<(int, int)>();
            for (int i = 0; i < verts.Count; i++)
                for (int j = i + 1; j < verts.Count; j++)
                {
                    var e = (Math.Min(verts[i], verts[j]), Math.Max(verts[i], verts[j]));
                    if (!edges.ContainsKey(e)) addable.Add(e);
                }
            var isolated = verts.Where(v => incident[v] == 0).ToList();
            double r = rng.NextDouble();
            if (r < 0.45 && addable.Count > 0)
            {
                var e = addable[rng.Next(addable.Count)];
                int eid = nextId++;
                edges[e] = eid;
                incident[e.Item1]++; incident[e.Item2]++;
                f.Add(eid, new[] { e.Item1, e.Item2 });
            }
            else if (r < 0.70 && edges.Count > 0)
            {
                var e = edges.Keys.ElementAt(rng.Next(edges.Count));
                f.Delete(edges[e]);
                edges.Remove(e);
                incident[e.Item1]--; incident[e.Item2]--;
            }
            else if (r < 0.85 && isolated.Count > 0)
            {
                int v = isolated[rng.Next(isolated.Count)];
                f.Delete(v);
                present.Remove(v);
                incident.Remove(v);
            }
            else
            {
                int v = nextId++;
                f.Add(v, new int[0]);
                present.Add(v);
                incident[v] = 0;
            }
        }
        return f;
    }

    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(7)] [InlineData(42)] [InlineData(99)]
    [InlineData(123)] [InlineData(2024)] [InlineData(31337)] [InlineData(8675309)] [InlineData(13)]
    public void RandomEdgeChurn_FastEqualsSlowAndOracles(int seed) => AssertParity(BuildEdgeChurn(seed));

    [Theory]
    [InlineData(1)] [InlineData(3)] [InlineData(8)] [InlineData(21)] [InlineData(55)]
    [InlineData(144)] [InlineData(7)] [InlineData(99)] [InlineData(2024)] [InlineData(512)]
    public void RandomVertexChurn_FastEqualsSlowAndOracles(int seed) => AssertParity(BuildVertexChurn(seed));
}
