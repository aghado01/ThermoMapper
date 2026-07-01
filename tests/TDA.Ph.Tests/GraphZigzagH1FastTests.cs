#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// Z5c Phase F — the near-linear H1 engine <see cref="GraphZigzagH1Fast"/> (dynamic-MSF bottleneck)
/// must agree with the slow Z5b oracle <see cref="GraphZigzagH1"/> and with Z1/Z2 on dimension 1,
/// across the hand zoo and random edge-/vertex-churn sweeps that birth and kill overlapping cycles.
/// </summary>
public sealed class GraphZigzagH1FastTests
{
    static IEnumerable<(double, double, int, int, int)> Sig1(Barcode bc) =>
        bc.Bars.Where(b => b.Dimension == 1)
              .Select(b => (b.Birth, b.Death, b.Dimension, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item1).ThenBy(x => x.Item2).ThenBy(x => x.Item4).ThenBy(x => x.Item5);

    static void AssertParity(ZigzagFiltration f)
    {
        var fast = Sig1(GraphZigzagH1Fast.Compute(f)).ToList();
        Assert.Equal(Sig1(GraphZigzagH1.Compute(f)), fast);           // vs slow Z5b
        Assert.Equal(Sig1(ZigzagBarcodeNaive.Compute(f, 1)), fast);   // vs Z1
        Assert.Equal(Sig1(FastZigzag.Compute(f, 1)), fast);          // vs Z2
    }

    [Fact]
    public void HandZoo()
    {
        var tri = new ZigzagFiltration();
        tri.Add(0, new int[0]); tri.Add(1, new int[0]); tri.Add(2, new int[0]);
        tri.Add(3, new[] { 0, 1 }); tri.Add(4, new[] { 1, 2 }); tri.Add(5, new[] { 0, 2 });
        tri.Delete(5);
        AssertParity(tri);

        var theta = new ZigzagFiltration();
        for (int v = 0; v < 4; v++) theta.Add(v, new int[0]);
        theta.Add(4, new[] { 0, 1 }); theta.Add(5, new[] { 1, 2 }); theta.Add(6, new[] { 0, 2 });
        theta.Add(7, new[] { 2, 3 }); theta.Add(8, new[] { 0, 3 });
        theta.Delete(6);   // earliest-birth pairing
        AssertParity(theta);
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
