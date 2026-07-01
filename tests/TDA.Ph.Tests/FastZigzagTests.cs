#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

using Maths.Topology;
namespace TDA.Ph.Tests;

/// <summary>FastZigzag (reduce-to-standard) must agree with the ZigzagBarcodeNaive oracle.</summary>
public sealed class FastZigzagTests
{
    static IEnumerable<(double, double, int, int, int)> Sig(Barcode bc) =>
        bc.Bars.Select(b => (b.Birth, b.Death, b.Dimension, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item3).ThenBy(x => x.Item1).ThenBy(x => x.Item2)
              .ThenBy(x => x.Item4).ThenBy(x => x.Item5);

    static void AssertParity(ZigzagFiltration f, int maxDim) =>
        Assert.Equal(Sig(ZigzagBarcodeNaive.Compute(f, maxDim)), Sig(FastZigzag.Compute(f, maxDim)));

    [Fact]
    public void AddThenDelete_cc()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Delete(0);
        AssertParity(f, 0);
    }

    [Fact]
    public void Taxonomy_co()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new[] { 0, 1 });
        AssertParity(f, 0);
    }

    [Fact]
    public void Taxonomy_oc()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new[] { 0, 1 });
        f.Delete(2); f.Delete(1);
        AssertParity(f, 0);
    }

    [Fact]
    public void Taxonomy_oo()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 });
        f.Delete(3); f.Add(5, new[] { 0, 1 });
        AssertParity(f, 0);
    }

    [Fact]
    public void TriangleFormsThenBreaks_H1()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 }); f.Add(5, new[] { 0, 2 });
        f.Delete(5);
        AssertParity(f, 1);
    }

    [Fact]
    public void DynamicGraphH0()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 });
        f.Delete(3);
        AssertParity(f, 0);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(7)]
    [InlineData(42)]
    [InlineData(99)]
    [InlineData(123)]
    public void RandomGraphZigzag_ParityWithOracle(int seed)
    {
        var rng = new Random(seed);
        int V = rng.Next(3, 7);
        var f = new ZigzagFiltration();

        // Vertices 0..V-1 added once (never deleted -> exercises non-empty endings / truncation).
        for (int v = 0; v < V; v++) f.Add(v, new int[0]);
        int nextId = V;

        // Dynamic edges with FRESH ids on each addition (non-repetitive).
        var edgePresent = new Dictionary<(int, int), int>();
        int steps = rng.Next(8, 22);
        for (int s = 0; s < steps; s++)
        {
            var addable = new List<(int, int)>();
            for (int u = 0; u < V; u++)
                for (int w = u + 1; w < V; w++)
                    if (!edgePresent.ContainsKey((u, w))) addable.Add((u, w));

            bool doAdd = edgePresent.Count == 0 || (addable.Count > 0 && rng.NextDouble() < 0.6);
            if (doAdd)
            {
                var e = addable[rng.Next(addable.Count)];
                int id = nextId++;
                edgePresent[e] = id;
                f.Add(id, new[] { e.Item1, e.Item2 });
            }
            else
            {
                var e = edgePresent.Keys.ElementAt(rng.Next(edgePresent.Count));
                f.Delete(edgePresent[e]);
                edgePresent.Remove(e);
            }
        }

        AssertParity(f, 1);
    }

    [Fact]
    public void Representatives_DoNotChangeBarcode()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new[] { 0, 1 });
        Assert.Equal(Sig(FastZigzag.Compute(f, 0)), Sig(FastZigzag.Compute(f, 0, representatives: true)));
    }

    [Fact]
    public void TriangleH1_LoopCycleRepresentative()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 }); f.Add(5, new[] { 0, 2 });
        f.Delete(5);

        var h1 = FastZigzag.Compute(f, 1, representatives: true).Bars.Single(b => b.Dimension == 1);

        Assert.NotNull(h1.Cycle);
        // The only H1 loop is the three edges e01=3, e12=4, e02=5 (original cell ids).
        Assert.Equal(new[] { 3, 4, 5 }, h1.Cycle!.OrderBy(x => x).ToArray());
    }
}
