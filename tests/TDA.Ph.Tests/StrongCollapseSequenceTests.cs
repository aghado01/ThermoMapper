#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

using Maths.Topology;
namespace TDA.Ph.Tests;

/// <summary>
/// §4 of 1809.10945, Theorem 2: strong-collapse core-assembly preserves the zigzag barcode.
/// StrongCollapseZigzag (collapse → general oracle) must equal the naive oracle on the original.
/// </summary>
public sealed class StrongCollapseSequenceTests
{
    static IEnumerable<(double, double, int, int, int)> Sig(Barcode bc) =>
        bc.Bars.Select(b => (b.Birth, b.Death, b.Dimension, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item3).ThenBy(x => x.Item1).ThenBy(x => x.Item2)
              .ThenBy(x => x.Item4).ThenBy(x => x.Item5);

    static void AssertCollapseInvariant(ZigzagFiltration f, int maxDim) =>
        Assert.Equal(Sig(ZigzagBarcodeNaive.Compute(f, maxDim)),
                     Sig(StrongCollapseZigzag.Compute(f, maxDim)));

    [Fact]
    public void AddThenDelete()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Delete(0);
        AssertCollapseInvariant(f, 0);
    }

    [Fact]
    public void Taxonomy_co()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new[] { 0, 1 });
        AssertCollapseInvariant(f, 0);
    }

    [Fact]
    public void Taxonomy_oc()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new[] { 0, 1 });
        f.Delete(2); f.Delete(1);
        AssertCollapseInvariant(f, 0);
    }

    [Fact]
    public void Taxonomy_oo()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 });
        f.Delete(3); f.Add(5, new[] { 0, 1 });
        AssertCollapseInvariant(f, 0);
    }

    [Fact]
    public void TriangleFormsThenBreaks_H1()
    {
        // K_5 is a path (collapses to a point); K_6 is the circle (minimal) — barcode must survive.
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 }); f.Add(5, new[] { 0, 2 });
        f.Delete(5);
        AssertCollapseInvariant(f, 1);
    }

    [Fact]
    public void DynamicGraphH0()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 });
        f.Delete(3);
        AssertCollapseInvariant(f, 0);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(7)]
    [InlineData(42)]
    [InlineData(99)]
    [InlineData(123)]
    public void RandomGraphZigzag_CollapseInvariant(int seed) => AssertCollapseInvariant(RandomZigzag(seed), 1);

    static ZigzagFiltration RandomZigzag(int seed)
    {
        var rng = new Random(seed);
        int V = rng.Next(3, 7);
        var f = new ZigzagFiltration();
        for (int v = 0; v < V; v++) f.Add(v, new int[0]);
        int nextId = V;
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
        return f;
    }
}
