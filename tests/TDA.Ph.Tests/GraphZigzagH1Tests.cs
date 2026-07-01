#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// Z5b — the slow-correct H1 graph-zigzag (<see cref="GraphZigzagH1"/>, Dey–Hou Algorithm 2) must agree
/// with both oracles (Z1 <see cref="ZigzagBarcodeNaive"/> and Z2 <see cref="FastZigzag"/>) on dimension 1.
/// The engine emits only dim-1 bars, so parity is taken over each oracle's dim-1 slice. Random edge
/// churn (few vertices, dense edges) is the real guard — it births and kills many overlapping cycles,
/// exercising the "earliest compatible birth" pairing (Remark 16: zigzag pairs the smallest, not the
/// youngest, positive index).
/// </summary>
public sealed class GraphZigzagH1Tests
{
    static IEnumerable<(double, double, int, int, int)> Sig1(Barcode bc) =>
        bc.Bars.Where(b => b.Dimension == 1)
              .Select(b => (b.Birth, b.Death, b.Dimension, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item1).ThenBy(x => x.Item2).ThenBy(x => x.Item4).ThenBy(x => x.Item5);

    static void AssertH1Parity(ZigzagFiltration f)
    {
        var graph = Sig1(GraphZigzagH1.Compute(f)).ToList();
        Assert.Equal(Sig1(ZigzagBarcodeNaive.Compute(f, 1)), graph);
        Assert.Equal(Sig1(FastZigzag.Compute(f, 1)), graph);
    }

    [Fact]
    public void NoEdges_NoH1()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Delete(0);
        Assert.Empty(GraphZigzagH1.Compute(f).Bars);
        AssertH1Parity(f);
    }

    [Fact]
    public void TriangleFormsThenBreaks_H1()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 }); f.Add(5, new[] { 0, 2 });
        f.Delete(5);
        AssertH1Parity(f);
        // Anchor (anti-circularity): the single loop is born when the closing edge (step 5) is added
        // and dies when it is removed (step 6) — a closed–closed bar [5, 6].
        var bar = Assert.Single(GraphZigzagH1.Compute(f).Bars);
        Assert.Equal((5.0, 6.0, 1, IntervalEnd.Closed, IntervalEnd.Closed),
            (bar.Birth, bar.Death, bar.Dimension, bar.BirthEnd, bar.DeathEnd));
    }

    [Fact]
    public void TrianglePersistsToEnd_H1_OpenDeath()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 }); f.Add(5, new[] { 0, 2 });
        AssertH1Parity(f);   // loop survives -> open death at m
    }

    [Fact]
    public void SquareFormsThenBreaks_H1()
    {
        var f = new ZigzagFiltration();
        for (int v = 0; v < 4; v++) f.Add(v, new int[0]);
        f.Add(4, new[] { 0, 1 }); f.Add(5, new[] { 1, 2 }); f.Add(6, new[] { 2, 3 });
        f.Add(7, new[] { 0, 3 });   // closes the 4-cycle
        f.Delete(5);                // break it
        AssertH1Parity(f);
    }

    [Fact]
    public void ThetaGraph_TwoCyclesShareEdge_H1()
    {
        // Vertices 0,1; three parallel paths 0-1 via a middle vertex each -> two independent cycles.
        var f = new ZigzagFiltration();
        for (int v = 0; v < 5; v++) f.Add(v, new int[0]);   // 0,1 endpoints; 2,3,4 midpoints
        f.Add(5, new[] { 0, 2 }); f.Add(6, new[] { 2, 1 }); // path A: 0-2-1
        f.Add(7, new[] { 0, 3 }); f.Add(8, new[] { 3, 1 }); // path B closes cycle 1
        f.Add(9, new[] { 0, 4 }); f.Add(10, new[] { 4, 1 }); // path C closes cycle 2
        f.Delete(6);   // remove a path-A edge: still two paths -> one cycle dies, one survives
        AssertH1Parity(f);
    }

    [Fact]
    public void LoopReform_EarliestBirthPairing_H1()
    {
        // A cycle is born, a second overlapping cycle is born later, then a shared edge is deleted:
        // the death must pair with the EARLIER birth (the surviving cycle keeps the later one).
        var f = new ZigzagFiltration();
        for (int v = 0; v < 4; v++) f.Add(v, new int[0]);
        f.Add(4, new[] { 0, 1 }); f.Add(5, new[] { 1, 2 }); f.Add(6, new[] { 0, 2 }); // triangle 0-1-2 (birth A)
        f.Add(7, new[] { 2, 3 }); f.Add(8, new[] { 0, 3 });                            // triangle 0-2-3 (birth B)
        f.Delete(6);   // delete shared chord 0-2: endpoints still connected -> a death
        AssertH1Parity(f);
    }

    // ---- random sweeps -------------------------------------------------------------------------
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
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(7)]
    [InlineData(42)]
    [InlineData(99)]
    [InlineData(123)]
    [InlineData(2024)]
    [InlineData(31337)]
    public void RandomEdgeChurn_H1ParityWithOracles(int seed) => AssertH1Parity(BuildEdgeChurn(seed));

    [Theory]
    [InlineData(1)]
    [InlineData(3)]
    [InlineData(8)]
    [InlineData(21)]
    [InlineData(55)]
    [InlineData(144)]
    public void RandomVertexChurn_H1ParityWithOracles(int seed) => AssertH1Parity(BuildVertexChurn(seed));
}
