#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// Z5a — the slow-correct H0 graph-zigzag (<see cref="GraphZigzag"/>) must agree with both oracles
/// (Z1 <see cref="ZigzagBarcodeNaive"/> and Z2 <see cref="FastZigzag"/>) on dimension 0. Random sweeps
/// are the real guard for the merge/departure calibration; the vertex-churn sweep is what exercises
/// the departure walk (the edge-only generator never deletes vertices).
/// </summary>
public sealed class GraphZigzagTests
{
    static IEnumerable<(double, double, int, int, int)> Sig(Barcode bc) =>
        bc.Bars.Select(b => (b.Birth, b.Death, b.Dimension, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item3).ThenBy(x => x.Item1).ThenBy(x => x.Item2)
              .ThenBy(x => x.Item4).ThenBy(x => x.Item5);

    // This class tests the Reference (forest/Kruskal) engine — the Z5a/Z5b bug guards live here, so pin it.
    static void AssertH0Parity(ZigzagFiltration f)
    {
        var graph = Sig(GraphZigzag.Compute(f, 0, algorithm: GraphZigzagAlgorithm.Reference));
        Assert.Equal(Sig(ZigzagBarcodeNaive.Compute(f, 0)), graph);
        Assert.Equal(Sig(FastZigzag.Compute(f, 0)), graph);
    }

    // The combined entry point: GraphZigzag.Compute(f, 1) returns H0 + the folded-in H1 (Z5b).
    static void AssertFullParity(ZigzagFiltration f)
    {
        var graph = Sig(GraphZigzag.Compute(f, 1, algorithm: GraphZigzagAlgorithm.Reference));
        Assert.Equal(Sig(ZigzagBarcodeNaive.Compute(f, 1)), graph);
        Assert.Equal(Sig(FastZigzag.Compute(f, 1)), graph);
    }

    [Fact]
    public void CombinedH0H1_TriangleFormsThenBreaks()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 }); f.Add(5, new[] { 0, 2 });
        f.Delete(5);
        AssertFullParity(f);   // H0 component + H1 loop, from one Compute(f, 1)
    }

    [Fact]
    public void CombinedH0H1_DynamicGraph()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]); f.Add(3, new int[0]);
        f.Add(4, new[] { 0, 1 }); f.Add(5, new[] { 1, 2 }); f.Add(6, new[] { 0, 2 }); // triangle (H1)
        f.Add(7, new[] { 2, 3 });   // attach a tail (H0 merge, no new cycle)
        f.Delete(5);                // break the triangle (H1 death), still connected
        f.Delete(7);                // detach vertex 3 (H0 split)
        AssertFullParity(f);
    }

    [Fact]
    public void AddThenDelete()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Delete(0);
        AssertH0Parity(f);
    }

    [Fact]
    public void Taxonomy_co()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new[] { 0, 1 });
        AssertH0Parity(f);
    }

    [Fact]
    public void Taxonomy_oc()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new[] { 0, 1 });
        f.Delete(2); f.Delete(1);
        AssertH0Parity(f);
    }

    [Fact]
    public void Taxonomy_oo()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 });
        f.Delete(3); f.Add(5, new[] { 0, 1 });
        AssertH0Parity(f);
    }

    [Fact]
    public void DynamicGraphH0()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 });
        f.Delete(3);
        AssertH0Parity(f);
    }

    [Fact]
    public void TriangleFormsThenBreaks_H0()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 }); f.Add(5, new[] { 0, 2 });
        f.Delete(5);
        AssertH0Parity(f);
    }

    /// <summary>
    /// Bug-1 guard (merge gluing / multi-level path identification). A component splits twice
    /// (nested split nodes), then a sequence of merges must glue a path spanning several levels so
    /// the deeper split node is carried into the merged component; a later same-tree merge then
    /// reads its HCA at the deep split level. The pre-fix top-down glue only identified the top
    /// level, dropping the deep split node — so the same-tree merge fell through to a lower HCA and
    /// the orphaned split leaked an end-interval. (Minimal failing shape, localized from edge-churn
    /// seed 7 by the forest tracer: the engine first diverged from the oracle after step 10.)
    /// </summary>
    [Fact]
    public void SplitThenMerge_MultiLevelGlue_H0()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]); f.Add(3, new int[0]);
        f.Add(4, new[] { 1, 2 });   // merge {1},{2}
        f.Add(5, new[] { 0, 2 });   // merge in {0}
        f.Delete(4);                // split off {1}
        f.Delete(5);                // split {0,2} -> {0} | {2}  (nested splits)
        f.Add(6, new[] { 2, 3 });   // merge {2}-side with {3}: multi-level glue
        f.Add(7, new[] { 1, 3 });   // same-tree merge (resolves the level-6 split)
        f.Add(8, new[] { 0, 1 });   // same-tree merge — HCA must be the deep (level-7) split
        AssertH0Parity(f);
    }

    /// <summary>
    /// Bug-2 guard (departure with nested splits). The path a-b-c-d is built, then broken at b-c
    /// and c-d so that the isolated vertex c has two splitting ancestors — the c|d split (higher
    /// level) and the {a,b}|{c,d} split (lower level). Deleting c must pair with the highest-LEVEL
    /// (first-encountered) splitting ancestor, leaving the lower split to survive to the end; the
    /// pre-fix code paired with the lowest-level split, mis-dating the bar.
    /// </summary>
    [Fact]
    public void NestedSplitThenDeparture_H0()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]); f.Add(3, new int[0]);
        f.Add(4, new[] { 0, 1 }); f.Add(5, new[] { 1, 2 }); f.Add(6, new[] { 2, 3 }); // path a-b-c-d
        f.Delete(5);   // -bc: split {a,b} | {c,d}
        f.Delete(6);   // -cd: split {c} | {d}  (c gains a second, higher-level split ancestor)
        f.Delete(2);   // depart isolated c — two splitting ancestors on its path
        AssertH0Parity(f);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(7)]
    [InlineData(42)]
    [InlineData(99)]
    [InlineData(123)]
    public void RandomEdgeChurn_H0ParityWithOracles(int seed)
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
        AssertH0Parity(f);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(3)]
    [InlineData(8)]
    [InlineData(21)]
    [InlineData(55)]
    [InlineData(144)]
    public void RandomVertexChurn_H0ParityWithOracles(int seed)
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
        AssertH0Parity(f);
    }
}
