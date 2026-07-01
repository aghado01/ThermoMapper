#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// Z5 representatives (H1) — <see cref="GraphZigzagH1.Compute(ZigzagFiltration, bool)"/> with
/// <c>representatives: true</c> must (a) leave the barcode unchanged vs the oracles, and (b) emit, for
/// every H1 bar, a representative 1-cycle that is genuinely sound: a Z/2 cycle (∂ = 0), present in
/// <i>every</i> intermediate graph across the bar's life, and containing both the birth and death edges.
/// These are the checkable necessary conditions of a valid zigzag representative (Dey–Hou Prop 17).
/// </summary>
public sealed class GraphZigzagH1RepresentativeTests
{
    static IEnumerable<(double, double, int, int, int)> Sig1(Barcode bc) =>
        bc.Bars.Where(b => b.Dimension == 1)
              .Select(b => (b.Birth, b.Death, b.Dimension, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item1).ThenBy(x => x.Item2).ThenBy(x => x.Item4).ThenBy(x => x.Item5);

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

    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(7)] [InlineData(42)] [InlineData(99)]
    [InlineData(123)] [InlineData(2024)] [InlineData(31337)] [InlineData(8675309)] [InlineData(13)]
    public void H1Representatives_AreSound(int seed) => AssertSound(BuildEdgeChurn(seed), GraphZigzagH1.Compute(BuildEdgeChurn(seed), representatives: true));

    // The near-linear engine carries the same sound H1 reps (MSF path from the dynamic MSF).
    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(7)] [InlineData(42)] [InlineData(99)]
    [InlineData(123)] [InlineData(2024)] [InlineData(31337)] [InlineData(8675309)] [InlineData(13)]
    public void Fast_H1Representatives_AreSound(int seed) => AssertSound(BuildEdgeChurn(seed), GraphZigzagH1Fast.Compute(BuildEdgeChurn(seed), representatives: true));

    static void AssertSound(ZigzagFiltration f, Barcode withReps)
    {
        int m = f.Count;

        // Cell -> endpoints, and the per-complex present-cell sets G_0..G_m.
        var ends = new Dictionary<int, (int U, int V)>();
        var present = new HashSet<int>();
        var K = new List<HashSet<int>> { new HashSet<int>() };
        foreach (var step in f)
        {
            if (step.Direction == ZigzagDirection.Add)
            {
                present.Add(step.GlobalCellId);
                if (step.BoundaryAtAdd!.Length == 2) ends[step.GlobalCellId] = (step.BoundaryAtAdd[0], step.BoundaryAtAdd[1]);
            }
            else present.Remove(step.GlobalCellId);
            K.Add(new HashSet<int>(present));
        }

        // (a) representatives must not perturb the barcode.
        Assert.Equal(Sig1(ZigzagBarcodeNaive.Compute(f, 1)), Sig1(withReps));
        Assert.Equal(Sig1(FastZigzag.Compute(f, 1)), Sig1(withReps));

        // (b) every H1 bar carries a sound representative cycle.
        foreach (var bar in withReps.Bars.Where(b => b.Dimension == 1))
        {
            Assert.NotNull(bar.Cycle);
            var z = bar.Cycle!;
            Assert.NotEmpty(z);

            // Z/2 cycle: every vertex has even degree among the cycle's edges.
            var deg = new Dictionary<int, int>();
            foreach (int cell in z)
            {
                var (u, v) = ends[cell];
                deg[u] = deg.GetValueOrDefault(u) + 1;
                deg[v] = deg.GetValueOrDefault(v) + 1;
            }
            Assert.All(deg.Values, d => Assert.True(d % 2 == 0, "cycle boundary nonzero"));

            int b = (int)bar.Birth + 1, d = (int)bar.Death;     // present in complexes G_b..G_d

            // Present in every intermediate graph across the bar.
            for (int k = b; k <= d; k++)
                foreach (int cell in z) Assert.Contains(cell, K[k]);

            // Contains the birth edge (added at step Birth) and, when finite, the death edge.
            Assert.Contains(f[(int)bar.Birth].GlobalCellId, z);
            if (d < m) Assert.Contains(f[d].GlobalCellId, z);
        }
    }
}
