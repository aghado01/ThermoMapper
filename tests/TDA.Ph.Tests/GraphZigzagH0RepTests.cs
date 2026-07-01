#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace TDA.Ph.Tests;

/// <summary>
/// Z5 representatives (H0) — <see cref="GraphZigzag.Compute(ZigzagFiltration, int, bool)"/> with
/// <c>representatives: true</c> must leave the H0 barcode unchanged, and every <i>emitted</i> H0 rep
/// 0-chain must be sound: present in every intermediate graph across its bar, and nonzero there. For a
/// fixed vertex set (edge churn) every bar gets a rep and, at each complex, the alive bars' reps form a
/// basis of H0(G_k). With vertex departures a rep that would reference a deleted vertex is left null
/// (deferred to the full evolving reduction) — so we check soundness of emitted reps, not completeness.
/// An entrance births {v}; a split births the kernel {minA, minB}; each bar's rep is its birth chain.
/// </summary>
public sealed class GraphZigzagH0RepTests
{
    static IEnumerable<(double, double, int, int)> Sig0(Barcode bc) =>
        bc.Bars.Where(b => b.Dimension == 0)
              .Select(b => (b.Birth, b.Death, (int)b.BirthEnd, (int)b.DeathEnd))
              .OrderBy(x => x.Item1).ThenBy(x => x.Item2).ThenBy(x => x.Item3).ThenBy(x => x.Item4);

    // Rank over Z/2 of a set of vectors, each given as the set of component-roots with odd parity.
    static int Rank(List<HashSet<int>> vectors)
    {
        var pivots = new Dictionary<int, HashSet<int>>();   // low component -> reduced vector
        int rank = 0;
        foreach (var v0 in vectors)
        {
            var v = new HashSet<int>(v0);
            while (v.Count > 0)
            {
                int low = v.Max();
                if (!pivots.TryGetValue(low, out var p)) { pivots[low] = v; rank++; break; }
                v.SymmetricExceptWith(p);
            }
        }
        return rank;
    }

    static void Check(ZigzagFiltration f) => Check(f, GraphZigzag.Compute(f, 0, representatives: true, algorithm: GraphZigzagAlgorithm.Reference));
    static void CheckFast(ZigzagFiltration f) => Check(f, GraphZigzagFast.Compute(f, 0, representatives: true));

    static void Check(ZigzagFiltration f, Barcode bc)
    {
        int m = f.Count;
        var isVertex = new Dictionary<int, bool>();
        var ends = new Dictionary<int, (int U, int V)>();
        var present = new HashSet<int>();
        var K = new List<HashSet<int>> { new() };
        foreach (var st in f)
        {
            if (st.Direction == ZigzagDirection.Add) { present.Add(st.GlobalCellId); bool v = st.BoundaryAtAdd!.Length == 0; isVertex[st.GlobalCellId] = v; if (!v) ends[st.GlobalCellId] = (st.BoundaryAtAdd[0], st.BoundaryAtAdd[1]); }
            else present.Remove(st.GlobalCellId);
            K.Add(new HashSet<int>(present));
        }

        Assert.Equal(Sig0(ZigzagBarcodeNaive.Compute(f, 0)), Sig0(bc));
        Assert.Equal(Sig0(FastZigzag.Compute(f, 0)), Sig0(bc));

        var h0 = bc.Bars.Where(b => b.Dimension == 0).ToList();

        for (int k = 0; k <= m; k++)
        {
            var verts = K[k].Where(c => isVertex[c]).ToHashSet();
            var parent = new Dictionary<int, int>();
            int Find(int x) { if (!parent.TryGetValue(x, out int p)) { parent[x] = x; return x; } while (p != x) { parent[x] = parent[p]; x = p; p = parent[x]; } return x; }
            foreach (int v in verts) Find(v);
            foreach (int cell in K[k]) if (!isVertex[cell]) { var (a, b) = ends[cell]; parent[Find(a)] = Find(b); }
            int compCount = verts.Select(Find).Distinct().Count();

            int alive = 0;
            var vectors = new List<HashSet<int>>();
            foreach (var bar in h0)
            {
                if (k < (int)bar.Birth + 1 || k > (int)bar.Death) continue;
                alive++;
                if (bar.Cycle == null) continue;             // departure-affected rep, left null (sound omission)
                var vec = new HashSet<int>();
                foreach (int x in bar.Cycle)
                {
                    Assert.Contains(x, verts);               // emitted rep present at G_k
                    int c = Find(x);
                    if (!vec.Add(c)) vec.Remove(c);          // Z/2 parity over components
                }
                Assert.NotEmpty(vec);                        // nonzero class at G_k
                vectors.Add(vec);
            }
            Assert.Equal(vectors.Count, Rank(vectors));      // emitted reps are independent
            if (alive == vectors.Count) Assert.Equal(compCount, vectors.Count);  // complete -> full H0 basis
        }
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
            for (int u = 0; u < V; u++) for (int w = u + 1; w < V; w++) if (!edges.ContainsKey((u, w))) addable.Add((u, w));
            bool doAdd = edges.Count == 0 || (addable.Count > 0 && rng.NextDouble() < 0.6);
            if (doAdd) { var e = addable[rng.Next(addable.Count)]; int eid = nextId++; edges[e] = eid; f.Add(eid, new[] { e.Item1, e.Item2 }); }
            else { var e = edges.Keys.ElementAt(rng.Next(edges.Count)); f.Delete(edges[e]); edges.Remove(e); }
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
            if (r < 0.45 && addable.Count > 0) { var e = addable[rng.Next(addable.Count)]; int eid = nextId++; edges[e] = eid; incident[e.Item1]++; incident[e.Item2]++; f.Add(eid, new[] { e.Item1, e.Item2 }); }
            else if (r < 0.70 && edges.Count > 0) { var e = edges.Keys.ElementAt(rng.Next(edges.Count)); f.Delete(edges[e]); edges.Remove(e); incident[e.Item1]--; incident[e.Item2]--; }
            else if (r < 0.85 && isolated.Count > 0) { int v = isolated[rng.Next(isolated.Count)]; f.Delete(v); present.Remove(v); incident.Remove(v); }
            else { int v = nextId++; f.Add(v, new int[0]); present.Add(v); incident[v] = 0; }
        }
        return f;
    }

    [Fact]
    public void StrictMode_FailsLoudInsteadOfSilentPartial()
    {
        // Lenient: never throws; barcode unchanged. Strict: completes (all H0 reps present) OR throws a
        // clear NotSupportedException — never a silent partial. And strict must actually fire on the gap.
        int strictThrows = 0;
        for (int seed = 1; seed <= 40; seed++)
        {
            var f = BuildVertexChurn(seed);
            var lenient = GraphZigzag.Compute(f, 0, representatives: true, algorithm: GraphZigzagAlgorithm.Reference);   // must not throw
            Assert.Equal(
                ZigzagBarcodeNaive.Compute(f, 0).Bars.Count(b => b.Dimension == 0),
                lenient.Bars.Count(b => b.Dimension == 0));
            try
            {
                var strict = GraphZigzag.Compute(f, 0, representatives: true, strictRepresentatives: true, algorithm: GraphZigzagAlgorithm.Reference);
                Assert.All(strict.Bars.Where(b => b.Dimension == 0), b => Assert.NotNull(b.Cycle));
            }
            catch (NotSupportedException) { strictThrows++; }
        }
        Assert.True(strictThrows > 0, "strict mode never exercised — expected some vertex-departure nulls");

        // Fixed vertex set (edge churn): strict must always succeed (every H0 bar is representable).
        for (int seed = 1; seed <= 10; seed++)
            GraphZigzag.Compute(BuildEdgeChurn(seed), 0, representatives: true, strictRepresentatives: true, algorithm: GraphZigzagAlgorithm.Reference);
    }

    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(7)] [InlineData(42)] [InlineData(99)]
    [InlineData(123)] [InlineData(2024)] [InlineData(31337)] [InlineData(8675309)] [InlineData(13)]
    public void EdgeChurn_H0Reps_FormBasisEachLevel(int seed) => Check(BuildEdgeChurn(seed));

    [Theory]
    [InlineData(1)] [InlineData(3)] [InlineData(8)] [InlineData(21)] [InlineData(55)]
    [InlineData(144)] [InlineData(7)] [InlineData(99)] [InlineData(2024)] [InlineData(512)]
    public void VertexChurn_H0Reps_FormBasisEachLevel(int seed) => Check(BuildVertexChurn(seed));

    // Near-linear engine carries the same sound H0 reps.
    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(7)] [InlineData(42)] [InlineData(99)]
    [InlineData(123)] [InlineData(2024)] [InlineData(31337)] [InlineData(8675309)] [InlineData(13)]
    public void Fast_EdgeChurn_H0Reps_FormBasisEachLevel(int seed) => CheckFast(BuildEdgeChurn(seed));

    [Theory]
    [InlineData(1)] [InlineData(3)] [InlineData(8)] [InlineData(21)] [InlineData(55)]
    [InlineData(144)] [InlineData(7)] [InlineData(99)] [InlineData(2024)] [InlineData(512)]
    public void Fast_VertexChurn_H0Reps_FormBasisEachLevel(int seed) => CheckFast(BuildVertexChurn(seed));
}
