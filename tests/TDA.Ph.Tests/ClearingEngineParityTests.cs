#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Topology;
using Xunit;
using Xunit.Abstractions;

namespace TDA.Ph.Tests;

/// <summary>
/// The clearing engines against <see cref="PersistentHomology"/> on filtrations that are not Rips: vertices born
/// apart, cells in any order the faces allow, dimensions up to tetrahedra. Rips parity has its own tests; what
/// differs here is that the elder rule matters. Checks: <c>issues/ph/ph-loose-ends-brief.md</c>, step 2.
/// </summary>
public sealed class ClearingEngineParityTests
{
    readonly ITestOutputHelper _output;

    public ClearingEngineParityTests(ITestOutputHelper output) => _output = output;

    // ── by hand ─────────────────────────────────────────────────────────────────

    [Fact]
    public void ElderSurvives_WhenTwoComponentsMerge()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0), new Simplex(1.0, 1), new Simplex(2.0, 0, 1),
        });

        var expected = new[] { (0, 0.0, double.PositiveInfinity), (0, 1.0, 2.0) };

        Assert.Equal(expected, Canon(PersistentHomology.Compute(filtration)));
        Assert.Equal(expected, Canon(PersistentCohomology.Compute(filtration)));
        Assert.Equal(expected, Canon(PersistentInvolutedHomology.Compute(filtration)));
    }

    [Fact]
    public void ElderSurvivesAcrossAChainOfMerges()
    {
        // The middle vertex is youngest: it dies first, then the component it joined keeps the elder's birth.
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0), new Simplex(1.0, 1), new Simplex(2.0, 2),
            new Simplex(3.0, 1, 2), new Simplex(4.0, 0, 1),
        });

        var expected = new[] { (0, 0.0, double.PositiveInfinity), (0, 1.0, 4.0), (0, 2.0, 3.0) };

        Assert.Equal(expected, Canon(PersistentHomology.Compute(filtration)));
        Assert.Equal(expected, Canon(PersistentCohomology.Compute(filtration)));
        Assert.Equal(expected, Canon(PersistentInvolutedHomology.Compute(filtration)));
    }

    // ── distinct vertex births ──────────────────────────────────────────────────

    [Fact]
    public void DistinctVertexBirths_BothEnginesEqualTheColumnReducer_InEveryDegree()
    {
        const int trials = 400;
        int invertedWouldDiffer = 0;
        var finite = new int[4];
        var essential = new int[4];

        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            // Half the trials draw cell births from a small grid, so edges and cofaces tie with each other and with
            // their faces, which gives zero-length pairs and tied merges.
            // A quarter of the trials run below zero.
            SimplicialFiltration filtration = RandomFiltration(rng, distinctVertexBirths: true, grid: seed % 2 == 0, shift: seed % 4 == 3 ? -20.0 : 0.0);

            foreach (int maxDimension in new[] { int.MaxValue, 0, 1, 2 })
            {
                var reference = Canon(PersistentHomology.Compute(filtration, maxDimension));
                Assert.True(reference.SequenceEqual(Canon(PersistentCohomology.Compute(filtration, maxDimension))),
                    $"cohomology, seed {seed}, maxDimension {maxDimension}");
                Assert.True(reference.SequenceEqual(Canon(PersistentInvolutedHomology.Compute(filtration, maxDimension))),
                    $"involuted, seed {seed}, maxDimension {maxDimension}");
            }

            var full = Canon(PersistentHomology.Compute(filtration));
            foreach (var bar in full)
            {
                if (double.IsPositiveInfinity(bar.Item3)) essential[bar.Item1]++;
                else finite[bar.Item1]++;
            }

            var h0 = full.Where(b => b.Item1 == 0).ToArray();
            if (!h0.SequenceEqual(H0(filtration, keepYounger: true)))
                invertedWouldDiffer++;
        }

        // The sweep has to reach what the fix is about, and the degrees it claims.
        _output.WriteLine($"trials where the inverted rule differs: {invertedWouldDiffer}/{trials}; " +
            $"finite H0-H3: {string.Join("/", finite)}; essential H0-H3: {string.Join("/", essential)}");
        Assert.True(invertedWouldDiffer > trials / 2);
        Assert.True(finite[0] > 0 && finite[1] > 0 && finite[2] > 0);
        Assert.True(essential[0] > 0 && essential[1] > 0 && essential[2] > 0);
    }

    // ── all vertices at 0: nothing may move ─────────────────────────────────────

    [Fact]
    public void VerticesAllAtZero_BarsDoNotDependOnTheRule()
    {
        const int trials = 200;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            SimplicialFiltration filtration = RandomFiltration(rng, distinctVertexBirths: false, grid: seed % 2 == 0);

            var reference = Canon(PersistentHomology.Compute(filtration));
            Assert.True(reference.SequenceEqual(Canon(PersistentCohomology.Compute(filtration))), $"cohomology, seed {seed}");
            Assert.True(reference.SequenceEqual(Canon(PersistentInvolutedHomology.Compute(filtration))), $"involuted, seed {seed}");

            // The rule the engines used to run gives these same bars, so no Rips result moved.
            Assert.True(H0(filtration, keepYounger: true).SequenceEqual(H0(filtration, keepYounger: false)), $"rule, seed {seed}");
            Assert.True(reference.Where(b => b.Item1 == 0).SequenceEqual(H0(filtration, keepYounger: false)), $"H0, seed {seed}");
        }
    }

    // ── representatives on non-Rips filtrations ─────────────────────────────────

    [Fact]
    public void DistinctVertexBirths_InvolutedCyclesCloseAndLieInTheBirthComplex()
    {
        const int trials = 200;
        int finite = 0, essential = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            SimplicialFiltration filtration = RandomFiltration(rng, distinctVertexBirths: true, grid: seed % 2 == 0, shift: seed % 4 == 3 ? -20.0 : 0.0);
            Barcode barcode = PersistentInvolutedHomology.Compute(filtration, 1, representatives: true);

            foreach (Bar bar in barcode.Bars.Where(b => b.Dimension == 1 && (b.IsInfinite || b.Death > b.Birth)))
            {
                int[] cycle = bar.Cycle!;
                var boundary = new HashSet<int>();
                foreach (int c in cycle)
                {
                    Assert.Equal(1, filtration.GetDimension(c));
                    Assert.True(filtration.GetBirth(c) <= bar.Birth, $"seed {seed}: a cycle cell born after its bar");
                    foreach (int f in filtration.GetBoundaryIndices(c))
                        if (!boundary.Remove(f)) boundary.Add(f);
                }

                Assert.True(boundary.Count == 0, $"seed {seed}: the cycle does not close");
                if (bar.IsInfinite) essential++;
                else finite++;
            }
        }

        _output.WriteLine($"cycles checked: {finite} finite, {essential} essential");
        Assert.True(finite > 0 && essential > 0);
    }

    // ── helpers ─────────────────────────────────────────────────────────────────

    // Any complex the faces allow, not a flag complex: each candidate coface enters with a chance, at or after its
    // latest face. Goes up to tetrahedra so H2 can die.
    static SimplicialFiltration RandomFiltration(Random rng, bool distinctVertexBirths, bool grid, double shift = 0.0)
    {
        int n = rng.Next(4, 10);
        double Later(double earliest) => grid ? earliest + rng.Next(0, 4) : earliest + rng.NextDouble() * 3.0;

        var births = new Dictionary<string, double>();
        var cells = new List<Simplex>();

        void Add(double birth, params int[] vertices)
        {
            births[Key(vertices)] = birth;
            cells.Add(new Simplex(birth + shift, vertices));
        }

        // A permutation makes the vertex births distinct whichever values the rest of the complex draws.
        int[] order = Enumerable.Range(0, n).OrderBy(_ => rng.Next()).ToArray();
        for (int v = 0; v < n; v++)
            Add(distinctVertexBirths ? order[v] : 0.0, v);

        var edges = new List<(int, int)>();
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
                if (rng.NextDouble() < 0.5)
                {
                    Add(Later(Math.Max(births[Key(a)], births[Key(b)])), a, b);
                    edges.Add((a, b));
                }

        var triangles = new List<(int, int, int)>();
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
                for (int c = b + 1; c < n; c++)
                {
                    if (!births.TryGetValue(Key(a, b), out double ab)
                        || !births.TryGetValue(Key(a, c), out double ac)
                        || !births.TryGetValue(Key(b, c), out double bc)
                        || rng.NextDouble() >= 0.6)
                        continue;

                    Add(Later(Math.Max(ab, Math.Max(ac, bc))), a, b, c);
                    triangles.Add((a, b, c));
                }

        foreach (var (a, b, c) in triangles)
            for (int d = c + 1; d < n; d++)
            {
                if (!births.TryGetValue(Key(a, b, d), out double abd)
                    || !births.TryGetValue(Key(a, c, d), out double acd)
                    || !births.TryGetValue(Key(b, c, d), out double bcd)
                    || rng.NextDouble() >= 0.6)
                    continue;

                Add(Later(Math.Max(births[Key(a, b, c)], Math.Max(abd, Math.Max(acd, bcd)))), a, b, c, d);
            }

        return new SimplicialFiltration(cells);
    }

    static string Key(params int[] vertices) => string.Join(",", vertices);

    // H0 by union-find over the filtration's edges, written apart from the engine, with either rule.
    static (int, double, double)[] H0(SimplicialFiltration filtration, bool keepYounger)
    {
        var parent = new Dictionary<int, int>();
        var birth = new Dictionary<int, double>();
        var bars = new List<(int, double, double)>();

        int Find(int x) => parent[x] == x ? x : parent[x] = Find(parent[x]);

        foreach (Simplex s in filtration.Simplices.Where(s => s.Dimension == 0))
        {
            parent[s.Vertices[0]] = s.Vertices[0];
            birth[s.Vertices[0]] = s.FiltrationValue;
        }

        foreach (Simplex s in filtration.Simplices.Where(s => s.Dimension == 1))
        {
            int a = Find(s.Vertices[0]), b = Find(s.Vertices[1]);
            if (a == b) continue;

            bool aDies = keepYounger ? birth[a] < birth[b] : birth[a] > birth[b];
            int dying = aDies ? a : b;
            int surviving = aDies ? b : a;
            if (s.FiltrationValue > birth[dying])
                bars.Add((0, birth[dying], s.FiltrationValue));
            parent[dying] = surviving;
        }

        foreach (int v in parent.Keys.Where(v => parent[v] == v))
            bars.Add((0, birth[v], double.PositiveInfinity));

        return bars.OrderBy(x => x).ToArray();
    }

    // Significant bars as a sorted multiset, compared exactly.
    static (int, double, double)[] Canon(Barcode barcode) =>
        barcode.Bars
            .Where(b => b.IsInfinite || b.Death > b.Birth)
            .Select(b => (b.Dimension, b.Birth, b.Death))
            .OrderBy(x => x)
            .ToArray();
}
