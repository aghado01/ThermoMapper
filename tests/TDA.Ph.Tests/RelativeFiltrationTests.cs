#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Graphs.Primitives;
using Maths.Topology;
using Xunit;
using Xunit.Abstractions;

namespace TDA.Ph.Tests;

/// <summary>
/// The relative reading <c>H(K_δ, K₀)</c>, against hand-worked complexes and two references that do not go
/// through the relative complex: the quotient identity read by the involuted engine on the union, and reduced
/// homology of <c>K ∪ cone(K₀)</c>. Checks and expectations: <c>issues/ph/ph-setup-brief.md</c> (R1–R11).
/// </summary>
public sealed class RelativeFiltrationTests
{
    readonly ITestOutputHelper _output;

    public RelativeFiltrationTests(ITestOutputHelper output) => _output = output;

    // ── R1 · K₀ = ∅ ─────────────────────────────────────────────────────────────

    [Fact]
    public void EmptyBackbone_RelativeEqualsAbsolute_InEveryDegree()
    {
        const int n = 6;
        var none = new List<(int, int)>();
        var content = new List<(int, int, double)>
        {
            (0, 1, 1.0), (1, 2, 2.0), (2, 3, 3.0), (3, 4, 4.0), (4, 5, 5.0), (5, 0, 6.0), (0, 2, 7.0),
        };

        CsrGraph graph = ConditionedFiltration.BuildGraph(n, none, content);
        var filtration = RipsFiltration.GraphRips(graph, FiltrationWeights.RawDistance, 2);
        Barcode absolute = PersistentHomology.Compute(filtration);

        Barcode relative = ConditionedFiltration.BuildRelative(n, none, content).ComputeBarcode();

        Assert.Equal(Canon(absolute), Canon(relative));
    }

    // ── R2 · an acyclic backbone ────────────────────────────────────────────────

    [Fact]
    public void PathBackbone_HasNoH0_AndAbsoluteH1()
    {
        const int n = 8;
        List<(int, int)> backbone = PathEdges(n);
        var content = new List<(int, int, double)> { (0, 4, 1.5), (2, 7, 2.5), (1, 6, 3.5) };

        Barcode absolute = ConditionedFiltration.ComputeBarcode(n, backbone, content);
        Barcode relative = ConditionedFiltration.ComputeRelativeBarcode(n, backbone, content);

        Assert.Empty(Significant(relative, 0));
        Assert.Equal(Canon(absolute, 1), Canon(relative, 1));
    }

    // ── R3 · a backbone with its own loop ───────────────────────────────────────

    [Fact]
    public void CyclicBackbone_LoopIsQuotiented_ChordStays()
    {
        const int n = 6;
        const double d = 2.5;
        List<(int, int)> backbone = CycleEdges(n);
        var content = new List<(int, int, double)> { (0, 3, d) };

        Barcode absolute = ConditionedFiltration.ComputeBarcode(n, backbone, content);
        Assert.Equal(new[] { (1, 0.0, double.PositiveInfinity), (1, d, double.PositiveInfinity) }, Canon(absolute, 1));

        RelativeFiltration relative = ConditionedFiltration.BuildRelative(n, backbone, content);
        Bar chord = Assert.Single(Significant(relative.ComputeBarcode(1), 1));

        Assert.True(chord.IsInfinite);
        Assert.Equal(d, chord.Birth);
        Assert.Equal(new[] { 0, 3 }, relative.Ambient.Simplices[chord.Generator!.Value].Vertices.OrderBy(v => v));
    }

    // ── R4 · a backbone loop that K fills ───────────────────────────────────────

    [Fact]
    public void FilledBackboneLoop_ReappearsInH2_BornWhereItIsFilled()
    {
        const int n = 4;
        const double d = 2.5;
        List<(int, int)> backbone = CycleEdges(n);
        var content = new List<(int, int, double)> { (0, 2, d) };

        Barcode absolute = ConditionedFiltration.ComputeBarcode(n, backbone, content);
        Assert.Equal(new[] { (1, 0.0, d) }, Canon(absolute, 1));

        Barcode relative = ConditionedFiltration.BuildRelative(n, backbone, content).ComputeBarcode();
        Assert.Empty(Significant(relative, 1));
        Assert.Equal(new[] { (2, d, double.PositiveInfinity) }, Canon(relative, 2));
    }

    // ── R5 · a backbone that does not span ──────────────────────────────────────

    [Fact]
    public void NonSpanningBackbone_H0AttachesToIt()
    {
        const int n = 6;
        const double a = 1.5, b = 2.5;
        var backbone = new List<(int, int)> { (0, 1), (1, 2) };
        var content = new List<(int, int, double)> { (2, 3, a), (4, 5, b) };

        Barcode relative = ConditionedFiltration.ComputeRelativeBarcode(n, backbone, content);

        Assert.Equal(
            new[] { (0, 0.0, a), (0, 0.0, b), (0, 0.0, double.PositiveInfinity) }.OrderBy(x => x).ToArray(),
            Canon(relative, 0));
    }

    // ── R6 · a predicted return stays ───────────────────────────────────────────

    [Fact]
    public void ContentAtResidualZero_StaysOutsideK0()
    {
        const int n = 6;
        List<(int, int)> backbone = PathEdges(n);
        var observations = new double[] { 0, 1, 2, 3, 4, 5 };
        var prior = new List<(int, int, double)> { (0, 5, 5.0) }; // predicts the gap exactly
        IReadOnlyList<(int i, int j, double r)> content = ResidualPrior.ResidualEdges(observations, prior);

        Barcode relative = ConditionedFiltration.ComputeRelativeBarcode(n, backbone, content);

        Assert.Equal(new[] { (1, 0.0, double.PositiveInfinity) }, Canon(relative, 1));
    }

    // ── R7 · the quotient identity, against the involuted engine on the union ───

    [Fact]
    public void ConnectedSpanningBackbone_RelativeH1_IsAbsoluteH1BornAfterK0()
    {
        const int trials = 200;
        int quotiented = 0, kept = 0, finite = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(5, 10);
            List<(int, int)> backbone = RandomConnectedSpanning(n, rng);
            List<(int, int, double)> content = RandomContent(n, rng, backbone, 0.35);

            Barcode absolute = ConditionedFiltration.ComputeBarcode(n, backbone, content);
            Barcode relative = ConditionedFiltration.ComputeRelativeBarcode(n, backbone, content);

            var absoluteH1 = Canon(absolute, 1);
            var expected = absoluteH1.Where(bar => bar.Item2 > 0.0).ToArray();
            var relativeH1 = Canon(relative, 1);
            Assert.True(expected.SequenceEqual(relativeH1), $"seed {seed}");

            quotiented += absoluteH1.Length - expected.Length;
            kept += relativeH1.Length;
            finite += relativeH1.Count(bar => !double.IsPositiveInfinity(bar.Item3));
        }

        // The sweep has to reach the cases it claims: backbone loops removed, returns kept, returns that die.
        _output.WriteLine($"backbone bars quotiented: {quotiented}, relative H1 bars: {kept}, finite: {finite}");
        Assert.True(quotiented > 0 && kept > 0 && finite > 0);
    }

    // ── R8 · reduced homology of K ∪ cone(K₀) ───────────────────────────────────

    [Fact]
    public void ArbitraryBackbone_MatchesReducedHomologyOfTheCone()
    {
        const int trials = 200;
        var byDegree = new int[3];
        int nonSpanning = 0, filledK0 = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(4, 9);
            List<(int, int)> backbone = RandomEdges(n, rng, 0.3);
            if (backbone.Count == 0) backbone.Add((0, 1));
            List<(int, int, double)> content = RandomContent(n, rng, backbone, 0.35);

            RelativeFiltration relative = ConditionedFiltration.BuildRelative(n, backbone, content);
            var relativeBars = Canon(relative.ComputeBarcode());

            HashSet<int> k0Vertices = backbone.SelectMany(e => new[] { e.Item1, e.Item2 }).ToHashSet();
            HashSet<(int, int)> k0Edges = backbone.Select(Ordered).ToHashSet();
            SimplicialFiltration coned = Cone(relative.Ambient, s => InK0(s, k0Vertices, k0Edges), apex: n);
            var conedBars = Canon(PersistentHomology.Compute(coned)).ToList();

            // Reduced homology drops the one H₀ class of the cone's own component.
            Assert.True(conedBars.Remove((0, 0.0, double.PositiveInfinity)), $"seed {seed}: no cone class");
            Assert.True(conedBars.SequenceEqual(relativeBars), $"seed {seed}");

            foreach (var bar in relativeBars) byDegree[bar.Item1]++;
            if (k0Vertices.Count < n) nonSpanning++;
            if (relative.Ambient.Simplices.Any(s => s.Dimension == 2 && InK0(s, k0Vertices, k0Edges))) filledK0++;
        }

        // The sweep has to reach every degree, backbones that do not span, and K₀ triangles.
        _output.WriteLine($"bars H0/H1/H2: {byDegree[0]}/{byDegree[1]}/{byDegree[2]}; non-spanning: {nonSpanning}; K0 triangles: {filledK0}");
        Assert.True(byDegree.All(c => c > 0) && nonSpanning > 0 && filledK0 > 0);
    }

    // ── R10 · a finite bar's representative, in the ambient's indices ───────────

    [Fact]
    public void FiniteBar_RepresentativeIsARelativeCycle_InAmbientIndices()
    {
        const int n = 5;
        List<(int, int)> backbone = PathEdges(n);
        var content = new List<(int, int, double)> { (0, 4, 1.0), (0, 2, 2.0), (2, 4, 3.0) };

        RelativeFiltration relative = ConditionedFiltration.BuildRelative(n, backbone, content);
        Bar bar = Assert.Single(Significant(relative.ComputeBarcode(1, representatives: true), 1));

        Assert.Equal((1.0, 3.0), (bar.Birth, bar.Death));
        UndirectedEdge chord = Assert.Single(BarCycleEdges.GetEdgePairs(bar, relative.Ambient));
        Assert.Equal((0, 4), (chord.Lo, chord.Hi));
    }

    // ── R11 · the convenience entry stops below the truncated top ───────────────

    [Fact]
    public void ConvenienceEntry_ReportsNoTruncatedTopDegree_AtAnySkeleton()
    {
        const int n = 4;
        List<(int, int)> backbone = CycleEdges(n);
        var content = new List<(int, int, double)> { (0, 2, 2.5) };

        Assert.DoesNotContain(ConditionedFiltration.ComputeRelativeBarcode(n, backbone, content, maxDimension: 3).Bars, b => b.Dimension == 2);
        Assert.DoesNotContain(ConditionedFiltration.ComputeRelativeBarcode(n, backbone, content, maxDimension: 1).Bars, b => b.Dimension >= 1);
    }

    // ── R9 · the guard ──────────────────────────────────────────────────────────

    [Fact]
    public void PredicateNotClosedUnderFaces_Throws()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0), new Simplex(0.0, 1), new Simplex(1.0, 0, 1),
        });

        Assert.Throws<ArgumentException>(() => new RelativeFiltration(filtration, s => s.Dimension == 1));
    }

    // ── helpers ─────────────────────────────────────────────────────────────────

    static List<(int, int)> PathEdges(int n) =>
        Enumerable.Range(0, n - 1).Select(i => (i, i + 1)).ToList();

    static List<(int, int)> CycleEdges(int n)
    {
        List<(int, int)> edges = PathEdges(n);
        edges.Add((0, n - 1));
        return edges;
    }

    // A random spanning tree plus a few extra edges, so K₀ is connected and usually holds cycles.
    static List<(int, int)> RandomConnectedSpanning(int n, Random rng)
    {
        var edges = new HashSet<(int, int)>();
        for (int v = 1; v < n; v++)
            edges.Add(Ordered((rng.Next(v), v)));

        for (int extra = rng.Next(0, n); extra > 0; extra--)
        {
            int a = rng.Next(n), b = rng.Next(n);
            if (a != b) edges.Add(Ordered((a, b)));
        }

        return edges.ToList();
    }

    static List<(int, int)> RandomEdges(int n, Random rng, double p)
    {
        var edges = new List<(int, int)>();
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
                if (rng.NextDouble() < p)
                    edges.Add((a, b));
        return edges;
    }

    // Content on the pairs outside the backbone, at positive values from a small grid so ties occur.
    static List<(int, int, double)> RandomContent(int n, Random rng, List<(int, int)> backbone, double p)
    {
        HashSet<(int, int)> taken = backbone.Select(Ordered).ToHashSet();
        var content = new List<(int, int, double)>();
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
                if (!taken.Contains((a, b)) && rng.NextDouble() < p)
                    content.Add((a, b, rng.Next(1, 12)));
        return content;
    }

    static (int, int) Ordered((int, int) e) => e.Item1 < e.Item2 ? e : (e.Item2, e.Item1);

    static bool InK0(Simplex s, HashSet<int> vertices, HashSet<(int, int)> edges)
    {
        int[] v = s.Vertices;
        if (v.Length == 1) return vertices.Contains(v[0]);
        for (int a = 0; a < v.Length; a++)
            for (int b = a + 1; b < v.Length; b++)
                if (!edges.Contains(Ordered((v[a], v[b]))))
                    return false;
        return true;
    }

    // K plus an apex joined to every simplex of K₀, all at 0, where K₀ already is.
    static SimplicialFiltration Cone(SimplicialFiltration k, Func<Simplex, bool> inK0, int apex)
    {
        var cells = new List<Simplex>(k.Simplices) { new Simplex(0.0, apex) };
        foreach (Simplex s in k.Simplices)
            if (inK0(s))
                cells.Add(new Simplex(0.0, s.Vertices.Append(apex).ToArray()));
        return new SimplicialFiltration(cells, "coned");
    }

    static IEnumerable<Bar> Significant(Barcode barcode, int dimension) =>
        barcode.Bars.Where(b => b.Dimension == dimension && (b.IsInfinite || b.Death > b.Birth));

    // Significant bars as a sorted multiset, compared exactly.
    static (int, double, double)[] Canon(Barcode barcode, int? dimension = null) =>
        barcode.Bars
            .Where(b => (dimension is null || b.Dimension == dimension) && (b.IsInfinite || b.Death > b.Birth))
            .Select(b => (b.Dimension, b.Birth, b.Death))
            .OrderBy(x => x)
            .ToArray();
}
