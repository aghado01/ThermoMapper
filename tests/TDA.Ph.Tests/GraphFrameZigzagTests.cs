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
/// The frame builder: frame sequences through the zigzag engines, read back on the frames. Hand decompositions,
/// nested sequences against ordinary persistence and union-find, four-engine parity on churn, the read-back against
/// a direct decomposition of the frame-level module, and the Diamond Principle between the two transitions. Checks
/// and expectations: <c>issues/ph/ph-setup-brief.md</c> (Z1–Z11, D1–D3).
/// </summary>
public sealed class GraphFrameZigzagTests
{
    readonly ITestOutputHelper _output;

    public GraphFrameZigzagTests(ITestOutputHelper output) => _output = output;

    // ── Z1 · a merge then a split ───────────────────────────────────────────────

    [Theory]
    [InlineData(FrameTransition.Union)]
    [InlineData(FrameTransition.Intersection)]
    public void MergeThenSplit_RecordsTheSplitAsItsOwnBar(FrameTransition transition)
    {
        var built = GraphFrameZigzag.FromEdgeSets(2, Sets(Array.Empty<(int, int)>(), new[] { (0, 1) }, Array.Empty<(int, int)>()), transition);

        Assert.Equal(new[] { (0, 0, 0), (0, 0, 2), (0, 2, 2) }, Frames(built, GraphZigzag.Compute(built.Filtration)));
    }

    // ── Z2 · re-entry ───────────────────────────────────────────────────────────

    [Theory]
    [InlineData(FrameTransition.Union)]
    [InlineData(FrameTransition.Intersection)]
    public void EdgeLeavesAndReturns_KeepsTheComponentBookkeeping(FrameTransition transition)
    {
        var built = GraphFrameZigzag.FromEdgeSets(2, Sets(new[] { (0, 1) }, Array.Empty<(int, int)>(), new[] { (0, 1) }), transition);

        Assert.Equal(new[] { (0, 0, 2), (0, 1, 1) }, Frames(built, GraphZigzag.Compute(built.Filtration)));
    }

    // ── Z3 · nested frames are ordinary persistence ─────────────────────────────

    [Fact]
    public void NestedFrames_MatchOrdinaryPersistence_AndUnionFindMerges()
    {
        const int trials = 100;
        int merges = 0, cycles = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(3, 9), frames = rng.Next(1, 7);

            // Each candidate edge takes a first frame, or never enters.
            var firstFrame = new Dictionary<(int, int), int>();
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    if (rng.NextDouble() < 0.5)
                        firstFrame[(a, b)] = rng.Next(frames);

            var sets = new IReadOnlyList<(int i, int j)>[frames];
            for (int t = 0; t < frames; t++)
                sets[t] = firstFrame.Where(e => e.Value <= t).Select(e => e.Key).ToList();

            var built = GraphFrameZigzag.FromEdgeSets(n, sets);
            var read = Frames(built, GraphZigzag.Compute(built.Filtration, maxDimension: 1));

            // Ordinary persistence on the 1-complex with each edge born at its first frame.
            var cells = new List<Simplex>();
            for (int v = 0; v < n; v++) cells.Add(new Simplex(0.0, v));
            foreach (var ((a, b), first) in firstFrame) cells.Add(new Simplex(first, a, b));
            var expected = PersistentHomology.Compute(new SimplicialFiltration(cells))
                .Bars
                .Where(b => b.IsInfinite || b.Death > b.Birth)
                .Select(b => (b.Dimension, (int)b.Birth, b.IsInfinite ? frames - 1 : (int)b.Death - 1))
                .OrderBy(x => x)
                .ToArray();
            Assert.True(expected.SequenceEqual(read), $"seed {seed}: ordinary persistence");

            // Kruskal on the first frames: each merge after frame 0 ends one H₀ bar at the frame before it.
            var uf = new UnionFind(n);
            var deaths = new List<int>();
            foreach (var ((a, b), first) in firstFrame.OrderBy(e => e.Value))
            {
                if (uf.Find(a) == uf.Find(b)) continue;
                uf.Union(a, b);
                if (first > 0) deaths.Add(first - 1);
            }

            var finiteH0 = read.Where(b => b.Item1 == 0 && b.Item3 < frames - 1).Select(b => b.Item3).OrderBy(x => x);
            var endingBeforeLast = deaths.Where(d => d < frames - 1).OrderBy(x => x);
            Assert.True(endingBeforeLast.SequenceEqual(finiteH0), $"seed {seed}: union-find");

            merges += deaths.Count;
            cycles += read.Count(b => b.Item1 == 1);
        }

        _output.WriteLine($"H0 merges read: {merges}, H1 bars: {cycles}");
        Assert.True(merges > 0 && cycles > 0);
    }

    // ── Z4 · four engines on edge churn with re-entry ───────────────────────────

    [Fact]
    public void EdgeChurn_FourEnginesAgree()
    {
        const int trials = 10; // the naive oracle's cost; ran green at 30 (issues/ph/ph-setup-brief.md)
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            var built = GraphFrameZigzag.FromEdgeSets(5, RandomEdgeSets(5, 5, rng), rng.Next(2) == 0 ? FrameTransition.Union : FrameTransition.Intersection);
            AssertEnginesAgree(built.Filtration, seed);
        }
    }

    // ── Z5 · four engines with vertices entering and leaving ────────────────────

    [Fact]
    public void VertexChurn_FourEnginesAgree()
    {
        const int trials = 10; // the naive oracle's cost; ran green at 30 (issues/ph/ph-setup-brief.md)
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            var built = GraphFrameZigzag.Build(6, RandomFrames(6, 5, rng), rng.Next(2) == 0 ? FrameTransition.Union : FrameTransition.Intersection);
            AssertEnginesAgree(built.Filtration, seed);
        }
    }

    // ── Z11 · the read-back against a direct decomposition of the frame-level module ───

    [Fact]
    public void ReadBack_MatchesDirectDecompositionOfTheCoarseModule()
    {
        const int trials = 40;
        int frameBars = 0, transients = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            var transition = seed % 2 == 0 ? FrameTransition.Union : FrameTransition.Intersection;
            GraphFrame[] frames = RandomFrames(5, 4, rng);

            var built = GraphFrameZigzag.Build(5, frames, transition);
            Barcode barcode = GraphZigzag.Compute(built.Filtration, maxDimension: 1);
            var readFrames = Frames(built, barcode);
            var readTransients = built.TransientBars(barcode).Select(b => (b.Dimension, b.AfterFrame)).OrderBy(x => x).ToArray();

            var (expectedFrames, expectedTransients) = CoarseDecomposition(5, frames, transition);

            Assert.True(expectedFrames.SequenceEqual(readFrames), $"seed {seed}: frame bars");
            Assert.True(expectedTransients.SequenceEqual(readTransients), $"seed {seed}: transients");

            frameBars += readFrames.Length;
            transients += readTransients.Length;
        }

        _output.WriteLine($"frame bars: {frameBars}, transients: {transients}");
        Assert.True(frameBars > 0 && transients > 0);
    }

    // The module on frames and between-frames complexes alone, F₀ – X₀ – F₁ – X₁ – …, decomposed with no single-cell
    // refinement and no read-back: even positions are frames, odd ones the union or intersection between them.
    static ((int, int, int)[] Frames, (int, int)[] Transients) CoarseDecomposition(int n, GraphFrame[] frames, FrameTransition transition)
    {
        var edgeIds = new Dictionary<(int, int), int>();
        foreach (GraphFrame frame in frames)
            foreach (var (i, j) in frame.Edges)
                edgeIds.TryAdd((i, j), n + edgeIds.Count);

        int cellCount = n + edgeIds.Count;
        var dimOf = new int[cellCount];
        var bndOf = new int[cellCount][];
        for (int c = 0; c < cellCount; c++) bndOf[c] = Array.Empty<int>();
        foreach (var ((i, j), id) in edgeIds)
        {
            dimOf[id] = 1;
            bndOf[id] = new[] { i, j };
        }

        HashSet<int> Cells(GraphFrame frame)
        {
            var cells = new HashSet<int>(frame.Vertices);
            foreach (var edge in frame.Edges) cells.Add(edgeIds[(edge.i, edge.j)]);
            return cells;
        }

        int[] identity = Enumerable.Range(0, cellCount).ToArray();
        bool viaUnion = transition == FrameTransition.Union;
        var complexes = new List<HashSet<int>>();
        var steps = new List<ZigzagMapBarcode.Step>();
        for (int t = 0; t < frames.Length; t++)
        {
            complexes.Add(Cells(frames[t]));
            if (t == frames.Length - 1) break;

            var between = Cells(frames[t]);
            if (viaUnion) between.UnionWith(Cells(frames[t + 1]));
            else between.IntersectWith(Cells(frames[t + 1]));
            complexes.Add(between);

            steps.Add(new ZigzagMapBarcode.Step(viaUnion, identity));
            steps.Add(new ZigzagMapBarcode.Step(!viaUnion, identity));
        }

        Barcode coarse = ZigzagMapBarcode.Compute(complexes, dimOf, bndOf, cellCount, steps, maxDimension: 1);

        var frameBars = new List<(int, int, int)>();
        var transients = new List<(int, int)>();
        foreach (Bar bar in coarse.Bars)
        {
            int first = (int)bar.Birth + 1, last = (int)bar.Death;
            int firstFrame = (first + 1) / 2, lastFrame = last / 2;
            if (firstFrame <= lastFrame) frameBars.Add((bar.Dimension, firstFrame, lastFrame));
            else transients.Add((bar.Dimension, (first - 1) / 2));
        }

        return (frameBars.OrderBy(x => x).ToArray(), transients.OrderBy(x => x).ToArray());
    }

    // ── Z6 · one frame ──────────────────────────────────────────────────────────

    [Fact]
    public void SingleFrame_EachClassIsOneBarAtThatFrame()
    {
        // A triangle's cycle plus an isolated vertex: two components, one independent cycle.
        var built = GraphFrameZigzag.FromEdgeSets(4, Sets(new[] { (0, 1), (1, 2), (0, 2) }));

        Assert.Equal(new[] { (0, 0, 0), (0, 0, 0), (1, 0, 0) }, Frames(built, GraphZigzag.Compute(built.Filtration, maxDimension: 1)));
    }

    [Fact]
    public void NoFrames_IsEmpty()
    {
        var built = GraphFrameZigzag.FromEdgeSets(3, Array.Empty<IReadOnlyList<(int i, int j)>>());

        Assert.Empty(built.Filtration);
        Assert.Equal(0, built.FrameCount);
    }

    // ── Z7 · a conditioned sweep ────────────────────────────────────────────────

    [Fact]
    public void ConditionedSweep_ReturnReadsInH1_WithOneH0Bar()
    {
        const int n = 6;
        var backbone = Enumerable.Range(0, n - 1).Select(i => (i, i + 1)).ToList();
        var chord = new List<(int i, int j, double d)> { (0, 5, 1.0) };
        var content = Enumerable.Repeat<IReadOnlyList<(int i, int j, double d)>>(chord, 5).ToList();
        var slack = new[] { 0.0, 2.0, 2.0, 0.0, 2.0 };

        var built = ConditionedFiltration.BuildSweep(n, backbone, content, slack);

        Assert.Equal(new[] { (0, 0, 4), (1, 1, 2), (1, 4, 4) }, Frames(built, GraphZigzag.Compute(built.Filtration, maxDimension: 1)));
    }

    [Fact]
    public void ConditionedSweep_MismatchedLengths_Throws()
    {
        var content = new List<IReadOnlyList<(int i, int j, double d)>> { new List<(int, int, double)>() };

        Assert.Throws<ArgumentException>(() => ConditionedFiltration.BuildSweep(2, new List<(int, int)>(), content, new[] { 1.0, 2.0 }));
    }

    // ── Z8 · a field on one substrate, thresholds along a path ──────────────────

    [Fact]
    public void EdgeField_AdmitsByThreshold_OnUpperSlotsOnly()
    {
        CsrGraph substrate = Triangle();

        // Frame 0 at θ = 0.5 holds (0,1); frame 1 at θ = 0.2 adds (1,2); frame 2 at θ = 0.9 holds (0,2) alone.
        double[] Field(double e01, double e12, double e02) => FieldOf(substrate, e01, e12, e02);
        var fields = new[] { Field(0.8, 0.1, 0.1), Field(0.8, 0.3, 0.1), Field(0.1, 0.1, 0.95) };
        var built = GraphFrameZigzag.FromEdgeField(substrate, fields, new[] { 0.5, 0.2, 0.9 });

        var expected = GraphFrameZigzag.FromEdgeSets(3, Sets(new[] { (0, 1) }, new[] { (0, 1), (1, 2) }, new[] { (0, 2) }));
        Assert.Equal(Steps(expected), Steps(built));
    }

    [Fact]
    public void EdgeField_WrongLength_Throws()
    {
        CsrGraph substrate = Triangle();

        Assert.Throws<ArgumentException>(() =>
            GraphFrameZigzag.FromEdgeField(substrate, new[] { new double[1] }, new[] { 0.5 }));
        Assert.Throws<ArgumentException>(() =>
            GraphFrameZigzag.FromEdgeField(substrate, new[] { new double[substrate.Targets.Length] }, new[] { 0.5, 0.6 }));
    }

    // ── Z9 · the CSR boundary ───────────────────────────────────────────────────

    [Fact]
    public void FromGraphs_MatchesFromEdgeSets()
    {
        var a = new[] { (0, 1), (2, 3) };
        var b = new[] { (1, 2), (2, 3) };

        var viaGraphs = GraphFrameZigzag.FromGraphs(new[] { Graph(4, a), Graph(4, b), Graph(4, a) });
        var viaSets = GraphFrameZigzag.FromEdgeSets(4, Sets(a, b, a));

        Assert.Equal(Steps(viaSets), Steps(viaGraphs));
    }

    [Fact]
    public void FromGraphs_DifferentVertexCounts_Throws()
    {
        Assert.Throws<ArgumentException>(() => GraphFrameZigzag.FromGraphs(new[] { Graph(3, new[] { (0, 1) }), Graph(4, new[] { (0, 1) }) }));
    }

    // ── Z10 · bad frames ────────────────────────────────────────────────────────

    [Fact]
    public void EdgeEndpointAbsentFromFrame_Throws()
    {
        var frame = new GraphFrame(new[] { 0, 1 }, new[] { (0, 2) });

        Assert.Throws<ArgumentException>(() => GraphFrameZigzag.Build(3, new[] { frame }));
    }

    [Fact]
    public void VertexOutOfRange_Throws()
    {
        var frame = new GraphFrame(new[] { 0, 7 }, Array.Empty<(int, int)>());

        Assert.Throws<ArgumentOutOfRangeException>(() => GraphFrameZigzag.Build(3, new[] { frame }));
    }

    [Fact]
    public void SelfLoopsAndDuplicates_Collapse()
    {
        var clean = GraphFrameZigzag.FromEdgeSets(3, Sets(new[] { (0, 1) }));
        var messy = GraphFrameZigzag.FromEdgeSets(3, Sets(new[] { (0, 1), (1, 0), (2, 2) }));

        Assert.Equal(Steps(clean), Steps(messy));
    }

    // ── D1, D2 · the Diamond Principle ──────────────────────────────────────────

    [Fact]
    public void Transitions_AgreeOnFrames_TransientsShiftOneDegree()
    {
        const int trials = 200;
        int totalTransients = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(3, 9), frames = rng.Next(2, 7);
            IReadOnlyList<(int i, int j)>[] sets = RandomEdgeSets(n, frames, rng);

            var viaUnion = GraphFrameZigzag.FromEdgeSets(n, sets, FrameTransition.Union);
            var viaIntersection = GraphFrameZigzag.FromEdgeSets(n, sets, FrameTransition.Intersection);
            Barcode union = GraphZigzag.Compute(viaUnion.Filtration, maxDimension: 1);
            Barcode intersection = GraphZigzag.Compute(viaIntersection.Filtration, maxDimension: 1);

            Assert.True(Frames(viaUnion, union).SequenceEqual(Frames(viaIntersection, intersection)), $"seed {seed}: frames");

            var transientsU = viaUnion.TransientBars(union);
            var transientsI = viaIntersection.TransientBars(intersection);
            Assert.True(transientsU.All(b => b.Dimension == 1), $"seed {seed}: an H₀ union transient");
            Assert.True(transientsI.All(b => b.Dimension == 0), $"seed {seed}: an H₁ intersection transient");
            Assert.True(
                transientsU.Select(b => b.AfterFrame).OrderBy(x => x).SequenceEqual(transientsI.Select(b => b.AfterFrame).OrderBy(x => x)),
                $"seed {seed}: transients");
            totalTransients += transientsU.Count;
        }

        _output.WriteLine($"transients: {totalTransients}");
        Assert.True(totalTransients > 0);
    }

    // ── D3 · one diamond by hand ────────────────────────────────────────────────

    [Fact]
    public void OneDiamond_UnionCycle_IntersectionComponent()
    {
        var sets = Sets(new[] { (0, 1), (1, 2) }, new[] { (0, 2), (1, 2) });

        var viaUnion = GraphFrameZigzag.FromEdgeSets(3, sets, FrameTransition.Union);
        var viaIntersection = GraphFrameZigzag.FromEdgeSets(3, sets, FrameTransition.Intersection);
        Barcode union = GraphZigzag.Compute(viaUnion.Filtration, maxDimension: 1);
        Barcode intersection = GraphZigzag.Compute(viaIntersection.Filtration, maxDimension: 1);

        Assert.Equal(new[] { (0, 0, 1) }, Frames(viaUnion, union));
        Assert.Equal(new[] { (0, 0, 1) }, Frames(viaIntersection, intersection));

        TransientBar cycle = Assert.Single(viaUnion.TransientBars(union));
        Assert.Equal((1, 0), (cycle.Dimension, cycle.AfterFrame));

        TransientBar component = Assert.Single(viaIntersection.TransientBars(intersection));
        Assert.Equal((0, 0), (component.Dimension, component.AfterFrame));
    }

    // ── helpers ─────────────────────────────────────────────────────────────────

    static IReadOnlyList<(int i, int j)>[] Sets(params (int, int)[][] frames) =>
        frames.Select(f => (IReadOnlyList<(int i, int j)>)f.Select(e => (e.Item1, e.Item2)).ToList()).ToArray();

    static IReadOnlyList<(int i, int j)>[] RandomEdgeSets(int n, int frames, Random rng)
    {
        var sets = new IReadOnlyList<(int i, int j)>[frames];
        for (int t = 0; t < frames; t++)
        {
            var edges = new List<(int i, int j)>();
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    if (rng.NextDouble() < 0.4)
                        edges.Add((a, b));
            sets[t] = edges;
        }
        return sets;
    }

    static GraphFrame[] RandomFrames(int n, int frames, Random rng)
    {
        var result = new GraphFrame[frames];
        for (int t = 0; t < frames; t++)
        {
            int[] vertices = Enumerable.Range(0, n).Where(_ => rng.NextDouble() < 0.7).ToArray();
            var edges = new List<(int i, int j)>();
            for (int a = 0; a < vertices.Length; a++)
                for (int b = a + 1; b < vertices.Length; b++)
                    if (rng.NextDouble() < 0.4)
                        edges.Add((vertices[a], vertices[b]));
            result[t] = new GraphFrame(vertices, edges);
        }
        return result;
    }

    static void AssertEnginesAgree(ZigzagFiltration f, int seed)
    {
        var reference = Step(GraphZigzag.Compute(f, maxDimension: 1, algorithm: GraphZigzagAlgorithm.Reference));
        var fast = Step(GraphZigzag.Compute(f, maxDimension: 1, algorithm: GraphZigzagAlgorithm.Fast));
        var cone = Step(FastZigzag.Compute(f, maxDimension: 1));
        var naive = Step(ZigzagBarcodeNaive.Compute(f, maxDimension: 1));

        Assert.True(naive.SequenceEqual(reference), $"seed {seed}: reference");
        Assert.True(naive.SequenceEqual(fast), $"seed {seed}: fast");
        Assert.True(naive.SequenceEqual(cone), $"seed {seed}: FastZigzag");
    }

    static (int, double, double)[] Step(Barcode barcode) =>
        barcode.Bars.Select(b => (b.Dimension, b.Birth, b.Death)).OrderBy(x => x).ToArray();

    static (int, int, int)[] Frames(GraphFrameFiltration built, Barcode barcode) =>
        built.FrameBars(barcode).Select(b => (b.Dimension, b.FirstFrame, b.LastFrame)).OrderBy(x => x).ToArray();

    static (int, ZigzagDirection, string)[] Steps(GraphFrameFiltration built) =>
        built.Filtration
            .Select(s => (s.GlobalCellId, s.Direction, s.BoundaryAtAdd is null ? "" : string.Join(",", s.BoundaryAtAdd)))
            .ToArray();

    static CsrGraph Graph(int n, (int, int)[] edges) =>
        CsrGraph.FromEdges(edges.Select(e => new Edge(e.Item1, e.Item2, 1.0)).ToArray(), n);

    static CsrGraph Triangle() => Graph(3, new[] { (0, 1), (1, 2), (0, 2) });

    // A field on the triangle with its values at j > i slots and an admissible-looking +∞ on every mirror slot.
    static double[] FieldOf(CsrGraph g, double e01, double e12, double e02)
    {
        var field = new double[g.Targets.Length];
        for (int u = 0; u < g.NodeCount; u++)
            for (int e = g.RowPointers[u]; e < g.RowPointers[u + 1]; e++)
            {
                int v = g.Targets[e];
                field[e] = v <= u ? double.PositiveInfinity : (u, v) switch
                {
                    (0, 1) => e01,
                    (1, 2) => e12,
                    (0, 2) => e02,
                    _ => throw new InvalidOperationException(),
                };
            }
        return field;
    }
}
