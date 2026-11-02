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
/// a direct decomposition of the frame-level module, the Diamond Principle between the two transitions, and flag frames
/// where a cycle can be filled. Checks and expectations: <c>issues/ph/ph-setup-brief.md</c> (Z1–Z11, D1–D3) and
/// <c>issues/ph/p1d-filled-h1-brief.md</c> (F1–F5, Z12, N1, B1, E1, C1, D1, D2, H1).
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

    // ── C1 · the read-back against a direct decomposition of the frame-level module ───

    [Theory]
    [InlineData(FrameComplex.Graph)]
    [InlineData(FrameComplex.Flag)]
    public void ReadBack_MatchesDirectDecompositionOfTheCoarseModule(FrameComplex complex)
    {
        const int trials = 80;
        int frameBars = 0, transients = 0, h1EndedByArrival = 0, h1StartedByDeparture = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            var transition = seed % 2 == 0 ? FrameTransition.Union : FrameTransition.Intersection;
            GraphFrame[] frames = seed % 4 < 2 ? RandomFrames(6, 4, rng, edgeProbability: 0.6) : EvolvingFrames(6, 5, rng, vertexFlip: 0.1);

            var built = GraphFrameZigzag.Build(6, frames, transition, complex);
            Barcode barcode = Read(built, maxDimension: 1);
            var readFrames = FramesWithEvents(built, barcode);
            var readTransients = built.TransientBars(barcode).Select(b => (b.Dimension, b.AfterFrame)).OrderBy(x => x).ToArray();

            var (expectedFrames, expectedTransients) = CoarseDecomposition(6, frames, transition, complex);

            Assert.True(expectedFrames.SequenceEqual(readFrames), $"seed {seed}: frame bars");
            Assert.True(expectedTransients.SequenceEqual(readTransients), $"seed {seed}: transients");

            frameBars += readFrames.Length;
            transients += readTransients.Length;
            h1EndedByArrival += readFrames.Count(b => b.Item1 == 1 && b.Item5 == FrameEvent.Arrival);
            h1StartedByDeparture += readFrames.Count(b => b.Item1 == 1 && b.Item4 == FrameEvent.Departure);
        }

        _output.WriteLine($"{complex}: frame bars {frameBars}, transients {transients}, H1 ended by arrival {h1EndedByArrival}, H1 started by departure {h1StartedByDeparture}");
        Assert.True(frameBars > 0 && transients > 0);

        // A fill and a hole opening exist only with triangles; an unfilled cycle cannot do either.
        if (complex == FrameComplex.Flag) Assert.True(h1EndedByArrival > 0 && h1StartedByDeparture > 0);
        else Assert.True(h1EndedByArrival == 0 && h1StartedByDeparture == 0);
    }

    // The module on frames and between-frames complexes alone, F₀ – X₀ – F₁ – X₁ – …, decomposed with no single-cell
    // refinement and no read-back: even positions are frames, odd ones the union or intersection between them.
    static ((int, int, int, FrameEvent, FrameEvent)[] Frames, (int, int)[] Transients) CoarseDecomposition(
        int n, GraphFrame[] frames, FrameTransition transition, FrameComplex complex)
    {
        var edgeIds = new Dictionary<(int, int), int>();
        foreach (GraphFrame frame in frames)
            foreach (var (i, j) in frame.Edges)
                edgeIds.TryAdd((i, j), n + edgeIds.Count);

        var triangleIds = new Dictionary<(int, int, int), int>();
        List<(int, int, int)> TrianglesOf(GraphFrame frame)
        {
            var found = new List<(int, int, int)>();
            if (complex != FrameComplex.Flag) return found;

            HashSet<(int, int)> edges = frame.Edges.Select(e => (e.i, e.j)).ToHashSet();
            int[] v = frame.Vertices.OrderBy(x => x).ToArray();
            for (int a = 0; a < v.Length; a++)
                for (int b = a + 1; b < v.Length; b++)
                    for (int c = b + 1; c < v.Length; c++)
                        if (edges.Contains((v[a], v[b])) && edges.Contains((v[a], v[c])) && edges.Contains((v[b], v[c])))
                            found.Add((v[a], v[b], v[c]));
            return found;
        }

        foreach (GraphFrame frame in frames)
            foreach (var triangle in TrianglesOf(frame))
                triangleIds.TryAdd(triangle, n + edgeIds.Count + triangleIds.Count);

        int cellCount = n + edgeIds.Count + triangleIds.Count;
        var dimOf = new int[cellCount];
        var bndOf = new int[cellCount][];
        for (int c = 0; c < cellCount; c++) bndOf[c] = Array.Empty<int>();
        foreach (var ((i, j), id) in edgeIds)
        {
            dimOf[id] = 1;
            bndOf[id] = new[] { i, j };
        }
        foreach (var ((a, b, c), id) in triangleIds)
        {
            dimOf[id] = 2;
            bndOf[id] = new[] { edgeIds[(a, b)], edgeIds[(a, c)], edgeIds[(b, c)] };
        }

        HashSet<int> Cells(GraphFrame frame)
        {
            var cells = new HashSet<int>(frame.Vertices);
            foreach (var edge in frame.Edges) cells.Add(edgeIds[(edge.i, edge.j)]);
            foreach (var triangle in TrianglesOf(frame)) cells.Add(triangleIds[triangle]);
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

        // A forward arrow is cells arriving and a backward one cells leaving, so the module's own arrow into a bar's
        // first position and out of its last name its start and end.
        var frameBars = new List<(int, int, int, FrameEvent, FrameEvent)>();
        var transients = new List<(int, int)>();
        foreach (Bar bar in coarse.Bars)
        {
            int first = (int)bar.Birth + 1, last = (int)bar.Death;
            int firstFrame = (first + 1) / 2, lastFrame = last / 2;
            if (firstFrame > lastFrame)
            {
                transients.Add((bar.Dimension, (first - 1) / 2));
                continue;
            }

            FrameEvent start = first == 0 ? FrameEvent.Unobserved
                : bar.BirthEnd == IntervalEnd.Closed ? FrameEvent.Arrival : FrameEvent.Departure;
            FrameEvent end = last == steps.Count ? FrameEvent.Unobserved
                : bar.DeathEnd == IntervalEnd.Closed ? FrameEvent.Departure : FrameEvent.Arrival;
            frameBars.Add((bar.Dimension, firstFrame, lastFrame, start, end));
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

            Assert.True(FramesWithEvents(viaUnion, union).SequenceEqual(FramesWithEvents(viaIntersection, intersection)), $"seed {seed}: frames");

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

    // ── F1–F5 · flag frames by hand ─────────────────────────────────────────────

    const FrameEvent U = FrameEvent.Unobserved, A = FrameEvent.Arrival, D = FrameEvent.Departure;

    static readonly (int, int)[] C4 = { (0, 1), (1, 2), (2, 3), (0, 3) };
    static readonly (int, int)[] C4Chord = { (0, 1), (1, 2), (2, 3), (0, 3), (0, 2) };
    static readonly (int, int)[] C4Broken = { (0, 1), (1, 2), (2, 3) };

    static (int, int, int, FrameEvent, FrameEvent)[] ReadFlag(FrameTransition transition, FrameComplex complex, params (int, int)[][] edgeSets)
    {
        var built = GraphFrameZigzag.FromEdgeSets(4, Sets(edgeSets), transition, complex);
        return FramesWithEvents(built, Read(built, maxDimension: 1));
    }

    static (int, int, int, FrameEvent, FrameEvent)[] Sorted(params (int, int, int, FrameEvent, FrameEvent)[] bars) =>
        bars.OrderBy(x => x).ToArray();

    [Theory]
    [InlineData(FrameTransition.Union)]
    [InlineData(FrameTransition.Intersection)]
    public void F1_CycleFilledByTriangles_EndsByArrival(FrameTransition transition)
    {
        Assert.Equal(Sorted((0, 0, 1, U, U), (1, 0, 0, U, A)), ReadFlag(transition, FrameComplex.Flag, C4, C4Chord));
    }

    [Theory]
    [InlineData(FrameTransition.Union)]
    [InlineData(FrameTransition.Intersection)]
    public void F2_CycleBrokenByEdgeLeaving_EndsByDeparture(FrameTransition transition)
    {
        Assert.Equal(Sorted((0, 0, 1, U, U), (1, 0, 0, U, D)), ReadFlag(transition, FrameComplex.Flag, C4, C4Broken));
    }

    [Theory]
    [InlineData(FrameTransition.Union)]
    [InlineData(FrameTransition.Intersection)]
    public void F3_TrianglesLeaving_OpenAHole_StartingByDeparture(FrameTransition transition)
    {
        Assert.Equal(Sorted((0, 0, 1, U, U), (1, 1, 1, D, U)), ReadFlag(transition, FrameComplex.Flag, C4Chord, C4));
    }

    [Theory]
    [InlineData(FrameTransition.Union)]
    [InlineData(FrameTransition.Intersection)]
    public void F4_FilledThenReopened_ReadsAsTwoBars_AndUnfilledAsOne(FrameTransition transition)
    {
        Assert.Equal(
            Sorted((0, 0, 2, U, U), (1, 0, 0, U, A), (1, 2, 2, D, U)),
            ReadFlag(transition, FrameComplex.Flag, C4, C4Chord, C4));

        Assert.Equal(
            Sorted((0, 0, 2, U, U), (1, 0, 2, U, U), (1, 1, 1, A, D)),
            ReadFlag(transition, FrameComplex.Graph, C4, C4Chord, C4));
    }

    [Theory]
    [InlineData(FrameTransition.Union)]
    [InlineData(FrameTransition.Intersection)]
    public void F5_ConditionedSweep_OneReturnDiesBothWays(FrameTransition transition)
    {
        var backbone = new List<(int, int)> { (0, 1), (1, 2), (2, 3) };
        var chord = (0, 3, 1.0);
        var filler = (0, 2, 1.0);
        var content = new List<IReadOnlyList<(int i, int j, double d)>>
        {
            new List<(int, int, double)> { chord },
            new List<(int, int, double)> { chord, filler },
            new List<(int, int, double)> { chord },
            new List<(int, int, double)>(),
        };

        var built = ConditionedFiltration.BuildSweep(4, backbone, content, new[] { 2.0, 2.0, 2.0, 2.0 }, transition, FrameComplex.Flag);

        Assert.Equal(
            Sorted((0, 0, 3, U, U), (1, 0, 0, U, A), (1, 2, 2, D, D)),
            FramesWithEvents(built, Read(built, maxDimension: 1)));
    }

    // ── Z12 · the graph engines refuse flag frames ──────────────────────────────

    [Theory]
    [InlineData(GraphZigzagAlgorithm.Reference)]
    [InlineData(GraphZigzagAlgorithm.Fast)]
    public void GraphEngines_OnFlagFrames_Throw(GraphZigzagAlgorithm algorithm)
    {
        var built = GraphFrameZigzag.FromEdgeSets(4, Sets(C4Chord), FrameTransition.Union, FrameComplex.Flag);

        Assert.ThrowsAny<Exception>(() => GraphZigzag.Compute(built.Filtration, maxDimension: 1, algorithm: algorithm));
    }

    // ── N1 · nested flag frames are ordinary persistence ────────────────────────

    [Fact]
    public void NestedFlagFrames_MatchOrdinaryPersistence()
    {
        const int trials = 100;
        int h1Finite = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(4, 9), frames = rng.Next(1, 7);

            var firstFrame = new Dictionary<(int, int), int>();
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    if (rng.NextDouble() < 0.6)
                        firstFrame[(a, b)] = rng.Next(frames);

            var sets = new IReadOnlyList<(int i, int j)>[frames];
            for (int t = 0; t < frames; t++)
                sets[t] = firstFrame.Where(e => e.Value <= t).Select(e => e.Key).ToList();

            var built = GraphFrameZigzag.FromEdgeSets(n, sets, FrameTransition.Union, FrameComplex.Flag);
            var read = Frames(built, FastZigzag.Compute(built.Filtration, maxDimension: 1));

            CsrGraph graph = CsrGraph.FromEdges(firstFrame.Select(e => new Edge(e.Key.Item1, e.Key.Item2, e.Value)).ToArray(), n);
            var filtration = RipsFiltration.GraphRips(graph, FiltrationWeights.RawDistance, 2);
            var expected = PersistentHomology.Compute(filtration, maxDimension: 1)
                .Bars
                .Where(b => b.IsInfinite || b.Death > b.Birth)
                .Select(b => (b.Dimension, (int)b.Birth, b.IsInfinite ? frames - 1 : (int)b.Death - 1))
                .OrderBy(x => x)
                .ToArray();

            Assert.True(expected.SequenceEqual(read), $"seed {seed}");
            h1Finite += read.Count(b => b.Item1 == 1 && b.Item3 < frames - 1);
        }

        _output.WriteLine($"H1 bars that ended before the last frame: {h1Finite}");
        Assert.True(h1Finite > 0);
    }

    // ── B1 · each frame's Betti numbers, each frame reduced alone ───────────────

    [Theory]
    [InlineData(FrameComplex.Graph)]
    [InlineData(FrameComplex.Flag)]
    public void BarsCoveringAFrame_EqualThatFramesBettiNumbers(FrameComplex complex)
    {
        const int trials = 100;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(3, 8), frames = rng.Next(1, 6);
            IReadOnlyList<(int i, int j)>[] sets = RandomEdgeSets(n, frames, rng, 0.5);

            var built = GraphFrameZigzag.FromEdgeSets(n, sets, FrameTransition.Union, complex);
            var read = Frames(built, Read(built, maxDimension: 1));

            for (int t = 0; t < frames; t++)
            {
                CsrGraph graph = CsrGraph.FromEdges(sets[t].Select(e => new Edge(e.i, e.j, 1.0)).ToArray(), n);
                var alone = RipsFiltration.GraphRips(graph, FiltrationWeights.RawDistance, complex == FrameComplex.Flag ? 2 : 1);
                Barcode betti = PersistentHomology.Compute(alone, maxDimension: 1);

                for (int dim = 0; dim <= 1; dim++)
                {
                    int expected = betti.Bars.Count(b => b.Dimension == dim && b.IsInfinite);
                    int covering = read.Count(b => b.Item1 == dim && b.Item2 <= t && t <= b.Item3);
                    Assert.True(expected == covering, $"seed {seed}, frame {t}, H{dim}: {expected} against {covering}");
                }
            }
        }
    }

    // ── E1 · three engines on flag frames ───────────────────────────────────────

    [Fact]
    public void FlagFrames_ThreeEnginesAgree()
    {
        const int trials = 10;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            var built = GraphFrameZigzag.FromEdgeSets(
                5, RandomEdgeSets(5, 3, rng, 0.6),
                rng.Next(2) == 0 ? FrameTransition.Union : FrameTransition.Intersection,
                FrameComplex.Flag);

            var cone = Step(FastZigzag.Compute(built.Filtration, maxDimension: 2));
            var naive = Step(ZigzagBarcodeNaive.Compute(built.Filtration, maxDimension: 2));
            var map = Step(ZigzagMapBarcode.ComputeFromZigzag(built.Filtration, maxDimension: 2));

            Assert.True(naive.SequenceEqual(cone), $"seed {seed}: FastZigzag");
            Assert.True(naive.SequenceEqual(map), $"seed {seed}: ZigzagMapBarcode");
        }
    }

    // ── D1, D2 · the Diamond Principle on flag frames ───────────────────────────

    [Fact]
    public void FlagFrames_TransitionsAgreeOnBarsAndEvents_TransientsShiftOneDegree()
    {
        const int trials = 200;
        var byDegree = new int[3];
        int h1EndedByArrival = 0, h1EndedByDeparture = 0, h1StartedByDeparture = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(4, 9), frames = rng.Next(2, 6);
            IReadOnlyList<(int i, int j)>[] sets = seed % 2 == 0
                ? RandomEdgeSets(n, frames, rng, 0.5)
                : EvolvingFrames(n, frames, rng, vertexFlip: 0.0).Select(f => f.Edges).ToArray();

            var viaUnion = GraphFrameZigzag.FromEdgeSets(n, sets, FrameTransition.Union, FrameComplex.Flag);
            var viaIntersection = GraphFrameZigzag.FromEdgeSets(n, sets, FrameTransition.Intersection, FrameComplex.Flag);
            Barcode union = FastZigzag.Compute(viaUnion.Filtration, maxDimension: 2);
            Barcode intersection = FastZigzag.Compute(viaIntersection.Filtration, maxDimension: 2);

            var read = FramesWithEvents(viaUnion, union);
            Assert.True(read.SequenceEqual(FramesWithEvents(viaIntersection, intersection)), $"seed {seed}: frames");

            h1EndedByArrival += read.Count(b => b.Item1 == 1 && b.Item5 == FrameEvent.Arrival);
            h1EndedByDeparture += read.Count(b => b.Item1 == 1 && b.Item5 == FrameEvent.Departure);
            h1StartedByDeparture += read.Count(b => b.Item1 == 1 && b.Item4 == FrameEvent.Departure);

            var transientsU = viaUnion.TransientBars(union);
            var transientsI = viaIntersection.TransientBars(intersection);
            Assert.True(transientsU.All(b => b.Dimension >= 1), $"seed {seed}: an H0 union transient");
            Assert.True(transientsI.All(b => b.Dimension <= 1), $"seed {seed}: an H2 intersection transient");
            Assert.True(
                transientsU.Select(b => (b.Dimension - 1, b.AfterFrame)).OrderBy(x => x)
                    .SequenceEqual(transientsI.Select(b => (b.Dimension, b.AfterFrame)).OrderBy(x => x)),
                $"seed {seed}: transients");

            foreach (var bar in transientsU) byDegree[bar.Dimension]++;
        }

        _output.WriteLine($"union transients in H1: {byDegree[1]}, in H2: {byDegree[2]}; H1 ended by arrival {h1EndedByArrival}, by departure {h1EndedByDeparture}; started by departure {h1StartedByDeparture}");
        Assert.True(byDegree[1] > 0 && byDegree[2] > 0);
        Assert.True(h1EndedByArrival > 0 && h1EndedByDeparture > 0 && h1StartedByDeparture > 0);
    }

    // ── H1 · triangles do not change connectivity ───────────────────────────────

    [Fact]
    public void FlagAndGraphFrames_AgreeInH0()
    {
        const int trials = 100;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(3, 8), frames = rng.Next(1, 6);
            IReadOnlyList<(int i, int j)>[] sets = RandomEdgeSets(n, frames, rng, 0.5);

            var flag = GraphFrameZigzag.FromEdgeSets(n, sets, FrameTransition.Union, FrameComplex.Flag);
            var graph = GraphFrameZigzag.FromEdgeSets(n, sets, FrameTransition.Union, FrameComplex.Graph);

            Assert.True(
                FramesWithEvents(graph, GraphZigzag.Compute(graph.Filtration))
                    .SequenceEqual(FramesWithEvents(flag, FastZigzag.Compute(flag.Filtration, maxDimension: 0))),
                $"seed {seed}");
        }
    }

    // ── helpers ─────────────────────────────────────────────────────────────────

    static IReadOnlyList<(int i, int j)>[] Sets(params (int, int)[][] frames) =>
        frames.Select(f => (IReadOnlyList<(int i, int j)>)f.Select(e => (e.Item1, e.Item2)).ToList()).ToArray();

    static IReadOnlyList<(int i, int j)>[] RandomEdgeSets(int n, int frames, Random rng, double edgeProbability = 0.4)
    {
        var sets = new IReadOnlyList<(int i, int j)>[frames];
        for (int t = 0; t < frames; t++)
        {
            var edges = new List<(int i, int j)>();
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    if (rng.NextDouble() < edgeProbability)
                        edges.Add((a, b));
            sets[t] = edges;
        }
        return sets;
    }

    static GraphFrame[] RandomFrames(int n, int frames, Random rng, double edgeProbability = 0.4)
    {
        var result = new GraphFrame[frames];
        for (int t = 0; t < frames; t++)
        {
            int[] vertices = Enumerable.Range(0, n).Where(_ => rng.NextDouble() < 0.7).ToArray();
            var edges = new List<(int i, int j)>();
            for (int a = 0; a < vertices.Length; a++)
                for (int b = a + 1; b < vertices.Length; b++)
                    if (rng.NextDouble() < edgeProbability)
                        edges.Add((vertices[a], vertices[b]));
            result[t] = new GraphFrame(vertices, edges);
        }
        return result;
    }

    // Each frame is the last with a few vertices and pairs flipped, so a cycle can outlive a frame and be filled later.
    // Independent frames share too few edges for that.
    internal static GraphFrame[] EvolvingFrames(int n, int frames, Random rng, double vertexFlip, double pairFlip = 0.2)
    {
        var present = Enumerable.Repeat(true, n).ToArray();
        var linked = new bool[n, n];
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
                linked[a, b] = rng.NextDouble() < 0.5;

        var result = new GraphFrame[frames];
        for (int t = 0; t < frames; t++)
        {
            if (t > 0)
            {
                for (int a = 0; a < n; a++)
                {
                    if (rng.NextDouble() < vertexFlip) present[a] = !present[a];
                    for (int b = a + 1; b < n; b++)
                        if (rng.NextDouble() < pairFlip) linked[a, b] = !linked[a, b];
                }
            }

            var edges = new List<(int i, int j)>();
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    if (linked[a, b] && present[a] && present[b])
                        edges.Add((a, b));

            result[t] = new GraphFrame(Enumerable.Range(0, n).Where(v => present[v]).ToArray(), edges);
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

    static Barcode Read(GraphFrameFiltration built, int maxDimension) =>
        built.Complex == FrameComplex.Flag
            ? FastZigzag.Compute(built.Filtration, maxDimension)
            : GraphZigzag.Compute(built.Filtration, maxDimension);

    static (int, int, int, FrameEvent, FrameEvent)[] FramesWithEvents(GraphFrameFiltration built, Barcode barcode) =>
        built.FrameBars(barcode).Select(b => (b.Dimension, b.FirstFrame, b.LastFrame, b.Start, b.End)).OrderBy(x => x).ToArray();

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
