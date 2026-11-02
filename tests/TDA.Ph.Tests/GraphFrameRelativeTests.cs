#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Xunit;
using Xunit.Abstractions;

namespace TDA.Ph.Tests;

/// <summary>
/// The relative reading over frames, <c>H(K_t, K₀)</c>, built on a cone over <c>K₀</c>. It answers to the monotone reading
/// (<see cref="RelativeFiltration"/>), to a decomposition of the relative chain complex that uses no cone at all, and to the
/// absolute reading where a tree or an empty <c>K₀</c> forces them equal. Checks: <c>issues/ph/ph-loose-ends-brief.md</c>, step 5.
/// </summary>
public sealed class GraphFrameRelativeTests
{
    readonly ITestOutputHelper _output;

    public GraphFrameRelativeTests(ITestOutputHelper output) => _output = output;

    static readonly FrameTransition[] Transitions = { FrameTransition.Union, FrameTransition.Intersection };
    static readonly FrameComplex[] Complexes = { FrameComplex.Graph, FrameComplex.Flag };

    // ── by hand ─────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(FrameComplex.Graph)]
    [InlineData(FrameComplex.Flag)]
    public void K0sOwnCycleVanishes_AChordStays(FrameComplex complex)
    {
        // A 6-cycle backbone; the chord (0,3) in frames 1 and 2 of 4. Absolutely H₁ also holds the backbone's loop across every frame.
        var backbone = new List<(int i, int j)> { (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5) };
        var content = Sweep(new[] { Array.Empty<(int, int, double)>(), new[] { (0, 3, 1.0) }, new[] { (0, 3, 1.0) }, Array.Empty<(int, int, double)>() });
        double[] slack = { 1.0, 1.0, 1.0, 1.0 };

        foreach (FrameTransition transition in Transitions)
        {
            GraphFrameFiltration absolute = ConditionedFiltration.BuildSweep(6, backbone, content, slack, transition, complex);
            Assert.Contains((1, 0, 3), Frames(absolute, 1));

            GraphFrameFiltration relative = ConditionedFiltration.BuildRelativeSweep(6, backbone, content, slack, transition, complex);
            Assert.True(relative.IsRelative);
            Assert.Equal(new[] { (1, 1, 2) }, Frames(relative, 1));
        }
    }

    [Fact]
    public void AK0CycleThatFills_MovesUpADegree()
    {
        // R4 over frames. K₀ is the 4-cycle; the chord and its two triangles arrive at frame 1 and leave at frame 2.
        var backbone = new List<(int i, int j)> { (0, 1), (1, 2), (2, 3), (0, 3) };
        var content = Sweep(new[] { Array.Empty<(int, int, double)>(), new[] { (0, 2, 1.0) }, Array.Empty<(int, int, double)>() });

        foreach (FrameTransition transition in Transitions)
        {
            GraphFrameFiltration relative = ConditionedFiltration.BuildRelativeSweep(4, backbone, content, new[] { 1.0, 1.0, 1.0 }, transition, FrameComplex.Flag);
            Barcode barcode = FastZigzag.Compute(relative.Filtration, 2);
            FrameBar[] bars = relative.FrameBars(barcode).ToArray();

            FrameBar moved = Assert.Single(bars);
            Assert.Equal((2, 1, 1, FrameEvent.Arrival, FrameEvent.Departure), (moved.Dimension, moved.FirstFrame, moved.LastFrame, moved.Start, moved.End));
        }
    }

    [Fact]
    public void ANonSpanningK0_H0AttachesToIt()
    {
        // R5 over frames: backbone on {0,1,2}, (2,3) from frame 1, (4,5) from frame 2.
        var backbone = new List<(int i, int j)> { (0, 1), (1, 2) };
        var content = Sweep(new[]
        {
            Array.Empty<(int, int, double)>(),
            new[] { (2, 3, 1.0) },
            new[] { (2, 3, 1.0), (4, 5, 1.0) },
        });

        GraphFrameFiltration relative = ConditionedFiltration.BuildRelativeSweep(6, backbone, content, new[] { 1.0, 1.0, 1.0 });

        Assert.Equal(new[] { (0, 0, 0), (0, 0, 1), (0, 0, 2) }, Frames(relative, 0));
    }

    [Fact]
    public void FramesWithoutK0_AreRefused()
    {
        var k0 = new GraphFrame(new[] { 0, 1 }, new[] { (0, 1) });
        var missingEdge = new GraphFrame(new[] { 0, 1 }, Array.Empty<(int, int)>());
        var missingVertex = new GraphFrame(new[] { 1, 2 }, new[] { (1, 2) });

        Assert.Throws<ArgumentException>(() => GraphFrameZigzag.BuildRelative(3, new[] { missingEdge }, k0));
        Assert.Throws<ArgumentException>(() => GraphFrameZigzag.BuildRelative(3, new[] { missingVertex }, new GraphFrame(new[] { 0 }, Array.Empty<(int, int)>())));
        Assert.Throws<ArgumentException>(() => GraphFrameZigzag.BuildRelative(3, Array.Empty<GraphFrame>(), new GraphFrame(new[] { 0 }, new[] { (0, 1) })));
    }

    [Fact]
    public void TheGraphEnginesRefuseAConeOverAnEdge()
    {
        var edge = new GraphFrame(new[] { 0, 1 }, new[] { (0, 1) });
        GraphFrameFiltration relative = GraphFrameZigzag.BuildRelative(2, new[] { edge }, edge);

        Assert.ThrowsAny<Exception>(() => GraphZigzag.Compute(relative.Filtration, 1, algorithm: GraphZigzagAlgorithm.Reference));
        Assert.ThrowsAny<Exception>(() => GraphZigzag.Compute(relative.Filtration, 1, algorithm: GraphZigzagAlgorithm.Fast));
    }

    [Fact]
    public void AbsoluteFiltration_HasNoRelativeCycle()
    {
        GraphFrameFiltration absolute = GraphFrameZigzag.Build(2, new[] { new GraphFrame(new[] { 0, 1 }, Array.Empty<(int, int)>()) });
        Assert.False(absolute.IsRelative);
        Assert.Throws<InvalidOperationException>(() => absolute.RelativeCycle(new Bar(0, 1, 0, Cycle: new[] { 0 })));
    }

    [Fact]
    public void AChordsRelativeCycle_IsTheChordAlone()
    {
        // Path backbone 0-1-2, chord (0,2) at frame 1: the absolute cycle runs through the backbone, the relative one is the chord.
        var backbone = new List<(int i, int j)> { (0, 1), (1, 2) };
        var content = Sweep(new[] { Array.Empty<(int, int, double)>(), new[] { (0, 2, 1.0) } });

        GraphFrameFiltration relative = ConditionedFiltration.BuildRelativeSweep(3, backbone, content, new[] { 1.0, 1.0 });
        Barcode barcode = FastZigzag.Compute(relative.Filtration, 1, representatives: true);
        FrameBar bar = relative.FrameBars(barcode).Single(b => b.Dimension == 1);

        int cell = Assert.Single(relative.RelativeCycle(bar.StepBar));
        Assert.Equal(new[] { 0, 2 }, relative.VerticesOf(cell));
    }

    // ── the bridge to the monotone reading ──────────────────────────────────────

    [Fact]
    public void MonotoneFrames_EqualTheRelativeFiltration()
    {
        const int trials = 100;
        int bars = 0, finite = 0, movedUpDegree2 = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(4, 8);
            List<(int i, int j)> backbone = RandomEdges(n, rng, 0.3);
            if (seed % 2 == 0) backbone = SpanningPlus(n, rng);
            HashSet<(int, int)> taken = backbone.Select(e => (e.i, e.j)).ToHashSet();

            // A small grid of values, 0 among them, so ties occur and content can sit at residual 0.
            var edges = new List<(int i, int j, double d)>();
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    if (!taken.Contains((a, b)) && rng.NextDouble() < 0.4)
                        edges.Add((a, b, rng.Next(0, 6)));

            double[] levels = edges.Select(e => e.d).Append(0.0).Distinct().OrderBy(x => x).ToArray();
            var content = Enumerable.Repeat<IReadOnlyList<(int i, int j, double d)>>(edges, levels.Length).ToArray();

            foreach (FrameComplex complex in Complexes)
            {
                int top = complex == FrameComplex.Flag ? 2 : 1;
                RelativeFiltration reference = ConditionedFiltration.BuildRelative(n, backbone, edges, maxDimension: top);

                // A class born at value b and killed at d covers the frames at values in [b, d).
                var expected = new List<(int, int, int)>();
                foreach (Bar bar in reference.ComputeBarcode().Bars.Where(b => b.IsInfinite || b.Death > b.Birth))
                    expected.Add((bar.Dimension, Array.IndexOf(levels, bar.Birth), bar.IsInfinite ? levels.Length - 1 : Array.IndexOf(levels, bar.Death) - 1));

                foreach (FrameTransition transition in Transitions)
                {
                    GraphFrameFiltration relative = ConditionedFiltration.BuildRelativeSweep(n, backbone, content, levels, transition, complex);
                    FrameBar[] read = relative.FrameBars(FastZigzag.Compute(relative.Filtration, top)).ToArray();
                    Assert.True(
                        expected.OrderBy(x => x).SequenceEqual(read.Select(b => (b.Dimension, b.FirstFrame, b.LastFrame)).OrderBy(x => x)),
                        $"seed {seed}, {complex}, {transition}");
                    Assert.DoesNotContain(read, b => b.Start == FrameEvent.Departure || b.End == FrameEvent.Departure);
                }

                bars += expected.Count;
                finite += expected.Count(b => b.Item3 < levels.Length - 1);
                if (complex == FrameComplex.Flag) movedUpDegree2 += expected.Count(b => b.Item1 == 2);
            }
        }

        _output.WriteLine($"bars {bars}, ended before the last frame {finite}, H2 under Flag {movedUpDegree2}");
        Assert.True(bars > 0 && finite > 0 && movedUpDegree2 > 0);
    }

    // ── against the relative chain complex, with no cone ────────────────────────

    [Fact]
    public void Churn_EqualsTheDecompositionOfTheRelativeChainComplex()
    {
        const int trials = 300;
        int bars = 0, h1 = 0, h2 = 0, quotiented = 0, opened = 0, filled = 0, transients = 0, relativeCycles = 0, conedTriangles = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(5, 8);
            GraphFrameRelativeCase c = RandomCase(n, rng, pairFlip: 0.1);
            FrameComplex complex = Complexes[seed % 2];
            FrameTransition transition = Transitions[(seed / 2) % 2];

            GraphFrameFiltration relative = GraphFrameZigzag.BuildRelative(n, c.Frames, c.K0, transition, complex);
            Barcode barcode = FastZigzag.Compute(relative.Filtration, 2, representatives: true);
            HashSet<int> quotient = SubcomplexCells(relative, c.K0, n);

            // The relative chain complex's own intervals, read back by the same view: same frames, same events, same transients.
            var direct = new List<Bar>();
            for (int p = 0; p <= 2; p++)
                foreach (var (birth, death) in new ZigzagModuleOracle(relative.Filtration, p, quotient).Intervals())
                    direct.Add(new Bar(birth, death, p));

            var directBarcode = new Barcode(direct);
            Assert.True(
                FrameKeys(relative, barcode).SequenceEqual(FrameKeys(relative, directBarcode)),
                $"seed {seed}, {complex}, {transition}: frame bars");
            Assert.True(
                TransientKeys(relative, barcode).SequenceEqual(TransientKeys(relative, directBarcode)),
                $"seed {seed}, {complex}, {transition}: transients");

            // Each frame bar's relative cycle is a representative in that complex.
            var oracles = new[] { new ZigzagModuleOracle(relative.Filtration, 0, quotient), new ZigzagModuleOracle(relative.Filtration, 1, quotient), new ZigzagModuleOracle(relative.Filtration, 2, quotient) };
            foreach (FrameBar bar in relative.FrameBars(barcode))
            {
                BigInteger chain = ZigzagModuleOracle.ChainOf(relative.RelativeCycle(bar.StepBar));
                string? failure = oracles[bar.Dimension].Check(bar.StepBar, chain);
                Assert.True(failure is null, $"seed {seed}, degree {bar.Dimension}: {failure}");
                relativeCycles++;
            }

            FrameBar[] read = relative.FrameBars(barcode).ToArray();
            bars += read.Length;
            h1 += read.Count(b => b.Dimension == 1);
            h2 += read.Count(b => b.Dimension == 2);
            conedTriangles += relative.Filtration.Count(step => relative.VerticesOf(step.GlobalCellId).Count == 4 && step.Direction == ZigzagDirection.Add);
            opened += read.Count(b => b.Dimension == 1 && b.Start == FrameEvent.Departure);
            filled += read.Count(b => b.Dimension == 1 && b.End == FrameEvent.Arrival);
            transients += relative.TransientBars(barcode).Count;

            GraphFrameFiltration absolute = GraphFrameZigzag.Build(n, c.Frames, transition, complex);
            quotiented += Math.Max(0, absolute.FrameBars(FastZigzag.Compute(absolute.Filtration, 1)).Count(b => b.Dimension == 1) - read.Count(b => b.Dimension == 1));
        }

        _output.WriteLine($"frame bars {bars}, H1 {h1} ({filled} filled, {opened} opened), H2 {h2}, transients {transients}; " +
            $"H1 bars the relative reading has fewer than the absolute: {quotiented}; relative cycles checked {relativeCycles}; K0 triangles coned: {conedTriangles}");
        Assert.True(h1 > 0 && h2 > 0 && filled > 0 && opened > 0 && transients > 0 && quotiented > 0 && conedTriangles > 0);
    }

    // ── three engines ───────────────────────────────────────────────────────────

    [Fact]
    public void ThreeEnginesAgree_AtStepLevel()
    {
        const int trials = 10;
        int bars = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            GraphFrameRelativeCase c = RandomCase(4, rng, frames: 3);
            FrameComplex complex = Complexes[seed % 2];
            GraphFrameFiltration relative = GraphFrameZigzag.BuildRelative(c.N, c.Frames, c.K0, FrameTransition.Union, complex);

            var fast = Step(FastZigzag.Compute(relative.Filtration, 2));
            Assert.True(fast.SequenceEqual(Step(ZigzagBarcodeNaive.Compute(relative.Filtration, 2))), $"naive, seed {seed}");
            Assert.True(fast.SequenceEqual(Step(ZigzagMapBarcode.ComputeFromZigzag(relative.Filtration, 2))), $"map, seed {seed}");
            bars += fast.Length;
        }

        _output.WriteLine($"step-level bars: {bars}");
        Assert.True(bars > 0);
    }

    // ── the Diamond Principle, relative ─────────────────────────────────────────

    [Fact]
    public void BothTransitions_GiveTheSameFrameBarsAndEvents_AndTransientsOneDegreeApart()
    {
        const int trials = 200;
        int bars = 0, withEvents = 0, transients = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(4, 7);
            GraphFrameRelativeCase c = RandomCase(n, rng);
            FrameComplex complex = Complexes[seed % 2];

            GraphFrameFiltration union = GraphFrameZigzag.BuildRelative(n, c.Frames, c.K0, FrameTransition.Union, complex);
            GraphFrameFiltration intersection = GraphFrameZigzag.BuildRelative(n, c.Frames, c.K0, FrameTransition.Intersection, complex);
            int read = complex == FrameComplex.Flag ? 2 : 1;
            Barcode unionBars = FastZigzag.Compute(union.Filtration, read);
            Barcode intersectionBars = FastZigzag.Compute(intersection.Filtration, read);

            var keys = FrameKeys(union, unionBars);
            Assert.True(keys.SequenceEqual(FrameKeys(intersection, intersectionBars)), $"seed {seed}: frame bars");

            // A union transient of degree d is an intersection transient of degree d−1, between the same frames.
            Assert.DoesNotContain(union.TransientBars(unionBars), t => t.Dimension == 0);
            var fromUnion = union.TransientBars(unionBars).Where(t => t.Dimension >= 1).Select(t => (t.Dimension - 1, t.AfterFrame)).OrderBy(x => x).ToArray();
            var fromIntersection = intersection.TransientBars(intersectionBars).Where(t => t.Dimension < read).Select(t => (t.Dimension, t.AfterFrame)).OrderBy(x => x).ToArray();
            Assert.True(fromUnion.SequenceEqual(fromIntersection), $"seed {seed}: transients");

            bars += keys.Length;
            withEvents += keys.Count(k => k.Item4 != FrameEvent.Unobserved || k.Item5 != FrameEvent.Unobserved);
            transients += fromUnion.Length;
        }

        _output.WriteLine($"frame bars {bars}, with an observed event {withEvents}, matched transients {transients}");
        Assert.True(bars > 0 && withEvents > 0 && transients > 0);
    }

    // ── where K₀ forces the two readings equal ──────────────────────────────────

    [Fact]
    public void ASpanningTree_LeavesTheAbsoluteH1_AndNoH0()
    {
        const int trials = 100;
        int h1 = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(4, 8);
            List<(int i, int j)> tree = new();
            for (int v = 1; v < n; v++) tree.Add((rng.Next(v), v));

            IReadOnlyList<(int i, int j, double d)>[] content = RandomContent(n, rng, tree, rng.Next(3, 6));
            double[] slack = Enumerable.Repeat(1.0, content.Length).ToArray();
            FrameComplex complex = Complexes[seed % 2];
            FrameTransition transition = Transitions[(seed / 2) % 2];

            GraphFrameFiltration absolute = ConditionedFiltration.BuildSweep(n, tree, content, slack, transition, complex);
            GraphFrameFiltration relative = ConditionedFiltration.BuildRelativeSweep(n, tree, content, slack, transition, complex);
            Barcode absoluteBars = FastZigzag.Compute(absolute.Filtration, 1);
            Barcode relativeBars = FastZigzag.Compute(relative.Filtration, 1);

            var absoluteH1 = FrameKeys(absolute, absoluteBars).Where(k => k.Item1 == 1).ToArray();
            Assert.True(absoluteH1.SequenceEqual(FrameKeys(relative, relativeBars).Where(k => k.Item1 == 1)), $"seed {seed}");
            Assert.DoesNotContain(FrameKeys(relative, relativeBars), k => k.Item1 == 0);
            h1 += absoluteH1.Length;
        }

        _output.WriteLine($"H1 frame bars: {h1}");
        Assert.True(h1 > 0);
    }

    [Fact]
    public void AnEmptyK0_IsTheAbsoluteReading()
    {
        const int trials = 100;
        int bars = 0;
        var none = new GraphFrame(Array.Empty<int>(), Array.Empty<(int, int)>());
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(4, 7);
            GraphFrame[] frames = GraphFrameZigzagTests.EvolvingFrames(n, rng.Next(3, 6), rng, 0.15);
            FrameComplex complex = Complexes[seed % 2];
            FrameTransition transition = Transitions[(seed / 2) % 2];

            GraphFrameFiltration absolute = GraphFrameZigzag.Build(n, frames, transition, complex);
            GraphFrameFiltration relative = GraphFrameZigzag.BuildRelative(n, frames, none, transition, complex);

            var keys = FrameKeys(absolute, FastZigzag.Compute(absolute.Filtration, 1));
            Assert.True(keys.SequenceEqual(FrameKeys(relative, FastZigzag.Compute(relative.Filtration, 1))), $"seed {seed}");
            bars += keys.Length;
        }

        Assert.True(bars > 0);
    }

    // ── what the cone costs against a quotient ──────────────────────────────────

    [Fact]
    public void TheCone_AddsOneCellPerK0Cell_AndAnApex()
    {
        var backbone = new List<(int i, int j)> { (0, 1), (1, 2), (2, 3), (3, 4), (0, 4) };
        var content = Sweep(new[] { new[] { (0, 2, 1.0) }, new[] { (1, 3, 1.0) } });

        foreach (FrameComplex complex in Complexes)
        {
            GraphFrameFiltration absolute = ConditionedFiltration.BuildSweep(5, backbone, content, new[] { 1.0, 1.0 }, complex: complex);
            GraphFrameFiltration relative = ConditionedFiltration.BuildRelativeSweep(5, backbone, content, new[] { 1.0, 1.0 }, complex: complex);

            // K₀ is 5 vertices and 5 edges, with no triangle among them.
            const int k0Cells = 10;
            Assert.Equal(absolute.Filtration.Count + 2 * (k0Cells + 1), relative.Filtration.Count);
            _output.WriteLine($"{complex}: absolute {absolute.Filtration.Count} steps, coned {relative.Filtration.Count}, a quotient would be {absolute.Filtration.Count - 2 * (k0Cells - 1)}");
        }
    }

    // ── helpers ─────────────────────────────────────────────────────────────────

    sealed record GraphFrameRelativeCase(int N, GraphFrame[] Frames, GraphFrame K0);

    // Evolving frames with a K₀ forced into every one: a random subgraph of the first frame, its vertices kept throughout.
    static GraphFrameRelativeCase RandomCase(int n, Random rng, int? frames = null, double pairFlip = 0.2)
    {
        GraphFrame[] raw = GraphFrameZigzagTests.EvolvingFrames(n, frames ?? rng.Next(3, 6), rng, rng.Next(2) == 0 ? 0.0 : 0.1, pairFlip);
        var edges = new List<(int i, int j)>();
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
                if (rng.NextDouble() < 0.3)
                    edges.Add((a, b));

        var vertices = edges.SelectMany(e => new[] { e.i, e.j }).ToHashSet();
        var withK0 = raw
            .Select(f => new GraphFrame(
                f.Vertices.Union(vertices).OrderBy(v => v).ToArray(),
                f.Edges.Union(edges).ToArray()))
            .ToArray();

        return new GraphFrameRelativeCase(n, withK0, new GraphFrame(vertices.OrderBy(v => v).ToArray(), edges));
    }

    // The ids of K₀'s cells, its cone's and the apex, found from the cells' vertices: the apex is vertex n, and K₀ under Flag
    // is the flag complex of its edges.
    static HashSet<int> SubcomplexCells(GraphFrameFiltration built, GraphFrame k0, int n)
    {
        HashSet<(int, int)> edges = k0.Edges.Select(e => e.i < e.j ? (e.i, e.j) : (e.j, e.i)).ToHashSet();
        HashSet<int> vertices = k0.Vertices.ToHashSet();

        bool InK0(IReadOnlyList<int> v)
        {
            if (v.Count == 1) return vertices.Contains(v[0]);
            for (int a = 0; a < v.Count; a++)
                for (int b = a + 1; b < v.Count; b++)
                    if (!edges.Contains((v[a], v[b])))
                        return false;
            return true;
        }

        var cells = new HashSet<int>();
        foreach (ZigzagStep step in built.Filtration)
        {
            IReadOnlyList<int> v = built.VerticesOf(step.GlobalCellId);
            bool coned = v[v.Count - 1] == n;
            if (coned || InK0(v))
                cells.Add(step.GlobalCellId);
        }

        return cells;
    }

    static List<(int i, int j)> SpanningPlus(int n, Random rng)
    {
        var edges = new HashSet<(int i, int j)>();
        for (int v = 1; v < n; v++) edges.Add((rng.Next(v), v));
        for (int extra = rng.Next(0, n); extra > 0; extra--)
        {
            int a = rng.Next(n), b = rng.Next(n);
            if (a < b) edges.Add((a, b));
            else if (b < a) edges.Add((b, a));
        }
        return edges.ToList();
    }

    static List<(int i, int j)> RandomEdges(int n, Random rng, double p)
    {
        var edges = new List<(int i, int j)>();
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
                if (rng.NextDouble() < p)
                    edges.Add((a, b));
        return edges;
    }

    // Per frame, content edges outside the backbone, each pair flipping now and then so cycles last and return.
    static IReadOnlyList<(int i, int j, double d)>[] RandomContent(int n, Random rng, List<(int i, int j)> backbone, int frames)
    {
        HashSet<(int, int)> taken = backbone.Select(e => e.i < e.j ? (e.i, e.j) : (e.j, e.i)).ToHashSet();
        var on = new Dictionary<(int, int), bool>();
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
                if (!taken.Contains((a, b)))
                    on[(a, b)] = rng.NextDouble() < 0.4;

        var result = new IReadOnlyList<(int i, int j, double d)>[frames];
        for (int t = 0; t < frames; t++)
        {
            if (t > 0)
                foreach (var key in on.Keys.ToArray())
                    if (rng.NextDouble() < 0.25)
                        on[key] = !on[key];

            result[t] = on.Where(kv => kv.Value).Select(kv => (kv.Key.Item1, kv.Key.Item2, 1.0)).ToArray();
        }

        return result;
    }

    static IReadOnlyList<(int i, int j, double d)>[] Sweep(IEnumerable<(int, int, double)>[] frames) =>
        frames.Select(f => (IReadOnlyList<(int i, int j, double d)>)f.Select(e => (e.Item1, e.Item2, e.Item3)).ToArray()).ToArray();

    static (int, int, int)[] Frames(GraphFrameFiltration built, int maxDegree) =>
        built.FrameBars(FastZigzag.Compute(built.Filtration, maxDegree))
            .Select(b => (b.Dimension, b.FirstFrame, b.LastFrame))
            .Where(b => b.Dimension <= maxDegree)
            .OrderBy(x => x)
            .ToArray();

    static (int, int, int, FrameEvent, FrameEvent)[] FrameKeys(GraphFrameFiltration built, Barcode barcode) =>
        built.FrameBars(barcode)
            .Select(b => (b.Dimension, b.FirstFrame, b.LastFrame, b.Start, b.End))
            .OrderBy(x => x)
            .ToArray();

    static (int, int)[] TransientKeys(GraphFrameFiltration built, Barcode barcode) =>
        built.TransientBars(barcode)
            .Select(t => (t.Dimension, t.AfterFrame))
            .OrderBy(x => x)
            .ToArray();

    static (int, double, double, IntervalEnd, IntervalEnd)[] Step(Barcode barcode) =>
        barcode.Bars.Select(b => (b.Dimension, b.Birth, b.Death, b.BirthEnd, b.DeathEnd)).OrderBy(x => x).ToArray();
}
