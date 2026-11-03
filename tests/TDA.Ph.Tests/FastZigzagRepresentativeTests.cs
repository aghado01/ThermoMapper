#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Xunit;
using Xunit.Abstractions;

namespace TDA.Ph.Tests;

/// <summary>
/// What <see cref="Bar.Cycle"/> means on a <see cref="FastZigzag"/> bar, against <see cref="ZigzagModuleOracle"/>, exact
/// linear algebra outside the engine: a cycle of the complex at the bar's first level whose class is the bar's own, for
/// births by addition and by deletion, with re-entering cells. Checks: <c>issues/ph/ph-loose-ends-brief.md</c>, step 4.
/// </summary>
public sealed class FastZigzagRepresentativeTests
{
    readonly ITestOutputHelper _output;

    public FastZigzagRepresentativeTests(ITestOutputHelper output) => _output = output;

    // ── the oracle answers to the barcode first ─────────────────────────────────

    [Fact]
    public void Oracle_ReadsTheSameBarcodeAsTheNaiveEngine()
    {
        int bars = 0;
        for (int seed = 0; seed < 12; seed++)
        {
            ZigzagFiltration f = RandomChurn(new Random(seed), steps: 30);
            Barcode naive = ZigzagBarcodeNaive.Compute(f, 2);
            for (int p = 0; p <= 2; p++)
            {
                var expected = naive.Bars.Where(b => b.Dimension == p).Select(b => (b.Birth, b.Death)).OrderBy(x => x).ToArray();
                Assert.True(expected.SequenceEqual(new ZigzagModuleOracle(f, p).Intervals()), $"seed {seed}, degree {p}");
                bars += expected.Length;
            }
        }

        _output.WriteLine($"bars: {bars}");
        Assert.True(bars > 0);
    }

    // ── by hand ─────────────────────────────────────────────────────────────────

    [Fact]
    public void AddBorn_CycleClosesAtBirth()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 }); f.Add(5, new[] { 0, 2 });
        f.Delete(5);

        Bar bar = FastZigzag.Compute(f, 1, representatives: true).Bars.Single(b => b.Dimension == 1);

        Assert.Equal(IntervalEnd.Closed, bar.BirthEnd);
        Assert.Equal(new[] { 3, 4, 5 }, bar.Cycle!.OrderBy(x => x));
    }

    [Fact]
    public void DeleteBorn_H0_IsTheEndpointsOfTheSplittingEdge()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new[] { 0, 1 });
        f.Delete(2);

        Bar split = FastZigzag.Compute(f, 0, representatives: true).Bars.Single(b => b.BirthEnd == IntervalEnd.Open);

        Assert.Equal(3.0, split.Birth);
        Assert.Equal(new[] { 0, 1 }, split.Cycle!.OrderBy(x => x));
    }

    [Fact]
    public void DeleteBorn_H1_IsTheBoundaryOfTheLeavingTriangle()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 }); f.Add(5, new[] { 0, 2 });
        f.Add(6, new[] { 3, 4, 5 });
        f.Delete(6);

        Bar hole = FastZigzag.Compute(f, 1, representatives: true).Bars.Single(b => b.Dimension == 1 && b.BirthEnd == IntervalEnd.Open);

        Assert.Equal(7.0, hole.Birth);
        Assert.Equal(new[] { 3, 4, 5 }, hole.Cycle!.OrderBy(x => x));
    }

    [Fact]
    public void ReenteringEdge_NamesTheCycleItClosesOnReturn()
    {
        // The edge (0,2) closes a loop, leaves, and closes it again. Two bars, each with the cell id of that edge.
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 });
        f.Add(5, new[] { 0, 2 });
        f.Delete(5);
        f.Add(5, new[] { 0, 2 });

        Barcode barcode = FastZigzag.Compute(f, 1, representatives: true);
        Bar[] loops = barcode.Bars.Where(b => b.Dimension == 1).OrderBy(b => b.Birth).ToArray();

        Assert.Equal(2, loops.Length);
        Assert.All(loops, b => Assert.Equal(new[] { 3, 4, 5 }, b.Cycle!.OrderBy(x => x)));
        AssertEveryBarRepresented(f, barcode, 1, new Reach());
    }

    // ── random churn, with re-entry ─────────────────────────────────────────────

    [Fact]
    public void RandomChurn_EveryBarNamesItsClass()
    {
        const int trials = 120;
        var reach = new Reach();
        for (int seed = 0; seed < trials; seed++)
        {
            ZigzagFiltration f = RandomChurn(new Random(seed), steps: 40 + seed % 60);
            Barcode barcode = FastZigzag.Compute(f, 2, representatives: true);
            Assert.Equal(Signature(FastZigzag.Compute(f, 2)), Signature(barcode));
            AssertEveryBarRepresented(f, barcode, 2, reach, $"seed {seed}");
        }

        _output.WriteLine(reach.ToString());
        foreach (int p in new[] { 0, 1, 2 })
        {
            Assert.True(reach.AddBorn[p] > 0, $"no add-born bar in degree {p}");
            Assert.True(reach.DeleteBorn[p] > 0, $"no delete-born bar in degree {p}");
        }

        Assert.True(reach.ReenteredAddBorn[1] > 0, "no add-born H1 cycle through a re-entered cell");
        Assert.True(reach.ReenteredDeleteBorn[1] > 0, "no delete-born H1 cycle through a re-entered cell");
        Assert.True(reach.ReenteredAddBorn[2] > 0, "no add-born H2 cycle through a re-entered cell");
    }

    // ── the strong condition is stronger than the birth step's ──────────────────

    [Fact]
    public void AddingAnOlderLongerClass_PassesAtBirth_AndFailsTheStrongCondition()
    {
        const int trials = 60;
        int passedAtBirth = 0, rejected = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            ZigzagFiltration f = RandomChurn(new Random(seed), steps: 40 + seed % 60);
            Barcode barcode = FastZigzag.Compute(f, 1, representatives: true);
            var oracle = new ZigzagModuleOracle(f, 1);
            foreach (Bar bar in barcode.Bars.Where(b => b.Dimension == 1 && b.BirthEnd == IntervalEnd.Closed))
            {
                BigInteger cycle = ZigzagModuleOracle.ChainOf(bar.Cycle!);
                foreach (BigInteger older in oracle.Z((int)bar.Birth).Basis)
                {
                    BigInteger other = cycle ^ older;
                    if (oracle.CheckAtBirth(bar, other) is not null) continue;
                    passedAtBirth++;
                    if (oracle.Check(bar, other) is not null) rejected++;
                }
            }
        }

        _output.WriteLine($"perturbed cycles passing every birth condition: {passedAtBirth}; of those, rejected: {rejected}");
        Assert.True(rejected > 0);
    }

    // ── through the frames ──────────────────────────────────────────────────────

    static readonly (int, int)[] C4 = { (0, 1), (0, 3), (1, 2), (2, 3) };

    [Theory]
    [InlineData(FrameTransition.Union)]
    [InlineData(FrameTransition.Intersection)]
    public void FilledBar_NamesTheCycleThatWasFilled_InFrameTerms(FrameTransition transition)
    {
        var chord = C4.Append((0, 2)).ToArray();
        GraphFrameFiltration built = GraphFrameZigzag.FromEdgeSets(4, new[] { (IReadOnlyList<(int i, int j)>)C4.Select(ToEdge).ToArray(), chord.Select(ToEdge).ToArray() }, transition, FrameComplex.Flag);

        FrameBar filled = built.FrameBars(FastZigzag.Compute(built.Filtration, 1, representatives: true)).Single(b => b.Dimension == 1);

        Assert.Equal(FrameEvent.Arrival, filled.End);
        Assert.Equal(C4, EdgesOf(built, filled.StepBar).OrderBy(x => x));
    }

    [Theory]
    [InlineData(FrameTransition.Union)]
    [InlineData(FrameTransition.Intersection)]
    public void OpenedHole_NamesATriangleBoundary_ThatTheFrameHoldsOnlyUpToHomology(FrameTransition transition)
    {
        // The class is the hole C4 at frame 1. Where it begins a triangle has just left, so its chain is that triangle's
        // boundary, which holds the chord, and the chord is gone by the frame.
        var chord = C4.Append((0, 2)).ToArray();
        GraphFrameFiltration built = GraphFrameZigzag.FromEdgeSets(4, new[] { (IReadOnlyList<(int i, int j)>)chord.Select(ToEdge).ToArray(), C4.Select(ToEdge).ToArray() }, transition, FrameComplex.Flag);

        FrameBar hole = built.FrameBars(FastZigzag.Compute(built.Filtration, 1, representatives: true)).Single(b => b.Dimension == 1);

        Assert.Equal(FrameEvent.Departure, hole.Start);
        (int, int)[] named = EdgesOf(built, hole.StepBar);
        Assert.Equal(3, named.Length);
        Assert.Contains((0, 2), named);

        int level = (int)hole.StepBar.Birth + 1;
        BigInteger c4 = ZigzagModuleOracle.ChainOf(hole.StepBar.Cycle!) ^ ChainOfEdges(built, level, C4);
        Assert.True(new ZigzagModuleOracle(built.Filtration, 1).B(level).Contains(c4), "the named chain is not homologous to C4 where the bar begins");
    }

    [Fact]
    public void OneReturnDyingBothWays_NamesEachCycle()
    {
        // F5 of the P1d brief: path 0-3, chord (0,3) in frames 0-2, chord (0,2) in frame 1.
        var path = new List<(int i, int j)> { (0, 1), (1, 2), (2, 3) };
        var content = new IReadOnlyList<(int i, int j, double d)>[]
        {
            new (int, int, double)[] { (0, 3, 1.0) },
            new (int, int, double)[] { (0, 3, 1.0), (0, 2, 1.0) },
            new (int, int, double)[] { (0, 3, 1.0) },
            Array.Empty<(int, int, double)>(),
        };

        GraphFrameFiltration built = ConditionedFiltration.BuildSweep(4, path, content, new[] { 1.0, 1.0, 1.0, 1.0 }, FrameTransition.Union, FrameComplex.Flag);
        Barcode barcode = FastZigzag.Compute(built.Filtration, 1, representatives: true);
        FrameBar[] bars = built.FrameBars(barcode).Where(b => b.Dimension == 1).OrderBy(b => b.FirstFrame).ToArray();

        Assert.Equal(2, bars.Length);
        Assert.Equal(new[] { (0, 1), (0, 3), (1, 2), (2, 3) }, EdgesOf(built, bars[0].StepBar).OrderBy(x => x));
        Assert.Contains((0, 2), EdgesOf(built, bars[1].StepBar));
        AssertEveryBarRepresented(built.Filtration, barcode, 1, new Reach());
    }

    [Fact]
    public void RandomChurn_EveryFrameBarNamesItsClass_UnderBothComplexes()
    {
        const int trials = 300;
        int filled = 0, opened = 0, broken = 0, filledOutsideFirstFrame = 0, openedOutsideFirstFrame = 0, otherOutsideFirstFrame = 0, graphBars = 0;
        var reach = new Reach();

        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(5, 8);
            GraphFrame[] frames = GraphFrameZigzagTests.EvolvingFrames(n, rng.Next(3, 6), rng, seed % 2 == 0 ? 0.0 : 0.1, pairFlip: 0.1);
            FrameTransition transition = seed % 4 < 2 ? FrameTransition.Union : FrameTransition.Intersection;
            HashSet<(int, int)>[] edges = frames.Select(fr => fr.Edges.Select(e => e.i < e.j ? (e.i, e.j) : (e.j, e.i)).ToHashSet()).ToArray();

            GraphFrameFiltration flag = GraphFrameZigzag.Build(n, frames, transition, FrameComplex.Flag);
            Barcode flagBars = FastZigzag.Compute(flag.Filtration, 1, representatives: true);
            AssertEveryBarRepresented(flag.Filtration, flagBars, 1, reach, $"seed {seed}");

            foreach (FrameBar bar in flag.FrameBars(flagBars).Where(b => b.Dimension == 1))
            {
                bool inFirstFrame = EdgesOf(flag, bar.StepBar).All(edges[bar.FirstFrame].Contains);
                if (bar.End == FrameEvent.Arrival) filled++;
                if (bar.End == FrameEvent.Departure) broken++;
                if (bar.Start == FrameEvent.Departure) opened++;

                if (inFirstFrame) continue;
                if (bar.Start == FrameEvent.Departure) openedOutsideFirstFrame++;
                else if (bar.End == FrameEvent.Arrival) filledOutsideFirstFrame++;
                else otherOutsideFirstFrame++;
            }

            // A graph has no boundaries in degree 1, so a class is its one chain, and that chain has to be in every frame.
            GraphFrameFiltration graph = GraphFrameZigzag.Build(n, frames, transition, FrameComplex.Graph);
            Barcode graphBarcode = FastZigzag.Compute(graph.Filtration, 1, representatives: true);
            foreach (FrameBar bar in graph.FrameBars(graphBarcode).Where(b => b.Dimension == 1))
            {
                graphBars++;
                for (int t = bar.FirstFrame; t <= bar.LastFrame; t++)
                    Assert.True(EdgesOf(graph, bar.StepBar).All(edges[t].Contains), $"seed {seed}: a graph cycle absent from frame {t}");
            }
        }

        _output.WriteLine($"H1 frame bars under Flag: {filled} filled, {broken} broken, {opened} opened; " +
            $"cycle not all in the first frame: {openedOutsideFirstFrame} opened, {filledOutsideFirstFrame} filled, {otherOutsideFirstFrame} other; " +
            $"under Graph: {graphBars} bars. {reach}");
        Assert.True(filled > 0 && opened > 0 && broken > 0 && graphBars > 0);
        // A triangle leaves a flag complex only when an edge of it does, so an opened hole's named chain never fits its first frame.
        Assert.Equal(opened, openedOutsideFirstFrame);
    }

    static (int, int) ToPair(IReadOnlyList<int> v) => v[0] < v[1] ? (v[0], v[1]) : (v[1], v[0]);

    static (int i, int j) ToEdge((int, int) e) => (e.Item1, e.Item2);

    static (int, int)[] EdgesOf(GraphFrameFiltration built, Bar bar) => bar.Cycle!.Select(c => ToPair(built.VerticesOf(c))).ToArray();

    // The chain of the given edges among the cells of the filtration, by their vertices.
    static BigInteger ChainOfEdges(GraphFrameFiltration built, int level, (int, int)[] edges)
    {
        var cells = new List<int>();
        for (int s = 0; s < level; s++)
        {
            int cell = built.Filtration[s].GlobalCellId;
            if (built.Filtration[s].Direction == ZigzagDirection.Add && built.VerticesOf(cell).Count == 2 && edges.Contains(ToPair(built.VerticesOf(cell))) && !cells.Contains(cell))
                cells.Add(cell);
        }

        return ZigzagModuleOracle.ChainOf(cells);
    }

    // ── helpers ─────────────────────────────────────────────────────────────────

    sealed class Reach
    {
        public readonly int[] AddBorn = new int[3], DeleteBorn = new int[3], ReenteredAddBorn = new int[3], ReenteredDeleteBorn = new int[3];

        public override string ToString() =>
            $"add-born H0/H1/H2 {string.Join("/", AddBorn)}, delete-born {string.Join("/", DeleteBorn)}; " +
            $"through a re-entered cell, add-born {string.Join("/", ReenteredAddBorn)}, delete-born {string.Join("/", ReenteredDeleteBorn)}";
    }

    static void AssertEveryBarRepresented(ZigzagFiltration f, Barcode barcode, int maxDegree, Reach reach, string context = "")
    {
        for (int p = 0; p <= maxDegree; p++)
        {
            var oracle = new ZigzagModuleOracle(f, p);
            foreach (Bar bar in barcode.Bars.Where(b => b.Dimension == p))
            {
                Assert.NotNull(bar.Cycle);
                Assert.Equal(bar.Cycle!.Length, bar.Cycle.Distinct().Count());
                string? failure = oracle.Check(bar, ZigzagModuleOracle.ChainOf(bar.Cycle));
                Assert.True(failure is null, $"{context}, degree {p}, bar [{bar.Birth}, {bar.Death}] {bar.BirthEnd}: {failure}");

                bool reentered = bar.Cycle.Any(cell => DeletedBefore(f, cell, (int)bar.Birth));
                if (bar.BirthEnd == IntervalEnd.Closed)
                {
                    reach.AddBorn[p]++;
                    if (reentered) reach.ReenteredAddBorn[p]++;
                }
                else
                {
                    reach.DeleteBorn[p]++;
                    if (reentered) reach.ReenteredDeleteBorn[p]++;
                }
            }
        }
    }

    // A cell in a bar's cycle that was deleted before the bar began is present again by a later incarnation: the case the
    // old docs left open, where the coned cycle might have held the earlier one.
    static bool DeletedBefore(ZigzagFiltration f, int cell, int step)
    {
        for (int s = 0; s < step; s++)
            if (f[s].GlobalCellId == cell && f[s].Direction == ZigzagDirection.Delete)
                return true;
        return false;
    }

    static (double, double, int, IntervalEnd, IntervalEnd)[] Signature(Barcode barcode) =>
        barcode.Bars.Select(b => (b.Birth, b.Death, b.Dimension, b.BirthEnd, b.DeathEnd)).OrderBy(x => x).ToArray();

    // Simplices through tetrahedra over a few vertices, each with one cell id across re-entries. Each step adds a cell whose
    // faces are present or deletes one with no coface present, so cells leave and return. Half the trials tear down to empty.
    internal static ZigzagFiltration RandomChurn(Random rng, int steps)
    {
        int n = rng.Next(4, 7);
        var cells = new List<int[]>();
        for (int v = 0; v < n; v++) cells.Add(new[] { v });
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) cells.Add(new[] { a, b });
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) for (int c = b + 1; c < n; c++) cells.Add(new[] { a, b, c });
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) for (int c = b + 1; c < n; c++) for (int d = c + 1; d < n; d++) cells.Add(new[] { a, b, c, d });

        var idOf = new Dictionary<string, int>();
        for (int id = 0; id < cells.Count; id++) idOf[string.Join(",", cells[id])] = id;

        int[] FacesOf(int id)
        {
            int[] v = cells[id];
            return v.Length == 1 ? Array.Empty<int>()
                : Enumerable.Range(0, v.Length).Select(k => idOf[string.Join(",", v.Where((_, i) => i != k))]).ToArray();
        }

        int[][] faces = Enumerable.Range(0, cells.Count).Select(FacesOf).ToArray();
        var cofaces = new List<int>[cells.Count];
        for (int id = 0; id < cells.Count; id++) cofaces[id] = new List<int>();
        for (int id = 0; id < cells.Count; id++) foreach (int face in faces[id]) cofaces[face].Add(id);

        var present = new bool[cells.Count];
        var f = new ZigzagFiltration();
        double addChance = 0.55 + rng.NextDouble() * 0.2;

        for (int s = 0; s < steps; s++)
        {
            var addable = Enumerable.Range(0, cells.Count).Where(id => !present[id] && faces[id].All(x => present[x])).ToList();
            var deletable = Enumerable.Range(0, cells.Count).Where(id => present[id] && cofaces[id].All(x => !present[x])).ToList();

            if (addable.Count > 0 && (deletable.Count == 0 || rng.NextDouble() < addChance))
            {
                int id = addable[rng.Next(addable.Count)];
                f.Add(id, faces[id]);
                present[id] = true;
            }
            else if (deletable.Count > 0)
            {
                int id = deletable[rng.Next(deletable.Count)];
                f.Delete(id);
                present[id] = false;
            }
        }

        if (rng.Next(2) == 0)
            for (int id = cells.Count - 1; id >= 0; id--)
                if (present[id])
                    f.Delete(id);

        return f;
    }
}
