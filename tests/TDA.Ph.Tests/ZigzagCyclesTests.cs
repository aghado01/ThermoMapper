#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Numerics;
using Xunit;
using Xunit.Abstractions;

namespace TDA.Ph.Tests;

/// <summary>
/// <see cref="ZigzagCycles"/> against <see cref="ZigzagModuleOracle"/>, exact linear algebra outside the engine: each bar's
/// cycles are a set of representatives, Dey and Hou's Definition 5, and at each level the alive bars' classes are a basis.
/// </summary>
public sealed class ZigzagCyclesTests
{
    readonly ITestOutputHelper _output;

    public ZigzagCyclesTests(ITestOutputHelper output) => _output = output;

    // ── by hand ─────────────────────────────────────────────────────────────────

    [Fact]
    public void NoBoundaries_TheCycleStaysPut()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]);
        f.Add(3, new[] { 0, 1 }); f.Add(4, new[] { 1, 2 }); f.Add(5, new[] { 0, 2 });
        f.Add(6, new int[0]);
        f.Delete(5);

        ZigzagCycles cycles = FastZigzag.ComputeCycles(f, 1);
        Bar loop = cycles.Barcode.Bars.Single(b => b.Dimension == 1);

        Assert.Equal(new[] { 3, 4, 5 }, cycles.CycleAt(loop, 6));
        Assert.Equal(new[] { 3, 4, 5 }, cycles.CycleAt(loop, 7));
    }

    [Fact]
    public void CycleMovesAcrossATriangle_BeforeAnEdgeOfItLeaves()
    {
        // The loop 01-12-23-03 closes, then chord 02 and triangle 012 arrive. The triangle leaves and edge 01 goes after it, so
        // the bar must already have moved off 01 and 12 by then, onto 02-23-03. It can only move while the triangle is there.
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]); f.Add(3, new int[0]);
        f.Add(4, new[] { 0, 1 }); f.Add(5, new[] { 1, 2 }); f.Add(6, new[] { 2, 3 });
        f.Add(7, new[] { 0, 3 });                                // step 7: closes the loop
        f.Add(8, new[] { 0, 2 });                                // step 8: the chord
        f.Add(9, new[] { 4, 5, 8 });                             // step 9: fills 01-12-02
        f.Delete(9);                                             // step 10: opens it again
        f.Delete(4);                                             // step 11: breaks it

        ZigzagCycles cycles = FastZigzag.ComputeCycles(f, 1);
        Bar loop = cycles.Barcode.Bars.Single(b => b.Dimension == 1 && b.Birth == 7.0);

        Assert.Equal(new[] { 4, 5, 6, 7 }, cycles.CycleAt(loop, 8));
        Assert.Equal(new[] { 4, 5, 6, 7 }, cycles.CycleAt(loop, 9));
        Assert.Equal(new[] { 6, 7, 8 }, cycles.CycleAt(loop, 10));
        Assert.Equal(new[] { 6, 7, 8 }, cycles.CycleAt(loop, 12));
        Assert.Equal(12.0, loop.Death);
        AssertEveryBarRepresented(f, cycles, 1, new Reach());
    }

    [Fact]
    public void OpenedHole_GrowsAsItsTrianglesLeave()
    {
        // Two triangles fill the square. The first to leave opens its own boundary; the second turns that into the square.
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]); f.Add(1, new int[0]); f.Add(2, new int[0]); f.Add(3, new int[0]);
        f.Add(4, new[] { 0, 1 }); f.Add(5, new[] { 1, 2 }); f.Add(6, new[] { 2, 3 }); f.Add(7, new[] { 0, 3 });
        f.Add(8, new[] { 0, 2 });
        f.Add(9, new[] { 4, 5, 8 }); f.Add(10, new[] { 6, 7, 8 });
        f.Delete(9);                                             // step 11
        f.Delete(10);                                            // step 12
        f.Delete(8);                                             // step 13

        ZigzagCycles cycles = FastZigzag.ComputeCycles(f, 1);
        Bar hole = cycles.Barcode.Bars.Single(b => b.Dimension == 1 && b.Birth == 11.0);

        Assert.Equal(IntervalEnd.Open, hole.BirthEnd);
        Assert.Equal(new[] { 4, 5, 8 }, cycles.CycleAt(hole, 12));
        Assert.Equal(new[] { 4, 5, 6, 7 }, cycles.CycleAt(hole, 13));
        Assert.Equal(new[] { 4, 5, 6, 7 }, cycles.CycleAt(hole, 14));
        AssertEveryBarRepresented(f, cycles, 1, new Reach());
    }

    [Fact]
    public void ComponentMovesOntoTheVertexThatOutlivesIt()
    {
        // Vertex 0's component takes in vertex 1, then 0 has to leave. The H₀ bar born at 0 must have moved onto 1 by the time
        // the edge is gone, or vertex 0's departure would end it.
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]);
        f.Add(1, new int[0]);
        f.Add(2, new[] { 0, 1 });                                // kills the younger component
        f.Delete(2);                                             // step 3: splits again
        f.Delete(0);                                             // step 4

        ZigzagCycles cycles = FastZigzag.ComputeCycles(f, 0);
        Bar elder = cycles.Barcode.Bars.Single(b => b.Birth == 0.0);

        Assert.Equal(5.0, elder.Death);
        Assert.Equal(new[] { 0 }, cycles.CycleAt(elder, 2));
        Assert.Equal(new[] { 1 }, cycles.CycleAt(elder, 3));
        Assert.Equal(new[] { 1 }, cycles.CycleAt(elder, 5));
        AssertEveryBarRepresented(f, cycles, 0, new Reach());
    }

    [Fact]
    public void Levels_OutsideTheBar_AreRefused_AndSoIsAnotherBarcodesBar()
    {
        var f = new ZigzagFiltration();
        f.Add(0, new int[0]);
        f.Delete(0);

        ZigzagCycles cycles = FastZigzag.ComputeCycles(f);
        Bar bar = cycles.Barcode.Bars.Single();

        Assert.Equal(new[] { 0 }, cycles.CycleAt(bar, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => cycles.CycleAt(bar, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => cycles.CycleAt(bar, 2));
        Assert.Throws<ArgumentException>(() => cycles.CycleAt(new Bar(5, 6, 0), 6));
        Assert.Throws<ArgumentException>(() => cycles.CycleAt(bar with { Dimension = 1 }, 1));
    }

    [Fact]
    public void EmptyFiltration_HasNoBars()
    {
        Assert.Empty(FastZigzag.ComputeCycles(new ZigzagFiltration()).Barcode.Bars);
    }

    // ── random churn, with re-entry ─────────────────────────────────────────────

    [Fact]
    public void RandomChurn_EveryBarsCyclesAreItsRepresentatives_AndABasisAtEachLevel()
    {
        const int trials = 120;
        var reach = new Reach();
        for (int seed = 0; seed < trials; seed++)
        {
            ZigzagFiltration f = FastZigzagRepresentativeTests.RandomChurn(new Random(seed), steps: 40 + seed % 60);
            ZigzagCycles cycles = FastZigzag.ComputeCycles(f, 2);
            Barcode plain = FastZigzag.Compute(f, 2, representatives: true);

            // The same bars as the other path, with the same cycle where each begins.
            Assert.Equal(Signature(plain), Signature(cycles.Barcode));
            foreach (Bar bar in cycles.Barcode.Bars)
            {
                Bar other = plain.Bars.Single(b => b.Birth == bar.Birth);
                Assert.Equal(other.Cycle!.OrderBy(x => x), cycles.CycleAt(bar, (int)bar.Birth + 1));
            }

            AssertEveryBarRepresented(f, cycles, 2, reach, $"seed {seed}");
        }

        _output.WriteLine(reach.ToString());
        // Every kind of bar, in every degree, and bars that actually move.
        for (int p = 0; p <= 2; p++)
            for (int kind = 0; kind < 4; kind++)
                Assert.True(reach.ByKind[p, kind] > 0, $"no bar of kind {kind} in degree {p}");
        Assert.True(reach.Moved[0] > 0 && reach.Moved[1] > 0 && reach.Moved[2] > 0);
        Assert.True(reach.ThroughReentry > 0);
    }

    // ── the check can fail ──────────────────────────────────────────────────────

    [Fact]
    public void WrongSequences_AreRefused()
    {
        const int trials = 60;
        int barsThatMove = 0, heldStill = 0, oneOff = 0, swapped = 0, swappedRefused = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            ZigzagFiltration f = FastZigzagRepresentativeTests.RandomChurn(new Random(seed), steps: 40 + seed % 60);
            ZigzagCycles cycles = FastZigzag.ComputeCycles(f, 1);
            for (int p = 0; p <= 1; p++)
            {
                var oracle = new ZigzagModuleOracle(f, p);
                Bar[] bars = cycles.Barcode.Bars.Where(b => b.Dimension == p).ToArray();
                foreach (Bar bar in bars)
                {
                    int first = (int)bar.Birth + 1, last = (int)bar.Death;
                    BigInteger At(int level) => ZigzagModuleOracle.ChainOf(cycles.CycleAt(bar, level));

                    // The cycle where the bar begins, held for its whole length, is what Bar.Cycle alone would give.
                    if (Enumerable.Range(first, last - first + 1).Any(level => At(level) != At(first)))
                    {
                        barsThatMove++;
                        if (oracle.CheckSequence(bar, _ => At(first)) is not null) heldStill++;
                        if (oracle.CheckSequence(bar, level => At(Math.Min(level + 1, last))) is not null) oneOff++;
                    }

                    // One interior level swapped for its sum with another bar alive there and at the next level. A bar alive
                    // at that level alone would leave a valid sequence: representatives aren't unique.
                    if (last - first < 2) continue;
                    int middle = (first + last) / 2;
                    Bar? other = bars.Cast<Bar?>().FirstOrDefault(b => b!.Value.Birth != bar.Birth && b.Value.Birth + 1 <= middle && middle + 1 <= b.Value.Death);
                    if (other is null) continue;

                    Bar neighbor = other.Value;
                    swapped++;
                    BigInteger Swapped(int level) => level == middle ? At(level) ^ ZigzagModuleOracle.ChainOf(cycles.CycleAt(neighbor, level)) : At(level);
                    if (oracle.CheckSequence(bar, Swapped) is not null) swappedRefused++;
                }
            }
        }

        _output.WriteLine($"bars that move: {barsThatMove}, refused when held still: {heldStill}, when one level off: {oneOff}; swapped {swapped}, refused {swappedRefused}");
        Assert.True(barsThatMove > 0 && swapped > 0);
        // A bar whose chain moves might still hold one cycle up to homology, so not all of these fail. Most must.
        Assert.True(heldStill > barsThatMove / 2, $"{heldStill} of {barsThatMove}");
        Assert.True(oneOff > barsThatMove / 2, $"{oneOff} of {barsThatMove}");
        Assert.Equal(swapped, swappedRefused);
    }

    // ── through the frames ──────────────────────────────────────────────────────

    static readonly (int, int)[] C4 = { (0, 1), (0, 3), (1, 2), (2, 3) };

    [Theory]
    [InlineData(FrameTransition.Union)]
    [InlineData(FrameTransition.Intersection)]
    public void OpenedHole_AtItsFrame_IsTheHoleItself(FrameTransition transition)
    {
        var chord = C4.Append((0, 2)).ToArray();
        GraphFrameFiltration built = GraphFrameZigzag.FromEdgeSets(4, new[] { (IReadOnlyList<(int i, int j)>)chord.Select(ToEdge).ToArray(), C4.Select(ToEdge).ToArray() }, transition, FrameComplex.Flag);

        ZigzagCycles cycles = FastZigzag.ComputeCycles(built.Filtration, 1);
        FrameBar hole = built.FrameBars(cycles.Barcode).Single(b => b.Dimension == 1);

        Assert.Equal(FrameEvent.Departure, hole.Start);
        Assert.Equal(C4, built.CycleAtFrame(cycles, hole.StepBar, 1).Select(c => ToPair(built.VerticesOf(c))).OrderBy(x => x));
    }

    [Fact]
    public void RandomFrames_EachBarsCycleAtAFrame_LiesInThatFrame_UnderBothComplexes()
    {
        const int trials = 200;
        int opened = 0, openedWhoseBirthChainMisses = 0, frameBars = 0, transients = 0;
        var reach = new Reach();

        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(5, 8);
            GraphFrame[] frames = GraphFrameZigzagTests.EvolvingFrames(n, rng.Next(3, 6), rng, seed % 2 == 0 ? 0.0 : 0.1, pairFlip: 0.1);
            FrameTransition transition = seed % 4 < 2 ? FrameTransition.Union : FrameTransition.Intersection;
            HashSet<(int, int)>[] edges = frames.Select(fr => fr.Edges.Select(e => e.i < e.j ? (e.i, e.j) : (e.j, e.i)).ToHashSet()).ToArray();

            foreach (FrameComplex complex in new[] { FrameComplex.Flag, FrameComplex.Graph })
            {
                GraphFrameFiltration built = GraphFrameZigzag.Build(n, frames, transition, complex);
                ZigzagCycles cycles = FastZigzag.ComputeCycles(built.Filtration, 1);
                AssertEveryBarRepresented(built.Filtration, cycles, 1, reach, $"seed {seed}, {complex}");

                foreach (FrameBar bar in built.FrameBars(cycles.Barcode))
                {
                    frameBars++;
                    for (int t = bar.FirstFrame; t <= bar.LastFrame; t++)
                    {
                        int[] cycle = built.CycleAtFrame(cycles, bar.StepBar, t);
                        foreach (int cell in cycle)
                        {
                            IReadOnlyList<int> vertices = built.VerticesOf(cell);
                            if (vertices.Count == 1) Assert.Contains(vertices[0], frames[t].Vertices);
                            else Assert.True(edges[t].Contains(ToPair(vertices)), $"seed {seed}: an edge outside frame {t}");
                        }
                    }

                    if (complex == FrameComplex.Flag && bar.Dimension == 1 && bar.Start == FrameEvent.Departure)
                    {
                        opened++;
                        if (bar.StepBar.Cycle!.Any(c => !edges[bar.FirstFrame].Contains(ToPair(built.VerticesOf(c)))))
                            openedWhoseBirthChainMisses++;
                    }
                }

                foreach (TransientBar transient in built.TransientBars(cycles.Barcode))
                {
                    transients++;
                    Assert.NotEmpty(cycles.CycleAt(transient.StepBar, built.LevelBetweenFrames(transient.AfterFrame)));
                }
            }
        }

        _output.WriteLine($"{frameBars} frame bars, {transients} transients; opened holes under Flag: {opened}, " +
            $"of those whose birth chain misses their first frame: {openedWhoseBirthChainMisses}. {reach}");
        Assert.True(opened > 0 && transients > 0);
        Assert.Equal(opened, openedWhoseBirthChainMisses);
    }

    [Fact]
    public void RelativeFrames_EachBarsRelativeCycle_RepresentsItInTheQuotient()
    {
        const int trials = 120;
        int bars = 0, levels = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(5, 8);
            GraphFrameRelativeTests.GraphFrameRelativeCase c = GraphFrameRelativeTests.RandomCase(n, rng, pairFlip: 0.1);
            FrameComplex complex = seed % 2 == 0 ? FrameComplex.Flag : FrameComplex.Graph;
            FrameTransition transition = seed % 4 < 2 ? FrameTransition.Union : FrameTransition.Intersection;
            GraphFrameFiltration relative = GraphFrameZigzag.BuildRelative(n, c.Frames, c.K0, transition, complex);
            HashSet<int> quotient = GraphFrameRelativeTests.SubcomplexCells(relative, c.K0, n);
            ZigzagCycles cycles = FastZigzag.ComputeCycles(relative.Filtration, 1);

            for (int p = 0; p <= 1; p++)
            {
                var oracle = new ZigzagModuleOracle(relative.Filtration, p, quotient);
                foreach (FrameBar bar in relative.FrameBars(cycles.Barcode).Where(b => b.Dimension == p))
                {
                    bars++;
                    for (int t = bar.FirstFrame; t <= bar.LastFrame; t++)
                    {
                        levels++;
                        int level = relative.LevelOfFrame(t);
                        BigInteger chain = ZigzagModuleOracle.ChainOf(relative.RelativeCycleAtFrame(cycles, bar.StepBar, t));
                        Assert.True(oracle.BoundaryOf(level, chain).IsZero, $"seed {seed}: no relative cycle at frame {t}");
                        Assert.False(oracle.B(level).Contains(chain), $"seed {seed}: a relative boundary at frame {t}");
                        Assert.Equal((int)bar.StepBar.Death, oracle.EndOf(level, chain));
                    }

                    // Through the levels between, the cone's cells gone, the chain is the same sequence read in the quotient.
                    Assert.Null(oracle.CheckSequence(bar.StepBar, level =>
                        ZigzagModuleOracle.ChainOf(cycles.CycleAt(bar.StepBar, level).Where(cell => !quotient.Contains(cell)))));
                }
            }
        }

        _output.WriteLine($"relative bars {bars}, frame cycles {levels}");
        Assert.True(bars > 0);
    }

    // ── what keeping V costs ────────────────────────────────────────────────────

    // Off by default: set ZIGZAG_CYCLES_COST to a file path and run in Release. It appends a line per case, so a long run loses
    // nothing, and asserts no speed. Allocated bytes barely depend on load and track the toggles of column entries, so they are
    // the figure to read on a busy machine; the times are there for a quiet one. A size over a budget is the last.
    [Fact]
    public void ComputeCycles_Cost()
    {
        if (Environment.GetEnvironmentVariable("ZIGZAG_CYCLES_COST") is not { Length: > 0 } path) return;

        const double budgetSeconds = 60.0;
        Process process = Process.GetCurrentProcess();
        process.PriorityClass = ProcessPriorityClass.BelowNormal;

        void Report(string line)
        {
            _output.WriteLine(line);
            File.AppendAllText(path, line + Environment.NewLine);
        }

        (double Cpu, double Megabytes, T Result) Measure<T>(Func<T> run)
        {
            process.Refresh();
            TimeSpan cpu = process.TotalProcessorTime;
            long allocated = GC.GetAllocatedBytesForCurrentThread();
            T result = run();
            allocated = GC.GetAllocatedBytesForCurrentThread() - allocated;
            process.Refresh();
            return ((process.TotalProcessorTime - cpu).TotalSeconds, allocated / 1048576.0, result);
        }

        foreach (int n in new[] { 20, 40, 80, 160 })
        {
            var rng = new Random(n);
            GraphFrame[] frames = GraphFrameZigzagTests.EvolvingFrames(n, 6, rng, 0.0, pairFlip: 0.05);
            ZigzagFiltration f = GraphFrameZigzag.Build(n, frames, FrameTransition.Union, FrameComplex.Flag).Filtration;

            var bars = Measure(() => FastZigzag.Compute(f, 1));
            var named = Measure(() => FastZigzag.Compute(f, 1, representatives: true));
            var kept = Measure(() => FastZigzag.ComputeCycles(f, 1));
            ZigzagCycles cycles = kept.Result;

            long levels = 0, cells = 0;
            var read = Measure(() =>
            {
                foreach (Bar bar in cycles.Barcode.Bars)
                    for (int level = (int)bar.Birth + 1; level <= (int)bar.Death; level += Math.Max(1, ((int)bar.Death - (int)bar.Birth) / 8))
                    {
                        levels++;
                        cells += cycles.CycleAt(bar, level).Length;
                    }
                return 0;
            });

            Report($"n {n}, steps {f.Count}, bars {cycles.Barcode.Bars.Count}: MB allocated: bars {bars.Megabytes:F1}, representatives {named.Megabytes:F1}, " +
                $"cycles {kept.Megabytes:F1}; cycles/representatives {kept.Megabytes / named.Megabytes:F2}; cpu s: {bars.Cpu:F3}, {named.Cpu:F3}, {kept.Cpu:F3}; " +
                $"reading {levels} cycles of mean size {(double)cells / Math.Max(1, levels):F1}: {read.Megabytes:F1} MB, {read.Cpu:F3} cpu s");

            if (kept.Cpu > budgetSeconds) break;
        }
    }

    // ── helpers ─────────────────────────────────────────────────────────────────

    sealed class Reach
    {
        // Degree, then the kind of bar: closed or open birth, closed or open death.
        public readonly int[,] ByKind = new int[3, 4];
        public readonly int[] Moved = new int[3];
        public int ThroughReentry;
        public int Levels;
        public int LongestChain;

        public override string ToString() =>
            $"closed-closed/closed-open/open-closed/open-open by degree: " +
            string.Join("; ", Enumerable.Range(0, 3).Select(p => string.Join("/", Enumerable.Range(0, 4).Select(k => ByKind[p, k])))) +
            $"; moved {string.Join("/", Moved)}; through a re-entered cell {ThroughReentry}; levels read {Levels}; longest cycle {LongestChain}";
    }

    static void AssertEveryBarRepresented(ZigzagFiltration f, ZigzagCycles cycles, int maxDegree, Reach reach, string context = "")
    {
        for (int p = 0; p <= maxDegree; p++)
        {
            var oracle = new ZigzagModuleOracle(f, p);
            Bar[] bars = cycles.Barcode.Bars.Where(b => b.Dimension == p).ToArray();
            foreach (Bar bar in bars)
            {
                int first = (int)bar.Birth + 1, last = (int)bar.Death;
                string? failure = oracle.CheckSequence(bar, level => ZigzagModuleOracle.ChainOf(cycles.CycleAt(bar, level)));
                Assert.True(failure is null, $"{context}, degree {p}, bar [{bar.Birth}, {bar.Death}] {bar.BirthEnd}-{bar.DeathEnd}: {failure}");

                reach.ByKind[p, (bar.BirthEnd == IntervalEnd.Closed ? 0 : 2) + (bar.DeathEnd == IntervalEnd.Closed ? 0 : 1)]++;

                int[] start = cycles.CycleAt(bar, first);
                bool moved = false, reentered = false;
                for (int level = first; level <= last; level++)
                {
                    int[] cycle = cycles.CycleAt(bar, level);
                    Assert.True(cycle.SequenceEqual(cycle.OrderBy(x => x).Distinct()), $"{context}: cell ids repeat or are unsorted");
                    moved |= !cycle.SequenceEqual(start);
                    reentered |= cycle.Any(cell => DeletedBefore(f, cell, level));
                    reach.Levels++;
                    reach.LongestChain = Math.Max(reach.LongestChain, cycle.Length);
                }

                if (moved) reach.Moved[p]++;
                if (reentered) reach.ThroughReentry++;
            }

            for (int level = 0; level <= f.Count; level++)
            {
                BigInteger[] alive = bars.Where(b => b.Birth + 1 <= level && level <= b.Death)
                    .Select(b => ZigzagModuleOracle.ChainOf(cycles.CycleAt(b, level))).ToArray();
                string? failure = oracle.CheckBasis(level, alive);
                Assert.True(failure is null, $"{context}, degree {p}: {failure}");
            }
        }
    }

    static bool DeletedBefore(ZigzagFiltration f, int cell, int level)
    {
        for (int s = 0; s < level; s++)
            if (f[s].GlobalCellId == cell && f[s].Direction == ZigzagDirection.Delete)
                return true;
        return false;
    }

    static (double, double, int, IntervalEnd, IntervalEnd)[] Signature(Barcode barcode) =>
        barcode.Bars.Select(b => (b.Birth, b.Death, b.Dimension, b.BirthEnd, b.DeathEnd)).OrderBy(x => x).ToArray();

    static (int, int) ToPair(IReadOnlyList<int> v) => v[0] < v[1] ? (v[0], v[1]) : (v[1], v[0]);

    static (int i, int j) ToEdge((int, int) e) => (e.Item1, e.Item2);
}
