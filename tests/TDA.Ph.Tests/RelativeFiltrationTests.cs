#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Graphs.Primitives;
using Maths.Topology;
using Xunit;
using Xunit.Abstractions;

namespace TDA.Ph.Tests;

/// <summary>
/// The relative reading <c>H(K_δ, K₀)</c>, against hand-worked complexes and two references that do not go
/// through the relative complex: the quotient identity read by the involuted engine on the union, and reduced
/// homology of <c>K ∪ cone(K₀)</c>. Checks and expectations: <c>issues/ph/ph-setup-brief.md</c> (R1–R12).
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

    // ── R12 · one absolute run with K₀ first determines the relative barcode ────

    [Fact]
    public void RelativeBarcode_IsARereadingOfTheAbsoluteRunWithK0First()
    {
        const int trials = 200;
        const double k0Birth = -1.0;
        int kept = 0, movedUp = 0, dropped = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(4, 9);
            List<(int, int)> backbone = RandomEdges(n, rng, 0.3);
            if (backbone.Count == 0) backbone.Add((0, 1));
            List<(int, int, double)> content = RandomContent(n, rng, backbone, 0.35);

            RelativeFiltration relative = ConditionedFiltration.BuildRelative(n, backbone, content);

            // The same cells, with K₀ strictly before everything else, so a bar's birth says whether it began in K₀.
            HashSet<int> k0Vertices = backbone.SelectMany(e => new[] { e.Item1, e.Item2 }).ToHashSet();
            HashSet<(int, int)> k0Edges = backbone.Select(Ordered).ToHashSet();
            IEnumerable<Simplex> cells = relative.Ambient.Simplices
                .Select(s => InK0(s, k0Vertices, k0Edges) ? new Simplex(k0Birth, s.Vertices) : s);
            Barcode absolute = PersistentHomology.Compute(new SimplicialFiltration(cells));

            // Born outside K₀: the same bar. Born in K₀ and killed outside it: one degree up, born at that death, and
            // essential. Born in K₀ and essential, or killed within it: gone.
            var expected = new List<(int, double, double)>();
            foreach (Bar bar in absolute.Bars)
            {
                if (bar.Birth > k0Birth)
                {
                    if (bar.IsInfinite || bar.Death > bar.Birth) { expected.Add((bar.Dimension, bar.Birth, bar.Death)); kept++; }
                }
                else if (!bar.IsInfinite && bar.Death > k0Birth)
                {
                    expected.Add((bar.Dimension + 1, bar.Death, double.PositiveInfinity));
                    movedUp++;
                }
                else
                {
                    dropped++;
                }
            }

            Assert.True(expected.OrderBy(x => x).SequenceEqual(Canon(relative.ComputeBarcode())), $"seed {seed}");
        }

        _output.WriteLine($"bars kept: {kept}, moved up a degree: {movedUp}, dropped: {dropped}");
        Assert.True(kept > 0 && movedUp > 0 && dropped > 0);
    }

    // ── R13–R16 · the reading itself: one clearing run with K₀ strictly first ───

    [Fact]
    public void BothReadings_HandCases_MatchTheReference()
    {
        var cases = new (int N, List<(int, int)> Backbone, List<(int, int, double)> Content)[]
        {
            (6, new(), new() { (0, 1, 1.0), (1, 2, 2.0), (2, 3, 3.0), (3, 4, 4.0), (4, 5, 5.0), (5, 0, 6.0), (0, 2, 7.0) }), // R1
            (8, PathEdges(8), new() { (0, 4, 1.5), (2, 7, 2.5), (1, 6, 3.5) }),                                              // R2
            (6, CycleEdges(6), new() { (0, 3, 2.5) }),                                                                        // R3
            (4, CycleEdges(4), new() { (0, 2, 2.5) }),                                                                        // R4
            (6, new() { (0, 1), (1, 2) }, new() { (2, 3, 1.5), (4, 5, 2.5) }),                                                // R5
            (6, PathEdges(6), new() { (0, 5, 0.0) }),                                                                         // R6
            (5, PathEdges(5), new() { (0, 4, 1.0), (0, 2, 2.0), (2, 4, 3.0) }),                                               // R10
        };

        foreach (var (n, backbone, content) in cases)
        {
            RelativeFiltration relative = ConditionedFiltration.BuildRelative(n, backbone, content);
            AssertReadingsMatchReference(relative);
            AssertAbsoluteIsAmbient(relative);
        }
    }

    [Fact]
    public void BothReadings_AVoidK0Fills_MovesUpToH3_AndNothingKillsIt()
    {
        // K₀ is the boundary of a tetrahedron, its 2-sphere; K fills it at d. Absolutely H₂ [0, d); relatively H₃ [d, ∞).
        const double d = 2.5;
        var sphere = SimplicialFiltrationFixtures.TetrahedronBoundary().Simplices.Select(s => new Simplex(0.0, s.Vertices));
        var ambient = new SimplicialFiltration(sphere.Append(new Simplex(d, 0, 1, 2, 3)));
        var relative = new RelativeFiltration(ambient, s => s.Dimension < 3);

        var (absolute, read) = relative.ComputeBothReadings();

        Assert.Equal(new[] { (0, 0.0, double.PositiveInfinity), (2, 0.0, d) }, Canon(absolute));
        Bar bar = Assert.Single(Significant(read, 3));
        Assert.Equal((d, double.PositiveInfinity), (bar.Birth, bar.Death));
        Assert.Equal(3, ambient.Simplices[bar.Generator!.Value].Dimension);
        Assert.Equal(Canon(relative.ComputeBarcode()), Canon(read));
    }

    [Fact]
    public void BothReadings_RandomBackbones_MatchTheReference_AndReachEveryRow()
    {
        const int trials = 200;
        var rows = new Rows();
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(4, 10);
            // The R7 family, then the R8 family: connected spanning K₀ with cycles, then any K₀.
            List<(int, int)> backbone = seed % 2 == 0 ? RandomConnectedSpanning(n, rng) : RandomEdges(n, rng, 0.3);
            if (backbone.Count == 0) backbone.Add((0, 1));
            List<(int, int, double)> content = RandomContent(n, rng, backbone, 0.35);

            RelativeFiltration relative = ConditionedFiltration.BuildRelative(n, backbone, content);
            AssertReadingsMatchReference(relative, $"seed {seed}");
            AssertAbsoluteIsAmbient(relative, $"seed {seed}");
            rows.Count(relative);
        }

        _output.WriteLine(rows.ToString());
        Assert.True(rows.KeptFinite > 0 && rows.KeptEssential > 0 && rows.DroppedEssential > 0 && rows.DroppedFinite > 0);
        Assert.True(rows.MovedUp[0] > 0 && rows.MovedUp[1] > 0);
    }

    [Fact]
    public void BothReadings_AnyFiltration_AnySubcomplex_MatchTheReference_AndReachEveryRow()
    {
        const int trials = 300;
        var rows = new Rows();
        int k0Late = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            // A dense third of the trials, whose K0 holds whole 2-spheres, so a tetrahedron outside it can fill one.
            bool dense = seed % 3 == 0;
            SimplicialFiltration ambient = RandomComplex(rng, grid: seed % 2 == 0, fill: dense ? 0.95 : 0.55);
            HashSet<int> k0Vertices = Enumerable.Range(0, 12).Where(_ => rng.NextDouble() < (dense ? 0.9 : 0.6)).ToHashSet();
            Func<Simplex, bool> inK0 = RandomSubcomplex(ambient, k0Vertices, dense ? 0.95 : rng.NextDouble(), rng);

            var relative = new RelativeFiltration(ambient, inK0);
            AssertReadingsMatchReference(relative, $"seed {seed}");
            rows.Count(relative);

            // K₀ need not sit at the base here. Its reading is K with K₀ moved there, a filtration built by hand.
            double baseValue = ambient.GetBirth(0);
            var atBase = new SimplicialFiltration(ambient.Simplices.Select(s => inK0(s) ? new Simplex(baseValue, s.Vertices) : s));
            Assert.True(
                Canon(PersistentHomology.Compute(atBase)).SequenceEqual(Canon(relative.ComputeBothReadings().Absolute)),
                $"absolute, seed {seed}");

            if (ambient.Simplices.Any(s => inK0(s) && s.FiltrationValue > baseValue)) k0Late++;
        }

        _output.WriteLine($"{rows}; trials with K0 above the base: {k0Late}");
        Assert.True(rows.KeptFinite > 0 && rows.KeptEssential > 0 && rows.DroppedEssential > 0 && rows.DroppedFinite > 0);
        Assert.True(rows.MovedUp[0] > 0 && rows.MovedUp[1] > 0 && rows.MovedUp[2] > 0);
        // A moved-up class is essential even with cells above it to kill it: the claim R12 could not reach.
        Assert.True(rows.MovedUpWithCellsAbove > 0 && k0Late > 0);
    }

    // The one timing, off by default: RELATIVE_READING_TIMING=1, in Release. It reports and asserts no speed, since the
    // machine is shared. The graph's edges enter by a uniform draw, K₀ is a spanning tree with a few extra edges.
    [Fact]
    public void BothReadings_Timing()
    {
        if (Environment.GetEnvironmentVariable("RELATIVE_READING_TIMING") is not "1") return;

        Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.BelowNormal;
        foreach (var (n, p) in new[] { (100, 0.30), (200, 0.20), (400, 0.15), (800, 0.10) })
        {
            var rng = new Random(n);
            List<(int, int)> backbone = RandomConnectedSpanning(n, rng);
            HashSet<(int, int)> taken = backbone.Select(Ordered).ToHashSet();
            var content = new List<(int, int, double)>();
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    if (!taken.Contains((a, b)) && rng.NextDouble() < p)
                        content.Add((a, b, rng.NextDouble()));

            var built = Stopwatch.StartNew();
            RelativeFiltration relative = ConditionedFiltration.BuildRelative(n, backbone, content);
            built.Stop();
            int triangles = relative.Ambient.Simplices.Count(s => s.Dimension == 2);

            foreach (int maxDimension in new[] { 1, int.MaxValue })
            {
                var referenceMs = new List<double>();
                var readingMs = new List<double>();
                (int, double, double, int?)[]? referenceBars = null, readingBars = null;
                for (int repeat = 0; repeat < 3; repeat++)
                {
                    var watch = Stopwatch.StartNew();
                    Barcode reference = relative.ComputeBarcode(maxDimension);
                    referenceMs.Add(watch.Elapsed.TotalMilliseconds);

                    watch.Restart();
                    Barcode reading = relative.ComputeBothReadings(maxDimension).Relative;
                    readingMs.Add(watch.Elapsed.TotalMilliseconds);

                    referenceBars = Pairs(reference);
                    readingBars = Pairs(reading);
                }

                Assert.True(referenceBars!.SequenceEqual(readingBars!), $"n {n}, maxDimension {maxDimension}");
                double r = referenceMs.Min(), f = readingMs.Min();
                _output.WriteLine($"n {n}, edges {backbone.Count + content.Count}, triangles {triangles}, " +
                    $"maxDimension {(maxDimension == int.MaxValue ? "all" : maxDimension)}: build {built.Elapsed.TotalMilliseconds:F0} ms, " +
                    $"reference {r:F1} ms, one run {f:F1} ms, reference/run {r / f:F2} (best of 3; all {string.Join(", ", referenceMs.Select(x => x.ToString("F0")))} vs {string.Join(", ", readingMs.Select(x => x.ToString("F0")))})");
            }
        }
    }

    // The same pairs as the reference, at every cutoff degree: every bar with its generator, zero-length ones included.
    static void AssertReadingsMatchReference(RelativeFiltration relative, string context = "")
    {
        foreach (int maxDimension in new[] { int.MaxValue, 0, 1, 2 })
        {
            var expected = Pairs(relative.ComputeBarcode(maxDimension));
            var actual = Pairs(relative.ComputeBothReadings(maxDimension).Relative);
            Assert.True(expected.SequenceEqual(actual), $"{context}, maxDimension {maxDimension}");
        }
    }

    // With K₀ at the base already, which a backbone at ε₀ = 0 is, the absolute reading is the ambient's own barcode.
    static void AssertAbsoluteIsAmbient(RelativeFiltration relative, string context = "")
    {
        var expected = Canon(PersistentHomology.Compute(relative.Ambient));
        Assert.True(expected.SequenceEqual(Canon(relative.ComputeBothReadings().Absolute)), $"{context}, absolute");
    }

    static (int, double, double, int?)[] Pairs(Barcode barcode) =>
        barcode.Bars
            .Select(b => (b.Dimension, b.Birth, b.Death, b.Generator))
            .OrderBy(x => x)
            .ToArray();

    sealed class Rows
    {
        public int KeptFinite, KeptEssential, DroppedFinite, DroppedEssential, MovedUpWithCellsAbove;
        public readonly int[] MovedUp = new int[4];

        public void Count(RelativeFiltration relative)
        {
            int top = relative.Ambient.Simplices.Max(s => s.Dimension);
            foreach (AbsolutePair pair in relative.AbsolutePairsWithK0First(int.MaxValue))
            {
                switch (pair.Row)
                {
                    case RelativeRow.Kept when pair.Death < 0: KeptEssential++; break;
                    case RelativeRow.Kept: KeptFinite++; break;
                    case RelativeRow.Dropped when pair.Death < 0: DroppedEssential++; break;
                    case RelativeRow.Dropped: DroppedFinite++; break;
                    case RelativeRow.MovedUp:
                        MovedUp[pair.Dimension]++;
                        if (pair.Dimension + 2 <= top) MovedUpWithCellsAbove++;
                        break;
                }
            }
        }

        public override string ToString() =>
            $"kept finite {KeptFinite}, kept essential {KeptEssential}, moved up from H0/H1/H2 {string.Join("/", MovedUp.Take(3))} " +
            $"({MovedUpWithCellsAbove} with cells above), dropped finite {DroppedFinite}, dropped essential {DroppedEssential}";
    }

    // Any complex the faces allow, with vertices born apart and cells up to tetrahedra.
    static SimplicialFiltration RandomComplex(Random rng, bool grid, double fill)
    {
        int n = rng.Next(4, 10);
        double Later(double earliest) => grid ? earliest + rng.Next(0, 4) : earliest + rng.NextDouble() * 3.0;

        var births = new Dictionary<string, double>();
        var cells = new List<Simplex>();
        foreach (int[] candidate in Faces(n))
        {
            string[] faces = candidate.Length == 1 ? Array.Empty<string>()
                : Enumerable.Range(0, candidate.Length).Select(k => Key(candidate.Where((_, i) => i != k))).ToArray();
            if (!faces.All(births.ContainsKey) || (candidate.Length > 1 && rng.NextDouble() >= fill))
                continue;

            double birth = candidate.Length == 1 ? rng.Next(0, 5) : Later(faces.Max(f => births[f]));
            births[Key(candidate)] = birth;
            cells.Add(new Simplex(birth, candidate));
        }

        return new SimplicialFiltration(cells);
    }

    // Candidates by increasing size, so faces are decided first.
    static IEnumerable<int[]> Faces(int n)
    {
        for (int a = 0; a < n; a++) yield return new[] { a };
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) yield return new[] { a, b };
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) for (int c = b + 1; c < n; c++) yield return new[] { a, b, c };
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) for (int c = b + 1; c < n; c++) for (int d = c + 1; d < n; d++)
            yield return new[] { a, b, c, d };
    }

    static string Key(IEnumerable<int> vertices) => string.Join(",", vertices);

    // A random subcomplex: a cell joins with chance p once all its faces have, on the given vertices, and a tetrahedron
    // with a smaller chance, since one outside K0 is what kills a 2-sphere inside it.
    static Func<Simplex, bool> RandomSubcomplex(SimplicialFiltration ambient, HashSet<int> vertices, double p, Random rng)
    {
        var members = new HashSet<string>();
        foreach (Simplex s in ambient.Simplices.OrderBy(s => s.Dimension))
        {
            int[] v = s.Vertices;
            bool joins = v.Length == 1
                ? vertices.Contains(v[0])
                : Enumerable.Range(0, v.Length).All(k => members.Contains(Key(v.Where((_, i) => i != k)))) && rng.NextDouble() < (v.Length == 4 ? p / 3 : p);
            if (joins) members.Add(Key(v));
        }

        return s => members.Contains(Key(s.Vertices));
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
