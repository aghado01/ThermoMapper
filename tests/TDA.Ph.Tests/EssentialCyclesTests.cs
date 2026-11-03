#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Numerics;
using Maths.Topology;
using Xunit;
using Xunit.Abstractions;

namespace TDA.Ph.Tests;

/// <summary>
/// Whether a barcode's cycles are a persistent basis, by exact linear algebra outside the engines: each cycle closes, lies
/// in the complex where its bar begins, and is a boundary where its bar dies; and at every value, the classes of the bars
/// alive there are a basis of the homology there.
/// </summary>
internal sealed class PersistenceBasisOracle
{
    readonly IFiltration _filtration;
    readonly BigInteger[] _boundary;
    readonly double[] _values;

    public PersistenceBasisOracle(IFiltration filtration)
    {
        _filtration = filtration;
        _boundary = new BigInteger[filtration.Count];
        for (int c = 0; c < filtration.Count; c++)
            _boundary[c] = ChainOf(filtration.GetBoundaryIndices(c));

        _values = Enumerable.Range(0, filtration.Count).Select(filtration.GetBirth).Distinct().OrderBy(v => v).ToArray();
    }

    public static BigInteger ChainOf(IEnumerable<int> cells)
    {
        BigInteger chain = BigInteger.Zero;
        foreach (int c in cells) chain ^= BigInteger.One << c;
        return chain;
    }

    /// <summary>Null when the bars' cycles are a persistent basis in the degrees up to <paramref name="topDegree"/>, else what is wrong.</summary>
    public string? Check(IReadOnlyList<Bar> bars, int topDegree)
    {
        foreach (Bar bar in bars)
        {
            if (bar.Dimension > topDegree) continue;

            string name = $"H{bar.Dimension} [{bar.Birth}, {bar.Death})";
            if (bar.Cycle is null) return $"{name} has no cycle";
            if (bar.Cycle.Any(c => _filtration.GetDimension(c) != bar.Dimension)) return $"{name} has a cell of another degree";
            if (bar.Cycle.Any(c => _filtration.GetBirth(c) > bar.Birth)) return $"{name} has a cell from after it begins";
            if (bar.Generator is int g && bar.Cycle.Max() != g) return $"{name}'s latest cell is not its generator";

            BigInteger chain = ChainOf(bar.Cycle);
            if (!BoundaryOf(chain).IsZero) return $"{name}'s cycle does not close";

            if (!bar.IsInfinite)
            {
                var boundaries = new Echelon();
                AddBoundaries(boundaries, bar.Dimension, bar.Death);
                if (boundaries.Insert(chain)) return $"{name}'s cycle is not a boundary where it dies";
            }
        }

        foreach (double value in _values)
            for (int degree = 0; degree <= topDegree; degree++)
            {
                var classes = new Echelon();
                AddBoundaries(classes, degree, value);
                int boundaryRank = classes.Rank;

                var chains = new Echelon();
                int cellCount = 0;
                for (int c = 0; c < _filtration.Count && _filtration.GetBirth(c) <= value; c++)
                    if (_filtration.GetDimension(c) == degree)
                    {
                        cellCount++;
                        chains.Insert(_boundary[c]);
                    }

                int betti = cellCount - chains.Rank - boundaryRank;

                int alive = 0;
                foreach (Bar bar in bars)
                {
                    if (bar.Dimension != degree || bar.Birth > value || bar.Death <= value) continue;

                    alive++;
                    if (!classes.Insert(ChainOf(bar.Cycle!)))
                        return $"at {value}, H{degree} [{bar.Birth}, {bar.Death})'s class depends on the others";
                }

                if (alive != betti) return $"at {value}, H{degree} has rank {betti} with {alive} bars alive";
            }

        return null;
    }

    BigInteger BoundaryOf(BigInteger chain)
    {
        BigInteger boundary = BigInteger.Zero;
        for (int c = 0; !chain.IsZero; c++, chain >>= 1)
            if (!chain.IsEven)
                boundary ^= _boundary[c];

        return boundary;
    }

    void AddBoundaries(Echelon echelon, int degree, double value)
    {
        for (int c = 0; c < _filtration.Count && _filtration.GetBirth(c) <= value; c++)
            if (_filtration.GetDimension(c) == degree + 1)
                echelon.Insert(_boundary[c]);
    }

    /// <summary>Vectors over Z/2 in echelon form.</summary>
    internal sealed class Echelon
    {
        readonly Dictionary<long, BigInteger> _byPivot = new();

        public int Rank => _byPivot.Count;

        /// <summary>Adds the vector, and says whether it raised the rank.</summary>
        public bool Insert(BigInteger vector)
        {
            while (!vector.IsZero)
            {
                long pivot = (long)vector.GetBitLength();
                if (!_byPivot.TryGetValue(pivot, out BigInteger owner))
                {
                    _byPivot[pivot] = vector;
                    return true;
                }

                vector ^= owner;
            }

            return false;
        }
    }
}

/// <summary>
/// A cycle on every bar, from <see cref="PersistentHomology"/> and <see cref="PersistentInvolutedHomology"/>, the essential
/// ones read from <c>V</c>: by hand on closed surfaces, and over random complexes against <see cref="PersistenceBasisOracle"/>.
/// </summary>
public sealed class EssentialCyclesTests
{
    readonly ITestOutputHelper _output;

    public EssentialCyclesTests(ITestOutputHelper output) => _output = output;

    // ── closed surfaces: over Z/2 the one class in the top degree is every triangle ──

    static readonly int[][] HollowTetrahedron = { new[] { 0, 1, 2 }, new[] { 0, 1, 3 }, new[] { 0, 2, 3 }, new[] { 1, 2, 3 } };

    static readonly int[][] ProjectivePlane =
    {
        new[] { 0, 1, 3 }, new[] { 0, 1, 5 }, new[] { 0, 2, 3 }, new[] { 0, 2, 4 }, new[] { 0, 4, 5 },
        new[] { 1, 2, 4 }, new[] { 1, 2, 5 }, new[] { 1, 3, 4 }, new[] { 2, 3, 5 }, new[] { 3, 4, 5 },
    };

    static readonly int[][] Torus = Enumerable.Range(0, 7)
        .SelectMany(i => new[] { new[] { i, (i + 1) % 7, (i + 3) % 7 }, new[] { i, (i + 2) % 7, (i + 3) % 7 } })
        .ToArray();

    public static IEnumerable<object[]> Surfaces()
    {
        yield return new object[] { "sphere", HollowTetrahedron, new[] { 1, 0, 1 } };
        yield return new object[] { "projective plane", ProjectivePlane, new[] { 1, 1, 1 } };
        yield return new object[] { "torus", Torus, new[] { 1, 2, 1 } };
    }

    [Theory]
    [MemberData(nameof(Surfaces))]
    public void ClosedSurface_TheVoidsCycleIsEveryTriangle(string name, int[][] triangles, int[] betti)
    {
        _output.WriteLine(name);
        SimplicialFiltration filtration = Surface(triangles);
        int[] everyTriangle = Enumerable.Range(0, filtration.Count).Where(c => filtration.GetDimension(c) == 2).ToArray();
        var oracle = new PersistenceBasisOracle(filtration);

        foreach (Barcode barcode in new[]
        {
            PersistentHomology.Compute(filtration, representatives: true),
            PersistentInvolutedHomology.Compute(filtration, representatives: true),
        })
        {
            for (int d = 0; d <= 2; d++)
                Assert.Equal(betti[d], barcode.Bars.Count(b => b.Dimension == d && b.IsInfinite));

            Bar voidBar = Assert.Single(barcode.Bars, b => b.Dimension == 2);
            Assert.Equal(everyTriangle, voidBar.Cycle);
            Assert.Null(oracle.Check(barcode.Bars, 2));
        }
    }

    [Fact]
    public void BoundaryOfAFourSimplex_TheCycleIsEveryTetrahedron()
    {
        var simplices = new List<Simplex>();
        int[] vertices = { 0, 1, 2, 3, 4 };
        foreach (int a in vertices)
        {
            simplices.Add(new Simplex(0.0, a));
            foreach (int b in vertices.Where(v => v > a))
            {
                simplices.Add(new Simplex(1.0, a, b));
                foreach (int c in vertices.Where(v => v > b))
                {
                    simplices.Add(new Simplex(2.0, a, b, c));
                    foreach (int d in vertices.Where(v => v > c))
                        simplices.Add(new Simplex(3.0 + a + b + c + d, a, b, c, d));
                }
            }
        }

        var filtration = new SimplicialFiltration(simplices);
        int[] everyTetrahedron = Enumerable.Range(0, filtration.Count).Where(c => filtration.GetDimension(c) == 3).ToArray();

        foreach (Barcode barcode in new[]
        {
            PersistentHomology.Compute(filtration, representatives: true),
            PersistentInvolutedHomology.Compute(filtration, representatives: true),
        })
        {
            Bar sphere = Assert.Single(barcode.Bars, b => b.Dimension == 3);
            Assert.True(sphere.IsInfinite);
            Assert.Equal(everyTetrahedron, sphere.Cycle);
            Assert.Null(new PersistenceBasisOracle(filtration).Check(barcode.Bars, 3));
        }
    }

    [Fact]
    public void CellsWithNoBoundary_AreTheirOwnCycles()
    {
        var torus = new CellFiltration(new (int, double, int[])[]
        {
            (0, 0, Array.Empty<int>()), (1, 1, Array.Empty<int>()), (1, 2, Array.Empty<int>()), (2, 3, Array.Empty<int>()),
        });

        Barcode barcode = PersistentHomology.Compute(torus, representatives: true);

        Assert.Equal(4, barcode.Bars.Count);
        Assert.All(barcode.Bars, b => Assert.Equal(new[] { b.Generator!.Value }, b.Cycle));
    }

    [Fact]
    public void H0_NamesItsRoots_WhereAColumnCanNameAVertexBelowOne()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0), new Simplex(1.0, 1), new Simplex(2.0, 2), new Simplex(3.0, 3),
            new Simplex(4.0, 0, 1), new Simplex(5.0, 2, 3), new Simplex(6.0, 1, 2),
        });

        List<Bar> H0(Barcode barcode) => barcode.Bars.Where(b => b.Dimension == 0).OrderBy(b => b.Birth).ToList();

        Barcode involuted = PersistentInvolutedHomology.Compute(filtration, representatives: true);
        Assert.Equal(new[] { 0 }, H0(involuted)[0].Cycle);
        Assert.Equal(new[] { 0, 1 }, H0(involuted)[1].Cycle);
        Assert.Equal(new[] { 0, 2 }, H0(involuted)[2].Cycle);
        Assert.Equal(new[] { 2, 3 }, H0(involuted)[3].Cycle);

        // The last edge joins 1 to 2, and the column stops at its pivot: vertex 1 stands for its root 0.
        Barcode column = PersistentHomology.Compute(filtration, representatives: true);
        Assert.Equal(new[] { 1, 2 }, H0(column)[2].Cycle);

        var oracle = new PersistenceBasisOracle(filtration);
        Assert.Null(oracle.Check(involuted.Bars, 0));
        Assert.Null(oracle.Check(column.Bars, 0));
    }

    // ── sweeps ──────────────────────────────────────────────────────────────────

    [Fact]
    public void RandomComplexes_EveryBarsCycleIsAPersistentBasis_InBothEngines()
    {
        const int trials = 300;
        var essential = new int[4];
        var finite = new int[4];
        int largest = 0;

        for (int seed = 0; seed < trials; seed++)
        {
            SimplicialFiltration filtration = RandomComplex(new Random(seed));
            var oracle = new PersistenceBasisOracle(filtration);

            Barcode column = PersistentHomology.Compute(filtration, representatives: true);
            Barcode involuted = PersistentInvolutedHomology.Compute(filtration, representatives: true);

            Assert.Equal((string?)null, oracle.Check(column.Bars, 3) is { } c ? $"column, seed {seed}: {c}" : null);
            Assert.Equal((string?)null, oracle.Check(involuted.Bars, 3) is { } i ? $"involuted, seed {seed}: {i}" : null);

            // The same reduction, so the same chains, but for H0's roots and an essential H1 bar's path, which only need the same class.
            Assert.True(BarsWithCycles(column).SequenceEqual(BarsWithCycles(involuted)), $"seed {seed}");

            // Stopped below the top degree, the degree it stops at still reads its chains.
            Assert.True(
                BarsWithCycles(PersistentHomology.Compute(filtration, 2, representatives: true))
                    .SequenceEqual(BarsWithCycles(PersistentInvolutedHomology.Compute(filtration, 2, representatives: true))),
                $"stopped at H2, seed {seed}");

            foreach (Bar bar in column.Bars)
            {
                (bar.IsInfinite ? essential : finite)[bar.Dimension]++;
                if (bar.IsInfinite) largest = Math.Max(largest, bar.Cycle!.Length);
            }
        }

        _output.WriteLine($"finite H0-H3: {string.Join("/", finite)}; essential H0-H3: {string.Join("/", essential)}; largest essential cycle: {largest}");
        Assert.All(essential, count => Assert.True(count > 0));
    }

    [Fact]
    public void RandomChainComplexes_EveryBarsCycleIsAPersistentBasis_AndH0ClearingNamesTheSame()
    {
        const int trials = 400;
        int essentialAbove0 = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            var cells = new CellFiltration(CellComplexClearingTests.RandomChainComplex(rng, graphLike: seed % 2 == 0));
            var oracle = new PersistenceBasisOracle(cells);

            Barcode column = PersistentHomology.Compute(cells, representatives: true);
            Assert.Equal((string?)null, oracle.Check(column.Bars, 3) is { } failure ? $"seed {seed}: {failure}" : null);
            essentialAbove0 += column.Bars.Count(b => b.IsInfinite && b.Dimension > 0);

            (List<Bar> h0, _) = PersistenceClearing.ComputeH0(cells, int.MaxValue, 0.0, cycles: true);
            Assert.True(
                BarsWithCycles(new Barcode(h0)).SequenceEqual(BarsWithCycles(new Barcode(column.Bars.Where(b => b.Dimension == 0).ToList()))),
                $"H0, seed {seed}");
            Assert.Equal((string?)null, oracle.Check(h0, 0) is { } h0Failure ? $"H0, seed {seed}: {h0Failure}" : null);
        }

        _output.WriteLine($"essential above H0: {essentialAbove0}");
        Assert.True(essentialAbove0 > trials);
    }

    [Fact]
    public void RandomRelativeComplexes_EssentialChainsAreABasis_AndFillWhatKFills()
    {
        const int trials = 300;
        int essential = 0, movedUp = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            SimplicialFiltration ambient = RandomComplex(rng);

            var inK0 = new bool[ambient.Count];
            var chosen = new HashSet<Simplex>();
            for (int c = 0; c < ambient.Count; c++)
            {
                inK0[c] = ambient.GetBoundaryIndices(c).All(f => inK0[f]) && rng.NextDouble() < 0.5;
                if (inK0[c]) chosen.Add(ambient.Simplices[c]);
            }

            var relative = new RelativeFiltration(ambient, chosen.Contains);
            Barcode barcode = relative.ComputeBarcode(representatives: true);

            var relativeOf = new Dictionary<int, int>();
            for (int r = 0; r < relative.Count; r++) relativeOf[relative.AmbientIndices[r]] = r;

            List<Bar> bars = barcode.Bars.Select(b => b with
            {
                Generator = relativeOf[b.Generator!.Value],
                Cycle = Array.ConvertAll(b.Cycle!, c => relativeOf[c]),
            }).ToList();

            var oracle = new PersistenceBasisOracle(relative);
            Assert.Equal((string?)null, oracle.Check(bars, 3) is { } failure ? $"seed {seed}: {failure}" : null);

            // The connecting map: the boundaries of the essential chains span what K fills of H(K₀), as many as moved up.
            var ambientBoundary = new BigInteger[ambient.Count];
            for (int c = 0; c < ambient.Count; c++) ambientBoundary[c] = PersistenceBasisOracle.ChainOf(ambient.GetBoundaryIndices(c));

            List<AbsolutePair> pairs = relative.AbsolutePairsWithK0First(int.MaxValue);
            for (int degree = 0; degree <= 2; degree++)
            {
                var classes = new PersistenceBasisOracle.Echelon();
                for (int c = 0; c < ambient.Count; c++)
                    if (inK0[c] && ambient.GetDimension(c) == degree + 1)
                        classes.Insert(ambientBoundary[c]);

                int filled = 0;
                foreach (Bar bar in barcode.Bars.Where(b => b.IsInfinite && b.Dimension == degree + 1))
                {
                    essential++;
                    BigInteger boundary = BigInteger.Zero;
                    foreach (int c in bar.Cycle!) boundary ^= ambientBoundary[c];

                    for (int c = 0; c < ambient.Count; c++)
                        if (!inK0[c])
                            Assert.True(((boundary >> c) & BigInteger.One).IsZero, $"seed {seed}: the boundary leaves K₀");

                    if (classes.Insert(boundary)) filled++;
                }

                int expected = pairs.Count(p => p.Row == RelativeRow.MovedUp && p.Dimension == degree);
                movedUp += expected;
                Assert.True(expected == filled, $"seed {seed}, degree {degree}: {filled} filled, {expected} moved up");
            }
        }

        _output.WriteLine($"essential relative bars above H0: {essential}; moved up: {movedUp}");
        Assert.True(movedUp > trials / 2);
    }

    [Theory]
    [InlineData(0.5)]
    [InlineData(1.0)]
    [InlineData(2.5)]
    public void Involuted_WithACutoff_IsTheColumnReducerFiltered(double cutoff)
    {
        const int trials = 200;
        int dropped = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            SimplicialFiltration filtration = RandomComplex(new Random(seed));

            Barcode column = PersistentHomology.Compute(filtration, representatives: true);
            var kept = column.Bars.Where(b => b.Persistence >= cutoff).ToList();
            dropped += column.Bars.Count - kept.Count;

            Barcode involuted = PersistentInvolutedHomology.Compute(filtration, cutoff: cutoff, representatives: true);
            Assert.True(BarsWithCycles(new Barcode(kept)).SequenceEqual(BarsWithCycles(involuted)), $"seed {seed}");
        }

        _output.WriteLine($"cutoff {cutoff}: {dropped} bars below it");
        Assert.True(dropped > trials);
    }

    [Fact]
    public void WrongCycles_AreRefused()
    {
        int birthCellAlone = 0, withALaterCycle = 0, repeated = 0, boundary = 0, neverDies = 0, tried = 0;
        for (int seed = 0; seed < 300; seed++)
        {
            SimplicialFiltration filtration = RandomComplex(new Random(seed));
            var oracle = new PersistenceBasisOracle(filtration);
            Barcode barcode = PersistentHomology.Compute(filtration, representatives: true);

            List<Bar> bars = barcode.Bars.ToList();
            string? With(int at, int[] cycle)
            {
                var changed = new List<Bar>(bars);
                changed[at] = changed[at] with { Cycle = cycle };
                return oracle.Check(changed, 3);
            }

            for (int b = 0; b < bars.Count; b++)
            {
                Bar bar = bars[b];

                // A finite bar named by a class that outlives it, born before it so the rest of the check still reads it.
                if (!bar.IsInfinite && bars.Where(o => o.IsInfinite && o.Dimension == bar.Dimension && o.Birth <= bar.Birth).Select(o => (Bar?)o).FirstOrDefault() is { } outlives)
                {
                    var changed = new List<Bar>(bars);
                    changed[b] = bar with { Cycle = outlives.Cycle, Generator = null };
                    string? failure = oracle.Check(changed, 3);
                    Assert.NotNull(failure);
                    Assert.Contains("not a boundary where it dies", failure);
                    neverDies++;
                }

                if (!bar.IsInfinite || bar.Dimension < 2) continue;
                tried++;

                // What the involuted engine gave H2 and up before.
                if (With(b, new[] { bar.Generator!.Value }) is not null) birthCellAlone++;

                // Another essential class of the degree, born later, or the bar's own twice over.
                foreach (Bar later in bars.Where(o => o.IsInfinite && o.Dimension == bar.Dimension && o.Generator > bar.Generator))
                {
                    int[] sum = bar.Cycle!.Except(later.Cycle!).Union(later.Cycle!.Except(bar.Cycle!)).OrderBy(c => c).ToArray();
                    Assert.NotNull(With(b, sum));
                    withALaterCycle++;

                    Assert.NotNull(With(b, later.Cycle!));
                    repeated++;
                }

                // The class that is zero, as the old H1 tie bug named: a boundary from the same degree, if the complex has one in time.
                Bar? killed = bars.Where(o => !o.IsInfinite && o.Dimension == bar.Dimension && o.Death <= bar.Birth).Select(o => (Bar?)o).FirstOrDefault();
                if (killed is { } dead)
                {
                    Assert.NotNull(With(b, dead.Cycle!));
                    boundary++;
                }
            }
        }

        _output.WriteLine($"tried {tried}: birth cell alone refused {birthCellAlone}; with a later cycle {withALaterCycle}, repeated {repeated}, a boundary {boundary}; a finite bar that never dies {neverDies}");
        Assert.Equal(tried, birthCellAlone);
        Assert.True(withALaterCycle > 0 && repeated > 0 && boundary > 0 && neverDies > 0);
    }

    [Fact]
    public void MaxCycleDimension_KeepsEveryBar_AndNamesCyclesOnlyUpToIt()
    {
        int named = 0, unnamed = 0;

        void Check(Barcode whole, Barcode limited, int limit)
        {
            Assert.Equal(whole.Bars.Count, limited.Bars.Count);
            for (int b = 0; b < whole.Bars.Count; b++)
            {
                Bar expected = whole.Bars[b], bar = limited.Bars[b];
                Assert.Equal((expected.Dimension, expected.Birth, expected.Death, expected.Generator), (bar.Dimension, bar.Birth, bar.Death, bar.Generator));

                if (bar.Dimension <= limit)
                {
                    Assert.Equal(expected.Cycle, bar.Cycle);
                    named++;
                }
                else
                {
                    Assert.Null(bar.Cycle);
                    unnamed++;
                }
            }
        }

        for (int seed = 0; seed < 200; seed++)
        {
            SimplicialFiltration filtration = RandomComplex(new Random(seed));
            var relative = new RelativeFiltration(filtration, s => s.Dimension == 0);

            Barcode column = PersistentHomology.Compute(filtration, representatives: true);
            Barcode involuted = PersistentInvolutedHomology.Compute(filtration, representatives: true);
            Barcode quotient = relative.ComputeBarcode(representatives: true);

            for (int limit = -1; limit <= 3; limit++)
            {
                Check(column, PersistentHomology.Compute(filtration, representatives: true, maxCycleDimension: limit), limit);
                Check(involuted, PersistentInvolutedHomology.Compute(filtration, representatives: true, maxCycleDimension: limit), limit);
                Check(quotient, relative.ComputeBarcode(representatives: true, maxCycleDimension: limit), limit);
            }
        }

        _output.WriteLine($"bars named {named}, left without a cycle {unnamed}");
        Assert.True(named > 1000 && unnamed > 1000);
    }

    [Fact]
    public void LargerSkeleton_BothEnginesNameTheSameCycles_AndTheyClose()
    {
        int voids = 0, largest = 0;
        foreach (int n in new[] { 20, 30, 40 })
        {
            SimplicialFiltration filtration = RandomTwoSkeleton(n, new Random(n));
            Barcode column = PersistentHomology.Compute(filtration, representatives: true);
            Barcode involuted = PersistentInvolutedHomology.Compute(filtration, representatives: true);

            Assert.True(BarsWithCycles(column).SequenceEqual(BarsWithCycles(involuted)), $"n {n}");

            foreach (Bar bar in column.Bars.Where(b => b.IsInfinite && b.Dimension == 2))
            {
                var open = new HashSet<int>();
                foreach (int c in bar.Cycle!)
                    foreach (int face in filtration.GetBoundaryIndices(c))
                        if (!open.Remove(face)) open.Add(face);

                Assert.Empty(open);
                voids++;
                largest = Math.Max(largest, bar.Cycle!.Length);
            }
        }

        _output.WriteLine($"essential H2: {voids}, the largest of {largest} triangles");
        Assert.True(voids > 1000 && largest > 4);
    }

    [Fact]
    public void ProtectSet_ReadsTheSameEdges_FromH1Alone()
    {
        for (int seed = 0; seed < 100; seed++)
        {
            SimplicialFiltration filtration = RandomComplex(new Random(seed));

            string[] Edges(Barcode barcode) => BarCycleEdges.H1Loops(barcode, filtration)
                .SelectMany(loop => loop.Edges).Select(e => $"{e.Lo}-{e.Hi}").Distinct().OrderBy(e => e, StringComparer.Ordinal).ToArray();

            Assert.Equal(
                Edges(PersistentInvolutedHomology.Compute(filtration, representatives: true)),
                Edges(PersistentInvolutedHomology.Compute(filtration, maxDimension: 1, representatives: true)));
        }
    }

    // ── cost: bytes, which a busy machine barely moves ──────────────────────────

    [Fact]
    public void EssentialCycles_Cost()
    {
        if (Environment.GetEnvironmentVariable("ESSENTIAL_CYCLES_COST") is not { Length: > 0 } path) return;

        const double budgetSeconds = 120.0;
        Process process = Process.GetCurrentProcess();
        process.PriorityClass = ProcessPriorityClass.BelowNormal;

        void Report(string line)
        {
            _output.WriteLine(line);
            File.AppendAllText(path, line + Environment.NewLine);
        }

        (double Cpu, double Megabytes, Barcode Result) Measure(Func<Barcode> run)
        {
            process.Refresh();
            TimeSpan cpu = process.TotalProcessorTime;
            long allocated = GC.GetAllocatedBytesForCurrentThread();
            Barcode result = run();
            allocated = GC.GetAllocatedBytesForCurrentThread() - allocated;
            process.Refresh();
            return ((process.TotalProcessorTime - cpu).TotalSeconds, allocated / 1048576.0, result);
        }

        foreach (int n in new[] { 20, 40, 80, 160 })
        {
            SimplicialFiltration filtration = RandomTwoSkeleton(n, new Random(n));

            var columnBars = Measure(() => PersistentHomology.Compute(filtration));
            var columnCycles = Measure(() => PersistentHomology.Compute(filtration, representatives: true));
            var involutedBars = Measure(() => PersistentInvolutedHomology.Compute(filtration));
            var involutedCycles = Measure(() => PersistentInvolutedHomology.Compute(filtration, representatives: true));
            var columnUpToH1 = Measure(() => PersistentHomology.Compute(filtration, representatives: true, maxCycleDimension: 1));
            var involutedUpToH1 = Measure(() => PersistentInvolutedHomology.Compute(filtration, representatives: true, maxCycleDimension: 1));

            var voids = columnCycles.Result.Bars.Where(b => b.IsInfinite && b.Dimension == 2).ToList();
            double meanSize = voids.Count == 0 ? 0.0 : voids.Average(b => b.Cycle!.Length);

            Report($"n {n}, cells {filtration.Count}, essential H2 {voids.Count} of mean size {meanSize:F1}, largest {(voids.Count == 0 ? 0 : voids.Max(b => b.Cycle!.Length))}: " +
                $"MB column {columnBars.Megabytes:F1} -> {columnCycles.Megabytes:F1} ({columnCycles.Megabytes / columnBars.Megabytes:F2}); " +
                $"involuted {involutedBars.Megabytes:F1} -> {involutedCycles.Megabytes:F1} ({involutedCycles.Megabytes / involutedBars.Megabytes:F2}); " +
                $"up to H1: column {columnUpToH1.Megabytes:F1} ({columnUpToH1.Megabytes / columnBars.Megabytes:F2}), involuted {involutedUpToH1.Megabytes:F1} ({involutedUpToH1.Megabytes / involutedBars.Megabytes:F2}); " +
                $"cpu s column {columnBars.Cpu:F3} -> {columnCycles.Cpu:F3} (up to H1 {columnUpToH1.Cpu:F3}), involuted {involutedBars.Cpu:F3} -> {involutedCycles.Cpu:F3} (up to H1 {involutedUpToH1.Cpu:F3})");

            if (columnCycles.Cpu + involutedCycles.Cpu > budgetSeconds) break;
        }
    }

    // ── helpers ─────────────────────────────────────────────────────────────────

    static SimplicialFiltration Surface(int[][] triangles)
    {
        var simplices = new List<Simplex>();
        var seen = new HashSet<string>();
        foreach (int[] triangle in triangles)
        {
            int[] t = triangle.OrderBy(v => v).ToArray();
            foreach (int v in t)
                if (seen.Add($"{v}")) simplices.Add(new Simplex(0.0, v));

            foreach ((int a, int b) in new[] { (t[0], t[1]), (t[0], t[2]), (t[1], t[2]) })
                if (seen.Add($"{a}-{b}")) simplices.Add(new Simplex(1.0, a, b));

            simplices.Add(new Simplex(2.0 + simplices.Count, t[0], t[1], t[2]));
        }

        return new SimplicialFiltration(simplices);
    }

    /// <summary>Every bar with its cycle in a canonical order, but for H0's finite roots and an essential H1 bar's path. Exact on the values, since both engines read them.</summary>
    static string[] BarsWithCycles(Barcode barcode) => barcode.Bars
        .Select(b => $"{b.Dimension}|{b.Birth:R}|{b.Death:R}|{(b.Dimension == 0 && !b.IsInfinite ? "roots" : b.Dimension == 1 && b.IsInfinite ? "path" : string.Join(",", b.Cycle!))}")
        .OrderBy(s => s, StringComparer.Ordinal)
        .ToArray();

    /// <summary>Up to tetrahedra on a few vertices, with values that tie or not and can be negative, and a third of them dense enough to close a 3-cycle.</summary>
    internal static SimplicialFiltration RandomComplex(Random rng)
    {
        int n = rng.Next(4, 9);
        bool grid = rng.Next(2) == 0;
        double shift = rng.Next(3) == 0 ? -5.0 : 0.0;
        bool dense = rng.Next(3) == 0;
        double Later(double low) => low + (grid ? rng.Next(0, 3) : rng.NextDouble() * 2.0);

        var simplices = new List<Simplex>();
        var vertex = new double[n];
        for (int v = 0; v < n; v++)
        {
            vertex[v] = Later(shift);
            simplices.Add(new Simplex(vertex[v], v));
        }

        var edge = new Dictionary<(int, int), double>();
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
                if (dense || rng.NextDouble() < 0.8)
                {
                    edge[(a, b)] = Later(Math.Max(vertex[a], vertex[b]));
                    simplices.Add(new Simplex(edge[(a, b)], a, b));
                }

        var triangle = new Dictionary<(int, int, int), double>();
        foreach (var ((a, b), ab) in edge)
            for (int c = b + 1; c < n; c++)
                if (edge.TryGetValue((a, c), out double ac) && edge.TryGetValue((b, c), out double bc) && (dense || rng.NextDouble() < 0.7))
                {
                    triangle[(a, b, c)] = Later(Math.Max(ab, Math.Max(ac, bc)));
                    simplices.Add(new Simplex(triangle[(a, b, c)], a, b, c));
                }

        foreach (var ((a, b, c), abc) in triangle)
            for (int d = c + 1; d < n; d++)
                if (triangle.TryGetValue((a, b, d), out double abd)
                    && triangle.TryGetValue((a, c, d), out double acd)
                    && triangle.TryGetValue((b, c, d), out double bcd)
                    && rng.NextDouble() < (dense ? 0.85 : 0.6))
                {
                    simplices.Add(new Simplex(Later(Math.Max(abc, Math.Max(abd, Math.Max(acd, bcd)))), a, b, c, d));
                }

        return new SimplicialFiltration(simplices);
    }

    /// <summary>A random graph's flag complex stopped at triangles, with random edge values: the regime with the most essential bars.</summary>
    static SimplicialFiltration RandomTwoSkeleton(int n, Random rng)
    {
        var simplices = new List<Simplex>();
        for (int v = 0; v < n; v++) simplices.Add(new Simplex(0.0, v));

        var edge = new double?[n, n];
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
                if (rng.NextDouble() < 0.5)
                {
                    edge[a, b] = rng.NextDouble();
                    simplices.Add(new Simplex(edge[a, b]!.Value, a, b));
                }

        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
                if (edge[a, b] is double ab)
                    for (int c = b + 1; c < n; c++)
                        if (edge[a, c] is double ac && edge[b, c] is double bc)
                            simplices.Add(new Simplex(Math.Max(ab, Math.Max(ac, bc)), a, b, c));

        return new SimplicialFiltration(simplices);
    }
}
