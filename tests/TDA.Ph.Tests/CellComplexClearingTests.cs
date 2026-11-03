#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Maths.Topology;
using Xunit;
using Xunit.Abstractions;

namespace TDA.Ph.Tests;

/// <summary>
/// <see cref="PersistentCohomology"/> on complexes that are not simplicial, against <see cref="PersistentHomology"/>: cells
/// with an explicit boundary, 1-cells with any number of faces, cells with none, relative quotients.
/// </summary>
public sealed class CellComplexClearingTests
{
    readonly ITestOutputHelper _output;

    public CellComplexClearingTests(ITestOutputHelper output) => _output = output;

    // ── by hand ─────────────────────────────────────────────────────────────────

    [Fact]
    public void Triangle_ReadsAsItsSimplicialCopy()
    {
        var cells = new CellFiltration(new (int, double, int[])[]
        {
            (0, 0, Array.Empty<int>()), (0, 1, Array.Empty<int>()), (0, 2, Array.Empty<int>()),
            (1, 3, new[] { 0, 1 }), (1, 4, new[] { 1, 2 }), (1, 5, new[] { 0, 2 }),
            (2, 6, new[] { 3, 4, 5 }),
        });

        var expected = new[] { (0, 0.0, double.PositiveInfinity), (0, 1.0, 3.0), (0, 2.0, 4.0), (1, 5.0, 6.0) };

        Assert.Equal(expected, Canon(PersistentHomology.Compute(cells)));
        Assert.Equal(expected, Canon(PersistentCohomology.Compute(cells)));
    }

    [Fact]
    public void Torus_OneCellPerDegree_EveryCellIsEssential()
    {
        // Over Z/2 every boundary of the minimal torus vanishes, so nothing here has a face or a coface to find it through.
        var cells = new CellFiltration(new (int, double, int[])[]
        {
            (0, 0, Array.Empty<int>()), (1, 1, Array.Empty<int>()), (1, 2, Array.Empty<int>()), (2, 3, Array.Empty<int>()),
        });

        var expected = new[]
        {
            (0, 0.0, double.PositiveInfinity), (1, 1.0, double.PositiveInfinity),
            (1, 2.0, double.PositiveInfinity), (2, 3.0, double.PositiveInfinity),
        };

        Assert.Equal(expected, Canon(PersistentHomology.Compute(cells)));
        Assert.Equal(expected, Canon(PersistentCohomology.Compute(cells)));
    }

    [Fact]
    public void ThreeSphere_AnEmptyDegreeBelowACell_StillReachesIt()
    {
        var cells = new CellFiltration(new (int, double, int[])[]
        {
            (0, 0, Array.Empty<int>()), (3, 1, Array.Empty<int>()),
        });

        var expected = new[] { (0, 0.0, double.PositiveInfinity), (3, 1.0, double.PositiveInfinity) };

        Assert.Equal(expected, Canon(PersistentHomology.Compute(cells)));
        Assert.Equal(expected, Canon(PersistentCohomology.Compute(cells)));
    }

    [Fact]
    public void OneFace_JoinsItsVertexToTheQuotient_WhichOutlivesTheElder()
    {
        // A relative edge: vertex 0 dies at 2 though it is the elder of everything that survives.
        var cells = new CellFiltration(new (int, double, int[])[]
        {
            (0, 0, Array.Empty<int>()), (0, 1, Array.Empty<int>()), (1, 2, new[] { 0 }), (1, 3, new[] { 0, 1 }),
        });

        var expected = new[] { (0, 0.0, 2.0), (0, 1.0, 3.0) };

        Assert.Equal(expected, Canon(PersistentHomology.Compute(cells)));
        Assert.Equal(expected, Canon(PersistentCohomology.Compute(cells)));
    }

    [Fact]
    public void ThreeFaces_AOneCellNoUnionFindReads_IsReducedByColumns()
    {
        var cells = new CellFiltration(new (int, double, int[])[]
        {
            (0, 0, Array.Empty<int>()), (0, 1, Array.Empty<int>()), (0, 2, Array.Empty<int>()), (0, 3, Array.Empty<int>()),
            (1, 4, new[] { 0, 1 }), (1, 5, new[] { 1, 2, 3 }), (1, 6, new[] { 0, 1, 2, 3 }), (1, 7, new[] { 0, 2 }), (1, 8, new[] { 1, 2 }),
        });

        // Column 6 reduces to vertex 0 alone, so the eldest dies, and column 8 reduces to nothing.
        var expected = new[] { (0, 0.0, 6.0), (0, 1.0, 4.0), (0, 2.0, 7.0), (0, 3.0, 5.0), (1, 8.0, double.PositiveInfinity) };

        Assert.Equal(expected, Canon(PersistentHomology.Compute(cells)));
        Assert.Equal(expected, Canon(PersistentCohomology.Compute(cells)));
    }

    // ── random chain complexes ──────────────────────────────────────────────────

    [Fact]
    public void RandomChainComplexes_EqualTheColumnReducer_InEveryDegree()
    {
        const int trials = 400;
        int cellsWithoutFaces = 0, oneFace = 0, manyFaces = 0;
        var finite = new int[4];
        var essential = new int[4];

        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            // Half the trials keep every 1-cell to a graph's shape, which the union-find reads, and half do not.
            var cells = RandomChainComplex(rng, graphLike: seed % 2 == 0);
            var filtration = new CellFiltration(cells);

            foreach (int maxDimension in new[] { int.MaxValue, 0, 1, 2 })
            {
                var reference = Canon(PersistentHomology.Compute(filtration, maxDimension));
                Assert.True(reference.SequenceEqual(Canon(PersistentCohomology.Compute(filtration, maxDimension))),
                    $"seed {seed}, maxDimension {maxDimension}");
            }

            var pairs = PersistentCohomology.ComputePairs(filtration, int.MaxValue);
            Assert.Equal(filtration.Count, pairs.Sum(p => p.Death < 0 ? 1 : 2));

            foreach (var (dimension, _, death) in Canon(PersistentHomology.Compute(filtration)))
            {
                if (double.IsPositiveInfinity(death)) essential[dimension]++;
                else finite[dimension]++;
            }

            foreach (var (dimension, _, boundary) in cells)
            {
                if (dimension == 0) continue;
                if (boundary.Length == 0) cellsWithoutFaces++;
                else if (dimension == 1 && boundary.Length == 1) oneFace++;
                else if (dimension == 1 && boundary.Length > 2) manyFaces++;
            }
        }

        _output.WriteLine($"cells above degree 0 with no faces: {cellsWithoutFaces}; 1-cells with one face {oneFace}, with over two {manyFaces}; " +
            $"finite H0-H3: {string.Join("/", finite)}; essential H0-H3: {string.Join("/", essential)}");
        Assert.True(cellsWithoutFaces > 0 && oneFace > 0 && manyFaces > 0);
        Assert.True(finite[0] > 0 && finite[1] > 0 && finite[2] > 0);
        Assert.True(essential[1] > 0 && essential[2] > 0 && essential[3] > 0);
    }

    // ── the relative quotient ───────────────────────────────────────────────────

    [Fact]
    public void RelativeQuotient_EqualsTheColumnReducer_InEveryDegree()
    {
        const int trials = 300;
        int emptyBoundaries = 0, bars = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(4, 9);
            var simplices = new List<Simplex>();
            var edges = new List<(int, int)>();
            for (int v = 0; v < n; v++) simplices.Add(new Simplex(rng.NextDouble(), v));
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    if (rng.NextDouble() < 0.6)
                    {
                        edges.Add((a, b));
                        simplices.Add(new Simplex(1.0 + rng.NextDouble(), a, b));
                    }

            var edgeSet = edges.ToHashSet();
            foreach (var (a, b) in edges)
                for (int c = b + 1; c < n; c++)
                    if (edgeSet.Contains((a, c)) && edgeSet.Contains((b, c)) && rng.NextDouble() < 0.7)
                        simplices.Add(new Simplex(2.0 + rng.NextDouble(), a, b, c));

            var ambient = new SimplicialFiltration(simplices);
            int k0Vertices = rng.Next(1, n);
            double keep = rng.NextDouble();
            var inK0 = new Dictionary<Simplex, bool>();
            foreach (Simplex s in ambient.Simplices.OrderBy(s => s.Dimension))
                inK0[s] = s.Dimension == 0
                    ? s.Vertices[0] < k0Vertices
                    : s.Vertices.All(v => v < k0Vertices) && rng.NextDouble() < keep
                        && ambient.Simplices.Where(f => f.Dimension == s.Dimension - 1 && f.Vertices.All(s.Vertices.Contains)).All(f => inK0[f]);

            var relative = new RelativeFiltration(ambient, s => inK0[s]);
            emptyBoundaries += Enumerable.Range(0, relative.Count).Count(c => relative.GetDimension(c) > 0 && relative.GetBoundaryIndices(c).Length == 0);

            foreach (int maxDimension in new[] { int.MaxValue, 0, 1 })
            {
                var reference = Canon(PersistentHomology.Compute(relative, maxDimension));
                Assert.True(reference.SequenceEqual(Canon(PersistentCohomology.Compute(relative, maxDimension))),
                    $"seed {seed}, maxDimension {maxDimension}");
                bars += reference.Length;
            }
        }

        _output.WriteLine($"relative cells with an empty boundary: {emptyBoundaries}; bars compared: {bars}");
        Assert.True(emptyBoundaries > 0 && bars > 0);
    }

    // ── helpers ─────────────────────────────────────────────────────────────────

    // A filtered Z/2 chain complex by construction: each cell's boundary is a sum of cycles one degree down among the cells
    // before it, so ∂∂ = 0, and nothing makes it a CW complex. Births are the cell indices, or a coarse grid that ties.
    static List<(int Dimension, double Birth, int[] Boundary)> RandomChainComplex(Random rng, bool graphLike)
    {
        const int top = 3;
        int count = rng.Next(8, 28);
        bool ties = rng.Next(2) == 0;
        var cells = new List<(int, double, int[])>();
        var cycles = new List<BigInteger>[top + 1];
        var pivots = new Dictionary<int, (BigInteger Boundary, BigInteger Chain)>[top + 1];
        for (int d = 0; d <= top; d++)
        {
            cycles[d] = new List<BigInteger>();
            pivots[d] = new Dictionary<int, (BigInteger, BigInteger)>();
        }

        double value = 0.0;
        for (int i = 0; i < count; i++)
        {
            int dimension = i < 3 ? 0 : rng.Next(0, top + 1);
            BigInteger boundary = BigInteger.Zero;
            if (dimension == 1 && graphLike)
            {
                // Two vertices, one or none.
                foreach (BigInteger vertex in cycles[0].OrderBy(_ => rng.Next()).Take(rng.Next(0, 3)))
                    boundary ^= vertex;
            }
            else if (dimension > 0)
            {
                foreach (BigInteger cycle in cycles[dimension - 1])
                    if (rng.Next(2) == 0)
                        boundary ^= cycle;
            }

            int[] faces = Enumerable.Range(0, i).Where(c => !((boundary >> c) & BigInteger.One).IsZero).ToArray();
            value += ties ? rng.Next(0, 2) : 1.0;
            cells.Add((dimension, value, faces));

            // This cell's own cycle, if its boundary is one of the boundaries before it.
            BigInteger b = boundary;
            BigInteger chain = BigInteger.One << i;
            Dictionary<int, (BigInteger Boundary, BigInteger Chain)> own = pivots[dimension];
            while (!b.IsZero && own.TryGetValue((int)b.GetBitLength() - 1, out var pivot))
            {
                b ^= pivot.Boundary;
                chain ^= pivot.Chain;
            }

            if (b.IsZero) cycles[dimension].Add(chain);
            else own[(int)b.GetBitLength() - 1] = (b, chain);
        }

        return cells;
    }

    static (int, double, double)[] Canon(Barcode barcode) =>
        barcode.Bars
            .Where(b => b.IsInfinite || b.Death > b.Birth)
            .Select(b => (b.Dimension, b.Birth, b.Death))
            .OrderBy(x => x)
            .ToArray();
}
