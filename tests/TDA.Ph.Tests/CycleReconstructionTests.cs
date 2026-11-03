#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Xunit;

using Maths.Topology;
namespace TDA.Ph.Tests;

public sealed class CycleReconstructionTests
{
    static void AssertH1CycleInvariants(Bar bar, SimplicialFiltration filtration)
    {
        Assert.NotNull(bar.Cycle);
        Assert.Contains(bar.Generator!.Value, bar.Cycle!);

        foreach (int idx in bar.Cycle!)
            Assert.True(filtration.Simplices[idx].FiltrationValue <= bar.Birth);

        var edges = BarCycleEdges.GetEdgePairs(bar, filtration);
        var degree = new Dictionary<int, int>();
        foreach (UndirectedEdge e in edges)
        {
            degree[e.Lo] = degree.GetValueOrDefault(e.Lo) + 1;
            degree[e.Hi] = degree.GetValueOrDefault(e.Hi) + 1;
        }

        Assert.All(degree.Values, d => Assert.Equal(2, d));
    }

    [Fact]
    public void Hexagon_InfiniteH1_FullSixEdgeCycle()
    {
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(0.0, 2),
            new Simplex(0.0, 3),
            new Simplex(0.0, 4),
            new Simplex(0.0, 5),
            new Simplex(1.0, 0, 1),
            new Simplex(1.0, 1, 2),
            new Simplex(1.0, 2, 3),
            new Simplex(1.0, 3, 4),
            new Simplex(1.0, 4, 5),
            new Simplex(1.0, 0, 5),
        }, "t");

        Barcode barcode = PersistentInvolutedHomology.Compute(filtration, representatives: true);
        var h1 = barcode.Bars.Single(b => b.Dimension == 1);

        var edges = BarCycleEdges.GetEdgePairs(h1, filtration);
        Assert.Equal(6, edges.Count);
        AssertH1CycleInvariants(h1, filtration);
    }

    [Fact]
    public void ShortestPath_PrefersLowWeightChord()
    {
        // Long boundary 0-1-2-3-4 (weight 1 each); chord 0-2 and 2-4 (weight 0.5); birth (0,4) at 2.
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(0.0, 2),
            new Simplex(0.0, 3),
            new Simplex(0.0, 4),
            new Simplex(1.0, 0, 1),
            new Simplex(1.0, 1, 2),
            new Simplex(1.0, 2, 3),
            new Simplex(1.0, 3, 4),
            new Simplex(0.5, 0, 2),
            new Simplex(0.5, 2, 4),
            new Simplex(2.0, 0, 4),
        }, "t");

        int birthIdx = filtration.Simplices.ToList().FindIndex(s =>
            s.Dimension == 1 && s.Vertices[0] == 0 && s.Vertices[1] == 4);

        int[] cycle = CycleReconstruction.ReconstructH1Cycle(filtration, birthIdx);
        var edgePairs = cycle
            .Where(i => filtration.Simplices[i].Dimension == 1)
            .Select(i => new UndirectedEdge(
                filtration.Simplices[i].Vertices[0],
                filtration.Simplices[i].Vertices[1]))
            .ToList();

        Assert.Contains(new UndirectedEdge(0, 2), edgePairs);
        Assert.Contains(new UndirectedEdge(2, 4), edgePairs);
        Assert.DoesNotContain(new UndirectedEdge(0, 1), edgePairs);
        Assert.DoesNotContain(new UndirectedEdge(1, 2), edgePairs);
    }

    [Fact]
    public void NegativeValues_ReadByEdgeCount_AndFinish()
    {
        // The long way round is three edges with a far smaller sum of values; the short way is two edges. The birth edge is (0,4).
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(-30.0, 0),
            new Simplex(-30.0, 1),
            new Simplex(-30.0, 2),
            new Simplex(-30.0, 3),
            new Simplex(-30.0, 4),
            new Simplex(-20.0, 0, 1),
            new Simplex(-20.0, 1, 2),
            new Simplex(-20.0, 2, 4),
            new Simplex(-1.0, 0, 3),
            new Simplex(-1.0, 3, 4),
            new Simplex(0.5, 0, 4),
        }, "t");

        Bar loop = PersistentInvolutedHomology.Compute(filtration, 1, representatives: true).Bars
            .Single(b => b.Dimension == 1 && b.Birth == 0.5);

        Assert.True(loop.IsInfinite);
        Assert.Equal(3, loop.Cycle!.Length);
        AssertH1CycleInvariants(loop, filtration);
        Assert.Contains(new UndirectedEdge(0, 3), BarCycleEdges.GetEdgePairs(loop, filtration));
    }

    [Fact]
    public void TiedEdgeAfterTheBirthEdge_IsNotOnThePath_SinceATiedTriangleCanBoundThem()
    {
        // (0,1) is born at 2 and essential; (0,2) ties with it and a tied triangle kills it at once. The path through (0,2) is
        // the shorter, and with (0,1) it is the triangle's boundary, a class that is zero.
        var filtration = new SimplicialFiltration(new[]
        {
            new Simplex(0.0, 0),
            new Simplex(0.0, 1),
            new Simplex(0.0, 2),
            new Simplex(0.0, 3),
            new Simplex(0.5, 1, 2),
            new Simplex(1.5, 0, 3),
            new Simplex(1.5, 1, 3),
            new Simplex(2.0, 0, 1),
            new Simplex(2.0, 0, 2),
            new Simplex(2.0, 0, 1, 2),
        }, "t");

        Bar essential = PersistentInvolutedHomology.Compute(filtration, 1, representatives: true).Bars
            .Single(b => b.Dimension == 1 && b.IsInfinite);

        Assert.Equal(new[] { new UndirectedEdge(0, 1), new UndirectedEdge(0, 3), new UndirectedEdge(1, 3) },
            BarCycleEdges.GetEdgePairs(essential, filtration).OrderBy(e => e.Lo).ThenBy(e => e.Hi));
    }

    [Theory]
    [InlineData(5e307, 1.0)]    // lengths whose sum overflows
    [InlineData(1.0, -1e308)]   // every value far below zero
    [InlineData(1.0, -20.0)]    // a shift below zero
    [InlineData(1.0, -1.5)]     // straddling zero
    public void RandomValues_EssentialCyclesAreTheBarsOwn_AcrossTheRange(double scale, double shift)
    {
        const int trials = 150;
        int essential = 0, tiedAfter = 0;
        for (int seed = 0; seed < trials; seed++)
        {
            var rng = new Random(seed);
            SimplicialFiltration filtration = RandomValued(rng, scale, shift, grid: seed % 2 == 0);
            Barcode barcode = PersistentInvolutedHomology.Compute(filtration, 1, representatives: true);

            foreach (Bar bar in barcode.Bars.Where(b => b.Dimension == 1 && b.IsInfinite))
            {
                essential++;
                int generator = bar.Generator!.Value;
                AssertH1CycleInvariants(bar, filtration);
                Assert.Equal(generator, bar.Cycle!.Max());
                if (Enumerable.Range(generator + 1, filtration.Count - generator - 1).Any(c => filtration.GetDimension(c) == 1 && filtration.GetBirth(c) == bar.Birth))
                    tiedAfter++;

                // The class never bounds, and nothing older can account for it: not a boundary of the whole complex plus cycles
                // from before the birth edge.
                var older = new ChainSpan();
                for (int c = 0; c < filtration.Count; c++)
                    if (filtration.GetDimension(c) == 2)
                        older.Add(ZigzagModuleOracle.ChainOf(filtration.GetBoundaryIndices(c)));
                foreach (BigInteger cycle in CyclesBefore(filtration, generator))
                    older.Add(cycle);

                Assert.False(older.Contains(ZigzagModuleOracle.ChainOf(bar.Cycle!)), $"seed {seed}: the cycle's class is an older one");
            }
        }

        Assert.True(tiedAfter > 0, "no essential bar had an edge tied with its birth edge and after it");
        Assert.True(essential > 0);
    }

    // The cycle space of the edges before one in the filtration, as chains over simplex indices.
    static IEnumerable<BigInteger> CyclesBefore(SimplicialFiltration filtration, int edge)
    {
        var pivots = new Dictionary<int, (BigInteger Boundary, BigInteger Chain)>();
        var cycles = new List<BigInteger>();
        for (int i = 0; i < edge; i++)
        {
            if (filtration.GetDimension(i) != 1) continue;

            BigInteger b = ZigzagModuleOracle.ChainOf(filtration.GetBoundaryIndices(i));
            BigInteger chain = BigInteger.One << i;
            while (!b.IsZero && pivots.TryGetValue((int)b.GetBitLength() - 1, out var pivot))
            {
                b ^= pivot.Boundary;
                chain ^= pivot.Chain;
            }

            if (b.IsZero) cycles.Add(chain);
            else pivots[(int)b.GetBitLength() - 1] = (b, chain);
        }

        return cycles;
    }

    // A graph with some triangles, values on a grid half the time so edges tie, scaled and shifted across the range.
    static SimplicialFiltration RandomValued(Random rng, double scale, double shift, bool grid)
    {
        int n = rng.Next(4, 9);
        double Value(double low) => shift + scale * (low + (grid ? rng.Next(0, 3) : rng.NextDouble() * 3.0));

        var simplices = new List<Simplex>();
        for (int v = 0; v < n; v++) simplices.Add(new Simplex(shift, v));

        var edges = new Dictionary<(int, int), double>();
        for (int a = 0; a < n; a++)
            for (int b = a + 1; b < n; b++)
                if (rng.NextDouble() < 0.6)
                {
                    double value = Value(0.0);
                    edges[(a, b)] = value;
                    simplices.Add(new Simplex(value, a, b));
                }

        foreach (var ((a, b), ab) in edges)
            for (int c = b + 1; c < n; c++)
                if (edges.TryGetValue((a, c), out double ac) && edges.TryGetValue((b, c), out double bc) && rng.NextDouble() < 0.15)
                    simplices.Add(new Simplex(Math.Max(ab, Math.Max(ac, bc)), a, b, c));

        return new SimplicialFiltration(simplices);
    }
}
