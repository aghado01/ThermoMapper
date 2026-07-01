#nullable enable
using System.Collections.Generic;
using System.Linq;
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
}
