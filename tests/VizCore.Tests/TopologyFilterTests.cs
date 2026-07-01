using System;
using Graphs;
using Graphs.Pipeline.Filters;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

public sealed class TopologyFilterTests
{
    [Fact]
    public void MutualKnnFilter_IsolatedNodes_ProduceLengthNArrays()
    {
        const int n = 5;
        // Directed 1-NN chain: 0→1, 1→2, 2→3, 3→4, 4→0 — no mutual pairs for 3,4
        var directed = new Neighbor[n][];
        var kth = new double[n];
        var nn = new double[n];
        for (int i = 0; i < n; i++)
        {
            directed[i] = new[] { new Neighbor { Index = (i + 1) % n, Distance = 1.0 } };
            kth[i] = 1.0;
            nn[i] = 1.0;
        }

        var directedSel = new NeighborSelection(directed, nn, kth);
        var filter = new MutualKnnFilter(MutualBandwidthSource.DirectedKth);
        NeighborSelection result = filter.Filter(directedSel, n, static (_, _) => 1.0);

        Assert.Equal(n, result.NearestNeighborDistances.Length);
        Assert.Equal(n, result.KthNeighborDistances.Length);
        Assert.Equal(double.PositiveInfinity, result.NearestNeighborDistances[3]);
        Assert.Equal(double.PositiveInfinity, result.NearestNeighborDistances[4]);
        Assert.Empty(result.AllNeighbors[3]);
        Assert.Empty(result.AllNeighbors[4]);
    }

    [Fact]
    public void PassThroughFilter_EmptyRow_ProducesLengthNNearestNeighborDistances()
    {
        const int n = 3;
        var directed = new Neighbor[n][]
        {
            new[] { new Neighbor { Index = 1, Distance = 1.0 } },
            new[] { new Neighbor { Index = 0, Distance = 1.0 } },
            Array.Empty<Neighbor>(),
        };
        var kth = new double[] { 1.0, 1.0, 0.5 };
        var nn = new double[] { 1.0, 1.0, double.PositiveInfinity };

        var directedSel = new NeighborSelection(directed, nn, kth);
        var filter = new PassThroughFilter();
        NeighborSelection result = filter.Filter(directedSel, n, static (_, _) => 1.0);

        Assert.Equal(n, result.NearestNeighborDistances.Length);
        Assert.Equal(double.PositiveInfinity, result.NearestNeighborDistances[2]);
    }
}
