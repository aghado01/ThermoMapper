using System;
using System.Collections.Generic;
using Graphs;
using Graphs.Primitives;
using Graphs.Primitives.Traversal;
using Graphs.Pipeline.Refinement;
using Xunit;

namespace VizCore.Tests;

public sealed class GraphPathNeighborTests
{
    [Fact]
    public void PathNeighborRefiner_RefinesEdgeDistancesUsingDijkstra()
    {
        var neighbors = new Neighbor[4][]
        {
            new[] { new Neighbor { Index = 1, Distance = 5.0 }, new Neighbor { Index = 2, Distance = 10.0 } },
            new[] { new Neighbor { Index = 2, Distance = 3.0 } },
            new[] { new Neighbor { Index = 3, Distance = 2.0 } },
            Array.Empty<Neighbor>()
        };

        var selection = new NeighborSelection(
            neighbors,
            new double[] { 5.0, 3.0, 2.0, double.PositiveInfinity },
            new double[] { 10.0, 3.0, 2.0, double.PositiveInfinity });

        var refined = new PathNeighborRefiner().Refine(selection, n: 4);

        Assert.Equal(5.0, refined.AllNeighbors[0][0].Distance);
        Assert.Equal(8.0, refined.AllNeighbors[0][1].Distance);
        Assert.Equal(3.0, refined.AllNeighbors[1][0].Distance);
        Assert.Equal(2.0, refined.AllNeighbors[2][0].Distance);
        Assert.Empty(refined.AllNeighbors[3]);

        Assert.Equal(5.0, refined.NearestNeighborDistances[0]);
        Assert.Equal(3.0, refined.NearestNeighborDistances[1]);
        Assert.Equal(2.0, refined.NearestNeighborDistances[2]);
        Assert.Equal(double.PositiveInfinity, refined.NearestNeighborDistances[3]);
    }

    [Fact]
    public void Dijkstra_ComputeBoundedDistances_MaskedTargetCorrectlySettlesTarget()
    {
        CsrGraph graph = CsrGraph.FromEdges(
            new[]
            {
                new Edge(0, 1, 1.0),
                new Edge(1, 2, 2.0),
                new Edge(2, 3, 3.0)
            },
            nodeCount: 5);

        var distances = new double[5];
        var hops = new int[5];
        var queue = new PriorityQueue<int, double>();
        var targetMask = new bool[5];
        targetMask[3] = true;

        Dijkstra.ComputeBoundedDistances(
            graph,
            sourceNode: 0,
            distances,
            hops,
            queue,
            targetMask,
            maxDistance: double.PositiveInfinity);

        Assert.Equal(6.0, distances[3]);
        Assert.Equal(1.0, distances[1]);
        Assert.Equal(3.0, distances[2]);
        Assert.Equal(double.PositiveInfinity, distances[4]);
    }
}
