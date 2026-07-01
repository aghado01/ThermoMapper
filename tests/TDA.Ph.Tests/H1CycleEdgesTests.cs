#nullable enable
using System.Linq;
using Graphs.Primitives;
using TDA.Primitives;
using Xunit;

using Maths.Topology;
namespace TDA.Ph.Tests;

public sealed class H1CycleEdgesTests
{
    [Fact]
    public void TwoLoopWedge_DistanceGraph_YieldsLoopClosingEdges()
    {
        CsrGraph distanceGraph = TwoLoopWedgeDistanceGraph();
        var edges = H1CycleEdges.FromDistanceGraph(distanceGraph);

        Assert.Contains((0, 2), edges);
        Assert.Contains((0, 4), edges);
        Assert.True(edges.Count >= 4);
    }

    static CsrGraph TwoLoopWedgeDistanceGraph()
    {
        Edge[] edges = new[]
        {
            new Edge(0, 1, 1.0),
            new Edge(1, 2, 1.0),
            new Edge(0, 2, 1.0),
            new Edge(0, 3, 1.0),
            new Edge(3, 4, 1.0),
            new Edge(0, 4, 1.0),
        };
        return CsrGraph.FromEdges(edges, 5);
    }
}
