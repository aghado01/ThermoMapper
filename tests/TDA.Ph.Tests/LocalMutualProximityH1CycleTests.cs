#nullable enable
using Graphs.Primitives;
using Graphs.Proximity;
using Xunit;

using Maths.Topology;
namespace TDA.Ph.Tests;

public sealed class LocalMutualProximityH1CycleTests
{
    [Fact]
    public void ProtectedH1Edges_PreserveInnerCoupling()
    {
        CsrGraph couplingGraph = SquareWithDanglingTailCouplings();
        CsrGraph distanceGraph = SquareWithDanglingTailDistances();

        var protectedEdges = H1CycleEdges.FromDistanceGraph(distanceGraph);
        Assert.Contains((0, 1), protectedEdges);

        CsrGraph demoted = LocalMutualProximity.ApplyLocalScaling(
            couplingGraph, weightsAreCouplings: true);
        CsrGraph preserved = LocalMutualProximity.ApplyLocalScaling(
            couplingGraph, weightsAreCouplings: true, protectedEdges: protectedEdges);

        Assert.Equal(0.0, GetUndirectedWeight(demoted, 3, 4));
        Assert.Equal(0.6, GetUndirectedWeight(preserved, 0, 1));
        Assert.True(GetUndirectedWeight(preserved, 3, 4) < GetUndirectedWeight(couplingGraph, 3, 4));
    }

    static CsrGraph SquareWithDanglingTailCouplings()
    {
        Edge[] edges = new[]
        {
            new Edge(0, 1, 0.6),
            new Edge(1, 2, 0.6),
            new Edge(2, 3, 0.6),
            new Edge(0, 3, 0.6),
            new Edge(3, 4, 0.6),
        };
        return CsrGraph.FromEdges(edges, 5);
    }

    static CsrGraph SquareWithDanglingTailDistances()
    {
        Edge[] edges = new[]
        {
            new Edge(0, 1, 1.0),
            new Edge(1, 2, 1.0),
            new Edge(2, 3, 1.0),
            new Edge(0, 3, 1.0),
            new Edge(3, 4, 1.0),
        };
        return CsrGraph.FromEdges(edges, 5);
    }

    static double GetUndirectedWeight(CsrGraph graph, int u, int v)
    {
        int lo = u < v ? u : v;
        int hi = u < v ? v : u;

        for (int e = graph.RowPointers[lo]; e < graph.RowPointers[lo + 1]; e++)
        {
            if (graph.Targets[e] == hi)
                return graph.Weights[e];
        }

        throw new Xunit.Sdk.XunitException($"Edge ({lo}, {hi}) not found in graph.");
    }
}
