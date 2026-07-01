using System;
using Graphs.Primitives;
using TDA.Mapper.Filters;
using Xunit;

namespace TDA.Mapper.Tests;

public sealed class GraphMapperTests
{
    [Fact]
    public void FiedlerVectorFilter_ConnectedPathGraph_MatchesHandComputedReference()
    {
        CsrGraph graph = BuildGraph(3,
            (0, 1, 1.0),
            (1, 2, 1.0));

        double[] actual = GraphFilters.FiedlerVector.Apply(graph);
        double invSqrt2 = 1.0 / Math.Sqrt(2.0);
        double[] expected = { invSqrt2, 0.0, -invSqrt2 };

        if (Dot(actual, expected) < 0.0)
        {
            for (int i = 0; i < actual.Length; i++)
                actual[i] = -actual[i];
        }

        for (int i = 0; i < expected.Length; i++)
            Assert.InRange(Math.Abs(actual[i] - expected[i]), 0.0, 1e-9);
    }

    [Fact]
    public void FiedlerVectorFilter_DisconnectedGraph_ThrowsConnectivityGuard()
    {
        CsrGraph graph = BuildGraph(4,
            (0, 1, 1.0),
            (2, 3, 1.0));

        InvalidOperationException ex = Assert.Throws<InvalidOperationException>(() => GraphFilters.FiedlerVector.Apply(graph));

        Assert.Contains("requires a connected graph", ex.Message, StringComparison.Ordinal);
    }

    private static double Dot(double[] left, double[] right)
    {
        double sum = 0.0;
        for (int i = 0; i < left.Length; i++)
            sum += left[i] * right[i];

        return sum;
    }

    private static CsrGraph BuildGraph(int nodeCount, params (int Source, int Target, double Weight)[] edges)
    {
        var graphEdges = new Edge[edges.Length];
        for (int i = 0; i < edges.Length; i++)
            graphEdges[i] = new Edge(edges[i].Source, edges[i].Target, edges[i].Weight);

        return CsrGraph.FromEdges(graphEdges, nodeCount);
    }
}
