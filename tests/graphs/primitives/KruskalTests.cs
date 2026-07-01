using System;
using Graphs.Primitives.Mst;
using Xunit;

namespace Graphs.Primitives.Tests;

public sealed class KruskalTests
{
    [Fact]
    public void BuildMinimumSpanningTree_ReturnsExpectedEdges()
    {
        var sortedEdges = new MstEdge[]
        {
            new MstEdge(0, 1, 1.0),
            new MstEdge(1, 2, 2.0),
            new MstEdge(2, 3, 3.0),
            new MstEdge(0, 2, 5.0),
            new MstEdge(0, 3, 10.0),
        };

        var output = new MstEdge[3];
        int written = Kruskal.BuildMinimumSpanningTree(sortedEdges, 4, output);

        Assert.Equal(3, written);
        Assert.Collection(output,
            edge => Assert.Equal(new MstEdge(0, 1, 1.0), edge),
            edge => Assert.Equal(new MstEdge(1, 2, 2.0), edge),
            edge => Assert.Equal(new MstEdge(2, 3, 3.0), edge));
    }

    [Fact]
    public void BuildMinimumSpanningTree_SkipsCycleEdges()
    {
        var sortedEdges = new MstEdge[]
        {
            new MstEdge(0, 1, 1.0),
            new MstEdge(1, 2, 2.0),
            new MstEdge(0, 2, 2.5),
            new MstEdge(2, 3, 3.0),
        };

        var output = new MstEdge[3];
        int written = Kruskal.BuildMinimumSpanningTree(sortedEdges, 4, output);

        Assert.Equal(3, written);
        var excluded = new MstEdge(0, 2, 2.5);
        for (int i = 0; i < written; i++)
        {
            Assert.False(output[i].Equals(excluded), $"Unexpected cycle edge found at index {i}.");
        }
    }
}
