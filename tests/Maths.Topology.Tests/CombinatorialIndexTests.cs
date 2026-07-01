#nullable enable
using Xunit;

namespace Maths.Topology.Tests;

public sealed class CombinatorialIndexTests
{
    [Fact]
    public void IndexAndVertices_RoundTrip()
    {
        int[] verts = { 1, 4, 7 };
        int idx = CombinatorialIndex.Index(verts);
        int[] round = CombinatorialIndex.Vertices(idx, dimension: 2);
        Assert.Equal(verts, round);
    }

    [Fact]
    public void PackKey_DistinctDimensions()
    {
        int[] edge = { 0, 2 };
        long k0 = CombinatorialIndex.PackKey(0, edge[..1]);
        long k1 = CombinatorialIndex.PackKey(1, edge);
        Assert.NotEqual(k0, k1);
    }
}
