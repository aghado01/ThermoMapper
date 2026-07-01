using System.IO;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// T3: per-node thermodynamic landscapes on the Swendsen–Wang sampler. A cold,
/// connected graph drives every node into one giant FK cluster each draw, so the
/// un-reduced 0-forms (mean cluster size, giant participation) must accumulate
/// sanely and survive checkpoint + serialization.
/// </summary>
public sealed class SwLandscapeTests
{
    // 4-cycle: connected, so at low T (bond prob ≈ 1) every draw fuses all four
    // nodes into a single FK cluster. Standard tier on purpose — landscapes are
    // runtime-gated, independent of the per-edge observable tier.
    private static SwRunSpec ColdSpec(bool clusterSize, bool order) => new()
    {
        Graph = BuildGraph(4, (0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 0, 1.0)),
        Temperature = 0.1,
        Q = 2,
        Accumulation = new AccumulationSpec { ClusterSizeLandscape = clusterSize, OrderLandscape = order },
        Seed = 1234,
        Budget = new RunBudget(10, 50),
    };

    [Fact]
    public void Landscapes_AccumulateSanely_OnColdConnectedGraph()
    {
        Accumulator acc = SwRunner.Run(ColdSpec(clusterSize: true, order: true)).Accumulator;

        // Both sinks materialized (gated arrays carried through GetCheckpoint).
        Assert.NotNull(acc.SumClusterSizePerNode);
        Assert.NotNull(acc.SumInGiantClusterPerNode);
        Assert.Equal(4, acc.SumClusterSizePerNode!.Length);
        Assert.Equal(4, acc.SumInGiantClusterPerNode!.Length);
        Assert.Equal(50, acc.DrawCount);

        double[] meanSize = SwLandscapes.MeanClusterSize(acc);
        double[] giant = SwLandscapes.GiantParticipation(acc);

        for (int i = 0; i < 4; i++)
        {
            // Cold + connected ⇒ nodes co-cluster ⇒ mean FK cluster size > 1
            // (and never exceeds N = 4).
            Assert.InRange(meanSize[i], 1.0 + 1e-9, 4.0);
            // Order parameter is a per-draw indicator average ⇒ always in [0,1];
            // cold ⇒ the node sits in the giant on most draws.
            Assert.InRange(giant[i], 0.5, 1.0);
        }
    }

    [Fact]
    public void Landscapes_AreIndependentlyGated()
    {
        Accumulator sizeOnly = SwRunner.Run(ColdSpec(clusterSize: true, order: false)).Accumulator;
        Assert.NotNull(sizeOnly.SumClusterSizePerNode);
        Assert.Null(sizeOnly.SumInGiantClusterPerNode);

        Accumulator orderOnly = SwRunner.Run(ColdSpec(clusterSize: false, order: true)).Accumulator;
        Assert.Null(orderOnly.SumClusterSizePerNode);
        Assert.NotNull(orderOnly.SumInGiantClusterPerNode);

        Accumulator neither = SwRunner.Run(ColdSpec(clusterSize: false, order: false)).Accumulator;
        Assert.Null(neither.SumClusterSizePerNode);
        Assert.Null(neither.SumInGiantClusterPerNode);
    }

    [Fact]
    public void Landscapes_SurviveSerializationRoundTrip()
    {
        Accumulator acc = SwRunner.Run(ColdSpec(clusterSize: true, order: true)).Accumulator;

        using var ms = new MemoryStream();
        AccumulatorSerializer.Instance.WriteTo(acc, ms);
        ms.Position = 0;
        Accumulator round = AccumulatorSerializer.Instance.ReadFrom(ms);

        Assert.Equal(acc.SumClusterSizePerNode, round.SumClusterSizePerNode);
        Assert.Equal(acc.SumInGiantClusterPerNode, round.SumInGiantClusterPerNode);
    }

    [Fact]
    public void Serialization_PreservesNullLandscapes_WhenUntracked()
    {
        Accumulator acc = SwRunner.Run(ColdSpec(clusterSize: false, order: false)).Accumulator;

        using var ms = new MemoryStream();
        AccumulatorSerializer.Instance.WriteTo(acc, ms);
        ms.Position = 0;
        Accumulator round = AccumulatorSerializer.Instance.ReadFrom(ms);

        Assert.Null(round.SumClusterSizePerNode);
        Assert.Null(round.SumInGiantClusterPerNode);
    }

    private static CsrGraph BuildGraph(int nodeCount, params (int Source, int Target, double Weight)[] edges)
    {
        var graphEdges = new Edge[edges.Length];
        for (int i = 0; i < edges.Length; i++)
            graphEdges[i] = new Edge(edges[i].Source, edges[i].Target, edges[i].Weight);

        return CsrGraph.FromEdges(graphEdges, nodeCount);
    }
}
