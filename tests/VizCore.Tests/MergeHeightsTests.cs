using System;
using System.Linq;
using Graphs;
using Graphs.Coupling;
using Graphs.Observables;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The single-linkage merge heights in log space. The oracles are hand-built graphs and the Gaussian kernel's closed form
/// on the gap's own distance.
/// </summary>
public sealed class MergeHeightsTests
{
    [Fact]
    public void APath_MergesStrongestFirst_WithTheComponentSizes()
    {
        // 0-1 at -1, 1-2 at -5, 2-3 at -2: the pairs {0,1} and {2,3} form first, then the weak middle edge joins them.
        var g = CsrGraph.FromEdges(new[]
        {
            new Edge(0, 1, Math.Exp(-1.0), -1.0), new Edge(1, 2, Math.Exp(-5.0), -5.0), new Edge(2, 3, Math.Exp(-2.0), -2.0),
        }, 4);

        var merges = MergeHeights.Compute(g);

        Assert.Equal(new[] { -1.0, -2.0, -5.0 }, merges.Select(m => m.LogCoupling));
        Assert.Equal(new[] { (1, 1), (1, 1), (2, 2) }, merges.Select(m => (m.SizeA, m.SizeB)));
    }

    [Fact]
    public void ACycleEdge_JoinsNothing_SoTheCountIsNodesMinusOne()
    {
        var g = CsrGraph.FromEdges(new[]
        {
            new Edge(0, 1, 0.9, Math.Log(0.9)), new Edge(1, 2, 0.8, Math.Log(0.8)), new Edge(0, 2, 0.1, Math.Log(0.1)),
        }, 3);

        var merges = MergeHeights.Compute(g);

        Assert.Equal(2, merges.Count);
        Assert.DoesNotContain(merges, m => m.LogCoupling == Math.Log(0.1));
    }

    [Fact]
    public void WithoutALog_ItReadsTheLogOfTheWeights_AndAZeroWeightJoinsNothing()
    {
        var g = CsrGraph.FromEdges(new[] { new Edge(0, 1, 0.5), new Edge(1, 2, 0.0), new Edge(2, 3, 0.25) }, 4);

        var merges = MergeHeights.Compute(g);

        Assert.Equal(2, merges.Count);                       // a forest: {0,1} and {2,3}
        Assert.Equal(Math.Log(0.5), merges[0].LogCoupling, 12);
        Assert.Equal(Math.Log(0.25), merges[1].LogCoupling, 12);
    }

    [Fact]
    public void AnEdgeAtMinusInfinity_JoinsNothing()
    {
        var g = CsrGraph.FromEdges(new[] { new Edge(0, 1, 0.5, Math.Log(0.5)), new Edge(1, 2, 0.5, double.NegativeInfinity) }, 3);

        Assert.Single(MergeHeights.Compute(g));
    }

    [Fact]
    public void EqualHeights_AreOrderedByTheirEndpoints()
    {
        var g = CsrGraph.FromEdges(new[]
        {
            new Edge(2, 3, 0.5, Math.Log(0.5)), new Edge(0, 1, 0.5, Math.Log(0.5)), new Edge(1, 2, 0.5, Math.Log(0.5)),
        }, 4);

        var merges = MergeHeights.Compute(g);

        Assert.Equal(new[] { (0, 1), (1, 2), (2, 3) }, merges.Select(m => (m.Source, m.Target)));
    }

    [Fact]
    public void TwoGaussianBlobs_MergeAtMinusDSquaredOver2ASquared_PastTheDoubleRange()
    {
        // Two lines 10 apart; the Gaussian bridge's coupling underflows, and its height is the closed form on its own distance.
        var points = new double[20][];
        for (int line = 0; line < 2; line++)
            for (int i = 0; i < 10; i++)
                points[line * 10 + i] = new[] { line * 10.0 + 0.1 * i + 0.01 * i * i, 0.0 };
        Func<int, int, double> dist = (i, j) => Math.Abs(points[i][0] - points[j][0]);

        GraphBuildResult r = GraphCompiler.Build(
            new GraphCompilerConfig
            {
                Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 3 },
                Filter = new FilterConfig { Kind = FilterKind.OrRule },
                Repair = new RepairConfig { Kind = RepairKind.MstMin },
                Projection = new CouplingProjection { Kernel = new Gaussian(), LmpRescale = false },
            },
            points.Length, new GraphMetric(dist));

        var merges = MergeHeights.Compute(r.Graph);

        Assert.Equal(points.Length - 1, merges.Count);
        LinkageMerge last = merges[^1];
        Assert.Equal((10, 10), (last.SizeA, last.SizeB));
        double d = dist(last.Source, last.Target);
        double a = r.SingleBandwidth!.Value;
        double expected = -(d * d) / (2.0 * a * a);
        Assert.True(expected < -1000.0);
        Assert.InRange(last.LogCoupling / expected, 1 - 1e-12, 1 + 1e-12);
        Assert.All(merges.Take(merges.Count - 1), m => Assert.True(m.LogCoupling > last.LogCoupling));
    }
}
