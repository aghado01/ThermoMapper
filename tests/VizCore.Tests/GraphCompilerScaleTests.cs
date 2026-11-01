using System;
using System.Collections.Generic;
using System.Linq;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The scale stage: which edges the bandwidth reads (graph-compiler-brief §4.1, decision 3) and
/// whether K moves it (§2.3, decision 7).
/// </summary>
public sealed class GraphCompilerScaleTests
{
    private static readonly double[][] Points = MakePoints();

    private static double[][] MakePoints()
    {
        var points = new double[36][];
        for (int i = 0; i < 18; i++)
        {
            points[i] = new[] { 0.1 * i + 0.01 * i * i, 0.0 };
            points[18 + i] = new[] { 12.0 + 0.13 * i + 0.011 * i * i, 0.5 };
        }
        return points;
    }

    private static GraphMetric Metric() => new((i, j) =>
    {
        double dx = Points[i][0] - Points[j][0], dy = Points[i][1] - Points[j][1];
        return Math.Sqrt(dx * dx + dy * dy);
    });

    private static GraphCompilerConfig Config(
        int k, IEdgeProjection projection, RepairKind repair = RepairKind.MstMin, TopologyKind topology = TopologyKind.Knn) =>
        new()
        {
            Topology = new TopologyConfig { Kind = topology, K = k, Epsilon = topology == TopologyKind.EpsilonBall ? 1.0 : null },
            Filter = new FilterConfig { Kind = FilterKind.OrRule },
            Repair = new RepairConfig { Kind = repair },
            Projection = projection,
        };

    private static CouplingProjection Mean(ScaleEdgeSet edges = ScaleEdgeSet.BeforeRepair, int? rank = null) => new()
    {
        Kernel = new Gaussian(),
        LmpRescale = false,
        BandwidthOverride = BandwidthStrategy.MeanEdgeDistance,
        ScaleEdges = edges,
        ScaleRank = rank,
    };

    private static double MeanWeight(CsrGraph g)
    {
        double sum = 0.0;
        int count = 0;
        for (int i = 0; i < g.NodeCount; i++)
            for (int slot = g.RowPointers[i]; slot < g.RowPointers[i + 1]; slot++)
                if (g.Targets[slot] > i) { sum += g.Weights[slot]; count++; }
        return sum / count;
    }

    [Fact]
    public void MeanEdgeBandwidth_BeforeRepair_IgnoresTheBridges_AfterRepairAveragesThem()
    {
        GraphBuildResult distances = GraphCompiler.Build(Config(4, new DistanceProjection()), Points.Length, Metric());
        double afterExpected = MeanWeight(distances.Graph);
        double beforeExpected = MeanWeight(distances.PreRepairGraph!.Value);

        double before = GraphCompiler.Build(Config(4, Mean(ScaleEdgeSet.BeforeRepair)), Points.Length, Metric()).SingleBandwidth!.Value;
        double after = GraphCompiler.Build(Config(4, Mean(ScaleEdgeSet.AfterRepair)), Points.Length, Metric()).SingleBandwidth!.Value;

        Assert.InRange(before / beforeExpected, 1 - 1e-12, 1 + 1e-12);
        Assert.InRange(after / afterExpected, 1 - 1e-12, 1 + 1e-12);
        Assert.True(after > before, "The bridge, the longest edge, raises the mean it joins.");
    }

    [Fact]
    public void MeanEdgeBandwidth_BeforeRepair_DoesNotDependOnWhetherRepairRan()
    {
        double repaired = GraphCompiler.Build(Config(4, Mean(), RepairKind.MstMin), Points.Length, Metric()).SingleBandwidth!.Value;
        double unrepaired = GraphCompiler.Build(Config(4, Mean(), RepairKind.NoRepair), Points.Length, Metric()).SingleBandwidth!.Value;

        Assert.Equal(unrepaired, repaired, 12);
    }

    [Fact]
    public void ScaleRank_FixesTheBandwidthWhileKGrows()
    {
        double k4 = GraphCompiler.Build(Config(4, Mean(rank: 3)), Points.Length, Metric()).SingleBandwidth!.Value;
        double k8 = GraphCompiler.Build(Config(8, Mean(rank: 3)), Points.Length, Metric()).SingleBandwidth!.Value;
        double plain4 = GraphCompiler.Build(Config(4, Mean()), Points.Length, Metric()).SingleBandwidth!.Value;
        double plain8 = GraphCompiler.Build(Config(8, Mean()), Points.Length, Metric()).SingleBandwidth!.Value;

        Assert.Equal(k4, k8, 12);
        Assert.True(plain8 > plain4, "Without a rank the bandwidth widens with K.");
    }

    [Fact]
    public void ScaleRank_ForTheMadFamily_ReadsTheRthNeighbour()
    {
        CouplingProjection mad(int? rank) => new() { Kernel = new Gaussian(), LmpRescale = false, ScaleRank = rank };

        double k4 = GraphCompiler.Build(Config(4, mad(3)), Points.Length, Metric()).SingleBandwidth!.Value;
        double k8 = GraphCompiler.Build(Config(8, mad(3)), Points.Length, Metric()).SingleBandwidth!.Value;
        double plain4 = GraphCompiler.Build(Config(4, mad(null)), Points.Length, Metric()).SingleBandwidth!.Value;
        double plain8 = GraphCompiler.Build(Config(8, mad(null)), Points.Length, Metric()).SingleBandwidth!.Value;

        Assert.Equal(k4, k8, 12);
        Assert.NotEqual(plain4, plain8);
    }

    [Fact]
    public void ScaleRank_Outside1ToK_OrWithoutKnn_Throws()
    {
        Assert.Throws<ArgumentException>(() => GraphCompiler.Build(Config(4, Mean(rank: 5)), Points.Length, Metric()));
        Assert.Throws<ArgumentException>(() => GraphCompiler.Build(Config(4, Mean(rank: 0)), Points.Length, Metric()));
        Assert.Throws<ArgumentException>(() =>
            GraphCompiler.Build(Config(4, Mean(rank: 2), topology: TopologyKind.EpsilonBall), Points.Length, Metric()));
    }
}
