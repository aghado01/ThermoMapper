using System;
using System.Collections.Generic;
using Graphs.Primitives;
using TDA.Mapper.Diagnostics;
using Xunit;

namespace TDA.Mapper.Tests;

public sealed class MapperDiagnosticsTests
{
    [Fact]
    public void NerveTopology_From_ReportsTreeCycleTriangleAndDisconnectedSignals()
    {
        var tree = CreateResult(
            nodes: CreateNodes((new[] { 0, 1 }, 0.0), (new[] { 1, 2 }, 1.0), (new[] { 2, 3 }, 2.0)),
            nerve: BuildGraph(3, (0, 1, 1.0), (1, 2, 1.0)));

        var treeReport = NerveTopology.From(tree);
        Assert.True(treeReport.IsTreeLike);
        Assert.Equal(0, treeReport.LoopCount);
        Assert.Equal(-1, treeReport.Girth);

        var cycle = CreateResult(
            nodes: CreateNodes((new[] { 0, 1 }, 0.0), (new[] { 1, 2 }, 1.0), (new[] { 2, 3 }, 2.0), (new[] { 3, 0 }, 3.0)),
            nerve: BuildGraph(4, (0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 0, 1.0)));

        var cycleReport = NerveTopology.From(cycle);
        Assert.False(cycleReport.IsTreeLike);
        Assert.Equal(1, cycleReport.LoopCount);
        Assert.Equal(4, cycleReport.Girth);

        var triangle = CreateResult(
            nodes: CreateNodes((new[] { 0, 1 }, 0.0), (new[] { 1, 2 }, 1.0), (new[] { 0, 2 }, 2.0)),
            nerve: BuildGraph(3, (0, 1, 1.0), (1, 2, 1.0), (0, 2, 1.0)));

        var triangleReport = NerveTopology.From(triangle);
        Assert.True(triangleReport.TriangleCount >= 1);
        Assert.Equal(3, triangleReport.Girth);

        var disconnected = CreateResult(
            nodes: CreateNodes((new[] { 0 }, 0.0), (new[] { 1 }, 1.0), (new[] { 2 }, 2.0), (new[] { 3 }, 3.0)),
            nerve: BuildGraph(4, (0, 1, 1.0), (2, 3, 1.0)));

        var disconnectedReport = NerveTopology.From(disconnected);
        Assert.True(disconnectedReport.ConnectedComponents > 1);
        Assert.True(disconnectedReport.LargestComponentSize < disconnectedReport.NerveNodeCount);
    }

    [Fact]
    public void NodeStats_From_ReportsSingleUniformMixedAndEmptyShapes()
    {
        var single = CreateResult(
            nodes: CreateNodes((new[] { 0, 1, 2 }, 2.5)),
            nerve: BuildGraph(1));

        var singleReport = NodeStats.From(single);
        Assert.Equal(1, singleReport.NodeCount);
        Assert.Equal(3, singleReport.MinSize);
        Assert.Equal(3, singleReport.MaxSize);
        Assert.Equal(3.0, singleReport.MedianSize);
        Assert.Equal(3.0, singleReport.MeanSize);

        var uniform = CreateResult(
            nodes: CreateNodes((new[] { 0, 1 }, -1.0), (new[] { 2, 3 }, 0.0), (new[] { 4, 5 }, 1.0)),
            nerve: BuildGraph(3, (0, 1, 1.0), (1, 2, 1.0)));

        var uniformReport = NodeStats.From(uniform);
        Assert.Equal(2, uniformReport.MinSize);
        Assert.Equal(2, uniformReport.MaxSize);
        Assert.Equal(2.0, uniformReport.MedianSize);
        Assert.Equal(2.0, uniformReport.MeanSize);
        Assert.Equal(new[] { -1.0, 0.0, 0.0, 1.0, 1.0 }, uniformReport.FilterValuePercentiles);

        var mixed = CreateResult(
            nodes: CreateNodes((new[] { 0 }, -3.0), (new[] { 1, 2, 3 }, -1.0), (new[] { 4, 5 }, 5.0), (new[] { 6, 7, 8, 9 }, 10.0)),
            nerve: BuildGraph(4, (0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0)));

        var mixedReport = NodeStats.From(mixed);
        Assert.Equal(1, mixedReport.MinSize);
        Assert.Equal(4, mixedReport.MaxSize);
        Assert.Equal(2.5, mixedReport.MedianSize);
        Assert.Equal(2.5, mixedReport.MeanSize);
        Assert.Equal(new[] { -3.0, -1.0, 5.0, 5.0, 10.0 }, mixedReport.FilterValuePercentiles);

        var empty = CreateResult(Array.Empty<MapperNode>(), BuildGraph(0));
        var emptyReport = NodeStats.From(empty);
        Assert.Equal(0, emptyReport.NodeCount);
        Assert.Empty(emptyReport.FilterValuePercentiles);
    }

    [Fact]
    public void Coverage_From_ReportsFullOverlapGapAndZeroPointCases()
    {
        var full = CreateResult(
            nodes: CreateNodes((new[] { 0, 1 }, 0.0), (new[] { 2, 3 }, 1.0)),
            nerve: BuildGraph(2, (0, 1, 1.0)),
            emptyBinCount: 1);

        var fullReport = Coverage.From(full, totalOriginalPoints: 4);
        Assert.Equal(0, fullReport.PointsUncovered);
        Assert.Equal(1.0, fullReport.CoverageFraction);
        Assert.Equal(1, fullReport.EmptyBinCount);

        var overlap = CreateResult(
            nodes: CreateNodes((new[] { 0, 1, 2 }, 0.0), (new[] { 2, 3, 4 }, 1.0), (new[] { 4, 5 }, 2.0)),
            nerve: BuildGraph(3, (0, 1, 1.0), (1, 2, 1.0)));

        var overlapReport = Coverage.From(overlap, totalOriginalPoints: 6);
        Assert.Equal(6, overlapReport.PointsCoveredAtLeastOnce);
        Assert.Equal(2, overlapReport.PointsCoveredMultipleTimes);
        Assert.True(overlapReport.MeanOverlapMultiplicity > 1.0);

        var gap = CreateResult(
            nodes: CreateNodes((new[] { 0, 1 }, 0.0), (new[] { 3 }, 1.0)),
            nerve: BuildGraph(2));

        var gapReport = Coverage.From(gap, totalOriginalPoints: 5);
        Assert.Equal(2, gapReport.PointsUncovered);
        Assert.True(gapReport.CoverageFraction < 1.0);

        var zeroReport = Coverage.From(CreateResult(Array.Empty<MapperNode>(), BuildGraph(0)), totalOriginalPoints: 0);
        Assert.Equal(0, zeroReport.TotalOriginalPoints);
        Assert.Equal(0.0, zeroReport.CoverageFraction);
        Assert.Equal(0.0, zeroReport.MeanOverlapMultiplicity);
    }

    [Fact]
    public void IntrinsicDimension_From_ReportsLineDiskBallAndSingletonNodes()
    {
        double[][] points =
        {
            new[] { -2.0, 0.0, 0.0 },
            new[] { -1.0, 0.0, 0.0 },
            new[] { 0.0, 0.0, 0.0 },
            new[] { 1.0, 0.0, 0.0 },
            new[] { 2.0, 0.0, 0.0 },

            new[] { -1.0, -1.0, 0.0 },
            new[] { -1.0, 1.0, 0.0 },
            new[] { 1.0, -1.0, 0.0 },
            new[] { 1.0, 1.0, 0.0 },
            new[] { 0.0, 0.0, 0.0 },

            new[] { 1.0, 0.0, 0.0 },
            new[] { -1.0, 0.0, 0.0 },
            new[] { 0.0, 1.0, 0.0 },
            new[] { 0.0, -1.0, 0.0 },
            new[] { 0.0, 0.0, 1.0 },
            new[] { 0.0, 0.0, -1.0 },

            new[] { 7.0, 7.0, 7.0 },
        };

        var result = CreateResult(
            nodes: CreateNodes(
                (new[] { 0, 1, 2, 3, 4 }, 0.0),
                (new[] { 5, 6, 7, 8, 9 }, 1.0),
                (new[] { 10, 11, 12, 13, 14, 15 }, 2.0),
                (new[] { 16 }, 3.0)),
            nerve: BuildGraph(4, (0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0)));

        var report = IntrinsicDimension.From(result, points);

        Assert.InRange(report.PerNodeIntrinsicDim[0], 0.99, 1.01);
        Assert.InRange(report.PerNodeIntrinsicDim[1], 1.9, 2.1);
        Assert.InRange(report.PerNodeIntrinsicDim[2], 2.9, 3.1);
        Assert.True(double.IsNaN(report.PerNodeIntrinsicDim[3]));
        Assert.Equal(1, report.EffectivelyOneDimensionalNodes);
        Assert.Equal(1, report.EffectivelyHigherDimNodes);
        Assert.InRange(report.MedianDim, 1.9, 2.1);
        Assert.InRange(report.MeanDim, 1.9, 2.1);

        var singletonOnly = CreateResult(
            nodes: CreateNodes((new[] { 0 }, 0.0)),
            nerve: BuildGraph(1));

        var singletonReport = IntrinsicDimension.From(singletonOnly, new[] { new[] { 1.0, 2.0, 3.0 } });
        Assert.True(double.IsNaN(singletonReport.PerNodeIntrinsicDim[0]));
        Assert.True(double.IsNaN(singletonReport.MedianDim));
        Assert.True(double.IsNaN(singletonReport.MeanDim));
    }

    [Fact]
    public void MapperWarnings_From_ReconstructsLegacyWarningShapes()
    {
        var empty = MapperWarnings.From(
            topology: new NerveTopologyReport(0, 0, 0, 0, 0, false, 0, -1),
            nodeStats: new MapperNodeReport(0, 0, 0, 0.0, 0.0, Array.Empty<double>()),
            coverage: new CoverageReport(0, 0, 0, 0, 0.0, 0.0, 0));

        Assert.Single(empty);
        Assert.Contains("produced no nodes", empty[0], StringComparison.Ordinal);

        var warnings = MapperWarnings.From(
            topology: new NerveTopologyReport(4, 3, 2, 3, 2, false, 0, 4),
            nodeStats: new MapperNodeReport(4, 2, 7, 4.5, 4.5, new[] { 0.0, 0.0, 0.0, 0.0, 0.0 }),
            coverage: new CoverageReport(10, 8, 3, 2, 0.8, 1.2, 5));

        Assert.Equal(4, warnings.Length);
        Assert.Contains("disconnected components", warnings[0], StringComparison.Ordinal);
        Assert.Contains("<3 members", warnings[1], StringComparison.Ordinal);
        Assert.Contains("5 cover bins are empty", warnings[2], StringComparison.Ordinal);
        Assert.Contains("2 loops", warnings[3], StringComparison.Ordinal);

        var none = MapperWarnings.From(
            topology: new NerveTopologyReport(4, 3, 1, 4, 0, true, 0, -1),
            nodeStats: new MapperNodeReport(4, 3, 6, 4.5, 4.5, new[] { 0.0, 0.0, 0.0, 0.0, 0.0 }),
            coverage: new CoverageReport(10, 10, 0, 0, 1.0, 1.0, 0));

        Assert.Empty(none);
    }

    private static MapperResult CreateResult(IReadOnlyList<MapperNode> nodes, CsrGraph nerve, int emptyBinCount = 0)
    {
        return new MapperResult
        {
            Nodes = nodes,
            Nerve = nerve,
            FilterName = "test-filter",
            CoverName = "test-cover",
            ClustererName = "test-clusterer",
            EmptyBinCount = emptyBinCount,
        };
    }

    private static MapperNode[] CreateNodes(params (int[] Members, double FilterMean)[] specs)
    {
        var nodes = new MapperNode[specs.Length];
        for (int i = 0; i < specs.Length; i++)
        {
            int[] members = specs[i].Members;
            nodes[i] = new MapperNode(
                BinId: i,
                LocalClusterId: 0,
                MemberIndices: members,
                FilterValueMean: specs[i].FilterMean,
                FilterValueMin: specs[i].FilterMean,
                FilterValueMax: specs[i].FilterMean);
        }

        return nodes;
    }

    private static CsrGraph BuildGraph(int nodeCount, params (int Source, int Target, double Weight)[] edges)
    {
        var graphEdges = new Edge[edges.Length];
        for (int i = 0; i < edges.Length; i++)
            graphEdges[i] = new Edge(edges[i].Source, edges[i].Target, edges[i].Weight);

        return CsrGraph.FromEdges(graphEdges, nodeCount);
    }
}
