using System;
using System.Linq;
using Graphs;
using Graphs.Diagnostics;
using Graphs.Primitives;
using Graphs.Proximity;
using Graphs.TestSupport;
using Xunit;

using Graphs.Observables;

namespace VizCore.Tests;

public sealed class GraphDiagnosticsTests
{
    [Fact]
    public void AlgebraicConnectivity_ConnectedGraph_MatchesHandComputedValue()
    {
        CsrGraph graph = BuildGraph(2,
            (0, 1, 1.0));

        AlgebraicConnectivityReport report = AlgebraicConnectivity.Compute(graph);

        Assert.True(report.Computed);
        Assert.Equal(2, report.NodeCount);
        Assert.InRange(Math.Abs(report.Lambda2 - 2.0), 0.0, 1e-9);
        Assert.False(report.LikelyWeaklyConnected);
    }

    [Fact]
    public void AlgebraicConnectivity_WeaklyConnectedGraph_HasSmallPositiveFiedlerValue()
    {
        CsrGraph graph = BuildGraph(2,
            (0, 1, 2.5e-7));

        AlgebraicConnectivityReport report = AlgebraicConnectivity.Compute(graph);

        Assert.True(report.Computed);
        Assert.Equal(2, report.NodeCount);
        Assert.True(report.Lambda2 > 0.0);
        Assert.InRange(Math.Abs(report.Lambda2 - 5.0e-7), 0.0, 1e-12);
        Assert.True(report.LikelyWeaklyConnected);
    }

    [Fact]
    public void AlgebraicConnectivity_DisconnectedGraph_HasZeroFiedlerValue()
    {
        CsrGraph graph = BuildGraph(4,
            (0, 1, 1.0),
            (2, 3, 1.0));

        AlgebraicConnectivityReport report = AlgebraicConnectivity.Compute(graph);

        Assert.True(report.Computed);
        Assert.Equal(4, report.NodeCount);
        Assert.InRange(Math.Abs(report.Lambda2), 0.0, 1e-9);
        Assert.True(report.LikelyWeaklyConnected);
    }

    [Fact]
    public void EdgeWeights_Summary_ComputesUndirectedStatistics()
    {
        CsrGraph graph = BuildGraph(3,
            (0, 1, 2.0),
            (1, 2, 0.5));

        EdgeWeightSummary summary = EdgeWeights.Summary(graph);

        Assert.Equal(2, summary.EdgeCount);
        Assert.Equal(0, summary.NearZeroBridges);
        Assert.InRange(Math.Abs(summary.MinWeight - 0.5), 0.0, 1e-12);
        Assert.InRange(Math.Abs(summary.MaxWeight - 2.0), 0.0, 1e-12);
        Assert.InRange(Math.Abs(summary.MedianWeight - 1.25), 0.0, 1e-12);
        Assert.InRange(Math.Abs(summary.MeanWeight - 1.25), 0.0, 1e-12);
    }

    [Fact]
    public void Degree_Distribution_ReportsIsolatedAndUndersampledNodes()
    {
        CsrGraph graph = BuildGraph(4,
            (0, 1, 1.0),
            (1, 2, 1.0));

        DegreeReport report = Degree.Distribution(graph);

        Assert.Equal(4, report.NodeCount);
        Assert.Equal(0, report.MinDegree);
        Assert.Equal(2, report.MaxDegree);
        Assert.InRange(Math.Abs(report.MeanDegree - 1.0), 0.0, 1e-12);
        Assert.Equal(1, report.IsolatedCount);
        Assert.Equal(2, report.UndersampledCount);
    }

    [Fact]
    public void NeighborhoodScale_Compute_ContrastsMutualAndDirectedNeighborhoods()
    {
        double[][] points =
        {
            new[] { 0.0 },
            new[] { 1.0 },
            new[] { 10.0 },
            new[] { 11.0 }
        };

        NeighborSelection directed = ProximityGraph.SelectKnn(points.Length, 2, (i, j) => EuclideanDistance(points[i], points[j]));

        NeighborhoodScaleReport report = NeighborhoodScale.Compute(directed, 2);

        Assert.Equal(2, report.K);
        Assert.InRange(Math.Abs(report.Median1NN - 1.0), 0.0, 1e-12);
        Assert.InRange(Math.Abs(report.MedianKthNN - 9.5), 0.0, 1e-12);
        Assert.InRange(Math.Abs(report.ScaleRatio - 9.5), 0.0, 1e-12);
    }

    [Fact]
    public void MstBridge_Compare_FindsInjectedBridgeEdges()
    {
        DisconnectedControl.Fixture control = DisconnectedControl.Generate(pointsPerComponent: 8, separation: 100.0);
        double[][] points = control.Points;

        CsrGraph withoutRepair = GraphCompilerTestPresets.BuildResult(
            n: points.Length,
            dist: (i, j) => EuclideanDistance(points[i], points[j]),
            topologyKind: TopologyKind.Knn,
            filterKind: FilterKind.OrRule,
            k: 1,
            kernel: KernelType.Gaussian,
            ensureConnected: false).Graph;

        CsrGraph withRepair = GraphCompilerTestPresets.BuildResult(
            n: points.Length,
            dist: (i, j) => EuclideanDistance(points[i], points[j]),
            topologyKind: TopologyKind.Knn,
            filterKind: FilterKind.OrRule,
            k: 1,
            kernel: KernelType.Gaussian,
            ensureConnected: true).Graph;

        MstBridgeReport report = MstBridge.Compare(withoutRepair, withRepair);

        Assert.True(report.BridgeCount >= 1);
        Assert.True(report.MinBridgeWeight > 0.0);
        Assert.True(report.MaxBridgeWeight >= report.MinBridgeWeight);
        Assert.True(report.MedianBridgeWeight >= report.MinBridgeWeight);
    }

    [Fact]
    public void Hubness_Analyze_StarLikeDirectedSelection_FindsSingleHub()
    {
        NeighborSelection directed = BuildSelection(
            new[] { 1 },
            new[] { 0 },
            new[] { 0 },
            new[] { 0 },
            new[] { 0 });

        HubnessReport report = Hubness.Analyze(directed, k: 1, hubMultiple: 3.0, antiHubMultiple: 0.0);

        Assert.Equal(1, report.K);
        Assert.Equal(5, report.NodeCount);
        Assert.Equal(4, report.MaxInDegree);
        Assert.Equal(1, report.HubCount);
        Assert.True(report.InDegreeSkewness > 1.0);
        Assert.InRange(report.TopHubCoverage, 0.8, 1.0);
    }

    [Fact]
    public void Hubness_Analyze_BalancedDirectedSelection_HasZeroSkew()
    {
        NeighborSelection directed = BuildSelection(
            new[] { 1, 2 },
            new[] { 2, 3 },
            new[] { 3, 4 },
            new[] { 4, 5 },
            new[] { 5, 6 },
            new[] { 6, 7 },
            new[] { 7, 8 },
            new[] { 8, 0 },
            new[] { 0, 1 });

        HubnessReport report = Hubness.Analyze(directed, k: 2, hubMultiple: 2.0, antiHubMultiple: 0.5);

        Assert.Equal(2, report.MeanInDegree);
        Assert.Equal(2, report.MaxInDegree);
        Assert.InRange(Math.Abs(report.InDegreeSkewness), 0.0, 1e-12);
        Assert.Equal(0, report.HubCount);
        Assert.Equal(0, report.AntiHubCount);
        Assert.InRange(Math.Abs(report.TopHubCoverage - (1.0 / 3.0)), 0.0, 1e-12);
    }

    [Fact]
    public void Hubness_Analyze_EmptySelection_ReturnsGuardValues()
    {
        HubnessReport report = Hubness.Analyze(new NeighborSelection(Array.Empty<Neighbor[]>(), Array.Empty<double>(), Array.Empty<double>()), k: 3);

        Assert.Equal(0, report.NodeCount);
        Assert.Equal(0, report.MaxInDegree);
        Assert.Equal(0.0, report.MeanInDegree);
        Assert.True(double.IsNaN(report.InDegreeSkewness));
        Assert.Equal(0, report.HubCount);
        Assert.Equal(0, report.AntiHubCount);
        Assert.Equal(0.0, report.TopHubCoverage);
    }

    [Fact]
    public void Cycles_Compute_Tree_IsAcyclic()
    {
        CsrGraph graph = BuildGraph(5,
            (0, 1, 1.0),
            (1, 2, 1.0),
            (1, 3, 1.0),
            (3, 4, 1.0));

        CycleReport report = Cycles.Compute(graph);

        Assert.Equal(0, report.CyclomaticComplexity);
        Assert.Equal(0, report.TriangleCount);
        Assert.Equal(-1, report.Girth);
        Assert.Equal(0.0, report.MeanCycleLength);
        Assert.Equal(0, report.MaxCycleLength);
    }

    [Fact]
    public void Cycles_Compute_C5_ReportsSingleLengthFiveCycle()
    {
        CsrGraph graph = BuildGraph(5,
            (0, 1, 1.0),
            (1, 2, 1.0),
            (2, 3, 1.0),
            (3, 4, 1.0),
            (4, 0, 1.0));

        CycleReport report = Cycles.Compute(graph);

        Assert.Equal(1, report.CyclomaticComplexity);
        Assert.Equal(0, report.TriangleCount);
        Assert.Equal(5, report.Girth);
        Assert.InRange(Math.Abs(report.MeanCycleLength - 5.0), 0.0, 1e-12);
        Assert.Equal(5, report.MaxCycleLength);
    }

    [Fact]
    public void Cycles_Compute_K3_IsTriangleOnly()
    {
        CsrGraph graph = BuildGraph(3,
            (0, 1, 1.0),
            (1, 2, 1.0),
            (0, 2, 1.0));

        CycleReport report = Cycles.Compute(graph);

        Assert.Equal(1, report.CyclomaticComplexity);
        Assert.Equal(1, report.TriangleCount);
        Assert.Equal(3, report.Girth);
        Assert.InRange(Math.Abs(report.TriangleSaturation - 1.0), 0.0, 1e-12);
        Assert.InRange(Math.Abs(report.MeanCycleLength - 3.0), 0.0, 1e-12);
        Assert.Equal(3, report.MaxCycleLength);
    }

    [Fact]
    public void Cycles_Compute_K4_ReportsTriangleDominatedMesh()
    {
        CsrGraph graph = BuildGraph(4,
            (0, 1, 1.0),
            (0, 2, 1.0),
            (0, 3, 1.0),
            (1, 2, 1.0),
            (1, 3, 1.0),
            (2, 3, 1.0));

        CycleReport report = Cycles.Compute(graph);

        Assert.Equal(3, report.CyclomaticComplexity);
        Assert.Equal(4, report.TriangleCount);
        Assert.Equal(3, report.Girth);
        Assert.InRange(Math.Abs(report.TriangleSaturation - (4.0 / 3.0)), 0.0, 1e-12);
        Assert.InRange(Math.Abs(report.MeanCycleLength - 3.0), 0.0, 1e-12);
        Assert.Equal(3, report.MaxCycleLength);
    }

    [Fact]
    public void Cycles_Compute_TwoDisconnectedTriangles_PreservesComponentAdjustedCyclomaticComplexity()
    {
        CsrGraph graph = BuildGraph(6,
            (0, 1, 1.0),
            (1, 2, 1.0),
            (0, 2, 1.0),
            (3, 4, 1.0),
            (4, 5, 1.0),
            (3, 5, 1.0));

        CycleReport report = Cycles.Compute(graph);

        Assert.Equal(2, report.CyclomaticComplexity);
        Assert.Equal(2, report.TriangleCount);
        Assert.Equal(3, report.Girth);
        Assert.InRange(Math.Abs(report.TriangleSaturation - 1.0), 0.0, 1e-12);
        Assert.InRange(Math.Abs(report.MeanCycleLength - 3.0), 0.0, 1e-12);
        Assert.Equal(3, report.MaxCycleLength);
    }

    [Fact]
    public void Cycles_Compute_ThetaGraph_CapturesDistinctLongerCycleLengths()
    {
        CsrGraph graph = BuildGraph(8,
            (0, 2, 1.0),
            (2, 1, 1.0),
            (0, 3, 1.0),
            (3, 4, 1.0),
            (4, 1, 1.0),
            (0, 5, 1.0),
            (5, 6, 1.0),
            (6, 7, 1.0),
            (7, 1, 1.0));

        CycleReport report = Cycles.Compute(graph);

        Assert.Equal(2, report.CyclomaticComplexity);
        Assert.Equal(0, report.TriangleCount);
        Assert.Equal(5, report.Girth);
        Assert.InRange(Math.Abs(report.MeanCycleLength - 6.0), 0.0, 1e-12);
        Assert.Equal(7, report.MaxCycleLength);
    }

    [Fact]
    public void Cycles_Compute_OverCapGraph_ReturnsSentinelCycleStatsWithoutTrianglePass()
    {
        CsrGraph graph = BuildGraph(5001,
            (0, 1, 1.0),
            (1, 2, 1.0),
            (0, 2, 1.0));

        CycleReport report = Cycles.Compute(graph, maxNodesForCycleStats: 5000);

        Assert.Equal(1, report.CyclomaticComplexity);
        Assert.Equal(-1, report.TriangleCount);
        Assert.Equal(-1, report.Girth);
        Assert.True(double.IsNaN(report.TriangleSaturation));
        Assert.True(double.IsNaN(report.MeanCycleLength));
        Assert.Equal(-1, report.MaxCycleLength);
    }

    private static CsrGraph BuildGraph(int nodeCount, params (int Source, int Target, double Weight)[] edges)
    {
        var graphEdges = new Edge[edges.Length];
        for (int i = 0; i < edges.Length; i++)
            graphEdges[i] = new Edge(edges[i].Source, edges[i].Target, edges[i].Weight);

        return CsrGraph.FromEdges(graphEdges, nodeCount);
    }

    private static NeighborSelection BuildSelection(params int[][] neighborIndices)
    {
        Neighbor[][] rows = neighborIndices
            .Select(indices => indices.Select(index => new Neighbor { Index = index, Distance = 1.0 }).ToArray())
            .ToArray();

        double[] nearestNeighborDistances = rows
            .Where(row => row.Length > 0)
            .Select(row => row[0].Distance)
            .ToArray();

        double[] kthNeighborDistances = rows
            .Where(row => row.Length > 0)
            .Select(row => row[row.Length - 1].Distance)
            .ToArray();

        return new NeighborSelection(rows, nearestNeighborDistances, kthNeighborDistances);
    }

    private static double EuclideanDistance(double[] left, double[] right)
    {
        double sum = 0.0;
        for (int i = 0; i < left.Length; i++)
        {
            double diff = left[i] - right[i];
            sum += diff * diff;
        }

        return Math.Sqrt(sum);
    }
}
