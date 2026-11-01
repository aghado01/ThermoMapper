using System;
using System.Linq;
using Graphs;
using Graphs.Coupling;
using Graphs.Diagnostics;
using Graphs.Observables;
using Graphs.Pipeline;
using Graphs.Pipeline.Scalers;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// What the compiler says about the graph it returns, and the coupling floor that
/// decides it (graph-compiler-brief §2.1, §2.2). The §2.1 pin reads today's
/// behaviour: the floor drops a repair bridge the repair stage had just counted.
/// </summary>
public sealed class GraphCompilerReturnedGraphTests
{
    /// <summary>Three short lines of points, far apart: a bridge between lines is the longest edge by construction.</summary>
    private static double[][] FarApartLines(int lines, int pointsPerLine, double separation)
    {
        var points = new double[lines * pointsPerLine][];
        for (int line = 0; line < lines; line++)
            for (int i = 0; i < pointsPerLine; i++)
                // Uneven gaps (0.1 + 0.02 i) keep the neighbour distances free of ties.
                points[line * pointsPerLine + i] = new[] { line * separation + 0.1 * i + 0.01 * i * i, 0.0 };
        return points;
    }

    private static GraphCompilerConfig Config(IKernelDescriptor kernel, int k, bool? lmp = false, PathologyInterruptConfig? interrupts = null) =>
        new()
        {
            Topology   = new TopologyConfig { Kind = TopologyKind.Knn, K = k },
            Filter     = new FilterConfig { Kind = FilterKind.OrRule },
            Repair     = new RepairConfig { Kind = RepairKind.MstMin },
            Projection = new CouplingProjection { Kernel = kernel, LmpRescale = lmp },
            Interrupts = interrupts ?? new PathologyInterruptConfig(),
        };

    private static GraphBuildResult Build(double[][] points, GraphCompilerConfig config) =>
        GraphCompiler.Build(config, points.Length, new GraphMetric((i, j) =>
        {
            double dx = points[i][0] - points[j][0], dy = points[i][1] - points[j][1];
            return Math.Sqrt(dx * dx + dy * dy);
        }));

    [Fact]
    public void GaussianRepairBridge_UnderTheFloor_StaysInTheReturnedGraph()
    {
        // §2.1: repair joins the two lines with a bridge whose Gaussian coupling e^(-d²/2a²)
        // (d ≈ 10, a ≈ 0.2) underflows under the floor. The floor no longer cuts repair edges.
        double[][] points = FarApartLines(lines: 2, pointsPerLine: 10, separation: 10.0);

        GraphBuildResult result = Build(points, Config(new Gaussian(), k: 3));

        Assert.True(result.EnsureConnectedApplied);
        Assert.Equal(1, Connectivity.Validate(result.Graph).ComponentCount);
        Assert.Equal(1, CoupledComponents(result.Graph));
        Assert.Contains(result.Diagnostics!.Messages, m =>
            m.Stage == "Connectivity" && m.Text.StartsWith("Returned graph: 1 component(s)", StringComparison.Ordinal));
        Assert.DoesNotContain(result.Diagnostics.Warnings, m => m.Stage == "Connectivity");

        // The bridge carries the clamp: positive, far under the floor.
        double bridge = BridgeWeight(result.Graph, points.Length / 2);
        Assert.InRange(bridge, double.Epsilon, 1e-10);
    }

    private static double WeightOf(CsrGraph graph, int from, int to)
    {
        for (int slot = graph.RowPointers[from]; slot < graph.RowPointers[from + 1]; slot++)
            if (graph.Targets[slot] == to) return graph.Weights[slot];
        throw new InvalidOperationException($"No edge {from}-{to}.");
    }

    /// <summary>The smallest weight on an edge that crosses from the first half of the nodes to the second.</summary>
    private static double BridgeWeight(CsrGraph graph, int half)
    {
        double smallest = double.PositiveInfinity;
        for (int i = 0; i < half; i++)
            for (int slot = graph.RowPointers[i]; slot < graph.RowPointers[i + 1]; slot++)
                if (graph.Targets[slot] >= half) smallest = Math.Min(smallest, graph.Weights[slot]);
        return smallest;
    }

    [Fact]
    public void CauchyRepairBridge_SurvivesTheFloor_AndTheLogReportsOneComponent()
    {
        // Same topology, power-law tail: the bridge keeps a coupling near 1e-4.
        double[][] points = FarApartLines(lines: 2, pointsPerLine: 10, separation: 10.0);

        GraphBuildResult result = Build(points, Config(new Cauchy(), k: 3));

        Assert.Equal(1, Connectivity.Validate(result.Graph).ComponentCount);
        Assert.Contains(result.Diagnostics!.Messages, m =>
            m.Stage == "Connectivity" && m.Text.StartsWith("Returned graph: 1 component(s)", StringComparison.Ordinal));
        Assert.DoesNotContain(result.Diagnostics.Warnings, m => m.Stage == "Connectivity");
    }

    [Fact]
    public void DistanceProjection_LogsTheReturnedGraphsConnectivity()
    {
        double[][] points = FarApartLines(lines: 2, pointsPerLine: 10, separation: 10.0);
        GraphCompilerConfig config = Config(new Gaussian(), k: 3) with { Projection = new DistanceProjection() };

        GraphBuildResult result = Build(points, config);

        Assert.Equal(1, Connectivity.Validate(result.Graph).ComponentCount);
        Assert.Contains(result.Diagnostics!.Messages, m =>
            m.Stage == "Connectivity" && m.Text.StartsWith("Returned graph: 1 component(s)", StringComparison.Ordinal));
    }

    private static Neighbor Near(int index, double distance) => new() { Index = index, Distance = distance };

    [Fact]
    public void GlobalBandwidthScaler_CountsTheEdgesTheFloorDrops()
    {
        // Two edges: d = 1 (coupling e^-0.5) and d = 10 (e^-50 ≈ 2e-22, under the floor).
        var rows = new[]
        {
            new[] { Near(1, 1.0), Near(2, 10.0) },
            new[] { Near(0, 1.0) },
            new[] { Near(0, 10.0) },
        };
        var selection = new NeighborSelection(rows, new[] { 1.0, 1.0, 10.0 }, new[] { 10.0, 1.0, 10.0 });

        ScalerResult scaled = new GlobalBandwidthScaler(KernelType.Gaussian, bandwidth: 1.0).Scale(selection, 3);

        Assert.Equal(1, scaled.CutEdges);
        Assert.Equal(2, scaled.Graph.Targets.Length);   // the d = 1 edge, in both CSR slots
    }

    /// <summary>A line of points and one outlier far off its end: the outlier's own nominations are the long edges.</summary>
    private static double[][] LineWithOutlier(int pointsPerLine, double outlierAt)
    {
        double[][] line = FarApartLines(lines: 1, pointsPerLine: pointsPerLine, separation: 0.0);
        return line.Append(new[] { outlierAt, 0.0 }).ToArray();
    }

    [Fact]
    public void FloorStillCutsATopologyEdge_AndTheLogSaysTheGraphIsLessConnectedThanRepairReported()
    {
        // The outlier is joined by its own nominations (no repair), and the floor drops all three:
        // only repair edges are exempt, so the outlier is isolated in the returned graph.
        double[][] points = LineWithOutlier(pointsPerLine: 10, outlierAt: 15.0);

        GraphBuildResult result = Build(points, Config(new Gaussian(), k: 3));

        ConnectivityReport returned = Connectivity.Validate(result.Graph);
        Assert.Equal(2, returned.ComponentCount);
        Assert.Contains(result.Diagnostics!.Messages, m =>
            m.Stage == "Scaling" && m.Text.StartsWith("3 edge(s) fell under the", StringComparison.Ordinal));
        Assert.Contains(result.Diagnostics.Warnings, m =>
            m.Text.Contains("repaired topology had 1 component(s); the returned graph has 2", StringComparison.Ordinal));
    }

    [Fact]
    public void NearZeroInterrupt_CountsTheEdgesTheFloorAlreadyDropped()
    {
        // Three dropped edges among ~20: the kept edges alone show a ratio of 0.
        double[][] points = LineWithOutlier(pointsPerLine: 10, outlierAt: 15.0);
        GraphCompilerConfig config = Config(
            new Gaussian(), k: 3, interrupts: new PathologyInterruptConfig { MaxNearZeroEdgeRatio = 0.01 });

        var thrown = Assert.Throws<GraphPathologyException>(() => Build(points, config));
        Assert.Contains("NearZeroEdgeRatio", thrown.Message, StringComparison.Ordinal);
    }

    [Fact]
    public void CouplingFloorOfZero_KeepsEveryEdge()
    {
        double[][] points = LineWithOutlier(pointsPerLine: 10, outlierAt: 15.0);
        GraphCompilerConfig config = Config(new Gaussian(), k: 3) with
        {
            Projection = new CouplingProjection { Kernel = new Gaussian(), LmpRescale = false, CouplingFloor = 0.0 },
        };

        GraphBuildResult result = Build(points, config);

        Assert.Equal(1, Connectivity.Validate(result.Graph).ComponentCount);
        Assert.DoesNotContain(result.Diagnostics!.Messages, m => m.Text.Contains("coupling floor", StringComparison.Ordinal));
    }

    [Fact]
    public void GlobalBandwidthScaler_KeepsARepairEdgeAtItsKernelCoupling_AndClampsAnUnderflow()
    {
        // Edges 0-1 (d = 1), 0-2 (d = 10: e^-50 ≈ 2e-22), 0-3 (d = 40: e^-800 underflows to 0).
        var rows = new[]
        {
            new[] { Near(1, 1.0), Near(2, 10.0), Near(3, 40.0) },
            new[] { Near(0, 1.0) },
            new[] { Near(0, 10.0) },
            new[] { Near(0, 40.0) },
        };
        var selection = new NeighborSelection(rows, new[] { 1.0, 1.0, 10.0, 40.0 }, new[] { 40.0, 1.0, 10.0, 40.0 });
        var repair = RepairEdgeSet.Between(
            new NeighborSelection(new[] { rows[0][..1], rows[1], Array.Empty<Neighbor>(), Array.Empty<Neighbor>() }, new double[4], new double[4]),
            selection, 4);
        Assert.Equal(2, repair.Count);

        ScalerResult scaled = new GlobalBandwidthScaler(
            KernelType.Gaussian, bandwidth: 1.0, context: new ScalingContext(repair)).Scale(selection, 4);

        Assert.Equal(0, scaled.CutEdges);
        Assert.Equal(Math.Exp(-0.5), WeightOf(scaled.Graph, 0, 1), 12);
        Assert.InRange(WeightOf(scaled.Graph, 0, 2) / Math.Exp(-50.0), 1 - 1e-9, 1 + 1e-9);
        Assert.Equal(ScalingContext.RepairEdgeMinimum, WeightOf(scaled.Graph, 0, 3));
    }

    [Fact]
    public void AutoLmp_ReadsTheKeptEdgesOnly_SoDroppedEdgesDoNotSwitchItOn()
    {
        // Three short lines, two bridges dropped: the dropped edges alone would exceed the 5%
        // threshold, but LMP fragments sparse graphs, so the trigger counts the kept edges only.
        double[][] points = FarApartLines(lines: 3, pointsPerLine: 5, separation: 10.0);

        GraphBuildResult result = Build(points, Config(new Gaussian(), k: 2, lmp: null));

        Assert.DoesNotContain(result.Diagnostics!.Warnings, m => m.Text.StartsWith("Auto-LMP enabled", StringComparison.Ordinal));
        Assert.Contains(result.Diagnostics.Infos, m => m.Text.Contains("LMP is not required", StringComparison.Ordinal));
    }

    /// <summary>Components over the edges of positive weight, by a union-find independent of the library's.</summary>
    private static int CoupledComponents(CsrGraph graph)
    {
        int[] parent = Enumerable.Range(0, graph.NodeCount).ToArray();
        int Find(int x) { while (parent[x] != x) x = parent[x] = parent[parent[x]]; return x; }
        for (int i = 0; i < graph.NodeCount; i++)
            for (int slot = graph.RowPointers[i]; slot < graph.RowPointers[i + 1]; slot++)
                if (graph.Weights[slot] > 0.0) parent[Find(i)] = Find(graph.Targets[slot]);
        return Enumerable.Range(0, graph.NodeCount).Count(x => Find(x) == x);
    }

    [Fact]
    public void Lmp_KeepsTheRepairBridgeAtItsKernelCoupling()
    {
        // The Cauchy bridge is the weakest edge at both ends, which the strict rank would write 0.
        // Repair edges join the protect set, so they keep their kernel coupling.
        double[][] points = FarApartLines(lines: 2, pointsPerLine: 10, separation: 10.0);

        GraphBuildResult withLmp = Build(points, Config(new Cauchy(), k: 3, lmp: true));
        GraphBuildResult withoutLmp = Build(points, Config(new Cauchy(), k: 3, lmp: false));

        double bridge = BridgeWeight(withLmp.Graph, points.Length / 2);
        Assert.True(bridge > 0.0);
        Assert.Equal(BridgeWeight(withoutLmp.Graph, points.Length / 2), bridge);
        Assert.Contains(withLmp.Diagnostics!.Infos, m => m.Text.StartsWith("LMP will keep the 1 repair edge(s)", StringComparison.Ordinal));
        Assert.DoesNotContain(withLmp.Diagnostics.Warnings, m =>
            m.Text.Contains("repair edge(s) carry zero coupling", StringComparison.Ordinal));
    }

    [Fact]
    public void StrictLmp_ZerosEdgesInsideTheLines_AndTheLogCountsOnlyTheCoupledEdges()
    {
        double[][] points = FarApartLines(lines: 2, pointsPerLine: 10, separation: 10.0);

        GraphBuildResult result = Build(points, Config(new Cauchy(), k: 3, lmp: true));

        int zeroEdges = Enumerable.Range(0, result.Graph.Weights.Length).Count(slot => result.Graph.Weights[slot] == 0.0) / 2;
        Assert.True(zeroEdges > 0, "The strict rank writes an exact 0 on each node's weakest edge.");
        int coupledComponents = CoupledComponents(result.Graph);
        Assert.Equal(coupledComponents, Connectivity.Validate(result.Graph, positiveWeightsOnly: true).ComponentCount);
        Assert.Contains(result.Diagnostics!.Messages, m =>
            m.Stage == "Connectivity"
            && m.Text.StartsWith($"{zeroEdges} of ", StringComparison.Ordinal)
            && m.Text.Contains($"the graph has {coupledComponents} component(s)", StringComparison.Ordinal));
    }

    [Fact]
    public void InclusiveLmp_WritesNoZero_AndKeepsTheCoupledGraphConnected()
    {
        double[][] points = FarApartLines(lines: 2, pointsPerLine: 10, separation: 10.0);
        GraphCompilerConfig config = Config(new Cauchy(), k: 3) with
        {
            Projection = new CouplingProjection { Kernel = new Cauchy(), LmpRescale = true, LmpRank = LmpRankRule.Inclusive },
        };

        GraphBuildResult result = Build(points, config);

        Assert.DoesNotContain(result.Graph.Weights, w => w == 0.0);
        Assert.Equal(1, CoupledComponents(result.Graph));
    }

    [Fact]
    public void WithoutLmp_EveryEdgeOfAConnectedCauchyGraphCarriesCoupling()
    {
        double[][] points = FarApartLines(lines: 2, pointsPerLine: 10, separation: 10.0);

        GraphBuildResult result = Build(points, Config(new Cauchy(), k: 3, lmp: false));

        Assert.Equal(1, CoupledComponents(result.Graph));
        Assert.DoesNotContain(result.Diagnostics!.Messages, m => m.Text.Contains("carry zero coupling", StringComparison.Ordinal));
    }

    [Fact]
    public void PreRepairGraph_IsTheReturnedGraphWithoutTheRepairEdges_InTheSameUnits()
    {
        double[][] points = FarApartLines(lines: 2, pointsPerLine: 10, separation: 10.0);

        GraphBuildResult result = Build(points, Config(new Cauchy(), k: 3));

        CsrGraph pre = result.PreRepairGraph!.Value;
        ConnectivityReport preReport = Connectivity.Validate(pre);
        Assert.Equal(2, preReport.ComponentCount);
        Assert.Equal(Connectivity.Validate(result.Graph).EdgeCount - 1, preReport.EdgeCount);
        for (int i = 0; i < pre.NodeCount; i++)
            for (int slot = pre.RowPointers[i]; slot < pre.RowPointers[i + 1]; slot++)
                Assert.Equal(WeightOf(result.Graph, i, pre.Targets[slot]), pre.Weights[slot]);

        MstBridgeReport bridges = MstBridge.Compare(pre, result.Graph);
        Assert.Equal(1, bridges.BridgeCount);
        Assert.Equal(BridgeWeight(result.Graph, points.Length / 2), bridges.MinBridgeWeight);
    }

    [Fact]
    public void PointsWithinTheFloor_ReportOneComponentAndNoWarning()
    {
        // Control: a single line is connected by its own topology, no repair, no cut.
        double[][] points = FarApartLines(lines: 1, pointsPerLine: 12, separation: 10.0);

        GraphBuildResult result = Build(points, Config(new Gaussian(), k: 3));

        Assert.Equal(1, Connectivity.Validate(result.Graph).ComponentCount);
        Assert.DoesNotContain(result.Diagnostics!.Warnings, m => m.Stage == "Connectivity");
    }
}
