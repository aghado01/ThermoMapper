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
    public void GaussianRepairBridge_UnderTheFloor_IsDroppedFromTheReturnedGraph()
    {
        // §2.1 as it stands: repair joins the two lines, then the 1e-10 floor drops the
        // bridge (Gaussian coupling e^(-d²/2a²) with d ≈ 10 and a ≈ 0.2 underflows).
        double[][] points = FarApartLines(lines: 2, pointsPerLine: 10, separation: 10.0);

        GraphBuildResult result = Build(points, Config(new Gaussian(), k: 3));

        Assert.True(result.EnsureConnectedApplied);
        ConnectivityReport returned = Connectivity.Validate(result.Graph);
        Assert.Equal(2, returned.ComponentCount);

        Assert.Contains(result.Diagnostics!.Messages, m =>
            m.Stage == "Connectivity" && m.Text.StartsWith("Returned graph: 2 component(s)", StringComparison.Ordinal));
        Assert.Contains(result.Diagnostics.Warnings, m =>
            m.Stage == "Connectivity" && m.Text.StartsWith("1 of 1 repair edge(s) are absent", StringComparison.Ordinal));
        Assert.Contains(result.Diagnostics.Warnings, m =>
            m.Text.Contains("repaired topology had 1 component(s); the returned graph has 2", StringComparison.Ordinal));
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

    [Fact]
    public void NearZeroInterrupt_CountsTheEdgesTheFloorAlreadyDropped()
    {
        // One bridge among ~20 edges is dropped by the floor: 1/21 of the produced edges
        // are collapsed. The kept edges alone show a ratio of 0, which the interrupt used to read.
        double[][] points = FarApartLines(lines: 2, pointsPerLine: 10, separation: 10.0);
        GraphCompilerConfig config = Config(
            new Gaussian(), k: 3, interrupts: new PathologyInterruptConfig { MaxNearZeroEdgeRatio = 0.01 });

        var thrown = Assert.Throws<GraphPathologyException>(() => Build(points, config));
        Assert.Contains("NearZeroEdgeRatio", thrown.Message, StringComparison.Ordinal);
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
    public void LmpWritesTheRepairBridgeAZero_AndTheLogCountsOnlyTheCoupledEdges()
    {
        // The Cauchy bridge survives the floor, then Local Mutual Proximity ranks it the weakest
        // edge at both ends and writes it exactly 0: in the graph, joining nothing.
        double[][] points = FarApartLines(lines: 2, pointsPerLine: 10, separation: 10.0);

        GraphBuildResult result = Build(points, Config(new Cauchy(), k: 3, lmp: true));

        Assert.Equal(1, Connectivity.Validate(result.Graph).ComponentCount);
        int coupledComponents = CoupledComponents(result.Graph);
        Assert.True(coupledComponents >= 2, $"Expected the zeroed bridge to split the coupled graph; got {coupledComponents}.");
        Assert.Equal(coupledComponents, Connectivity.Validate(result.Graph, positiveWeightsOnly: true).ComponentCount);

        Assert.Contains(result.Diagnostics!.Warnings, m =>
            m.Stage == "Connectivity"
            && m.Text.Contains("carry zero coupling", StringComparison.Ordinal)
            && m.Text.Contains($"the graph has {coupledComponents} component(s)", StringComparison.Ordinal));
        Assert.Contains(result.Diagnostics.Warnings, m =>
            m.Text.StartsWith("1 of 1 repair edge(s) carry zero coupling", StringComparison.Ordinal));
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
    public void PointsWithinTheFloor_ReportOneComponentAndNoWarning()
    {
        // Control: a single line is connected by its own topology, no repair, no cut.
        double[][] points = FarApartLines(lines: 1, pointsPerLine: 12, separation: 10.0);

        GraphBuildResult result = Build(points, Config(new Gaussian(), k: 3));

        Assert.Equal(1, Connectivity.Validate(result.Graph).ComponentCount);
        Assert.DoesNotContain(result.Diagnostics!.Warnings, m => m.Stage == "Connectivity");
    }
}
