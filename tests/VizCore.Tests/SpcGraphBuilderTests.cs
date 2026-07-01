using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Graphical.SPC;
using Clustering.Graphical.SPC.Runtime.Execution;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Distance.Geodesic;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

public sealed class SpcGraphBuilderTests
{
    private readonly struct ConstantMetric : IDistanceMetric
    {
        private static readonly MetricProperties _props = new(
            IsBounded:           true,
            MaxValue:            1.0,
            RequiresProbability: false,
            RequiresUnitNorm:    false,
            FixedDimension:      null,
            Geometry:            SpaceGeometry.Euclidean,
            BandwidthStrategy:   BandwidthStrategy.MadConsistencyFactor,
            Name:                "Constant");

        public MetricProperties Properties => _props;

        public double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b) => 1.0;
    }

    [Fact]
    public void Build_FromFeatures_CreatesGraphWithExpectedNodes()
    {
        double[][] features =
        {
            new[] { 0.0, 0.0 },
            new[] { 1.0, 0.0 },
            new[] { 0.0, 1.0 },
        };

        var config = new GraphCompilerConfig
        {
            Topology   = new TopologyConfig { Kind = TopologyKind.Knn, K = 2 },
            Filter     = new FilterConfig { Kind = FilterKind.OrRule, MutualBandwidthSource = MutualBandwidthSource.DirectedKth },
            Repair     = new RepairConfig { Kind = RepairKind.MstMin },
            Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
            Projection = new CouplingProjection { Kernel = new Gaussian(1.0), LmpRescale = false },
        };

        CsrGraph graph = SpcGraphBuilder.BuildResult(features, config).Graph;

        Assert.Equal(3, graph.NodeCount);
        Assert.True(graph.Targets.Length >= 3);
    }

    [Fact]
    public void Build_UsesCustomDistanceMetricAndLmpPostProcessing()
    {
        double[][] features =
        {
            new[] { 0.0, 0.0 },
            new[] { 1.0, 0.0 },
            new[] { 0.0, 1.0 },
        };

        var metric = new ConstantMetric();
        var config = new GraphCompilerConfig
        {
            Topology   = new TopologyConfig { Kind = TopologyKind.Knn, K = 2 },
            Filter     = new FilterConfig { Kind = FilterKind.OrRule, MutualBandwidthSource = MutualBandwidthSource.DirectedKth },
            Repair     = new RepairConfig { Kind = RepairKind.NoRepair },
            Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
            Projection = new CouplingProjection { Kernel = new Gaussian(1.0), LmpRescale = true },
        };

        GraphBuildResult build = SpcGraphBuilder.BuildResult(features, config, metric);
        CsrGraph graph = build.Graph;

        Assert.Equal(3, graph.NodeCount);
        Assert.NotEmpty(graph.Weights);
        Assert.All(graph.Weights, weight => Assert.InRange(weight, 0.0, 1.0));
        Assert.Equal(metric.Properties, build.Metric);
    }

    [Fact]
    public void Build_WithSphericalGeodesicMetric_UsesMetricBandwidthStrategy()
    {
        double[][] features =
        {
            new[] { 1.0, 0.0 },
            new[] { 0.0, 1.0 },
            new[] { -1.0, 0.0 },
        };

        var metric = new SphericalGeodesicMetric();
        var config = new GraphCompilerConfig
        {
            Topology   = new TopologyConfig { Kind = TopologyKind.Knn, K = 2 },
            Filter     = new FilterConfig { Kind = FilterKind.OrRule, MutualBandwidthSource = MutualBandwidthSource.DirectedKth },
            Repair     = new RepairConfig { Kind = RepairKind.NoRepair },
            Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
            Projection = new CouplingProjection { Kernel = new Gaussian(1.0), LmpRescale = false },
        };

        GraphBuildResult build = SpcGraphBuilder.BuildResult(features, config, metric);
        CsrGraph graph = build.Graph;

        Assert.Equal(3, graph.NodeCount);
        Assert.NotEmpty(graph.Targets);
        Assert.NotEmpty(graph.Weights);
        Assert.Equal(metric.Properties, build.Metric);
        IEnumerable<Graphs.Diagnostics.DiagnosticMessage> diagnostics =
            (IEnumerable<Graphs.Diagnostics.DiagnosticMessage>?)build.Diagnostics?.Messages
            ?? Array.Empty<Graphs.Diagnostics.DiagnosticMessage>();
        Assert.Contains(
            diagnostics,
            message => message.Stage == "Scaling"
                && message.Text.Contains($"Strategy={metric.Properties.BandwidthStrategy}", StringComparison.Ordinal));
    }

    [Fact]
    public void Build_WithWasserstein1Metric_CreatesGraph()
    {
        double[][] features =
        {
            new[] { 0.2, 0.8 },
            new[] { 0.5, 0.5 },
            new[] { 0.1, 0.9 },
        };

        var metric = new Wasserstein1Metric();
        var config = new GraphCompilerConfig
        {
            Topology   = new TopologyConfig { Kind = TopologyKind.Knn, K = 1 },
            Filter     = new FilterConfig { Kind = FilterKind.OrRule, MutualBandwidthSource = MutualBandwidthSource.DirectedKth },
            Repair     = new RepairConfig { Kind = RepairKind.NoRepair },
            Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
            Projection = new CouplingProjection { Kernel = new Gaussian(1.0), LmpRescale = false },
        };

        GraphBuildResult build = SpcGraphBuilder.BuildResult(features, config, metric);
        CsrGraph graph = build.Graph;

        Assert.Equal(3, graph.NodeCount);
        Assert.Equal(3, graph.RowPointers.Length - 1);
        Assert.NotEmpty(graph.Weights);
        Assert.Equal(metric.Properties, build.Metric);
    }
}
