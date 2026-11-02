using System;
using System.IO;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Observables;
using UserRepl.Commands;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The resolved construction a build returns beside the config that requested it, and the
/// <see cref="AutoThresholds"/> every <c>Auto</c> decides by (graph-compiler-brief, decision 6).
/// </summary>
public sealed class GraphCompilerResolvedTests
{
    private static readonly double[][] Points = MakePoints();

    private static double[][] MakePoints()
    {
        var points = new double[24][];
        for (int i = 0; i < 12; i++)
        {
            points[i] = new[] { 0.1 * i + 0.01 * i * i, 0.0 };
            points[12 + i] = new[] { 8.0 + 0.12 * i + 0.012 * i * i, 0.5 };
        }
        return points;
    }

    private static GraphMetric Metric() => new((i, j) =>
    {
        double dx = Points[i][0] - Points[j][0], dy = Points[i][1] - Points[j][1];
        return Math.Sqrt(dx * dx + dy * dy);
    });

    private static GraphCompilerConfig Config(
        IEdgeProjection projection, FilterKind filter = FilterKind.Auto, RepairKind repair = RepairKind.Auto, AutoThresholds? auto = null) =>
        new()
        {
            Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 4 },
            Filter = new FilterConfig { Kind = filter },
            Repair = new RepairConfig { Kind = repair },
            Projection = projection,
            Auto = auto ?? new AutoThresholds(),
        };

    private static CouplingProjection Coupling(IKernelDescriptor? kernel = null, bool? lmp = false) =>
        new() { Kernel = kernel ?? new Cauchy(), LmpRescale = lmp };

    private static GraphBuildResult Build(GraphCompilerConfig config) =>
        GraphCompiler.Build(config, Points.Length, Metric());

    [Fact]
    public void AutoFilter_ResolvesBySkewnessAgainstTheConfiguredThreshold()
    {
        GraphBuildResult byDefault = Build(Config(Coupling()));
        GraphBuildResult lowered = Build(Config(Coupling(), auto: new AutoThresholds { FilterSkewness = -1000.0 }));

        Assert.Equal(FilterKind.OrRule, byDefault.Resolved!.Filter);
        Assert.Equal(FilterKind.MutualKnn, lowered.Resolved!.Filter);
        Assert.NotNull(byDefault.Resolved.InDegreeSkewness);
    }

    [Fact]
    public void AutoRepair_ResolvesByCoverageAgainstTheConfiguredThreshold()
    {
        // Two lines far apart: coverage after the OR-rule filter is 0.5.
        GraphBuildResult byDefault = Build(Config(Coupling()));
        GraphBuildResult relaxed = Build(Config(Coupling(), auto: new AutoThresholds { RepairCoverage = 0.4 }));

        Assert.Equal(RepairKind.MstMin, byDefault.Resolved!.Repair);
        Assert.Equal(RepairKind.NoRepair, relaxed.Resolved!.Repair);
        Assert.Equal(0.5, byDefault.Resolved.CoverageBeforeRepair!.Value, 12);
        Assert.Equal(1, byDefault.Resolved.RepairEdges);
        Assert.Equal(2, relaxed.Resolved.Components);
    }

    [Fact]
    public void AutoLmp_RecordsTheRatioItReadAndWhatItDecided()
    {
        GraphBuildResult byDefault = Build(Config(Coupling(lmp: null)));
        // Every coupling is under 2, so the whole graph counts as collapsed.
        GraphBuildResult lowered = Build(Config(Coupling(lmp: null), auto: new AutoThresholds { NearZeroCoupling = 2.0 }));

        Assert.False(byDefault.Resolved!.Coupling!.Lmp);
        Assert.NotNull(byDefault.Resolved.Coupling.LmpNearZeroRatio);
        Assert.True(lowered.Resolved!.Coupling!.Lmp);
        Assert.Equal(1.0, lowered.Resolved.Coupling.LmpNearZeroRatio!.Value, 12);
    }

    [Fact]
    public void ExplicitLmp_RecordsNoRatio()
    {
        Assert.Null(Build(Config(Coupling(lmp: true))).Resolved!.Coupling!.LmpNearZeroRatio);
        Assert.Null(Build(Config(Coupling(lmp: false))).Resolved!.Coupling!.LmpNearZeroRatio);
    }

    [Fact]
    public void Resolved_RecordsTheUnitsAndTheConnectivityOfTheReturnedGraph()
    {
        var projection = new CouplingProjection
        {
            Kernel = new Gaussian(),
            LmpRescale = false,
            BandwidthOverride = BandwidthStrategy.MeanEdgeDistance,
            Gauge = GaugeKind.InverseMeanDegree,
            CouplingFloor = 1e-12,
        };

        GraphBuildResult result = Build(Config(projection));

        CouplingResolution c = result.Resolved!.Coupling!;
        Assert.Equal("Gaussian", c.Kernel);
        Assert.Equal(BandwidthStrategy.MeanEdgeDistance, c.BandwidthStrategy);
        Assert.Equal(result.SingleBandwidth!.Value, c.Bandwidth!.Value, 12);
        Assert.Equal(1e-12, c.CouplingFloor);
        Assert.Equal(GaugeKind.InverseMeanDegree, c.Gauge);
        Assert.Equal(result.Graph.NodeCount / (double)result.Graph.Targets.Length, c.GaugeFactor, 12);
        Assert.True(c.HasBwdUnits);
        Assert.Equal(Connectivity.Validate(result.Graph).ComponentCount, result.Resolved.Components);
        Assert.Equal(1, c.CoupledComponents);
        Assert.Equal(0, c.RepairEdgesDropped);
    }

    [Fact]
    public void HasBwdUnits_NeedsGaussian_MeanEdgeBandwidth_TheGauge_AndNoLmp()
    {
        CouplingProjection classic = new()
        {
            Kernel = new Gaussian(), LmpRescale = false,
            BandwidthOverride = BandwidthStrategy.MeanEdgeDistance, Gauge = GaugeKind.InverseMeanDegree,
        };

        Assert.True(Build(Config(classic)).Resolved!.Coupling!.HasBwdUnits);
        Assert.False(Build(Config(classic with { Kernel = new Cauchy() })).Resolved!.Coupling!.HasBwdUnits);
        Assert.False(Build(Config(classic with { Gauge = GaugeKind.None })).Resolved!.Coupling!.HasBwdUnits);
        Assert.False(Build(Config(classic with { BandwidthOverride = BandwidthStrategy.MadConsistencyFactor })).Resolved!.Coupling!.HasBwdUnits);
        Assert.False(Build(Config(classic with { LmpRescale = true })).Resolved!.Coupling!.HasBwdUnits);
    }

    [Fact]
    public void DistanceBuild_ResolvesNoCoupling()
    {
        GraphBuildResult result = Build(Config(new DistanceProjection()));

        Assert.Null(result.Resolved!.Coupling);
        Assert.Equal(1, result.Resolved.Components);
    }

    [Fact]
    public void Manifest_CarriesTheResolvedConstructionBesideTheRequestedConfig()
    {
        GraphCompilerConfig config = Config(new CouplingProjection
        {
            Kernel = new Gaussian(), LmpRescale = null,
            BandwidthOverride = BandwidthStrategy.MeanEdgeDistance, Gauge = GaugeKind.InverseMeanDegree,
        });
        GraphBuildResult result = Build(config);
        GraphConstructionManifest manifest = GraphConstructionPersistence.MaterializeManifest(config, result, "fingerprint");

        string directory = Path.Combine(Path.GetTempPath(), "graph-manifest-" + Guid.NewGuid().ToString("N"));
        try
        {
            manifest.WriteTo(directory);
            GraphConstructionManifest read = GraphConstructionManifest.ReadFrom(directory);

            Assert.Equal(result.Resolved, read.Resolved);
            Assert.Equal(FilterKind.Auto, read.Config.Filter.Kind);      // requested
            Assert.Equal(FilterKind.OrRule, read.Resolved!.Filter);      // resolved
            Assert.Null(((CouplingProjection)read.Config.Projection).LmpRescale);
            Assert.False(read.Resolved.Coupling!.Lmp);
        }
        finally
        {
            if (Directory.Exists(directory)) Directory.Delete(directory, recursive: true);
        }
    }
}
