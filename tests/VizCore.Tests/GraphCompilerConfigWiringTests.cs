using System;
using System.IO;
using System.Text.Json;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using UserRepl.Commands;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The REPL's side of the compiler revision: the BWD presets keep their kernel by declaring the gauge
/// and the post-repair bandwidth, the new flags reach the projection, and the projection survives the
/// manifest's JSON.
/// </summary>
public sealed class GraphCompilerConfigWiringTests
{
    private static string PresetPath(string name)
    {
        for (DirectoryInfo? dir = new(AppContext.BaseDirectory); dir is not null; dir = dir.Parent)
        {
            string candidate = Path.Combine(dir.FullName, "presets", "spc", name + ".json");
            if (File.Exists(candidate)) return candidate;
        }
        throw new FileNotFoundException($"presets/spc/{name}.json not found above the test binary.");
    }

    [Theory]
    [InlineData("bwd1995-toy")]
    [InlineData("bwd1996-iris")]
    [InlineData("domany1999-iris")]
    public void BwdPresets_DeclareTheMeanEdgeBandwidth_TheGauge_AndTheAfterRepairEdges(string name)
    {
        var options = new SpcCommand.Options();
        SpcPreset.Load(PresetPath(name)).ApplyTo(options);

        var projection = (CouplingProjection)GraphRun.BuildGraphConfig(options).Projection;

        Assert.IsType<Gaussian>(projection.Kernel);
        Assert.Equal(BandwidthStrategy.MeanEdgeDistance, projection.BandwidthOverride);
        Assert.Equal(GaugeKind.InverseMeanDegree, projection.Gauge);
        Assert.Equal(ScaleEdgeSet.AfterRepair, projection.ScaleEdges);
        Assert.False(projection.LmpRescale);
    }

    [Fact]
    public void TheNewFlags_ReachTheProjection()
    {
        string[] args =
        {
            "--gauge", "inverse-mean-degree", "--scale-edges", "after-repair", "--scale-rank", "3",
            "--coupling-floor", "1e-12", "--lmp-rank", "inclusive",
        };
        var options = new SpcCommand.Options();
        for (int i = 0; i < args.Length; i++) Assert.True(GraphRun.TryParse(options, args, ref i));

        var projection = (CouplingProjection)GraphRun.BuildGraphConfig(options).Projection;

        Assert.Equal(GaugeKind.InverseMeanDegree, projection.Gauge);
        Assert.Equal(ScaleEdgeSet.AfterRepair, projection.ScaleEdges);
        Assert.Equal(3, projection.ScaleRank);
        Assert.Equal(1e-12, projection.CouplingFloor);
        Assert.Equal(LmpRankRule.Inclusive, projection.LmpRank);
    }

    [Fact]
    public void TheDefaults_AreNoGauge_BeforeRepair_AndTheHistoricalFloorAndRank()
    {
        var projection = (CouplingProjection)GraphRun.BuildGraphConfig(new SpcCommand.Options()).Projection;

        Assert.Equal(GaugeKind.None, projection.Gauge);
        Assert.Equal(ScaleEdgeSet.BeforeRepair, projection.ScaleEdges);
        Assert.Null(projection.ScaleRank);
        Assert.Equal(1e-10, projection.CouplingFloor);
        Assert.Equal(LmpRankRule.Strict, projection.LmpRank);
    }

    [Fact]
    public void TheProjectionAndTheAutoThresholds_RoundTripThroughTheManifestJson()
    {
        var config = new GraphCompilerConfig
        {
            Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 7 },
            Projection = new CouplingProjection
            {
                Kernel = new Cauchy(0.25),
                LmpRescale = null,
                BandwidthOverride = BandwidthStrategy.MeanEdgeDistance,
                Gauge = GaugeKind.InverseMeanDegree,
                ScaleEdges = ScaleEdgeSet.AfterRepair,
                ScaleRank = 4,
                CouplingFloor = 1e-12,
                LmpRank = LmpRankRule.Inclusive,
            },
            Auto = new AutoThresholds { FilterSkewness = 2.5, RepairCoverage = 0.9, NearZeroCoupling = 1e-9, LmpNearZeroRatio = 0.1 },
        };

        string json = JsonSerializer.Serialize(config, GraphConstructionManifest.JsonOptions);
        GraphCompilerConfig read = JsonSerializer.Deserialize<GraphCompilerConfig>(json, GraphConstructionManifest.JsonOptions)!;

        Assert.Equal(config.Auto, read.Auto);
        Assert.Equal(config.Projection, read.Projection);
    }
}
