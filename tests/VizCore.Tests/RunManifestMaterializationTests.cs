using System;
using System.Collections.Generic;
using System.IO;
using System.Text.Json;
using Clustering.Graphical.SPC.Runtime.Execution;
using Graphs;
using Graphs.Coupling;
using Graphs.Diagnostics;
using Graphs.Distance;
using Graphs.Distance.Geodesic;
using UserRepl.Commands;
using Xunit;

namespace VizCore.Tests;

public sealed class RunManifestMaterializationTests
{
    [Fact]
    public void Materialize_SyntheticDataset_RoundTripsStringGeneratorParametersFromManifest()
    {
        var manifest = new RunManifest(
            SchemaVersion: RunManifest.CurrentSchemaVersion,
            CreatedUtc: DateTime.UtcNow,
            Algorithm: "spc",
            CommandLine: "spc --dataset TwoMoons --param pointsPerMoon=150",
            Dataset: new DatasetSpec(
                Source: "synthetic",
                GeneratorName: "TwoMoons",
                GeneratorParameters: new Dictionary<string, object?>
                {
                    ["pointsPerMoon"] = "150",
                    ["noise"] = "0.1",
                },
                Seed: 42,
                CsvPath: null,
                LabelColumn: null,
                HasHeader: null,
                Delimiter: null),
            Graph: null,
            Sweep: null,
            Hdbscan: null,
            Output: new OutputSpec(
                RunDirectory: "artifacts/test",
                CheckpointDirectory: null));

        string json = JsonSerializer.Serialize(manifest, RunManifest.JsonOptions);
        RunManifest roundTripped = JsonSerializer.Deserialize<RunManifest>(json, RunManifest.JsonOptions)
            ?? throw new InvalidOperationException("Manifest deserialized to null.");

        var dataset = roundTripped.Dataset.Materialize();

        Assert.Equal(300, dataset.Features.Length);
        Assert.Equal(2, dataset.ClusterCount);
        Assert.Equal(42, dataset.Metadata["seed"]);
        Assert.Equal("TwoMoons", dataset.Metadata["Generator"]);
    }

    [Fact]
    public void GraphConstructionManifest_WriteRead_RoundTripsMetricProperties()
    {
        string graphDirectory = Path.Combine(
            Path.GetTempPath(),
            $"graph-manifest-{Guid.NewGuid():N}");

        try
        {
            double[][] features =
            {
                new[] { 0.00, 0.00 },
                new[] { 0.18, 0.05 },
                new[] { 0.34, 0.11 },
                new[] { 0.49, 0.16 },
            };
            var config = new GraphCompilerConfig
            {
                Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 2 },
                Filter = new FilterConfig { Kind = FilterKind.OrRule },
                Repair = new RepairConfig { Kind = RepairKind.NoRepair },
                Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
                Projection = new CouplingProjection { Kernel = new Gaussian(0.0), LmpRescale = false },
            };
            GraphBuildResult build = SpcGraphBuilder.BuildResult(features, config, new PoincareMetric());
            var manifest = new GraphConstructionManifest(
                SchemaVersion: GraphConstructionManifest.CurrentSchemaVersion,
                CreatedUtc: new DateTime(2026, 6, 2, 0, 0, 0, DateTimeKind.Utc),
                DatasetFingerprint: "abc123",
                Config: config,
                Metric: build.Metric,
                Diagnostics: new DiagnosticsLog());

            manifest.WriteTo(graphDirectory);
            string json = File.ReadAllText(GraphConstructionManifest.PathFor(graphDirectory));
            using JsonDocument document = JsonDocument.Parse(json);

            JsonElement metricElement = document.RootElement.GetProperty("metric");
            Assert.True(metricElement.TryGetProperty("geometry", out JsonElement geometryElement));
            if (geometryElement.ValueKind == JsonValueKind.String)
                Assert.Equal("hyperbolic", geometryElement.GetString(), ignoreCase: true);
            else
                Assert.Equal((int)SpaceGeometry.Hyperbolic, geometryElement.GetInt32());

            GraphConstructionManifest roundTripped = GraphConstructionManifest.ReadFrom(graphDirectory);

            Assert.Equal(manifest.SchemaVersion, roundTripped.SchemaVersion);
            Assert.Equal(manifest.DatasetFingerprint, roundTripped.DatasetFingerprint);
            Assert.Equal(build.Metric, roundTripped.Metric);
            Assert.Equal(SpaceGeometry.Hyperbolic, roundTripped.Metric?.Geometry);
        }
        finally
        {
            if (Directory.Exists(graphDirectory))
                Directory.Delete(graphDirectory, recursive: true);
        }
    }
}
