using System;
using System.Collections.Generic;
using System.IO;
using Clustering.Evaluation.External;
using Clustering.Graphical.SPC;
using Clustering.Graphical.SPC.Export;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs;
using Graphs.Coupling;
using Graphs.Primitives;
using Repo.TestHarness;
using Synthetic.Euclidean;
using Xunit;

namespace VizCore.Tests;

public sealed class SpcSharedRunDirectoryFacts
{
    [Fact]
    public void WriteSpcArtifacts_ToHarnessSharedRunDirectory()
    {
        var dataset = BlattThreeCluster.Generate(pointsPerCluster: 20, stdDev: 1.0, seed: 42);

        var graphConfig = new GraphCompilerConfig
        {
            Topology   = new TopologyConfig { Kind = TopologyKind.Knn, K = 10 },
            Filter     = new FilterConfig { Kind = FilterKind.OrRule, MutualBandwidthSource = MutualBandwidthSource.DirectedKth },
            Repair     = new RepairConfig { Kind = RepairKind.MstMin },
            Refinement = new RefinementConfig { Kind = RefinementKind.Auto },
            Projection = new CouplingProjection { Kernel = new Gaussian(2.0), LmpRescale = false },
        };

        CsrGraph graph = SpcGraphBuilder.BuildResult(dataset.Features, graphConfig).Graph;
        SpcSessionResult result = SpcClusteringSession.Run(
            graph,
            sweepStrategy: new FixedGridSweepStrategy(AutoGridFixedSweep.BuildConfig(graph)),
            externalEvaluators: new[] { new Purity() },
            referenceLabels: dataset.Labels);

        ArtifactRun run = HarnessArtifacts.Create(
            runKind: "test-runs",
            suiteName: nameof(SpcSharedRunDirectoryFacts),
            runName: nameof(WriteSpcArtifacts_ToHarnessSharedRunDirectory),
            metadata: new Dictionary<string, object?>
            {
                ["Dataset"] = "BlattThreeCluster",
                ["PointsPerCluster"] = 20,
                ["GraphK"] = 10,
            });

        string runDirectory = SpcOutputPathHelper.CreateRunDirectory(run.RunDirectory, nameof(WriteSpcArtifacts_ToHarnessSharedRunDirectory));
        string sweepPath = SpcOutputPathHelper.GetSweepCsvPath(runDirectory);
        string partitionPath = SpcOutputPathHelper.GetPartitionCsvPath(runDirectory);
        string criteriaPath = SpcOutputPathHelper.GetCriteriaCsvPath(runDirectory);
        string sessionPath = SpcOutputPathHelper.GetSessionCsvPath(runDirectory);

        SpcCsvWriter.WriteSweepProfile(result.Profile, sweepPath);
        SpcCsvWriter.WritePartition(result.Partition, partitionPath, features: dataset.Features, trueLabels: dataset.Labels);
        SpcCsvWriter.WriteCriteria(result.ProfileCriteria, criteriaPath);
        SpcCsvWriter.WriteSessionSummary(result, sessionPath);

        run.WriteRunText("summary", "SPC CSV artifacts generated for shared harness directory.");

        Console.WriteLine($"RunRoot\t{run.RunDirectory}");
        Console.WriteLine($"Manifest\t{run.ManifestPath}");
        Console.WriteLine($"Sweep\t{sweepPath}");
        Console.WriteLine($"Partition\t{partitionPath}");
        Console.WriteLine($"Criteria\t{criteriaPath}");
        Console.WriteLine($"Session\t{sessionPath}");

        Assert.True(File.Exists(sweepPath));
        Assert.True(File.Exists(partitionPath));
        Assert.True(File.Exists(criteriaPath));
        Assert.True(File.Exists(sessionPath));
        Assert.InRange(result.Purity.GetValueOrDefault(), 0.0, 1.0);
    }
}
