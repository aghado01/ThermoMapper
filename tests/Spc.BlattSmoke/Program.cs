using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Clustering.Graphical.SPC;
using Clustering.Graphical.SPC.Export;
using Clustering.Graphical.SPC.Partitions;
using Clustering.Graphical.SPC.Partitions.Strategies;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Clustering.Evaluation.External;
using Maths.Distance.Euclidean;
using Graphs.Primitives;
using Repo.TestHarness;
using Synthetic;
using Synthetic.Euclidean;

namespace Spc.BlattSmoke;

/// <summary>
/// End-to-end smoke harness: Blatt 1997 Euclidean hierarchy → mutual-kNN
/// CsrGraph → SpcScheduler dump. Produces one checkpoint file per
/// (Temperature, Replica) task in an artifacts directory ready for
/// offline χ(T) / purity analysis.
/// </summary>
/// <remarks>
/// <para>Graph construction is inlined here rather than going through
/// <c>GraphBuilder.Build</c> from <c>Graphs.Proximity</c> — that project
/// has its own pending rename ripple. Distance computation itself comes
/// from <c>Graphs.Distance.Euclidean</c> (the real assembly), so the
/// only inline math is the kNN selection + Gaussian kernel application,
/// not the metric. When Graphs.Proximity stabilizes, swap the helpers
/// below for the canonical builder; no other code in this harness
/// changes.</para>
/// </remarks>
internal static class Program
{
    // ── Blatt 1997 canonical configuration ────────────────────────────────
    private const int                 K              = 10;
    private const int                 Q              = 20;
    private const int                 CyclesPerTask  = 1000;
    private const int                 BurnInCycles   = 200;
    private const int                 NumReplicas    = 3;
    private const int                        BaseSeed       = 42;
    private static readonly AccumulationSpec   AccSpec        = AccumulationSpec.Currencies;

    private static int Main(string[] args)
    {
        // Top-down T schedule per SPC convention: T_max → T_min.
        // 0.20 → 0.005 step 0.005 brackets Blatt's reported T_fs ≈ 0.075.
        double[] temperatures = Enumerable.Range(0, 40)
            .Select(i => 0.20 - i * 0.005)
            .Where(t => t > 0)
            .ToArray();

        ArtifactRun run = args.Length > 0
            ? HarnessArtifacts.Attach(
                runKind: "spc",
                suiteName: "Spc.BlattSmoke",
                runName: "Main",
                runDirectory: args[0],
                metadata: new Dictionary<string, object?>
                {
                    ["OutputDirectoryMode"] = "user-specified",
                })
            : HarnessArtifacts.Create(
                runKind: "spc",
                suiteName: "Spc.BlattSmoke",
                runName: "Main",
                metadata: new Dictionary<string, object?>
                {
                    ["OutputDirectoryMode"] = "default",
                });
        string outputDir = run.RunDirectory;
        string checkpointDir = SpcOutputPathHelper.GetCheckpointDirectory(outputDir);
        string csvDir = SpcOutputPathHelper.GetCsvDirectory(outputDir);

        // ── Stage 1: synthetic hierarchy ──────────────────────────────────
        var dataset = BlattHierarchy.Generate(seed: BaseSeed);
        int n = dataset.Features.Length;
        int levels = dataset.LabelsByLevel?.Length ?? 0;

        // ── Stage 2: pairwise Euclidean distances ─────────────────────────
        var distances = ComputePairwiseDistances(dataset.Features);

        // ── Stage 3: mutual-kNN selection ─────────────────────────────────
        var pairs = SelectMutualKnn(distances, K);

        // ── Stage 4: bandwidth — median of 1-NN distances ─────────────────
        double bandwidth = MedianOneNnDistance(distances);

        // ── Stage 5: Gaussian-weighted Edge[] → CsrGraph ──────────────────
        var edges = new Edge[pairs.Count];
        double twoBwSq = 2.0 * bandwidth * bandwidth;
        for (int e = 0; e < pairs.Count; e++)
        {
            var (i, j) = pairs[e];
            double d = distances[i, j];
            double coupling = Math.Exp(-(d * d) / twoBwSq);
            edges[e] = new Edge(i, j, coupling);
        }
        var graph = CsrGraph.FromEdges(edges, n);

        // ── Stage 6: flat task list ───────────────────────────────────────
        var executor = new SpcExecutor();
        var tasks = executor.BuildTaskList(
            temperatures:        temperatures,
            numReplicas:         NumReplicas,
            q:                   Q,
            accumulation:        AccSpec,
            budget:              new RunBudget(BurnInCycles, CyclesPerTask),
            checkpointDirectory: checkpointDir,
            baseSeed:            BaseSeed);

        // ── Stage 7: run ──────────────────────────────────────────────────
        var sw = Stopwatch.StartNew();
        executor.Run(graph, tasks);
        sw.Stop();

        int written = Directory.GetFiles(checkpointDir, "*.spcx").Length;

        // ── Stage 8: SPC partition cut + per-level purity per T ───────────
        var groundTruth = dataset.LabelsByLevel;
        string? schedulePath = null;
        int scheduleRowCount = 0;
        if (groundTruth is not null && groundTruth.Length > 0)
            (schedulePath, scheduleRowCount) = WritePartitionScheduleCsv(checkpointDir, csvDir, graph, groundTruth, theta: 0.5);

        string summaryPath = run.WriteJson(
            "summary.json",
            new
            {
                Dataset = new
                {
                    PointCount = n,
                    dataset.ClusterCount,
                    GroundTruthLevels = levels,
                },
                Graph = new
                {
                    MutualK = K,
                    UndirectedEdgeCount = pairs.Count,
                    Bandwidth = bandwidth,
                    CsrEntries = graph.Targets.Length,
                },
                Schedule = new
                {
                    Temperatures = temperatures,
                    NumReplicas,
                    TaskCount   = tasks.Count,
                    SweepBurnIn = BurnInCycles,
                    SweepCycles = CyclesPerTask,
                    Accumulation = "Currencies",
                },
                Run = new
                {
                    CheckpointsWritten = written,
                    ExpectedCheckpoints = tasks.Count,
                    WallSeconds = sw.Elapsed.TotalSeconds,
                },
                PartitionSchedule = new
                {
                    Path = schedulePath,
                    RowCount = scheduleRowCount,
                },
            });

        Console.WriteLine($"RunRoot\t{run.RunDirectory}");
        Console.WriteLine($"Manifest\t{run.ManifestPath}");
        Console.WriteLine($"Summary\t{summaryPath}");
        if (schedulePath is not null)
            Console.WriteLine($"PartitionScheduleCsv\t{schedulePath}");

        return 0;
    }

    private static readonly string[] LevelNames = { "coarse", "medium", "fine" };

    private static (string? Path, int RowCount) WritePartitionScheduleCsv(
        string spceDirectory,
        string outputDir,
        CsrGraph graph,
        int[][] levels,
        double theta)
    {
        var spcxFiles = Directory.GetFiles(spceDirectory, "*.spcx");
        if (spcxFiles.Length == 0)
            return (null, 0);

        var observables = spcxFiles.Select(AccumulatorSerializer.Instance.ReadFromFile);
        var rollups = SpcScheduleHelpers.BuildPartitionScheduleRollups(
            graph,
            observables,
            levels,
            theta,
            LevelNames);

        if (rollups.Count == 0)
            return (null, 0);

        string schedulePath = Path.Combine(outputDir, SpcCsvWriter.PartitionScheduleFileName);
        SpcCsvWriter.WritePartitionScheduleRollups(rollups, schedulePath);
        return (schedulePath, rollups.Count);
    }

    // ── Inline graph construction (TODO: swap for GraphBuilder.Build) ─────

    private static double[,] ComputePairwiseDistances(double[][] points)
    {
        int n = points.Length;
        var d = new double[n, n];
        Parallel.For(0, n, i =>
        {
            for (int j = i + 1; j < n; j++)
            {
                double dist = Minkowski.Distance(points[i], points[j], 2.0);
                d[i, j] = dist;
                d[j, i] = dist;
            }
        });
        return d;
    }

    private static List<(int i, int j)> SelectMutualKnn(double[,] distances, int k)
    {
        int n = distances.GetLength(0);
        var knn = new int[n][];

        Parallel.For(0, n, i =>
        {
            var pairs = new (double d, int j)[n - 1];
            int idx = 0;
            for (int j = 0; j < n; j++)
            {
                if (i == j) continue;
                pairs[idx++] = (distances[i, j], j);
            }
            Array.Sort(pairs, (a, b) => a.d.CompareTo(b.d));
            knn[i] = pairs.Take(k).Select(p => p.j).ToArray();
        });

        var sets = knn.Select(row => new HashSet<int>(row)).ToArray();
        var edges = new List<(int, int)>(n * k / 2);
        for (int i = 0; i < n; i++)
        {
            foreach (int j in knn[i])
            {
                if (j > i && sets[j].Contains(i))
                    edges.Add((i, j));
            }
        }
        return edges;
    }

    private static double MedianOneNnDistance(double[,] distances)
    {
        int n = distances.GetLength(0);
        var oneNn = new double[n];
        Parallel.For(0, n, i =>
        {
            double min = double.PositiveInfinity;
            for (int j = 0; j < n; j++)
            {
                if (i == j) continue;
                double d = distances[i, j];
                if (d < min) min = d;
            }
            oneNn[i] = min;
        });
        Array.Sort(oneNn);
        return n % 2 == 0
            ? 0.5 * (oneNn[n / 2 - 1] + oneNn[n / 2])
            : oneNn[n / 2];
    }
}
