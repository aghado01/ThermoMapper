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
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Clustering.Evaluation.External;
using Maths.Distance.Geodesic;
using Graphs.Primitives;
using Repo.TestHarness;
using Synthetic;
using Synthetic.Manifolds;

namespace Spc.HyperbolicSmoke;

/// <summary>
/// End-to-end smoke harness: hyperbolic Blatt hierarchy in the Poincaré ball
/// → Cauchy-weighted mutual-kNN CsrGraph → SpcScheduler dump. Mirrors the
/// structure of <c>Spc.BlattSmoke</c> with four targeted swaps for the
/// hyperbolic regime per the mapper-spc kernel-decay analysis.
/// </summary>
/// <remarks>
/// <para>Differences from <c>Spc.BlattSmoke</c> (driver-level, not PottsModel-level):</para>
/// <list type="bullet">
///   <item>Generator: <c>GenerateHyperbolicBlattHierarchy</c> (points emerge already
///     in the Poincaré ball; no projection step needed).</item>
///   <item>Distance: <see cref="Poincare.Distance"/> — arcosh formula in the unit ball.</item>
///   <item>Kernel: Cauchy <c>1 / (1 + d²/δ²)</c> instead of Gaussian. Hyperbolic
///     distance distributions are heavy-tailed; Gaussian's <c>exp(-d²/2δ²)</c>
///     over-decays and crushes the multi-scale hierarchical couplings that the
///     SPC sweep is meant to discover. At <c>d/δ ≈ 5</c>, Cauchy gives a coupling
///     <c>~10⁴×</c> larger than Gaussian.</item>
///   <item>Q = 25 (vs Blatt's 20). Cauchy's polynomial tail populates the
///     intermediate-density region of the correlation distribution at SP-phase;
///     larger Q pushes <c>1/q</c> farther from that region for cleaner cut
///     discrimination.</item>
///   <item>T schedule wider, shifted up: <c>1.5 → 0.05</c> instead of <c>0.20 → 0.005</c>.
///     The heavier J distribution shifts <c>T_fs</c> upward (expect ≈0.3 vs Blatt's
///     0.075).</item>
/// </list>
///
/// <para>Inline graph construction (same as Spc.BlattSmoke) bypasses the
/// in-flux Graphs.Proximity project. The bandwidth diagnostic prints the
/// 1-NN <c>max/median</c> ratio so heavy-tail pathologies are visible at
/// run-time — if &gt; 10, the median-based δ may underestimate the
/// distribution spread.</para>
/// </remarks>
internal static class Program
{
    private const int                 K              = 10;
    private const int                 Q              = 25;
    private const int                 CyclesPerTask  = 1000;
    private const int                 BurnInCycles   = 200;
    private const int                 NumReplicas    = 3;
    private const int                        BaseSeed       = 42;
    private static readonly AccumulationSpec   AccSpec        = AccumulationSpec.Currencies;

    private static int Main(string[] args)
    {
        // Top-down per SPC convention, wider than Euclidean Blatt.
        double[] temperatures = Enumerable.Range(0, 30)
            .Select(i => 1.50 - i * 0.05)
            .Where(t => t >= 0.05)
            .ToArray();

        ArtifactRun run = args.Length > 0
            ? HarnessArtifacts.Attach(
                runKind: "spc",
                suiteName: "Spc.HyperbolicSmoke",
                runName: "Main",
                runDirectory: args[0],
                metadata: new Dictionary<string, object?>
                {
                    ["OutputDirectoryMode"] = "user-specified",
                })
            : HarnessArtifacts.Create(
                runKind: "spc",
                suiteName: "Spc.HyperbolicSmoke",
                runName: "Main",
                metadata: new Dictionary<string, object?>
                {
                    ["OutputDirectoryMode"] = "default",
                });
        string outputDir = run.RunDirectory;
        string checkpointDir = SpcOutputPathHelper.GetCheckpointDirectory(outputDir);
        string csvDir = SpcOutputPathHelper.GetCsvDirectory(outputDir);

        // ── Stage 1: synthetic Poincaré-ball hierarchy ────────────────────
        var dataset = HyperbolicHierarchy.Generate(seed: BaseSeed);
        int n = dataset.Features.Length;
        int levels = dataset.LabelsByLevel?.Length ?? 0;

        // ── Stage 2: pairwise hyperbolic distances ────────────────────────
        var distances = ComputePairwiseDistances(dataset.Features);

        // ── Stage 3: mutual-kNN selection ─────────────────────────────────
        var pairs = SelectMutualKnn(distances, K);

        // ── Stage 4: bandwidth + heavy-tail diagnostic ────────────────────
        var (bandwidth, maxOverMedian) = EstimateBandwidthAndDiagnostic(distances);

        // ── Stage 5: Cauchy-weighted edges → CsrGraph ─────────────────────
        var edges = new Edge[pairs.Count];
        double bwSq = bandwidth * bandwidth;
        for (int e = 0; e < pairs.Count; e++)
        {
            var (i, j) = pairs[e];
            double d = distances[i, j];
            double coupling = 1.0 / (1.0 + (d * d) / bwSq);
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
                    OneNearestNeighborMaxOverMedian = maxOverMedian,
                    CsrEntries = graph.Targets.Length,
                    Kernel = "Cauchy",
                },
                Schedule = new
                {
                    Temperatures = temperatures,
                    NumReplicas,
                    TaskCount   = tasks.Count,
                    SweepBurnIn = BurnInCycles,
                    SweepCycles = CyclesPerTask,
                    Q,
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

    // ── Inline graph construction ─────────────────────────────────────────

    private static double[,] ComputePairwiseDistances(double[][] points)
    {
        int n = points.Length;
        var d = new double[n, n];
        Parallel.For(0, n, i =>
        {
            for (int j = i + 1; j < n; j++)
            {
                double dist = Poincare.Distance(points[i], points[j]);
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

    private static (double bandwidth, double maxOverMedian) EstimateBandwidthAndDiagnostic(double[,] distances)
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
        double median = n % 2 == 0
            ? 0.5 * (oneNn[n / 2 - 1] + oneNn[n / 2])
            : oneNn[n / 2];
        double max = oneNn[n - 1];
        double ratio = median > 0 ? max / median : 0;
        return (median, ratio);
    }
}
