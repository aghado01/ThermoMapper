using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Threading.Tasks;
using Clustering.Dendrograms;
using Clustering.Evaluation.External;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs.Models.Potts;
using Graphs.Primitives;
using Maths.Rng;
using UserRepl;
using UserRepl.Commands;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The closure gap on data, under the protocol pre-registered in <c>issues/spc/bayesian-clustering-design.md</c>
/// ("Measurement of the closure gap"). A measurement harness, inert unless <c>CLOSURE_GAP_MEASURE=&lt;output dir&gt;</c>
/// is set; it asserts nothing. <c>CLOSURE_GAP_THREADS</c> caps the replicas run at once (default 4). It runs at
/// below-normal priority, appends each temperature's row to <c>closure-gap.tsv</c> as it finishes, and skips finished
/// rows on a rerun.
/// </summary>
public sealed class ClosureGapHarness
{
    private const string EnvOut = "CLOSURE_GAP_MEASURE";
    private const string EnvThreads = "CLOSURE_GAP_THREADS";
    private const int Q = 20, Replicas = 8, BurnIn = 1000, Cycles = 2000, KeepEvery = 10, GridPoints = 65, MasterSeed = 20261102;
    private const double Theta = 0.5;

    private static readonly string[] Sets = ["blobs", "dumbbell", "iris", "toy", "landsat", "isolet"];

    private const string Header =
        "set\tt\tn\tedges\tari\te_clusters\te_loss\te_rule_loss\te_gap\te_cross\te_merges\te_integrated\t" +
        "p_clusters\tp_loss\tp_rule_loss\tp_gap\tp_cross\tp_merges\tp_merge_cost\tp_splits\tp_split_cost\tp_integrated\t" +
        "contradicted\tlabels_side_with_rule\tsame_pairs\tclosure_joined\tclosure_joined_same\trule_linked\trule_linked_same\t" +
        "ms_sample\tms_reduce";

    [Theory]
    [InlineData("blobs")]
    [InlineData("dumbbell")]
    [InlineData("iris")]
    [InlineData("toy")]
    [InlineData("landsat")]
    [InlineData("isolet")]
    public void Measure(string set)
    {
        string? outDir = Environment.GetEnvironmentVariable(EnvOut);
        if (string.IsNullOrEmpty(outDir)) return;

        Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.BelowNormal;
        int threads = int.TryParse(Environment.GetEnvironmentVariable(EnvThreads), out int t) ? t : 4;

        Directory.CreateDirectory(outDir);
        string path = Path.Combine(outDir, "closure-gap.tsv");
        lock (Sets)
            if (!File.Exists(path)) File.WriteAllText(path, Header + Environment.NewLine);

        (double lo, double hi) = BwdPottsCriticalEstimate.Bracket(Q);
        double[] grid = SpcScheduleHelpers.LogSpaceGrid(lo, hi, GridPoints);
        HashSet<string> finished = File.ReadLines(path).Skip(1)
            .Select(line => line.Split('\t'))
            .Where(cells => cells[0] == set)
            .Select(cells => cells[1])
            .ToHashSet();
        if (grid.All(temperature => finished.Contains(Text(temperature)))) return;

        (double[][] features, int[] labels, int k) = Load(set);
        CsrGraph graph = SpcGraphBuilder.BuildResult(features, BwdParityTests.BuildGraphConfig(k, 0), null).Graph;
        long edges = graph.Targets.Length / 2;

        int setSeed = SeedTree.Derive(MasterSeed, Sets.Length)[Array.IndexOf(Sets, set)];
        int[] seeds = SeedTree.Derive(setSeed, grid.Length * Replicas);

        for (int g = 0; g < grid.Length; g++)
        {
            double temperature = grid[g];
            if (finished.Contains(Text(temperature))) continue;

            var clock = Stopwatch.StartNew();
            var frames = new Accumulator[Replicas];
            int start = g * Replicas;
            Parallel.For(0, Replicas, new ParallelOptions { MaxDegreeOfParallelism = threads }, r =>
            {
                frames[r] = SwRunner.Run(new SwRunSpec
                {
                    Graph = graph, Temperature = temperature, Q = Q, Seed = seeds[start + r], ReplicaIndex = r,
                    Budget = new RunBudget(BurnIn, Cycles),
                    Accumulation = new AccumulationSpec { CoMembership = true, PartitionEvery = KeepEvery },
                }).Accumulator;
            });
            long msSample = clock.ElapsedMilliseconds;

            clock.Restart();
            string row = Reduce(set, temperature, graph, edges, frames, labels, msSample, clock);
            File.AppendAllText(path, row + Environment.NewLine);
        }
    }

    private static string Reduce(
        string set, double temperature, CsrGraph graph, long edges, Accumulator[] frames, int[] labels, long msSample, Stopwatch clock)
    {
        ClosureGapCurve onEdges = SweepClosureGap.OverEdges(graph, frames, Theta);
        double cut = onEdges.Theta;
        ClosureGapResult e = onEdges.Results[0];

        double[] column = Column(frames);
        double ari = new AdjustedRandIndex().Evaluate(AffinityThreshold.Connect(graph, column, cut).Labels, labels);

        int n = graph.NodeCount;
        DenseCoMembership first = DenseCoMembership.FromPartitions(n, frames.Where((_, r) => r % 2 == 0).Select(f => f.KeptPartitions!).ToArray());
        DenseCoMembership second = DenseCoMembership.FromPartitions(n, frames.Where((_, r) => r % 2 == 1).Select(f => f.KeptPartitions!).ToArray());
        DenseCoMembership all = DenseCoMembership.Pool(first, second);

        ClosureGapResult p = ClosureGap.OverPairs(graph, all, cut);
        double cross = 0.5 * (ClosureGap.OverPairs(graph, first, cut, second).Gap + ClosureGap.OverPairs(graph, second, cut, first).Gap);
        double integrated = ClosureGap.IntegratedOverPairs(graph, all);

        // Which side the labels take, pair by pair.
        int[] closure = AffinityThreshold.Connect(graph, all.EdgeColumn(graph), cut).Labels;
        long contradicted = 0, sideWithRule = 0, samePairs = 0, joinedCount = 0, joinedSame = 0, linkedCount = 0, linkedSame = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
            {
                bool joined = closure[i] == closure[j];
                bool linked = all.At(i, j) > cut;
                bool same = labels[i] == labels[j];
                if (same) samePairs++;
                if (joined)
                {
                    joinedCount++;
                    if (same) joinedSame++;
                }
                if (linked)
                {
                    linkedCount++;
                    if (same) linkedSame++;
                }
                if (joined != linked)
                {
                    contradicted++;
                    if (same == linked) sideWithRule++;
                }
            }

        return string.Join('\t',
            set, Text(temperature), Text(n), Text(edges), Text(ari),
            Text(e.ClusterCount), Text(e.ClosureLoss), Text(e.PairRuleLoss), Text(e.Gap), Text(onEdges.CrossFitGap[0]), Text(e.FalseMerges), Text(onEdges.Integrated[0]),
            Text(p.ClusterCount), Text(p.ClosureLoss), Text(p.PairRuleLoss), Text(p.Gap), Text(cross), Text(p.FalseMerges), Text(p.FalseMergeCost), Text(p.FalseSplits), Text(p.FalseSplitCost), Text(integrated),
            Text(contradicted), Text(sideWithRule), Text(samePairs), Text(joinedCount), Text(joinedSame), Text(linkedCount), Text(linkedSame),
            Text(msSample), Text(clock.ElapsedMilliseconds));
    }

    private static double[] Column(Accumulator[] frames)
    {
        var column = new double[frames[0].CoMembershipCount!.Length];
        long draws = frames.Sum(f => (long)f.DrawCount);
        foreach (Accumulator frame in frames)
            for (int e = 0; e < column.Length; e++)
                column[e] += (double)frame.CoMembershipCount![e] / draws;
        return column;
    }

    // ── Data ────────────────────────────────────────────────────────────────

    private static (double[][] Features, int[] Labels, int K) Load(string set)
    {
        switch (set)
        {
            case "blobs": return Blobs(bridge: false);
            case "dumbbell": return Blobs(bridge: true);
            case "iris":
            {
                var data = SpcUserSession.FromCsv(Datasets.Path("iris.csv"), null, false, ',');
                return (data.Features, data.Labels!, 5);
            }
            case "landsat":
            {
                var data = SpcUserSession.FromCsv(Datasets.Path("landsat.csv"), null, false, ',');
                return (data.Features, data.Labels!, 5);
            }
            case "toy":
            {
                var data = Synthetic.Euclidean.Bwd1995Toy.Generate(seed: 42);
                return (data.Features, data.Labels, 10);
            }
            case "isolet": return Isolet();
            default: throw new ArgumentException($"Unknown set '{set}'.");
        }
    }

    private static (double[][] Features, int[] Labels, int K) Blobs(bool bridge)
    {
        const int perBlob = 200, bridgePoints = 24;
        const double centre = 12.0;
        var rng = new Xoshiro256PlusPlus(MasterSeed);
        var features = new List<double[]>();
        var labels = new List<int>();
        for (int blob = 0; blob < 2; blob++)
            for (int i = 0; i < perBlob; i++)
            {
                features.Add([(blob == 0 ? -centre : centre) + Gaussian(rng), Gaussian(rng)]);
                labels.Add(blob);
            }
        if (bridge)
            for (int i = 0; i < bridgePoints; i++)
            {
                double x = -9.0 + 18.0 * i / (bridgePoints - 1);
                features.Add([x, 0.0]);
                labels.Add(x < 0.0 ? 0 : 1);
            }
        return (features.ToArray(), labels.ToArray(), 10);
    }

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static (double[][] Features, int[] Labels, int K) Isolet()
    {
        using var file = File.OpenRead(Datasets.Path("isolet.csv.gz"));
        using var gz = new GZipStream(file, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var features = new List<double[]>(8000);
        var labels = new List<int>(8000);
        while (reader.ReadLine() is { } line)
        {
            if (string.IsNullOrWhiteSpace(line)) continue;
            string[] cells = line.Split(',');
            var row = new double[cells.Length - 1];
            for (int d = 0; d < row.Length; d++)
                row[d] = double.Parse(cells[d], CultureInfo.InvariantCulture);
            features.Add(row);
            labels.Add((int)Math.Round(double.Parse(cells[^1], CultureInfo.InvariantCulture)));
        }
        return (features.ToArray(), labels.ToArray(), 10);
    }

    private static string Text(double value) => value.ToString("R", CultureInfo.InvariantCulture);

    private static string Text(long value) => value.ToString(CultureInfo.InvariantCulture);
}
