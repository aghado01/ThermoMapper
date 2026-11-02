using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Clustering.Dendrograms;
using Clustering.Evaluation.External;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs;
using Graphs.Primitives;
using Maths.Rng;
using UserRepl;
using UserRepl.Commands;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// Landsat's closure gap probed: how much of it a partition reaches, on more than one seed and more than one graph,
/// under the protocol pre-registered in <c>issues/spc/bayesian-clustering-design.md</c> ("Probe of Landsat's case"). A
/// measurement harness, inert unless <c>CLOSURE_GAP_PROBE=&lt;output dir&gt;</c> is set; it asserts nothing.
/// <c>CLOSURE_GAP_THREADS</c> caps the replicas run at once (default 4). It runs at below-normal priority, works from
/// the hot end down, appends each temperature's row and partitions as it finishes, and skips finished rows on a rerun.
/// </summary>
public sealed class ClosureGapProbeHarness
{
    private const string EnvOut = "CLOSURE_GAP_PROBE";
    private const string EnvThreads = "CLOSURE_GAP_THREADS";
    private const int Q = 20, Replicas = 8, BurnIn = 1000, Cycles = 2000, KeepEvery = 10, GridPoints = 65;
    private const int SearchOrders = 4, SetCount = 6, LandsatIndex = 4;
    private const double Theta = 0.5;

    // Where to stop going colder: once the closure's ARI has reached the first, three levels in a row under the second.
    private const double SignalAri = 0.05, FadedAri = 0.02;
    private const int FadedLevels = 3;

    private const string Header =
        "arm\tmaster\tt\tn\tedges\tari\t" +
        "c_clusters\tc_ari\tc_loss\trule_loss\tgap\tcross_gap\tfalse_merges\tfalse_splits\t" +
        "s_clusters\ts_ari\ts_start\tsaving\tcross_saving\t" +
        "same_pairs\tc_joined\tc_joined_same\tr_linked\tr_linked_same\ts_joined\ts_joined_same\t" +
        "ms_sample\tms_reduce\tms_search";

    [Theory]
    [InlineData("k05", 20261102)]
    [InlineData("k05", 20261103)]
    [InlineData("k05", 20261104)]
    [InlineData("k10", 20261102)]
    [InlineData("k20", 20261102)]
    [InlineData("cauchy", 20261102)]
    public void Measure(string arm, int master)
    {
        string? outDir = Environment.GetEnvironmentVariable(EnvOut);
        if (string.IsNullOrEmpty(outDir)) return;

        Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.BelowNormal;
        int threads = int.TryParse(Environment.GetEnvironmentVariable(EnvThreads), out int t) ? t : 4;

        Directory.CreateDirectory(outDir);
        string rowsPath = Path.Combine(outDir, "landsat-probe.tsv");
        string labelsPath = Path.Combine(outDir, "landsat-probe-labels.tsv");
        if (!File.Exists(rowsPath)) File.WriteAllText(rowsPath, Header + Environment.NewLine);

        Dictionary<string, double> finished = File.ReadLines(rowsPath).Skip(1)
            .Select(line => line.Split('\t'))
            .Where(cells => cells[0] == arm && cells[1] == Text(master))
            .ToDictionary(cells => cells[2], cells => double.Parse(cells[5], CultureInfo.InvariantCulture));

        (double[][] features, int[] labels, _) = ClosureGapHarness.Load("landsat");
        GraphBuildResult build = SpcGraphBuilder.BuildResult(features, Config(arm), null);
        CsrGraph graph = build.Graph;
        long edges = graph.Targets.Length / 2;

        (double lo, double hi, _) = SpcScheduleHelpers.AutoBracket(graph, Q, build.Resolved);
        double[] grid = SpcScheduleHelpers.LogSpaceGrid(lo, hi, GridPoints);

        // The same derivation as ClosureGapHarness, so k05 on its master seed draws the measurement's own chains.
        int setSeed = SeedTree.Derive(master, SetCount)[LandsatIndex];
        int[] seeds = SeedTree.Derive(setSeed, grid.Length * Replicas);

        bool signal = false;
        int faded = 0;
        for (int g = grid.Length - 1; g >= 0 && faded < FadedLevels; g--)
        {
            double temperature = grid[g];
            if (!finished.TryGetValue(Text(temperature), out double ari))
            {
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

                string row = Reduce(arm, master, temperature, graph, edges, frames, labels, seeds[start], labelsPath, msSample, out ari);
                File.AppendAllText(rowsPath, row + Environment.NewLine);
            }

            signal |= ari >= SignalAri;
            faded = signal && ari < FadedAri ? faded + 1 : 0;
        }
    }

    /// <summary>Every graph behind the measurement and the probe, as the compiler resolved it. No sampling.</summary>
    [Fact]
    public void ResolvedGraphs()
    {
        string? outDir = Environment.GetEnvironmentVariable(EnvOut);
        if (string.IsNullOrEmpty(outDir)) return;
        Directory.CreateDirectory(outDir);

        var lines = new List<string>();
        foreach (string set in new[] { "blobs", "dumbbell", "iris", "toy", "landsat", "isolet" })
        {
            (double[][] features, _, int k) = ClosureGapHarness.Load(set);
            lines.Add(Describe($"{set} K={k}", SpcGraphBuilder.BuildResult(features, BwdParityTests.BuildGraphConfig(k, 0), null)));
        }
        (double[][] landsat, _, _) = ClosureGapHarness.Load("landsat");
        foreach (string arm in new[] { "k05", "k10", "k20", "cauchy" })
            lines.Add(Describe($"arm {arm}", SpcGraphBuilder.BuildResult(landsat, Config(arm), null)));
        File.WriteAllLines(Path.Combine(outDir, "graphs.txt"), lines);
    }

    private static string Describe(string name, GraphBuildResult build)
    {
        (double lo, double hi, string basis) = SpcScheduleHelpers.AutoBracket(build.Graph, Q, build.Resolved);
        return $"{name}: n={build.Graph.NodeCount} edges={build.Graph.Targets.Length / 2} q={Q} bracket=[{Text(lo)}, {Text(hi)}] basis={basis}"
            + Environment.NewLine + $"  {build.Resolved}";
    }

    private static GraphCompilerConfig Config(string arm)
    {
        switch (arm)
        {
            case "k05": return BwdParityTests.BuildGraphConfig(5, 0);
            case "k10": return BwdParityTests.BuildGraphConfig(10, 0);
            case "k20": return BwdParityTests.BuildGraphConfig(20, 0);
            case "cauchy":
            {
                GraphCompilerConfig parity = BwdParityTests.BuildGraphConfig(5, 0);
                return parity with { Projection = (CouplingProjection)parity.Projection with { Kernel = new Graphs.Coupling.Cauchy(0) } };
            }
            default: throw new ArgumentException($"Unknown arm '{arm}'.");
        }
    }

    private static string Reduce(
        string arm, int master, double temperature, CsrGraph graph, long edges, Accumulator[] frames, int[] labels,
        int seed, string labelsPath, long msSample, out double ari)
    {
        var clock = Stopwatch.StartNew();
        double cut = SweepClosureGap.CoMembershipTheta(Theta, Q);
        var scorer = new AdjustedRandIndex();
        ari = scorer.Evaluate(AffinityThreshold.Connect(graph, Column(frames), cut).Labels, labels);

        int n = graph.NodeCount;
        DenseCoMembership first = DenseCoMembership.FromPartitions(n, frames.Where((_, r) => r % 2 == 0).Select(f => f.KeptPartitions!).ToArray());
        DenseCoMembership second = DenseCoMembership.FromPartitions(n, frames.Where((_, r) => r % 2 == 1).Select(f => f.KeptPartitions!).ToArray());
        DenseCoMembership all = DenseCoMembership.Pool(first, second);

        int[] closure = AffinityThreshold.Connect(graph, all.EdgeColumn(graph), cut).Labels;
        ClosureGapResult p = ClosureGap.OverPairs(graph, all, cut);
        double crossGap = 0.5 * (ClosureGap.OverPairs(graph, first, cut, second).Gap + ClosureGap.OverPairs(graph, second, cut, first).Gap);
        long msReduce = clock.ElapsedMilliseconds;

        // Each half decides, and the other half evaluates, as the cross-fit gap does.
        clock.Restart();
        int[] searchSeeds = SeedTree.Derive(seed, 3);
        int[] closureFirst = AffinityThreshold.Connect(graph, first.EdgeColumn(graph), cut).Labels;
        int[] closureSecond = AffinityThreshold.Connect(graph, second.EdgeColumn(graph), cut).Labels;
        PairLossSearch.Found onAll = PairLossSearch.Best(all, cut, closure, SearchOrders, searchSeeds[0]);
        PairLossSearch.Found onFirst = PairLossSearch.Best(first, cut, closureFirst, SearchOrders, searchSeeds[1]);
        PairLossSearch.Found onSecond = PairLossSearch.Best(second, cut, closureSecond, SearchOrders, searchSeeds[2]);
        long msSearch = clock.ElapsedMilliseconds;

        clock.Restart();
        double saving = PairLossSearch.Score(all, onAll.Labels, cut) - PairLossSearch.Score(all, closure, cut);
        double crossSaving = 0.5 * (PairLossSearch.Score(second, onFirst.Labels, cut) - PairLossSearch.Score(second, closureFirst, cut)
                                  + PairLossSearch.Score(first, onSecond.Labels, cut) - PairLossSearch.Score(first, closureSecond, cut));

        long same = 0, cJoined = 0, cSame = 0, rLinked = 0, rSame = 0, sJoined = 0, sSame = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
            {
                bool together = labels[i] == labels[j];
                if (together) same++;
                if (closure[i] == closure[j])
                {
                    cJoined++;
                    if (together) cSame++;
                }
                if (all.At(i, j) > cut)
                {
                    rLinked++;
                    if (together) rSame++;
                }
                if (onAll.Labels[i] == onAll.Labels[j])
                {
                    sJoined++;
                    if (together) sSame++;
                }
            }

        string prefix = $"{arm}\t{Text(master)}\t{Text(temperature)}";
        File.AppendAllLines(labelsPath,
        [
            $"{prefix}\tclosure\t{string.Join(',', closure)}",
            $"{prefix}\tsearch\t{string.Join(',', onAll.Labels)}",
        ]);

        return string.Join('\t',
            arm, Text(master), Text(temperature), Text(n), Text(edges), Text(ari),
            Text(p.ClusterCount), Text(scorer.Evaluate(closure, labels)), Text(p.ClosureLoss), Text(p.PairRuleLoss), Text(p.Gap), Text(crossGap), Text(p.FalseMerges), Text(p.FalseSplits),
            Text(onAll.Labels.Max() + 1), Text(scorer.Evaluate(onAll.Labels, labels)), onAll.Start, Text(saving), Text(crossSaving),
            Text(same), Text(cJoined), Text(cSame), Text(rLinked), Text(rSame), Text(sJoined), Text(sSame),
            Text(msSample), Text(msReduce + clock.ElapsedMilliseconds), Text(msSearch));
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

    private static string Text(double value) => value.ToString("R", CultureInfo.InvariantCulture);

    private static string Text(long value) => value.ToString(CultureInfo.InvariantCulture);
}
