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
using Clustering.Primitives;
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

    private static readonly double[] Thetas = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99];

    private const string ThetaHeader =
        "arm\tmaster\tt\ttheta\tclusters\tsingletons\tari\tjoined\tjoined_same\tlargest\tlargest_purity\t" +
        "c50_clusters\tc50_covered\tc50_purity\tc20_clusters\tc20_covered\tc20_purity\tms_sample";

    // The thin-cut test's settings: no cut of at most K edges between two sides of at least M points.
    private static readonly int[] ThinCutEdges = [1, 2, 3, 4, 6], ThinCutSizes = [1, 20, 50, 100, 200, 400];
    private const int MaxThinCut = 6;
    private static readonly string ThinCutHeader = ThetaHeader.Replace("ms_sample", "ms_tree") + "\tk\tm\tbroken\tc50_classes";

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

        (GraphBuildResult build, int[] labels, double[] grid, int[] seeds) = Prepare(arm, master);
        CsrGraph graph = build.Graph;
        long edges = graph.Targets.Length / 2;

        bool signal = false;
        int faded = 0;
        for (int g = grid.Length - 1; g >= 0 && faded < FadedLevels; g--)
        {
            double temperature = grid[g];
            if (!finished.TryGetValue(Text(temperature), out double ari))
            {
                var clock = Stopwatch.StartNew();
                Accumulator[] frames = Sample(graph, temperature, g * Replicas, seeds, threads);
                long msSample = clock.ElapsedMilliseconds;

                string row = Reduce(arm, master, temperature, graph, edges, frames, labels, seeds[g * Replicas], labelsPath, msSample, out ari);
                File.AppendAllText(rowsPath, row + Environment.NewLine);
            }

            signal |= ari >= SignalAri;
            faded = signal && ari < FadedAri ? faded + 1 : 0;
        }
    }

    /// <summary>
    /// Closure alone at each θ in <see cref="Thetas"/>, on the probe's own chains and over the whole grid, under the
    /// protocol in the design note ("Closure at a higher θ"). Each level's edge column is kept beside the edge list, so a
    /// later question about edges or θ needs no new sampling.
    /// </summary>
    [Theory]
    [InlineData("k05", 20261102)]
    [InlineData("k05", 20261103)]
    [InlineData("k05", 20261104)]
    [InlineData("k10", 20261102)]
    [InlineData("k20", 20261102)]
    [InlineData("cauchy", 20261102)]
    public void MeasureTheta(string arm, int master)
    {
        string? outDir = Environment.GetEnvironmentVariable(EnvOut);
        if (string.IsNullOrEmpty(outDir)) return;

        Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.BelowNormal;
        int threads = int.TryParse(Environment.GetEnvironmentVariable(EnvThreads), out int t) ? t : 4;

        string columnsDir = Path.Combine(outDir, "landsat-theta-columns");
        Directory.CreateDirectory(columnsDir);
        string rowsPath = Path.Combine(outDir, "landsat-theta.tsv");
        if (!File.Exists(rowsPath)) File.WriteAllText(rowsPath, ThetaHeader + Environment.NewLine);

        HashSet<string> finished = File.ReadLines(rowsPath).Skip(1)
            .Select(line => line.Split('\t'))
            .Where(cells => cells[0] == arm && cells[1] == Text(master))
            .Select(cells => cells[2])
            .ToHashSet();

        (GraphBuildResult build, int[] labels, double[] grid, int[] seeds) = Prepare(arm, master);
        CsrGraph graph = build.Graph;

        var slots = new List<int>();
        var edgeLines = new List<string>();
        HashSet<(int, int)> beforeRepair = EdgeSet(build.PreRepairGraph ?? graph);
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
        {
            slots.Add(edge.Slot);
            edgeLines.Add($"{edge.Source}\t{edge.Target}\t{(beforeRepair.Contains((edge.Source, edge.Target)) ? 0 : 1)}");
        }
        File.WriteAllLines(Path.Combine(columnsDir, $"{arm}-edges.tsv"), edgeLines);

        var scorer = new AdjustedRandIndex();
        for (int g = grid.Length - 1; g >= 0; g--)
        {
            double temperature = grid[g];
            if (finished.Contains(Text(temperature))) continue;

            var clock = Stopwatch.StartNew();
            double[] column = Column(Sample(graph, temperature, g * Replicas, seeds, threads));
            long msSample = clock.ElapsedMilliseconds;

            using (var writer = new BinaryWriter(File.Open(Path.Combine(columnsDir, $"{arm}-{Text(master)}.bin"), FileMode.Append)))
            {
                writer.Write(temperature);
                foreach (int slot in slots) writer.Write(column[slot]);
            }

            File.AppendAllLines(rowsPath, Thetas.Select(theta =>
                ThetaRow(arm, master, temperature, theta, AffinityThreshold.Connect(graph, column, SweepClosureGap.CoMembershipTheta(theta, Q)).Labels, labels, scorer, msSample)));
        }
    }

    /// <summary>
    /// The thin-cut closure at every setting, on the columns <see cref="MeasureTheta"/> kept, under the protocol in the
    /// design note ("A closure that won't cross a thin cut"). No sampling. Going from the hot end down, a θ's rows are
    /// appended once all its levels finish, and a rerun skips finished θ.
    /// </summary>
    [Theory]
    [InlineData("k05", 20261102)]
    [InlineData("k05", 20261103)]
    [InlineData("k05", 20261104)]
    [InlineData("k10", 20261102)]
    [InlineData("k20", 20261102)]
    [InlineData("cauchy", 20261102)]
    public void MeasureThinCut(string arm, int master)
    {
        string? outDir = Environment.GetEnvironmentVariable(EnvOut);
        if (string.IsNullOrEmpty(outDir)) return;

        Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.BelowNormal;
        int threads = int.TryParse(Environment.GetEnvironmentVariable(EnvThreads), out int t) ? t : 4;

        string columnsDir = Path.Combine(outDir, "landsat-theta-columns");
        string[][] edges = File.ReadLines(Path.Combine(columnsDir, $"{arm}-edges.tsv")).Select(line => line.Split('\t')).ToArray();
        int[] sources = edges.Select(cells => int.Parse(cells[0], CultureInfo.InvariantCulture)).ToArray();
        int[] targets = edges.Select(cells => int.Parse(cells[1], CultureInfo.InvariantCulture)).ToArray();
        (double Temperature, double[] Column)[] levels = KeptColumns(Path.Combine(columnsDir, $"{arm}-{Text(master)}.bin"), edges.Length);
        (_, int[] labels, _) = ClosureGapHarness.Load("landsat");

        string rowsPath = Path.Combine(outDir, "landsat-thincut.tsv");
        if (!File.Exists(rowsPath)) File.WriteAllText(rowsPath, ThinCutHeader + Environment.NewLine);
        HashSet<string> finished = File.ReadLines(rowsPath).Skip(1)
            .Select(line => line.Split('\t'))
            .Where(cells => cells[0] == arm && cells[1] == Text(master))
            .Select(cells => cells[3])
            .ToHashSet();

        // K = 10 and 20 have no search target and the densest trees, so they read θ = ½ alone.
        double[] thetas = arm is "k10" or "k20" ? [Theta] : Thetas;
        (int K, int M)[] settings = [(0, 1), .. ThinCutEdges.SelectMany(k => ThinCutSizes.Select(m => (k, m)))];
        var writing = new object();

        Parallel.ForEach(thetas.Where(theta => !finished.Contains(Text(theta))), new ParallelOptions { MaxDegreeOfParallelism = threads }, theta =>
        {
            var scorer = new AdjustedRandIndex();
            var rows = new List<string>();
            var hotter = new Dictionary<(int, int), int[]>();
            double cut = SweepClosureGap.CoMembershipTheta(theta, Q);

            foreach ((double temperature, double[] column) in levels)
            {
                var clock = Stopwatch.StartNew();
                ThinCutClosure tree = ThinCutClosure.Above(labels.Length, sources, targets, column, cut, MaxThinCut);
                long msTree = clock.ElapsedMilliseconds;

                foreach ((int k, int m) in settings)
                {
                    int[] partition = tree.Partition(k, m);
                    long broken = hotter.TryGetValue((k, m), out int[]? above) ? Broken(above, partition) : 0;
                    hotter[(k, m)] = partition;
                    rows.Add(string.Join('\t', ThetaRow(arm, master, temperature, theta, partition, labels, scorer, msTree),
                        Text(k), Text(m), Text(broken), Text(ClassesHeld(partition, labels))));
                }
            }

            lock (writing) File.AppendAllLines(rowsPath, rows);
        });
    }

    private static (double Temperature, double[] Column)[] KeptColumns(string path, int edges)
    {
        var levels = new List<(double, double[])>();
        using var reader = new BinaryReader(File.OpenRead(path));
        while (reader.BaseStream.Position < reader.BaseStream.Length)
        {
            double temperature = reader.ReadDouble();
            var column = new double[edges];
            for (int e = 0; e < edges; e++) column[e] = reader.ReadDouble();
            levels.Add((temperature, column));
        }
        return levels.OrderByDescending(level => level.Item1).ToArray();
    }

    // Pairs joined at the hotter level and apart at the colder, which exact closure can't have.
    private static long Broken(int[] hotter, int[] colder)
    {
        long Joined(long size) => size * (size - 1) / 2;
        return hotter.GroupBy(label => label).Sum(c => Joined(c.Count()))
             - Enumerable.Range(0, hotter.Length).GroupBy(i => (hotter[i], colder[i])).Sum(c => Joined(c.Count()));
    }

    // How many classes lead some cluster of at least 50.
    private static long ClassesHeld(int[] partition, int[] labels) =>
        Enumerable.Range(0, partition.Length).GroupBy(i => partition[i]).Where(c => c.Count() >= 50)
            .Select(c => c.GroupBy(i => labels[i]).MaxBy(g => g.Count())!.Key).Distinct().Count();

    private static HashSet<(int, int)> EdgeSet(CsrGraph graph)
    {
        var set = new HashSet<(int, int)>();
        foreach (UndirectedEdge edge in graph.UndirectedEdges()) set.Add((edge.Source, edge.Target));
        return set;
    }

    private static string ThetaRow(string arm, int master, double temperature, double theta, int[] closure, int[] labels, AdjustedRandIndex scorer, long msSample)
    {
        Dictionary<int, long> sizes = closure.GroupBy(label => label).ToDictionary(c => c.Key, c => (long)c.Count());
        long Joined(long size) => size * (size - 1) / 2;

        // Purity scores assigned points only, so what isn't read goes unassigned.
        double PurityOver(Func<int, bool> read) =>
            Purity.Compute(closure.Select(label => read(label) ? label : Assignment.Unassigned).ToArray(), labels);
        string Covered(long minimum) => string.Join('\t',
            Text((long)sizes.Values.Count(size => size >= minimum)),
            Text(sizes.Values.Where(size => size >= minimum).Sum()),
            Text(PurityOver(label => sizes[label] >= minimum)));

        int largest = sizes.MaxBy(pair => pair.Value).Key;
        long joinedSame = Enumerable.Range(0, closure.Length).GroupBy(i => (closure[i], labels[i])).Sum(pair => Joined(pair.Count()));
        return string.Join('\t',
            arm, Text(master), Text(temperature), Text(theta),
            Text((long)sizes.Count), Text((long)sizes.Values.Count(size => size == 1)), Text(scorer.Evaluate(closure, labels)),
            Text(sizes.Values.Sum(Joined)), Text(joinedSame),
            Text(sizes[largest]), Text(PurityOver(label => label == largest)), Covered(50), Covered(20), Text(msSample));
    }

    private static (GraphBuildResult Build, int[] Labels, double[] Grid, int[] Seeds) Prepare(string arm, int master)
    {
        (double[][] features, int[] labels, _) = ClosureGapHarness.Load("landsat");
        GraphBuildResult build = SpcGraphBuilder.BuildResult(features, Config(arm), null);

        (double lo, double hi, _) = SpcScheduleHelpers.AutoBracket(build.Graph, Q, build.Resolved);
        double[] grid = SpcScheduleHelpers.LogSpaceGrid(lo, hi, GridPoints);

        // The same derivation as ClosureGapHarness, so k05 on its master seed draws the measurement's own chains.
        int setSeed = SeedTree.Derive(master, SetCount)[LandsatIndex];
        return (build, labels, grid, SeedTree.Derive(setSeed, grid.Length * Replicas));
    }

    private static Accumulator[] Sample(CsrGraph graph, double temperature, int start, int[] seeds, int threads)
    {
        var frames = new Accumulator[Replicas];
        Parallel.For(0, Replicas, new ParallelOptions { MaxDegreeOfParallelism = threads }, r =>
        {
            frames[r] = SwRunner.Run(new SwRunSpec
            {
                Graph = graph, Temperature = temperature, Q = Q, Seed = seeds[start + r], ReplicaIndex = r,
                Budget = new RunBudget(BurnIn, Cycles),
                Accumulation = new AccumulationSpec { CoMembership = true, PartitionEvery = KeepEvery },
            }).Accumulator;
        });
        return frames;
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
