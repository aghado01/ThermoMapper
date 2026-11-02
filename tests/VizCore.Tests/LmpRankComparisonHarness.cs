using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using Clustering.Dendrograms;
using Clustering.Graphical.SPC.Partitions.Hierarchical;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Profiling.Signals;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Clustering.Primitives;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Models.Potts;
using Graphs.Pipeline;
using Graphs.Primitives;
using Maths.Rng;
using UserRepl;
using Xunit;
using Xunit.Abstractions;

namespace VizCore.Tests;

/// <summary>
/// The strict-versus-inclusive Local Mutual Proximity comparison (graph-compiler-brief, "LMP rank comparison").
/// A measurement harness, inert unless <c>LMP_RANK_COMPARISON=&lt;output dir&gt;</c> is set; it asserts nothing.
/// <c>LMP_RANK_ONLY=&lt;substring&gt;</c> keeps the cells whose id contains it. Parts, each its own fact:
/// <c>Graphs</c> (compile only: connectivity, bridge units, merge-order inversion, chaining edges),
/// <c>Grid</c> (one seed, SW sweep + stage scoring + HierarchyEom), <c>Seeds</c> (the deciding cells on six
/// seeds), <c>ColdEdge</c> (BlattHierarchy 0.15, K = 10 with the grid's cold edge at 1e-4).
/// Rows go to <c>graphs.tsv</c> and <c>runs.tsv</c>. The protocol follows stopping-theory-primer §12.29: q = 20,
/// 141 temperatures over [1e-6, 3], 4 replicas x 1,000 cycles (200 burn-in), mutual kNN + MST-min, Cauchy
/// inner kernel with the mean-edge bandwidth, no gauge, no injected protect set.
/// </summary>
public sealed class LmpRankComparisonHarness
{
    private const string EnvOut = "LMP_RANK_COMPARISON";
    private const string EnvOnly = "LMP_RANK_ONLY";
    private const int Q = 20;
    private const int Points = 141;
    private const double ColdEdge = 1e-6, HotEdge = 3.0;
    private const int GridSeed = 42;

    private readonly ITestOutputHelper _out;
    public LmpRankComparisonHarness(ITestOutputHelper output) => _out = output;

    // ── Data ────────────────────────────────────────────────────────────────
    private sealed record Data(string Name, double[][] X, int[][] Levels, string[] LevelNames);

    private static readonly Dictionary<string, Data> DataCache = new();
    private static readonly Dictionary<string, GraphBuildResult> GraphCache = new();

    private static Data Load(string name)
    {
        lock (DataCache)
        {
            if (DataCache.TryGetValue(name, out Data? cached)) return cached;
            Data d = name switch
            {
                "iris" => FromCsv(name, "iris.csv"),
                "landsat" => FromCsv(name, "landsat.csv"),
                "blatt015" => FromBlatt(name, 0.15),
                "blatt030" => FromBlatt(name, 0.30),
                "aniso" => FromAniso(name),
                _ => throw new ArgumentException($"Unknown dataset '{name}'."),
            };
            return DataCache[name] = d;
        }
    }

    private static Data FromCsv(string name, string file)
    {
        var ds = SpcUserSession.FromCsv(Datasets.Path(file), null, false, ',');
        return new Data(name, ds.Features, new[] { ds.Labels! }, new[] { "k" });
    }

    private static Data FromBlatt(string name, double spread)
    {
        var ds = Synthetic.Euclidean.BlattHierarchy.Generate(pointsPerFine: 50, leafSpread: spread, seed: 42);
        return new Data(name, ds.Features, ds.LabelsByLevel!, new[] { "coarse", "medium", "fine" });
    }

    private static Data FromAniso(string name)
    {
        var ds = Synthetic.Euclidean.AnisotropicGaussian.Generate(
            pointsPerCluster: 400, meanSeparation: 4, anisotropyRatio: 1, seed: 55);
        return new Data(name, ds.Features, new[] { ds.Labels }, new[] { "k" });
    }

    // ── Graph ────────────────────────────────────────────────────────────────
    private sealed record Cell(string Dataset, int K, LmpRankRule Rank, string Kernel)
    {
        public string Id => $"{Dataset}|{Kernel}|K{K}|{Rank}";
    }

    private static GraphBuildResult Graph(Cell c)
    {
        lock (GraphCache)
        {
            if (GraphCache.TryGetValue(c.Id, out GraphBuildResult? cached)) return cached;
            Data d = Load(c.Dataset);
            var config = new GraphCompilerConfig
            {
                Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = c.K },
                Filter = new FilterConfig { Kind = FilterKind.MutualKnn },
                Repair = new RepairConfig { Kind = RepairKind.MstMin },
                Projection = new CouplingProjection
                {
                    Kernel = c.Kernel == "Gaussian" ? new Gaussian() : new Cauchy(),
                    LmpRescale = true,
                    BandwidthOverride = BandwidthStrategy.MeanEdgeDistance,
                    LmpRank = c.Rank,
                },
            };
            return GraphCache[c.Id] = SpcGraphBuilder.BuildResult(d.X, config, null);
        }
    }

    private static IEnumerable<Cell> GridCells()
    {
        foreach (string ds in new[] { "iris", "blatt015", "blatt030", "aniso", "landsat" })
            foreach (int k in new[] { 5, 10 })
                foreach (LmpRankRule rank in new[] { LmpRankRule.Strict, LmpRankRule.Inclusive })
                    yield return new Cell(ds, k, rank, "Cauchy");
        foreach (int k in new[] { 5, 10 })
            foreach (LmpRankRule rank in new[] { LmpRankRule.Strict, LmpRankRule.Inclusive })
                yield return new Cell("blatt015", k, rank, "Gaussian");
    }

    private static IEnumerable<Cell> DecidingCells()
    {
        foreach (LmpRankRule rank in new[] { LmpRankRule.Strict, LmpRankRule.Inclusive })
        {
            yield return new Cell("aniso", 5, rank, "Cauchy");
            yield return new Cell("blatt015", 10, rank, "Cauchy");
        }
    }

    private static IEnumerable<Cell> Only(IEnumerable<Cell> cells)
    {
        string? only = Environment.GetEnvironmentVariable(EnvOnly);
        return string.IsNullOrWhiteSpace(only) ? cells : cells.Where(c => c.Id.Contains(only, StringComparison.OrdinalIgnoreCase));
    }

    // ── Output ───────────────────────────────────────────────────────────────
    private static string? OutDir() => Environment.GetEnvironmentVariable(EnvOut) is { Length: > 0 } dir ? dir : null;

    private static string F(double v) => v.ToString("G6", CultureInfo.InvariantCulture);

    private void Emit(string dir, string file, List<(string Name, string Value)> row)
    {
        Directory.CreateDirectory(dir);
        string path = Path.Combine(dir, file);
        bool fresh = !File.Exists(path);
        using var w = new StreamWriter(path, append: true);
        if (fresh) w.WriteLine(string.Join('\t', row.Select(r => r.Name)));
        w.WriteLine(string.Join('\t', row.Select(r => r.Value)));
        _out.WriteLine(string.Join("  ", row.Where(r => r.Name is "part" or "cell" or "seed").Select(r => $"{r.Name}={r.Value}")));
    }

    // ── Graph-only measurements ───────────────────────────────────────────────
    [Fact]
    public void Graphs()
    {
        if (OutDir() is not string dir) return;
        foreach (Cell cell in Only(GridCells()))
            Emit(dir, "graphs.tsv", GraphRow(cell));
    }

    private static List<(string, string)> GraphRow(Cell cell)
    {
        Data d = Load(cell.Dataset);
        GraphBuildResult build = Graph(cell);
        CsrGraph g = build.Graph;
        int n = g.NodeCount;
        var row = new List<(string, string)>
        {
            ("cell", cell.Id), ("dataset", cell.Dataset), ("kernel", cell.Kernel), ("K", cell.K.ToString()), ("rank", cell.Rank.ToString()),
            ("nodes", n.ToString()),
        };

        // Edges, with weight, once each.
        var edges = new List<(int U, int V, double W)>();
        foreach (UndirectedEdge e in g.UndirectedEdges()) edges.Add((e.Source, e.Target, g.Weights[e.Slot]));
        int zero = edges.Count(e => !(e.W > 0.0));
        var coupledDegree = new int[n];
        foreach (var (u, v, w) in edges)
            if (w > 0.0) { coupledDegree[u]++; coupledDegree[v]++; }
        ResolvedConstruction r = build.Resolved!;
        row.Add(("edges", edges.Count.ToString()));
        row.Add(("zeroEdges", zero.ToString()));
        row.Add(("structComps", r.Components.ToString()));
        row.Add(("coupledComps", r.Coupling!.CoupledComponents.ToString()));
        row.Add(("repairEdges", r.RepairEdges.ToString()));
        row.Add(("degLe2Struct", Enumerable.Range(0, n).Count(i => g.Degree(i) <= 2).ToString()));
        row.Add(("degLe2Coupled", coupledDegree.Count(c => c <= 2).ToString()));
        row.Add(("isolatedCoupled", coupledDegree.Count(c => c == 0).ToString()));

        // Bridge units: repair edges against the interior (non-repair) edges at their endpoints.
        var repair = new HashSet<long>();
        if (build.PreRepairGraph is CsrGraph pre)
        {
            var preKeys = new HashSet<long>();
            foreach (UndirectedEdge e in pre.UndirectedEdges()) preKeys.Add(RepairEdgeSet.KeyOf(e.Source, e.Target, n));
            foreach (var (u, v, _) in edges)
                if (!preKeys.Contains(RepairEdgeSet.KeyOf(u, v, n))) repair.Add(RepairEdgeSet.KeyOf(u, v, n));
        }
        var interior = new List<double>[n];
        for (int i = 0; i < n; i++) interior[i] = new List<double>();
        foreach (var (u, v, w) in edges)
            if (!repair.Contains(RepairEdgeSet.KeyOf(u, v, n))) { interior[u].Add(w); interior[v].Add(w); }
        int[] fine = d.Levels[^1];
        int bridges = 0, minEither = 0, minBoth = 0, medEither = 0, medBoth = 0, noInterior = 0;
        int crossFine = 0, crossFineMinEither = 0, crossFineMedEither = 0;
        foreach (var (u, v, w) in edges)
        {
            if (!repair.Contains(RepairEdgeSet.KeyOf(u, v, n))) continue;
            bridges++;
            bool cross = fine[u] != fine[v];
            if (cross) crossFine++;
            if (interior[u].Count == 0 || interior[v].Count == 0) { noInterior++; continue; }
            bool gtMinU = w > interior[u].Min(), gtMinV = w > interior[v].Min();
            bool gtMedU = w > Median(interior[u]), gtMedV = w > Median(interior[v]);
            if (gtMinU || gtMinV) { minEither++; if (cross) crossFineMinEither++; }
            if (gtMinU && gtMinV) minBoth++;
            if (gtMedU || gtMedV) { medEither++; if (cross) crossFineMedEither++; }
            if (gtMedU && gtMedV) medBoth++;
        }
        row.Add(("bridges", bridges.ToString()));
        row.Add(("bridgesNoInterior", noInterior.ToString()));
        row.Add(("bridgeGtWeakestEither", minEither.ToString()));
        row.Add(("bridgeGtWeakestBoth", minBoth.ToString()));
        row.Add(("bridgeGtMedianEither", medEither.ToString()));
        row.Add(("bridgeGtMedianBoth", medBoth.ToString()));
        row.Add(("crossFineBridges", crossFine.ToString()));
        row.Add(("crossFineGtWeakestEither", crossFineMinEither.ToString()));
        row.Add(("crossFineGtMedianEither", crossFineMedEither.ToString()));

        // Single-linkage order per label level: does the strongest cross-class edge beat the weakest
        // within-class bottleneck? (an inversion: a link between classes before a class is whole).
        for (int li = 0; li < d.Levels.Length; li++)
        {
            string tag = d.Levels.Length == 1 ? "L" : d.LevelNames[li][..1].ToUpperInvariant();
            int[] y = Dense(d.Levels[li], out int k);
            (double wInt, double wCross, int chain, int crossEdges, double intraMedian, int crossLowDeg, int crossAboveIntraMedian) =
                MergeOrder(edges, y, k, n, coupledDegree);
            row.Add(($"{tag}_weakestClassBottleneck", F(wInt)));
            row.Add(($"{tag}_strongestCrossEdge", F(wCross)));
            row.Add(($"{tag}_inverted", (wCross > wInt).ToString()));
            row.Add(($"{tag}_crossEdges", crossEdges.ToString()));
            row.Add(($"{tag}_crossAboveBottleneck", chain.ToString()));
            row.Add(($"{tag}_intraMedian", F(intraMedian)));
            row.Add(($"{tag}_crossAboveIntraMedian", crossAboveIntraMedian.ToString()));
            row.Add(($"{tag}_crossAtDegLe2", crossLowDeg.ToString()));
        }
        return row;
    }

    private static double Median(List<double> xs)
    {
        var s = xs.OrderBy(x => x).ToArray();
        return s.Length % 2 == 1 ? s[s.Length / 2] : 0.5 * (s[s.Length / 2 - 1] + s[s.Length / 2]);
    }

    private static (double WInt, double WCross, int CrossAboveInt, int CrossEdges, double IntraMedian, int CrossLowDeg, int CrossAboveMedian)
        MergeOrder(List<(int U, int V, double W)> edges, int[] y, int k, int n, int[] coupledDegree)
    {
        var parent = Enumerable.Range(0, n).ToArray();
        int Find(int x) { while (parent[x] != x) x = parent[x] = parent[parent[x]]; return x; }
        var comps = new int[k];
        foreach (int c in y) comps[c]++;
        var bottleneck = new double[k];
        for (int c = 0; c < k; c++) bottleneck[c] = comps[c] <= 1 ? double.PositiveInfinity : 0.0;

        foreach (var (u, v, w) in edges.Where(e => e.W > 0.0 && y[e.U] == y[e.V]).OrderByDescending(e => e.W))
        {
            int ru = Find(u), rv = Find(v);
            if (ru == rv) continue;
            parent[ru] = rv;
            if (--comps[y[u]] == 1) bottleneck[y[u]] = w;
        }

        double wInt = bottleneck.Min();
        var intra = edges.Where(e => y[e.U] == y[e.V]).Select(e => e.W).ToList();
        double intraMedian = intra.Count == 0 ? 0.0 : Median(intra);
        double wCross = 0.0;
        int crossEdges = 0, aboveInt = 0, aboveMedian = 0, lowDeg = 0;
        foreach (var (u, v, w) in edges)
        {
            if (y[u] == y[v]) continue;
            crossEdges++;
            wCross = Math.Max(wCross, w);
            if (w > wInt) aboveInt++;
            if (w > intraMedian) aboveMedian++;
            if (w > 0.0 && (coupledDegree[u] <= 2 || coupledDegree[v] <= 2)) lowDeg++;
        }
        return (wInt, wCross, aboveInt, crossEdges, intraMedian, lowDeg, aboveMedian);
    }

    // ── SW sweeps ─────────────────────────────────────────────────────────────
    [Fact]
    public void Grid()
    {
        if (OutDir() is not string dir) return;
        foreach (Cell cell in Only(GridCells()))
            Emit(dir, "runs.tsv", RunRow("grid", cell, GridSeed, ColdEdge));
    }

    [Fact]
    public void Seeds()
    {
        if (OutDir() is not string dir) return;
        int[] seeds = SeedTree.Derive(GridSeed, 6);
        foreach (Cell cell in Only(DecidingCells()))
            foreach (int seed in seeds)
                Emit(dir, "runs.tsv", RunRow("seeds", cell, seed, ColdEdge));
    }

    [Fact]
    public void ColdEdgeRerun()
    {
        if (OutDir() is not string dir) return;
        foreach (LmpRankRule rank in new[] { LmpRankRule.Strict, LmpRankRule.Inclusive })
        {
            var cell = new Cell("blatt015", 10, rank, "Cauchy");
            if (Only(new[] { cell }).Any())
                Emit(dir, "runs.tsv", RunRow("cold1e-4", cell, GridSeed, 1e-4));
        }
    }

    private List<(string, string)> RunRow(string part, Cell cell, int seed, double tLo)
    {
        var clock = Stopwatch.StartNew();
        Data d = Load(cell.Dataset);
        GraphBuildResult build = Graph(cell);
        CsrGraph graph = build.Graph;
        int n = graph.NodeCount;
        double[] temps = SpcScheduleHelpers.LogSpaceGrid(tLo, HotEdge, Points);

        var cfg = new FixedGridSweepConfig
        {
            Temperatures = temps,
            Replicas = 4,
            SweepBudget = new RunBudget(200, 1000),
            Sampler = new PottsModelConfig { Q = Q },
            EquilibriumBudget = new RunBudget(10, 50),
            Accumulation = new AccumulationSpec { Affinities = true, CoMembership = true, ClusterSizeLandscape = true },
            SusceptibilityKind = SusceptibilityKind.FkReduced,
            BaseSeed = seed,
        };
        SweepResult sweep = new FixedGridSweepStrategy(cfg).Run(graph, new ChiPeakSignalAnalyzer());
        Accumulator[] frames = sweep.SweepRuns.Select(r => r.Accumulator).ToArray();
        PartitionHierarchy stack = SweepDenseTStack.Build(graph, frames, 0.5);
        HierarchyEomResult eom = SweepHierarchyEom.Resolve(graph, frames, 0.5, minClusterSize: n / 50);

        ResolvedConstruction r = build.Resolved!;
        var row = new List<(string, string)>
        {
            ("part", part), ("cell", cell.Id), ("dataset", cell.Dataset), ("kernel", cell.Kernel), ("K", cell.K.ToString()),
            ("rank", cell.Rank.ToString()), ("seed", seed.ToString()), ("tLo", F(tLo)),
            ("nodes", n.ToString()), ("edges", (graph.Targets.Length / 2).ToString()),
            ("zeroEdges", r.Coupling!.ZeroCouplingEdges.ToString()), ("structComps", r.Components.ToString()),
            ("coupledComps", r.Coupling.CoupledComponents.ToString()),
        };

        // The coldest partition against the coupled components: do the zeros, and not temperature, separate the clusters?
        int[] cold = stack.Levels[0].Partition.Labels;
        int[] coupledLabels = CoupledComponentLabels(graph);
        row.Add(("coldClusters", cold.Where(l => l >= 0).Distinct().Count().ToString()));
        row.Add(("coldEqualsCoupledComps", SamePartition(cold, coupledLabels).ToString()));

        double step = Math.Log(temps[1] / temps[0]);
        for (int li = 0; li < d.Levels.Length; li++)
        {
            string tag = d.Levels.Length == 1 ? "L" : new[] { "C", "M", "F" }[li];
            int[] y = Dense(d.Levels[li], out int k);
            LevelScore s = ScoreStack(stack, y, k, temps);
            LevelEval ev = Evaluate(eom.Assignment.Labels, y, k, s.MinReal);
            row.Add(($"{tag}_k", k.ToString()));
            row.Add(($"{tag}_stageFound", (s.Width >= 1.0).ToString()));
            row.Add(($"{tag}_stagePts", s.Points.ToString()));
            row.Add(($"{tag}_stageWidth", F(s.Width)));
            row.Add(($"{tag}_stageTlo", s.Points > 0 ? F(s.TLo) : ""));
            row.Add(($"{tag}_stageThi", s.Points > 0 ? F(s.THi) : ""));
            row.Add(($"{tag}_stageFromColdEnd", (s.Points > 0 && s.StartIndex == 0).ToString()));
            row.Add(($"{tag}_stageToHotEnd", (s.Points > 0 && s.StartIndex + s.Points == temps.Length).ToString()));
            row.Add(($"{tag}_bestMatched", F(s.BestMatched)));
            row.Add(($"{tag}_bestMatchedT", F(s.BestMatchedT)));
            row.Add(($"{tag}_eomReal", ev.Real.ToString()));
            row.Add(($"{tag}_eomMatches", ev.Stage.ToString()));
            row.Add(($"{tag}_eomMatched", F(ev.Matched)));
        }

        int[] eomLabels = eom.Assignment.Labels;
        int[] finest = Dense(d.Levels[^1], out int kFine);
        var (impure, selected) = Impurity(eomLabels, finest, kFine);
        row.Add(("eomClusters", eom.Assignment.Count.ToString()));
        row.Add(("eomAbstained", eomLabels.Count(l => l < 0).ToString()));
        row.Add(("eomImpureAtFinest", impure.ToString()));
        row.Add(("eomSelected", selected.ToString()));
        row.Add(("rawNested", eom.RawNestingHeld.ToString()));
        row.Add(("seconds", F(clock.Elapsed.TotalSeconds)));
        _out.WriteLine($"{cell.Id} seed {seed}: {clock.Elapsed.TotalSeconds:F0}s");
        return row;
    }

    // ── Scoring (primer §12.29 by BWD's rule) ─────────────────────────────────
    private readonly record struct LevelEval(bool Stage, double Matched, int Real);

    private sealed record LevelScore(
        int MinReal, int Points, int StartIndex, double TLo, double THi, double Width, double BestMatched, double BestMatchedT);

    /// <summary>
    /// A stage: exactly k real clusters (each holding at least 0.4 of the level's smallest class), with k distinct
    /// majority labels. The matched share credits each class once, by the largest real cluster it leads.
    /// </summary>
    private static LevelEval Evaluate(int[] labels, int[] y, int k, int minReal)
    {
        int clusters = labels.Length == 0 ? 0 : labels.Max() + 1;
        if (clusters <= 0) return new LevelEval(false, 0.0, 0);
        var counts = new int[clusters * k];
        var size = new int[clusters];
        for (int i = 0; i < labels.Length; i++)
        {
            if (labels[i] < 0) continue;
            counts[labels[i] * k + y[i]]++;
            size[labels[i]]++;
        }

        var real = new List<int>();
        for (int c = 0; c < clusters; c++) if (size[c] >= minReal) real.Add(c);
        // Each real cluster belongs to the class it holds most of; a class is credited once, by its largest such cluster.
        // (A giant cluster of several classes therefore credits one class only.)
        var majorities = new HashSet<int>();
        var credit = new int[k];
        foreach (int c in real)
        {
            int best = 0;
            for (int cls = 1; cls < k; cls++) if (counts[c * k + cls] > counts[c * k + best]) best = cls;
            majorities.Add(best);
            credit[best] = Math.Max(credit[best], counts[c * k + best]);
        }
        double matched = credit.Sum();
        return new LevelEval(real.Count == k && majorities.Count == k, matched / labels.Length, real.Count);
    }

    private static LevelScore ScoreStack(PartitionHierarchy stack, int[] y, int k, double[] temps)
    {
        int[] sizes = new int[k];
        foreach (int c in y) sizes[c]++;
        int minReal = (int)Math.Ceiling(0.4 * sizes.Min());

        int bestStart = 0, bestLen = 0, runStart = -1;
        double bestMatched = 0.0, bestMatchedT = double.NaN;
        for (int t = 0; t < stack.Levels.Count; t++)
        {
            LevelEval ev = Evaluate(stack.Levels[t].Partition.Labels, y, k, minReal);
            if (ev.Matched > bestMatched) { bestMatched = ev.Matched; bestMatchedT = stack.Levels[t].Temperature; }
            if (ev.Stage)
            {
                if (runStart < 0) runStart = t;
                int len = t - runStart + 1;
                if (len > bestLen) { bestLen = len; bestStart = runStart; }
            }
            else runStart = -1;
        }

        double lo = bestLen > 0 ? stack.Levels[bestStart].Temperature : double.NaN;
        double hi = bestLen > 0 ? stack.Levels[bestStart + bestLen - 1].Temperature : double.NaN;
        double width = bestLen > 1 ? Math.Log(hi / lo) : 0.0;
        return new LevelScore(minReal, bestLen, bestStart, lo, hi, width, bestMatched, bestMatchedT);
    }

    /// <summary>Selected clusters at the finest level whose majority class holds under 90% of them.</summary>
    private static (int Impure, int Selected) Impurity(int[] labels, int[] y, int k)
    {
        int clusters = labels.Length == 0 ? 0 : labels.Max() + 1;
        if (clusters <= 0) return (0, 0);
        var counts = new int[clusters * k];
        var size = new int[clusters];
        for (int i = 0; i < labels.Length; i++)
            if (labels[i] >= 0) { counts[labels[i] * k + y[i]]++; size[labels[i]]++; }
        int impure = 0, selected = 0;
        for (int c = 0; c < clusters; c++)
        {
            if (size[c] == 0) continue;
            selected++;
            int best = 0;
            for (int cls = 0; cls < k; cls++) best = Math.Max(best, counts[c * k + cls]);
            if (best < 0.9 * size[c]) impure++;
        }
        return (impure, selected);
    }

    private static int[] Dense(int[] labels, out int k)
    {
        var map = labels.Distinct().OrderBy(x => x).Select((v, i) => (v, i)).ToDictionary(p => p.v, p => p.i);
        k = map.Count;
        return labels.Select(l => map[l]).ToArray();
    }

    private static int[] CoupledComponentLabels(CsrGraph g)
    {
        var parent = Enumerable.Range(0, g.NodeCount).ToArray();
        int Find(int x) { while (parent[x] != x) x = parent[x] = parent[parent[x]]; return x; }
        for (int i = 0; i < g.NodeCount; i++)
            for (int s = g.RowPointers[i]; s < g.RowPointers[i + 1]; s++)
                if (g.Weights[s] > 0.0) parent[Find(i)] = Find(g.Targets[s]);
        return Dense(Enumerable.Range(0, g.NodeCount).Select(Find).ToArray(), out _);
    }

    /// <summary>Two labelings are the same partition when they differ only by renaming (unassigned points break it).</summary>
    private static bool SamePartition(int[] a, int[] b)
    {
        var ab = new Dictionary<int, int>();
        var ba = new Dictionary<int, int>();
        for (int i = 0; i < a.Length; i++)
        {
            if (a[i] < 0) return false;
            if (!ab.TryAdd(a[i], b[i]) && ab[a[i]] != b[i]) return false;
            if (!ba.TryAdd(b[i], a[i]) && ba[b[i]] != a[i]) return false;
        }
        return true;
    }
}
