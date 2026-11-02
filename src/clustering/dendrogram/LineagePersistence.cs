using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Dendrograms;

/// <summary>Which of a lineage's levels stands for it. These keep different information, so each is
/// selectable; <see cref="ClusterLineage.Levels"/> keeps every level, so any other reduction is a view.</summary>
public enum LineageRepresentative
{
    /// <summary>The largest level: the cluster when it formed, so capture. Covers the most, at some
    /// cost in purity at a low sampling budget.</summary>
    Largest,

    /// <summary>The level nearest the lineage's midpoint in ln T: a core, with the periphery that erodes
    /// in the first half of its life gone. Purer and steadier than <see cref="Largest"/>, and covers less.</summary>
    Midpoint,

    /// <summary>The least-decided level, the resolver's first representative (Azriel's leaf-call): the
    /// hottest and most eroded, so the smallest.</summary>
    LeastDecided,
}

/// <summary>What a lineage-persistence resolution reads and how: the whole record, echoed on the result.</summary>
public sealed record LineagePersistenceConfig
{
    /// <summary>The default <see cref="SplitShare"/>. Fitted to five graphs (iris, the BWD toy, Landsat and
    /// BlattHierarchy under two kernels, each at one K) and provisional: on Landsat the best share moves with K.
    /// Stopping primer §13.15–§13.16.</summary>
    public const double DefaultSplitShare = 0.10;

    /// <summary>Clusters below this size are ignored when building lineages: the overclustering-tail filter.</summary>
    public int MinClusterSize { get; init; } = 3;

    /// <summary>A step ends the lineage when its second successor holds at least this share of the cluster;
    /// below it, the sliver is erosion. A tiny positive share gives the size rule, where any second successor
    /// of at least <see cref="MinClusterSize"/> ends it.</summary>
    public double SplitShare { get; init; } = DefaultSplitShare;

    /// <summary>Read only the levels in [Lo, Hi], such as the SP-plateau from <c>SpcProfileAnalysis.SpPlateau</c>.
    /// Null reads the full grid.</summary>
    public (double Lo, double Hi)? TemperatureWindow { get; init; }

    /// <summary>Which level stands for each lineage. The default is read on the graphs of stopping primer §13.16.</summary>
    public LineageRepresentative Representative { get; init; } = LineageRepresentative.Largest;

    /// <summary>Report the lineages alive at the level nearest this temperature in ln T, which must lie within
    /// the levels read. Null takes the longest-lived lineage on each branch instead. Which level of a hierarchy
    /// to report is a modelling choice: this one asks for it by number, in the units of this graph, so a T* read on
    /// one graph doesn't carry over to another.</summary>
    public double? SliceTemperature { get; init; }

    /// <summary>Keep at most this many lineages, highest score first. Null keeps all that the reading gives.</summary>
    public int? SelectTopK { get; init; }
}

/// <summary>One level of a lineage: one cluster at one temperature. A stack that doesn't nest can put two
/// clusters of a lineage at one temperature, each with its own entry.</summary>
/// <param name="Temperature">The grid temperature.</param>
/// <param name="Members">The cluster's points.</param>
/// <param name="Decidedness">Mean δ̄ over the cluster's internal edges.</param>
/// <param name="CellWidth">The level's cell in ln T: the half-gaps to its neighbours among the levels read.</param>
public sealed record LineageLevel(double Temperature, int[] Members, double Decidedness, double CellWidth);

/// <summary>One cluster lineage tracked across the T-stack: the chain of clusters at successive temperatures
/// that are the same cluster, eroding and not splitting.</summary>
/// <param name="Members">The representative level's points, per <see cref="LineagePersistenceConfig.Representative"/>.</param>
/// <param name="TBirth">Coldest temperature the lineage exists at.</param>
/// <param name="TDeath">Hottest temperature the lineage exists at.</param>
/// <param name="Lifetime">Its levels' cells summed, ∫ d ln T in e-folds. Invariant to the grid's density, and
/// not to its extent: see <paramref name="ColdCensored"/>.</param>
/// <param name="Decidedness">δ̄ averaged over the lineage's ln T: <paramref name="Score"/> / <paramref name="Lifetime"/>.</param>
/// <param name="Score">∫ δ̄ d ln T over the lineage, with δ̄ averaged within a level holding two of its clusters.</param>
/// <param name="LevelCount">Distinct temperatures the lineage survives.</param>
/// <param name="ColdCensored">The lineage reaches the coldest level read, so it may extend colder than the grid does
/// and its lifetime is a lower bound. A longest-lived reading that selects such a lineage depends on the cold edge.</param>
/// <param name="Levels">Every level of the lineage, cold to hot; reductions are views of these.</param>
public sealed record ClusterLineage(
    int[]  Members,
    double TBirth,
    double TDeath,
    double Lifetime,
    double Decidedness,
    double Score,
    int    LevelCount,
    bool   ColdCensored,
    IReadOnlyList<LineageLevel> Levels)
{
    /// <summary>The level that stands for this lineage under <paramref name="kind"/>; ties go to the coldest.</summary>
    public LineageLevel Representative(LineageRepresentative kind) => LineagePersistence.Representative(Levels, kind);
}

/// <summary>Everything the lineage-persistence resolution produced, whole for inspection.</summary>
/// <param name="Stack">The stack read.</param>
/// <param name="AllLineages">Every lineage tracked, the menu <see cref="LineagePersistenceConfig.SliceTemperature"/> reads from.</param>
/// <param name="Selected">The lineages reported, highest score first; <c>Selected[k]</c> is label k.</param>
/// <param name="Config">The configuration resolved.</param>
/// <param name="SliceLevel">The grid temperature a slice read; null on the longest-lived reading.</param>
/// <param name="Assignment">Each point's label, or <see cref="Assignment.Unassigned"/> if no selected lineage holds it.</param>
public sealed record LineagePersistenceResult(
    PartitionHierarchy Stack,
    IReadOnlyList<ClusterLineage> AllLineages,
    IReadOnlyList<ClusterLineage> Selected,
    LineagePersistenceConfig Config,
    double? SliceLevel,
    Assignment Assignment);

/// <summary>
/// Track 2, the <b>lineage-persistence</b> resolver. Named for what it does (it is only conceptually wave_clus,
/// ChaureReyQuiroga2018): wave_clus's insight lifted to one principle, <b>read cluster lineages by persistence over
/// the T-stack</b>. wave_clus's three mechanisms collapse here: peak-select → persistence, inclusion/overlap-dedup →
/// overlap-linked lineages, regime-border → the SP-plateau temperature bound. A distinct member from the Blatt/Domany
/// dendrogram resolver: it never builds a dendrogram, it tracks clusters across the (non-uniform) stack and scores
/// their survival.
/// </summary>
/// <remarks>
/// <para><b>Lineages.</b> Each cluster at T_n is linked to its continuation at T_{n+1}, the successor holding the
/// plurality of its members, unless it splits: when a second successor holds at least
/// <see cref="LineagePersistenceConfig.SplitShare"/> of the cluster, the lineage ends and both children begin fresh
/// ones. The share separates erosion from a split. It began as a gap in the observed shares (Azriel's leaf-call);
/// measured, those shares have no gap to find, so it is a fixed, exposed value.</para>
///
/// <para><b>Persistence.</b> A lineage scores ∫ δ̄ d ln T, each level owning the half-gaps to its neighbours. That
/// holds across a grid's density, where a count of levels does not. It doesn't hold across its extent, since the cold
/// edge censors any lineage reaching it, which <see cref="ClusterLineage.ColdCensored"/> marks.</para>
///
/// <para><b>The reading.</b> By default, the longest-lived lineage on each branch: rank by score, and drop a lineage
/// whose representative shares at least half of the smaller set with a better one (wave_clus's inclusion criterion,
/// lifted off its k_O = 0.9). Every lineage left is kept, since a lineage disjoint from the others has no reason to be
/// dropped for living less long. It takes the widest stage the grid shows. Given
/// <see cref="LineagePersistenceConfig.SliceTemperature"/>, the lineages alive at that level instead, each
/// represented whole, which reaches every level the grid holds. The evidence and the trade-offs are in
/// <c>docs/Clustering/SpcAnalysis.md</c> and stopping primer §13.15–§13.16.</para>
///
/// <para><b>Output.</b> A point in several selected representatives goes to the higher-scoring lineage, and a point in
/// none is an honest <see cref="Assignment.Unassigned"/>. A lineage left with no point isn't selected. Validated against
/// the BWD/Domany oracles (published numbers + true labels), not against wave_clus's exact outputs.</para>
/// </remarks>
public static class LineagePersistence
{
    /// <summary>
    /// Producer-agnostic core: resolve a dense stack + its δ̄ columns. SW supplies pooled sampled columns
    /// (SPC's <c>SweepLineagePersistence</c>); a solver or BARS-monotonized columns feed this directly.
    /// </summary>
    public static LineagePersistenceResult Resolve(
        CsrGraph graph,
        PartitionHierarchy stack,
        IReadOnlyList<double[]> deltaByGridPoint,
        int q,
        LineagePersistenceConfig? config = null)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(stack);
        ArgumentNullException.ThrowIfNull(deltaByGridPoint);
        if (deltaByGridPoint.Count != stack.Count)
            throw new ArgumentException(
                $"One δ̄ column per stack level: {deltaByGridPoint.Count} vs {stack.Count}.",
                nameof(deltaByGridPoint));

        config ??= new LineagePersistenceConfig();
        if (!(config.SplitShare > 0.0 && config.SplitShare <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(config), config.SplitShare, "SplitShare must lie in (0, 1].");
        if (config.SliceTemperature is double requested && !(requested > 0.0))
            throw new ArgumentOutOfRangeException(nameof(config), requested, "SliceTemperature must be positive.");

        int n = graph.NodeCount;
        double noiseFloor = 1.0 / q;

        // Levels inside the window (regime border).
        var levelIdx = new List<int>();
        for (int li = 0; li < stack.Count; li++)
        {
            double t = stack.Levels[li].Temperature;
            if (config.TemperatureWindow is { } w && (t < w.Lo || t > w.Hi)) continue;
            levelIdx.Add(li);
        }
        if (levelIdx.Count == 0)
            return Abstain(stack, config, n);

        double[] temperatures = levelIdx.Select(li => stack.Levels[li].Temperature).ToArray();
        double[] cellWidths = CellWidths(temperatures);
        double? sliceLevel = config.SliceTemperature is double tStar ? NearestLevel(temperatures, tStar) : null;

        // ── Enumerate clusters (level → label → members), size-filtered. ──
        // clusterIndexByLevelLabel[w] maps the dense label → global cluster index for the
        // w-th level read; clusters[gi] = (windowPos, label, members). Cold levels come first.
        var clusters = new List<(int WindowPos, int Label, int[] Members)>();
        var clusterIndexByLevelLabel = new List<Dictionary<int, int>>(levelIdx.Count);
        var decidedness = new List<double>();
        for (int w = 0; w < levelIdx.Count; w++)
        {
            int li = levelIdx[w];
            int[] labels = stack.Levels[li].Partition.Labels;
            double[] col = deltaByGridPoint[li];

            var members = new Dictionary<int, List<int>>();
            for (int p = 0; p < n; p++)
            {
                int lab = labels[p];
                if (lab == Assignment.Unassigned) continue;
                if (!members.TryGetValue(lab, out var list)) { list = new List<int>(); members[lab] = list; }
                list.Add(p);
            }

            // Intra-cluster δ̄ sums for decidedness.
            var dSum   = new Dictionary<int, double>();
            var dCount = new Dictionary<int, int>();
            foreach (UndirectedEdge edge in graph.UndirectedEdges())
            {
                int la = labels[edge.Source];
                if (la == Assignment.Unassigned || la != labels[edge.Target]) continue;
                dSum[la]   = (dSum.TryGetValue(la, out double s) ? s : 0.0) + col[edge.Slot];
                dCount[la] = (dCount.TryGetValue(la, out int c) ? c : 0) + 1;
            }

            var map = new Dictionary<int, int>();
            foreach (var (lab, list) in members)
            {
                if (list.Count < config.MinClusterSize) continue;
                int gi = clusters.Count;
                clusters.Add((w, lab, list.ToArray()));
                map[lab] = gi;
                double dec = dCount.TryGetValue(lab, out int cnt) && cnt > 0 ? dSum[lab] / cnt : noiseFloor;
                decidedness.Add(dec);
            }
            clusterIndexByLevelLabel.Add(map);
        }
        if (clusters.Count == 0)
            return Abstain(stack, config, n);

        // ── Link continuations (not splits) into lineages. ──
        // Distribute each cluster's members over the clusters one level hotter: the plurality is its
        // continuation, and a second successor holding the split share ends the lineage.
        var uf = new UnionFind(clusters.Count);
        for (int gi = 0; gi < clusters.Count; gi++)
        {
            var (w, _, members) = clusters[gi];
            if (w + 1 >= levelIdx.Count) continue;
            int[] nextLabels = stack.Levels[levelIdx[w + 1]].Partition.Labels;
            Dictionary<int, int> nextMap = clusterIndexByLevelLabel[w + 1];

            var tally = new Dictionary<int, int>();   // successor gi' → shared count
            foreach (int p in members)
            {
                int lab = nextLabels[p];
                if (lab == Assignment.Unassigned || !nextMap.TryGetValue(lab, out int gj)) continue;
                tally[gj] = (tally.TryGetValue(gj, out int c) ? c : 0) + 1;
            }
            if (tally.Count == 0) continue;

            var ordered = tally.OrderByDescending(kv => kv.Value).ToArray();
            double secondShare = ordered.Length > 1 ? (double)ordered[1].Value / members.Length : 0.0;
            if (secondShare >= config.SplitShare) continue;   // a split: the lineage ends
            uf.Union(gi, ordered[0].Key);
        }

        // ── Reduce each lineage, keeping its levels. ──
        var lineageClusters = new Dictionary<int, List<int>>();
        for (int gi = 0; gi < clusters.Count; gi++)
        {
            int root = uf.Find(gi);
            if (!lineageClusters.TryGetValue(root, out var list)) { list = new List<int>(); lineageClusters[root] = list; }
            list.Add(gi);
        }

        var lineages = new List<ClusterLineage>(lineageClusters.Count);
        foreach (var group in lineageClusters.Values)
        {
            var levels = new List<LineageLevel>(group.Count);
            foreach (int gi in group)   // ascending, so cold to hot
            {
                var (w, _, members) = clusters[gi];
                levels.Add(new LineageLevel(temperatures[w], members, decidedness[gi], cellWidths[w]));
            }

            double lifetime = 0.0, score = 0.0;
            int levelCount = 0;
            foreach (var atTemperature in levels.GroupBy(l => l.Temperature))
            {
                double width = atTemperature.First().CellWidth;
                lifetime += width;
                score += width * atTemperature.Average(l => l.Decidedness);
                levelCount++;
            }

            lineages.Add(new ClusterLineage(
                Members:      Representative(levels, config.Representative).Members,
                TBirth:       levels[0].Temperature,
                TDeath:       levels[^1].Temperature,
                Lifetime:     lifetime,
                Decidedness:  score / lifetime,
                Score:        score,
                LevelCount:   levelCount,
                ColdCensored: levels[0].Temperature == temperatures[0],
                Levels:       levels));
        }

        // ── The reading: longest-lived on each branch, or the lineages alive at the slice's level. ──
        var ranked = lineages.Where(l => l.Score > 0.0).OrderByDescending(l => l.Score).ToList();
        IEnumerable<ClusterLineage> read = sliceLevel is double level
            ? ranked.Where(l => l.Levels.Any(x => x.Temperature == level))
            : DedupByOverlap(ranked, overlapThreshold: 0.5);
        if (config.SelectTopK is int k) read = read.Take(Math.Max(0, k));

        // ── Resolve to Assignment: higher score claims first. ──
        var labelsOut = new int[n];
        Array.Fill(labelsOut, Assignment.Unassigned);
        var selected = new List<ClusterLineage>();
        foreach (ClusterLineage lineage in read)
        {
            bool claimed = false;
            foreach (int p in lineage.Members)
            {
                if (labelsOut[p] != Assignment.Unassigned) continue;
                labelsOut[p] = selected.Count;
                claimed = true;
            }
            if (claimed) selected.Add(lineage);
        }

        var assignment = new Assignment { Labels = labelsOut, Count = selected.Count };
        return new LineagePersistenceResult(stack, lineages, selected, config, sliceLevel, assignment);
    }

    internal static LineageLevel Representative(IReadOnlyList<LineageLevel> levels, LineageRepresentative kind)
    {
        LineageLevel best = levels[0];
        switch (kind)
        {
            case LineageRepresentative.Largest:
                foreach (LineageLevel level in levels)
                    if (level.Members.Length > best.Members.Length) best = level;
                return best;

            case LineageRepresentative.Midpoint:
                double mid = 0.5 * (Math.Log(levels[0].Temperature) + Math.Log(levels[^1].Temperature));
                double nearest = double.PositiveInfinity;
                foreach (LineageLevel level in levels)
                {
                    double d = Math.Abs(Math.Log(level.Temperature) - mid);
                    if (d < nearest) { nearest = d; best = level; }
                }
                return best;

            case LineageRepresentative.LeastDecided:
                foreach (LineageLevel level in levels)
                    if (level.Decidedness < best.Decidedness) best = level;
                return best;

            default:
                throw new ArgumentOutOfRangeException(nameof(kind), kind, "Unknown representative.");
        }
    }

    /// <summary>Each level's cell in ln T: the half-gaps to its neighbours, so the edges own one half-gap and the
    /// cells sum to the span. A lone level has no span, and gets a unit cell so that its score stays positive.</summary>
    private static double[] CellWidths(double[] temperatures)
    {
        if (temperatures[0] <= 0.0)
            throw new ArgumentException($"Lifetimes read ln T, so temperatures must be positive; got {temperatures[0]}.", nameof(temperatures));
        for (int i = 1; i < temperatures.Length; i++)
            if (temperatures[i] <= temperatures[i - 1])
                throw new ArgumentException("Stack temperatures must be strictly ascending.", nameof(temperatures));

        var widths = new double[temperatures.Length];
        if (temperatures.Length == 1) { widths[0] = 1.0; return widths; }
        for (int i = 0; i < temperatures.Length; i++)
        {
            double below = i > 0 ? Math.Log(temperatures[i] / temperatures[i - 1]) : 0.0;
            double above = i + 1 < temperatures.Length ? Math.Log(temperatures[i + 1] / temperatures[i]) : 0.0;
            widths[i] = 0.5 * (below + above);
        }
        return widths;
    }

    /// <summary>The level nearest <paramref name="tStar"/> in ln T. A temperature outside the levels read has no
    /// level to report, and snapping to the edge would silently read a different stage.</summary>
    private static double NearestLevel(double[] temperatures, double tStar)
    {
        if (tStar < temperatures[0] || tStar > temperatures[^1])
            throw new ArgumentOutOfRangeException(nameof(tStar), tStar,
                $"SliceTemperature lies outside the levels read, [{temperatures[0]:G4}, {temperatures[^1]:G4}].");
        double best = temperatures[0], nearest = double.PositiveInfinity;
        foreach (double t in temperatures)
        {
            double d = Math.Abs(Math.Log(t / tStar));
            if (d < nearest) { nearest = d; best = t; }
        }
        return best;
    }

    /// <summary>The inclusion criterion: walk the score-ranked lineages keeping only those whose representative is
    /// not already covered by a higher-ranked keeper (overlap coefficient |∩|/min ≥ <paramref name="overlapThreshold"/>).
    /// Drops the same cluster tracked as several near-duplicate lineages, and a parent against its longer-lived
    /// children or the reverse; keeps genuinely distinct ones.</summary>
    private static List<ClusterLineage> DedupByOverlap(List<ClusterLineage> ranked, double overlapThreshold)
    {
        var kept = new List<ClusterLineage>();
        var keptSets = new List<HashSet<int>>();
        foreach (var lineage in ranked)
        {
            var set = new HashSet<int>(lineage.Members);
            bool redundant = false;
            for (int k = 0; k < keptSets.Count; k++)
            {
                int inter = 0;
                foreach (int p in set) if (keptSets[k].Contains(p)) inter++;
                int min = Math.Min(set.Count, keptSets[k].Count);
                if (min > 0 && (double)inter / min >= overlapThreshold) { redundant = true; break; }
            }
            if (!redundant) { kept.Add(lineage); keptSets.Add(set); }
        }
        return kept;
    }

    private static LineagePersistenceResult Abstain(PartitionHierarchy stack, LineagePersistenceConfig config, int n)
    {
        var labels = new int[n];
        Array.Fill(labels, Assignment.Unassigned);
        return new LineagePersistenceResult(stack, Array.Empty<ClusterLineage>(), Array.Empty<ClusterLineage>(),
            config, null, new Assignment { Labels = labels, Count = 0 });
    }
}
