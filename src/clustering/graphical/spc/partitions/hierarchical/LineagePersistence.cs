using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Hierarchical;

/// <summary>One cluster lineage tracked across the T-stack: the chain of
/// clusters at successive temperatures that are the SAME cluster (eroding
/// periphery, not splitting), reduced to a representative member set and a
/// persistence score.</summary>
/// <param name="Members">The representative cut's member points (the cluster at
/// the lineage's least-decided / most-inclusive level — the leaf-call default
/// representative-T).</param>
/// <param name="TBirth">Coldest temperature the lineage exists at.</param>
/// <param name="TDeath">Hottest temperature the lineage exists at.</param>
/// <param name="TSpan">T-range survival = TDeath − TBirth (the actual,
/// possibly non-uniform, temperature span — reported for inspection).</param>
/// <param name="Decidedness">Mean intra-cluster δ̄ over the lineage — crispness;
/// a merge of two species reads low here, a real cluster high.</param>
/// <param name="Score">LevelCount × Decidedness — the selection score.
/// Persistence is measured in RESOLUTION STEPS survived, not absolute ΔT: the
/// only grid-agnostic reading (absolute span over-rewards the hot end of a log
/// grid; log-span over-rewards the cold end of a linear one — step-count is
/// neutral, the schedule-agnostic contract).</param>
/// <param name="LevelCount">Number of T-levels the lineage survives — the
/// persistence axis.</param>
public sealed record ClusterLineage(
    int[]  Members,
    double TBirth,
    double TDeath,
    double TSpan,
    double Decidedness,
    double Score,
    int    LevelCount);

/// <summary>Everything the lineage-persistence resolution produced, whole for
/// inspection.</summary>
public sealed record LineagePersistenceResult(
    PartitionHierarchy Stack,
    IReadOnlyList<ClusterLineage> AllLineages,
    IReadOnlyList<ClusterLineage> Selected,
    double SplitShare,
    Assignment Assignment);

/// <summary>
/// Track 2 — the <b>lineage-persistence</b> resolver. Named for what it does
/// (it is only CONCEPTUALLY wave_clus, ChaureReyQuiroga2018): wave_clus's
/// insight lifted to ONE principle — <b>select cluster lineages by persistence
/// over the T-stack</b>. wave_clus's three mechanisms collapse here — peak-select → persistence,
/// inclusion/overlap-dedup → overlap-linked lineages, regime-border → the
/// SP-plateau temperature bound — without the elbow grease (no N_inc growth,
/// no Thr_border=0.4, no k_O=0.9). A DISTINCT member from the Blatt/Domany
/// dendrogram resolver: it never builds a dendrogram, it tracks clusters across
/// the (non-uniform) stack and scores their survival.
/// </summary>
/// <remarks>
/// <para><b>Lineages.</b> Each cluster at temperature T_n is linked to its
/// continuation at T_{n+1} — the successor holding the plurality of its members
/// — UNLESS it splits: when a second successor takes a real share
/// (<see cref="LineagePersistenceResult.SplitShare"/> of the cluster), the lineage ends
/// and both children begin fresh lineages. The split share is the only knob
/// (Azriel's leaf-call: a GAP in the second-successor-share distribution, else
/// an exposed threshold) — it separates periphery erosion (a sliver leaks
/// elsewhere) from a genuine split (a half does), grid-spacing-robust where a
/// fixed overlap coefficient is not.</para>
///
/// <para><b>Persistence.</b> Each lineage scores T-range survival × decidedness:
/// the temperature span it survives times its mean intra-cluster δ̄ (the
/// fuzzy/co-membership crispness). The decidedness factor is what demotes a
/// merge-of-two-species (weakly coupled, low δ̄) below the real species; span
/// alone would not. Selection takes the lineages above the largest gap in the
/// sorted score distribution (the persistence elbow), else an exposed top-K.</para>
///
/// <para><b>Output.</b> Each selected lineage's representative member set
/// becomes a cluster; a point in several selected reps goes to the higher-score
/// lineage; a point in none is an honest <see cref="Assignment.Unassigned"/>.
/// Validated ADVERSARIALLY against the BWD/Domany oracles (published numbers +
/// true labels), not against wave_clus's exact outputs.</para>
/// </remarks>
public static class LineagePersistence
{
    /// <summary>
    /// Resolve a rich sweep's frames by lineage persistence. Requires
    /// <c>AccumulationSpec.CoMembership</c> (the per-edge co-membership counts —
    /// the partition columns AND the decidedness field). No landscape needed.
    /// </summary>
    /// <param name="minClusterSize">Clusters below this size are ignored when
    /// building lineages (the overclustering-tail filter; default 3).</param>
    /// <param name="splitShare">Second-successor share above which a step is a
    /// split, not erosion. Null ⇒ gap-based on the observed distribution
    /// (fallback 0.25).</param>
    /// <param name="temperatureWindow">Restrict to [Lo, Hi] — the SP-plateau
    /// bound (regime border: T_fs/T_ps from <see cref="SpcProfileAnalysis.SpPlateau"/>).
    /// Null ⇒ the full grid.</param>
    /// <param name="selectTopK">Force exactly this many lineages instead of the
    /// gap rule. Null ⇒ gap-based selection.</param>
    public static LineagePersistenceResult Resolve(
        CsrGraph graph,
        IReadOnlyList<Accumulator> frames,
        double theta = 0.5,
        int minClusterSize = 3,
        double? splitShare = null,
        (double Lo, double Hi)? temperatureWindow = null,
        int? selectTopK = null)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(frames);

        int q = frames[0].Q;
        var (temperatures, delta) = SweepEdgeCurves.CoMembershipDelta(frames);
        PartitionHierarchy stack = DenseTStack.Build(graph, temperatures, delta, theta);

        return Resolve(graph, stack, delta, q, minClusterSize, splitShare, temperatureWindow, selectTopK);
    }

    /// <summary>
    /// Producer-agnostic core: resolve a dense stack + its δ̄ columns. SW supplies
    /// pooled sampled columns (the overload above); a solver / BARS-monotonized
    /// columns feed this directly.
    /// </summary>
    public static LineagePersistenceResult Resolve(
        CsrGraph graph,
        PartitionHierarchy stack,
        IReadOnlyList<double[]> deltaByGridPoint,
        int q,
        int minClusterSize = 3,
        double? splitShare = null,
        (double Lo, double Hi)? temperatureWindow = null,
        int? selectTopK = null)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(stack);
        ArgumentNullException.ThrowIfNull(deltaByGridPoint);
        if (deltaByGridPoint.Count != stack.Count)
            throw new ArgumentException(
                $"One δ̄ column per stack level: {deltaByGridPoint.Count} vs {stack.Count}.",
                nameof(deltaByGridPoint));

        int n = graph.NodeCount;
        double noiseFloor = 1.0 / q;

        // Levels inside the SP-plateau window (regime border).
        var levelIdx = new List<int>();
        for (int li = 0; li < stack.Count; li++)
        {
            double t = stack.Levels[li].Temperature;
            if (temperatureWindow is { } w && (t < w.Lo || t > w.Hi)) continue;
            levelIdx.Add(li);
        }
        if (levelIdx.Count == 0)
            return new LineagePersistenceResult(stack, Array.Empty<ClusterLineage>(),
                Array.Empty<ClusterLineage>(), splitShare ?? 0.25, AllAbstain(n));

        // ── Enumerate clusters (level → label → members), size-filtered. ──
        // clusterByLevel[w] maps the dense label → global cluster index for the
        // w-th in-window level; clusters[gi] = (windowPos, label, members).
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
                if (list.Count < minClusterSize) continue;
                int gi = clusters.Count;
                clusters.Add((w, lab, list.ToArray()));
                map[lab] = gi;
                double dec = dCount.TryGetValue(lab, out int cnt) && cnt > 0 ? dSum[lab] / cnt : noiseFloor;
                decidedness.Add(dec);
            }
            clusterIndexByLevelLabel.Add(map);
        }
        if (clusters.Count == 0)
            return new LineagePersistenceResult(stack, Array.Empty<ClusterLineage>(),
                Array.Empty<ClusterLineage>(), splitShare ?? 0.25, AllAbstain(n));

        // ── Continuation analysis: each cluster's successor shares at w+1. ──
        // For cluster gi at window-level w, distribute its members over the
        // clusters at w+1; record (best successor gi', s1, s2 as fractions).
        var bestSucc = new int[clusters.Count];
        var s2Share  = new double[clusters.Count];
        Array.Fill(bestSucc, -1);
        var s2Observed = new List<double>();
        for (int gi = 0; gi < clusters.Count; gi++)
        {
            var (w, _, members) = clusters[gi];
            if (w + 1 >= levelIdx.Count) continue;
            int liNext = levelIdx[w + 1];
            int[] nextLabels = stack.Levels[liNext].Partition.Labels;
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
            double inv = 1.0 / members.Length;
            bestSucc[gi] = ordered[0].Key;
            double share2 = ordered.Length > 1 ? ordered[1].Value * inv : 0.0;
            s2Share[gi] = share2;
            s2Observed.Add(share2);
        }

        // A real split sends a real share (≥ ~a fifth) to a second child; below
        // that is periphery erosion. Gap-based within that band, else 0.3.
        double splitThreshold = splitShare ?? GapThreshold(s2Observed, lo: 0.2, hi: 0.5, fallback: 0.3);

        // ── Link continuations (not splits) into lineages. ──
        var uf = new UnionFind(clusters.Count);
        for (int gi = 0; gi < clusters.Count; gi++)
        {
            if (bestSucc[gi] < 0) continue;
            if (s2Share[gi] >= splitThreshold) continue;   // a split — end the lineage
            uf.Union(gi, bestSucc[gi]);
        }

        // ── Reduce each lineage to a scored representative. ──
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
            double tBirth = double.PositiveInfinity, tDeath = double.NegativeInfinity;
            double decSum = 0.0;
            int repCluster = group[0];
            double repDec = double.PositiveInfinity;
            foreach (int gi in group)
            {
                int li = levelIdx[clusters[gi].WindowPos];
                double t = stack.Levels[li].Temperature;
                if (t < tBirth) tBirth = t;
                if (t > tDeath) tDeath = t;
                decSum += decidedness[gi];
                // Representative = the least-decided level (Azriel's leaf-call).
                // The crisp eroded core keeps DISTINCT lineages cleanly separated
                // for the inclusion dedup; the more-inclusive pre-split blob would
                // overlap its own children and collapse them. (Leaf-call to refine
                // — e.g. an inclusive rep for capture, an eroded rep for dedup.)
                if (decidedness[gi] < repDec) { repDec = decidedness[gi]; repCluster = gi; }
            }
            double meanDec = decSum / group.Count;
            double span = tDeath - tBirth;
            lineages.Add(new ClusterLineage(
                Members:      clusters[repCluster].Members,
                TBirth:       tBirth,
                TDeath:       tDeath,
                TSpan:        span,
                Decidedness:  meanDec,
                Score:        group.Count * meanDec,
                LevelCount:   group.Count));
        }

        // ── Inclusion criterion (dedup) then the persistence gap. ──
        // Rank by persistence; drop a lineage whose representative is essentially
        // a higher-ranked one seen again (overlap ≥ ½ of the smaller) — wave_clus's
        // inclusion criterion, lifted off the k_O=0.9 constant; then cut at the
        // persistence elbow (or a forced top-K).
        var ranked = lineages.Where(l => l.Score > 0.0).OrderByDescending(l => l.Score).ToList();
        var deduped = DedupByOverlap(ranked, overlapThreshold: 0.5);
        var selected = SelectByGap(deduped, selectTopK);

        // ── Resolve to Assignment: representative reps, higher score wins ties. ──
        var labelsOut = new int[n];
        Array.Fill(labelsOut, Assignment.Unassigned);
        for (int k = 0; k < selected.Count; k++)   // selected already score-descending
            foreach (int p in selected[k].Members)
                if (labelsOut[p] == Assignment.Unassigned) labelsOut[p] = k;

        var assignment = new Assignment { Labels = labelsOut, Count = selected.Count };
        return new LineagePersistenceResult(stack, lineages, selected, splitThreshold, assignment);
    }

    /// <summary>The inclusion criterion: walk the score-ranked lineages keeping
    /// only those whose representative is NOT already covered by a higher-ranked
    /// keeper (overlap coefficient |∩|/min ≥ <paramref name="overlapThreshold"/>)
    /// — drops the same cluster tracked as several near-duplicate lineages, keeps
    /// genuinely-distinct ones.</summary>
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

    private static Assignment AllAbstain(int n)
    {
        var labels = new int[n];
        Array.Fill(labels, Assignment.Unassigned);
        return new Assignment { Labels = labels, Count = 0 };
    }

    /// <summary>Take the lineages above the largest relative gap in the sorted
    /// score sequence (the persistence elbow); a forced <paramref name="topK"/>
    /// overrides. Requires the gap ratio ≥ 1.5 to fire, else falls back to
    /// scores ≥ half the maximum (the exposed default).</summary>
    private static List<ClusterLineage> SelectByGap(List<ClusterLineage> ranked, int? topK)
    {
        if (ranked.Count == 0) return ranked;
        if (topK is int k) return ranked.Take(Math.Max(0, k)).ToList();
        if (ranked.Count == 1) return ranked;

        int cut = -1;
        double bestRatio = 1.0;
        for (int i = 0; i < ranked.Count - 1; i++)
        {
            double a = ranked[i].Score, b = ranked[i + 1].Score;
            if (b <= 0.0) { cut = i; break; }
            double ratio = a / b;
            if (ratio > bestRatio) { bestRatio = ratio; cut = i; }
        }

        if (cut >= 0 && bestRatio >= 1.5)
            return ranked.Take(cut + 1).ToList();

        double half = 0.5 * ranked[0].Score;
        return ranked.Where(l => l.Score >= half).ToList();
    }

    /// <summary>Midpoint of the largest gap in the sorted <paramref name="values"/>
    /// restricted to (<paramref name="lo"/>, <paramref name="hi"/>); the
    /// data-driven split-share threshold. Falls back to
    /// <paramref name="fallback"/> when no values land in the band.</summary>
    private static double GapThreshold(IReadOnlyList<double> values, double lo, double hi, double fallback)
    {
        var band = values.Where(v => v > lo && v < hi).OrderBy(v => v).ToArray();
        if (band.Length == 0) return fallback;
        if (band.Length == 1) return Math.Clamp(band[0], lo, hi);

        double bestGap = 0.0, threshold = fallback;
        for (int i = 0; i < band.Length - 1; i++)
        {
            double gap = band[i + 1] - band[i];
            if (gap > bestGap) { bestGap = gap; threshold = 0.5 * (band[i] + band[i + 1]); }
        }
        // No meaningful gap ⇒ the exposed default.
        return bestGap > 0.05 ? threshold : fallback;
    }
}
