using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Clustering.Dendrograms;
using Clustering.Graphical.SPC.Partitions.Hierarchical;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;
using Maths.Rng;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The closure gap against ground truth from outside it: co-membership by enumerating every bond configuration of the
/// Fortuin–Kasteleyn measure, the best partition by searching every partition, and the join level by a widest-path
/// recursion. A dense field is the FK law's own partitions, each repeated in proportion to its weight, so it is a true
/// mixture of partitions, and each check reads its truth off that field. Tolerances were fixed before the first run:
/// 1e-10 on exact identities, 1e-3 on the θ integral's midpoint rule (15 pairs × a jump of at most 1 × half of h =
/// 5e-5, so 7.5e-4) and on the rounded fields (32 configurations × ½/10⁶ × 15 pairs, so 2.4e-4, and likewise below it
/// for the two poles), 0.02 per pair for Swendsen–Wang against exact, and for the noise check a plug-in mean above 0.2
/// (against about 0.26 from 7.5 false merges × E|N(0, 0.044)|) and a cross-fit mean within four standard errors of 0.
/// </summary>
public sealed class ClosureGapTests
{
    private const double Exact = 1e-10;

    private static readonly double[] Thetas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
    private static readonly int[] Qs = [2, 20];
    private static readonly double[] Betas = [0.3, 0.6, 0.9, 1.2, 1.8, 2.4, 3.5, 4.8, 6.5, 9.6, 14.0, 20.0];

    private static readonly string[] AllGraphs = ["path", "bridged", "mixed", "strong-weak", "chorded-ring", "two-poles"];

    public static TheoryData<string> GraphNames
    {
        get
        {
            var names = new TheoryData<string>();
            foreach (string name in AllGraphs) names.Add(name);
            return names;
        }
    }

    private static CsrGraph Fixture(string name) => name switch
    {
        "path" => Graph(6, (0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 4, 1.0), (4, 5, 1.0)),
        "bridged" => Graph(6, (0, 1, 1.0), (0, 2, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 4, 1.0), (3, 5, 1.0), (4, 5, 1.0)),
        "mixed" => Graph(6, (0, 1, 0.9), (1, 2, 0.35), (0, 2, 0.7), (2, 3, 0.2), (3, 4, 0.8), (4, 5, 0.6), (3, 5, 0.4), (1, 4, 0.25)),
        "strong-weak" => Graph(6, (0, 1, 0.5), (0, 2, 0.5), (1, 2, 0.5), (1, 3, 0.15), (2, 3, 1.0), (3, 4, 0.3), (4, 5, 0.95), (3, 5, 0.45), (0, 5, 0.2)),
        "chorded-ring" => Graph(6, (0, 1, 0.6), (1, 2, 0.6), (2, 3, 0.6), (3, 4, 0.6), (4, 5, 0.6), (0, 5, 0.6), (0, 3, 0.16), (1, 4, 0.33), (2, 5, 0.81), (0, 2, 0.22)),
        "two-poles" => Graph(6, (0, 2, 1.0), (0, 3, 1.0), (0, 4, 1.0), (0, 5, 1.0), (1, 2, 1.0), (1, 3, 1.0), (1, 4, 1.0), (1, 5, 1.0)),
        _ => throw new ArgumentException(name),
    };

    // ── The sandwich, against an exhaustive search ──────────────────────────────────────────────

    [Theory]
    [MemberData(nameof(GraphNames))]
    public void OverEdges_SandwichesTheBestPartition(string name)
    {
        CsrGraph graph = Fixture(name);
        foreach ((double[,] tau, double theta) in Cases(graph))
        {
            ClosureGapResult r = ClosureGap.OverEdges(graph, Column(graph, tau), theta);
            List<(int I, int J)> edges = Edges(graph);

            double best = Partitions(graph.NodeCount).Min(p => Loss(p, edges, tau, theta));
            double free = edges.Sum(e => Math.Min(theta * (1 - tau[e.I, e.J]), (1 - theta) * tau[e.I, e.J]));

            Assert.Equal(free, r.PairRuleLoss, Exact);
            Assert.True(r.PairRuleLoss <= best + Exact && best <= r.ClosureLoss + Exact);
            Assert.Equal(r.ClosureLoss - r.PairRuleLoss, r.Gap, Exact);
            Assert.Equal(0, r.FalseSplits);    // every edge above θ is inside its closure
            Assert.Equal(edges.Count, r.Pairs);
            if (r.Gap <= Exact) Assert.Equal(best, r.ClosureLoss, Exact);
        }
    }

    [Theory]
    [MemberData(nameof(GraphNames))]
    public void OverPairs_SandwichesTheBestPartition(string name)
    {
        CsrGraph graph = Fixture(name);
        List<(int I, int J)> pairs = AllPairs(graph.NodeCount);
        List<int[]> partitions = Partitions(graph.NodeCount);
        foreach (int q in Qs)
            foreach (double beta in Betas)
            {
                DenseCoMembership dense = Dense(graph, beta, q, 20000);
                double[,] tau = Tau(dense);
                foreach (double theta in Thetas)
                {
                    ClosureGapResult r = ClosureGap.OverPairs(graph, dense, theta);

                    double best = partitions.Min(p => Loss(p, pairs, tau, theta));
                    double free = pairs.Sum(e => Math.Min(theta * (1 - tau[e.I, e.J]), (1 - theta) * tau[e.I, e.J]));

                    Assert.Equal(free, r.PairRuleLoss, Exact);
                    Assert.True(r.PairRuleLoss <= best + Exact && best <= r.ClosureLoss + Exact);
                    Assert.Equal(r.ClosureLoss - r.PairRuleLoss, r.Gap, Exact);
                    Assert.Equal(pairs.Count, r.Pairs);
                    if (r.Gap <= Exact) Assert.Equal(best, r.ClosureLoss, Exact);
                }
            }
    }

    // ── The worked chain: link loss can't see it, pair loss can ─────────────────────────────────

    [Fact]
    public void Chain_IsInvisibleOverEdges_AndCountedOverPairs()
    {
        CsrGraph graph = Fixture("path");
        double[,] tau = ExactFk(graph, Math.Log(20 * 1.5 + 1), 20);    // each bond's odds w/q = 1.5, so ⟨n_e⟩ = 0.6
        for (int j = 1; j < 6; j++)
            Assert.Equal(Math.Pow(0.6, j), tau[0, j], Exact);

        ClosureGapResult edges = ClosureGap.OverEdges(graph, Column(graph, tau), 0.5);
        Assert.Equal(1, edges.ClusterCount);
        Assert.Equal(0.0, edges.Gap, Exact);
        Assert.Equal(1.0, edges.ClosureLoss, Exact);    // 5 edges × ½ × 0.4

        // By hand, per pair distance d with 6 − d pairs at 0.6^d, θ = ½.
        double closure = 0.5 * (5 * 0.4 + 4 * 0.64 + 3 * 0.784 + 2 * 0.8704 + 0.92224);
        double gap = 4 * 0.14 + 3 * 0.284 + 2 * 0.3704 + 0.42224;
        ClosureGapResult all = ClosureGap.OverPairs(graph, Dense(graph, Math.Log(20 * 1.5 + 1), 20, 1_000_000), 0.5);
        Assert.Equal(closure, all.ClosureLoss, 1e-3);
        Assert.Equal(gap, all.Gap, 1e-3);
        Assert.Equal(10, all.FalseMerges);
        Assert.Equal(0, all.FalseSplits);

        // Searched on the exact field: three pairs, against the closure's whole path.
        Assert.Equal(2.41248, Partitions(6).Min(p => Loss(p, AllPairs(6), tau, 0.5)), Exact);
    }

    // ── The poles: a pair with no edge can pass where no edge does ──────────────────────────────

    [Fact]
    public void Poles_AreAFalseSplit_ThatTheEdgesCantSee()
    {
        CsrGraph graph = Fixture("two-poles");
        double[,] tau = ExactFk(graph, 1.0, 2);
        double edge = tau[0, 2], poles = tau[0, 1];
        Assert.True(poles > edge + 0.05);    // 0.700 through four parallel routes, against 0.647 on each edge

        const double theta = 0.67;
        ClosureGapResult onEdges = ClosureGap.OverEdges(graph, Column(graph, tau), theta);
        Assert.Equal(6, onEdges.ClusterCount);
        Assert.Equal(0.0, onEdges.Gap, Exact);

        ClosureGapResult onPairs = ClosureGap.OverPairs(graph, Dense(graph, 1.0, 2, 1_000_000), theta);
        Assert.Equal(6, onPairs.ClusterCount);
        Assert.Equal(1, onPairs.FalseSplits);
        Assert.Equal(0, onPairs.FalseMerges);
        Assert.Equal(poles - theta, onPairs.Gap, 1e-3);
    }

    [Fact]
    public void Fixtures_ReachEveryKindOfContradiction()
    {
        // The sandwich checks above would pass vacuously on fields where the closure is always right.
        bool edgeGapWithRoomToSearch = false, pairMerges = false, pairSplits = false;
        foreach (string name in AllGraphs)
        {
            CsrGraph graph = Fixture(name);
            List<(int I, int J)> edges = Edges(graph);
            List<int[]> partitions = Partitions(graph.NodeCount);
            foreach ((double[,] tau, double theta) in Cases(graph))
            {
                ClosureGapResult r = ClosureGap.OverEdges(graph, Column(graph, tau), theta);
                if (r.Gap > 1e-3 && partitions.Min(p => Loss(p, edges, tau, theta)) < r.ClosureLoss - 1e-3)
                    edgeGapWithRoomToSearch = true;
            }
            foreach (int q in Qs)
                foreach (double beta in Betas)
                {
                    DenseCoMembership dense = Dense(graph, beta, q, 20000);
                    foreach (double theta in Thetas)
                    {
                        ClosureGapResult r = ClosureGap.OverPairs(graph, dense, theta);
                        pairMerges |= r.FalseMerges > 0;
                        pairSplits |= r.FalseSplits > 0;
                    }
                }
        }
        Assert.True(edgeGapWithRoomToSearch, "no edge case where a search beats the closure");
        Assert.True(pairMerges, "no pair case with a false merge");
        Assert.True(pairSplits, "no pair case with a false split");
    }

    // ── Integrated over θ ───────────────────────────────────────────────────────────────────────

    [Theory]
    [MemberData(nameof(GraphNames))]
    public void Integrated_EqualsDistortionAgainstWidestPaths_AndTheGapsIntegral(string name)
    {
        CsrGraph graph = Fixture(name);
        foreach (int q in Qs)
            for (int b = 0; b < Betas.Length; b++)
            {
                double beta = Betas[b];
                bool integrate = q == 20 && b % 2 == 0;    // the midpoint rule is ten thousand gaps a case

                double[,] tau = ExactFk(graph, beta, q);
                double[,] m = WidestPaths(graph, tau);
                double[] column = Column(graph, tau);
                double overEdges = 0.5 * Edges(graph).Sum(e => Square(m[e.I, e.J] - tau[e.I, e.J]));
                Assert.Equal(overEdges, ClosureGap.IntegratedOverEdges(graph, column), Exact);
                if (integrate) Assert.Equal(overEdges, Midpoint(t => ClosureGap.OverEdges(graph, column, t).Gap), 1e-3);

                DenseCoMembership dense = Dense(graph, beta, q, 20000);
                double[,] denseTau = Tau(dense);
                double[,] denseM = WidestPaths(graph, denseTau);
                double overPairs = 0.5 * AllPairs(6).Sum(e => Square(denseM[e.I, e.J] - denseTau[e.I, e.J]));
                Assert.Equal(overPairs, ClosureGap.IntegratedOverPairs(graph, dense), Exact);
                if (integrate) Assert.Equal(overPairs, Midpoint(t => ClosureGap.OverPairs(graph, dense, t).Gap), 1e-3);
            }
    }

    [Fact]
    public void Integrated_ReadsPairsNoLevelJoins_AsLevelZero()
    {
        CsrGraph graph = Graph(4, (0, 1, 1.0), (2, 3, 1.0));
        var column = new double[graph.Targets.Length];
        foreach (UndirectedEdge edge in graph.UndirectedEdges()) column[edge.Slot] = 0.5;
        Assert.Equal(0.0, ClosureGap.IntegratedOverEdges(graph, column), Exact);

        // Two kept partitions: {0,1,2,3} once and all apart once, so each pair reads ½. Cross-component pairs have m = 0.
        DenseCoMembership dense = DenseCoMembership.FromPartitions(4, [[0, 0, 0, 0, 0, 1, 2, 3]]);
        Assert.Equal(0.5 * 4 * 0.25, ClosureGap.IntegratedOverPairs(graph, dense), Exact);
    }

    // ── Deciding on one estimate, evaluating on another ─────────────────────────────────────────

    [Fact]
    public void Evaluate_ReadsCostsWhileDecideSetsTheClosure()
    {
        CsrGraph graph = Graph(3, (0, 1, 1.0), (1, 2, 1.0), (0, 2, 1.0));
        double[] decide = Slots(graph, ((0, 1), 0.52), ((1, 2), 0.52), ((0, 2), 0.48));
        double[] evaluate = Slots(graph, ((0, 1), 0.48), ((1, 2), 0.48), ((0, 2), 0.52));

        Assert.Equal(0.02, ClosureGap.OverEdges(graph, decide, 0.5).Gap, Exact);

        // The same closure and the same contradicted edge, but on the other estimate that edge sits above θ.
        ClosureGapResult r = ClosureGap.OverEdges(graph, decide, 0.5, evaluate);
        Assert.Equal(1, r.ClusterCount);
        Assert.Equal(1, r.FalseMerges);
        Assert.Equal(-0.02, r.Gap, Exact);
    }

    [Fact]
    public void Noise_InflatesTheSameEstimate_AndIndependentEvaluationCancelsIt()
    {
        // K6 with every true co-membership exactly on the cut, so the true gap is zero. Each half estimates each edge
        // from 64 fair coins.
        var edges = new List<(int, int, double)>();
        foreach ((int i, int j) in AllPairs(6)) edges.Add((i, j, 1.0));
        CsrGraph graph = Graph(6, edges.ToArray());

        var rng = new Xoshiro256PlusPlus(20261101);
        const int trials = 2000;
        var pooled = new double[trials];
        var crossFit = new double[trials];
        for (int t = 0; t < trials; t++)
        {
            double[] first = Coins(graph, rng), second = Coins(graph, rng);
            double[] both = first.Zip(second, (a, b) => 0.5 * (a + b)).ToArray();
            pooled[t] = ClosureGap.OverEdges(graph, both, 0.5).Gap;
            crossFit[t] = 0.5 * (ClosureGap.OverEdges(graph, first, 0.5, second).Gap
                               + ClosureGap.OverEdges(graph, second, 0.5, first).Gap);
        }

        Assert.True(pooled.Average() > 0.2, $"plug-in mean {pooled.Average()}");
        double mean = crossFit.Average();
        double se = Math.Sqrt(crossFit.Sum(g => Square(g - mean)) / (trials - 1) / trials);
        Assert.True(Math.Abs(mean) <= 4 * se, $"cross-fit mean {mean}, standard error {se}");
    }

    private static double[] Coins(CsrGraph graph, Xoshiro256PlusPlus rng)
    {
        var column = new double[graph.Targets.Length];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
            column[edge.Slot] = BitOperations.PopCount(rng.NextUInt64()) / 64.0;
        return column;
    }

    // ── DenseCoMembership ───────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(1)]
    [InlineData(7)]
    public void Dense_ReproducesTheEdgeCounts_FromASweepsKeptDraws(int every)
    {
        CsrGraph graph = SwAccumulatorProbeTests.Random16();
        Accumulator acc = Run(graph, 0.6, 20, new AccumulationSpec { CoMembership = true, PartitionEvery = every }, 3000, 11, 0);
        DenseCoMembership dense = DenseCoMembership.FromPartitions(16, [acc.KeptPartitions!]);

        Assert.Equal((acc.DrawCount + every - 1) / every, dense.PartitionCount);
        if (every == 1)
        {
            double[] column = dense.EdgeColumn(graph);
            foreach (UndirectedEdge edge in graph.UndirectedEdges())
                Assert.Equal((double)acc.CoMembershipCount![edge.Slot] / acc.DrawCount, column[edge.Slot], Exact);
        }

        // Every pair against a direct count over the kept partitions.
        int[] labels = acc.KeptPartitions!;
        for (int i = 0; i < 16; i++)
            for (int j = 0; j < 16; j++)
            {
                int together = 0;
                for (int p = 0; p < dense.PartitionCount; p++)
                    if (labels[p * 16 + i] == labels[p * 16 + j]) together++;
                Assert.Equal((double)together / dense.PartitionCount, dense.At(i, j), Exact);
            }
    }

    [Fact]
    public void Dense_Pool_EqualsAccumulatingTogether()
    {
        int[] a = [0, 0, 2, 2, 0, 1, 1, 3];
        int[] b = [0, 1, 1, 3];
        DenseCoMembership pooled = DenseCoMembership.Pool(DenseCoMembership.FromPartitions(4, [a]), DenseCoMembership.FromPartitions(4, [b]));
        DenseCoMembership together = DenseCoMembership.FromPartitions(4, [a, b]);

        Assert.Equal(3, pooled.PartitionCount);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                Assert.Equal(together.At(i, j), pooled.At(i, j), Exact);
    }

    [Fact]
    public void Dense_RefusesWhatIsNotWholePartitionsInRange()
    {
        Assert.Throws<ArgumentException>(() => DenseCoMembership.FromPartitions(4, [[0, 0, 0]]));
        Assert.Throws<ArgumentException>(() => DenseCoMembership.FromPartitions(4, [[0, 0, 0, 4]]));
        Assert.Throws<ArgumentException>(() => DenseCoMembership.FromPartitions(4, [[0, 0, 0, -1]]));
        Assert.Throws<ArgumentException>(() => DenseCoMembership.FromPartitions(4, Array.Empty<int[]>()));
    }

    // ── The probe's search, against every partition ─────────────────────────────────────────────

    [Fact]
    public void Search_MatchesEveryPartition_OnTheFixtures()
    {
        // Fixed before the first run: within 1e-10 of the exhaustive best on at least 99% of the cells.
        int cells = 0, misses = 0;
        foreach (string name in AllGraphs)
        {
            CsrGraph graph = Fixture(name);
            List<(int I, int J)> pairs = AllPairs(graph.NodeCount);
            List<int[]> partitions = Partitions(graph.NodeCount);
            foreach (int q in Qs)
                foreach (double beta in Betas)
                {
                    DenseCoMembership dense = Dense(graph, beta, q, 20000);
                    double[,] tau = Tau(dense);
                    foreach (double theta in Thetas)
                    {
                        int[] closure = AffinityThreshold.Connect(graph, dense.EdgeColumn(graph), theta).Labels;
                        PairLossSearch.Found found = PairLossSearch.Best(dense, theta, closure, 4, cells);

                        double foundLoss = Loss(found.Labels, pairs, tau, theta), closureLoss = Loss(closure, pairs, tau, theta);
                        double best = partitions.Min(p => Loss(p, pairs, tau, theta));

                        Assert.Equal(PairLossSearch.Score(dense, found.Labels, theta), found.Score, Exact);
                        Assert.Equal(closureLoss - foundLoss, found.Score - PairLossSearch.Score(dense, closure, theta), Exact);
                        Assert.True(foundLoss <= closureLoss + Exact);
                        Assert.True(foundLoss >= best - Exact);

                        cells++;
                        if (foundLoss > best + Exact) misses++;
                    }
                }
        }
        Assert.True(misses <= 0.01 * cells, $"{misses} of {cells} cells short of the best partition");
    }

    [Fact]
    public void Search_MatchesEveryPartition_OnRandomMixtures()
    {
        // Fixed before the first run: within 1e-10 of the exhaustive best on at least 90% of 200 fields of nine nodes,
        // and no field short by over 5% of its best score.
        const int nodes = 9, fields = 200;
        List<int[]> partitions = Partitions(nodes);
        var rng = new Xoshiro256PlusPlus(20261103);
        int misses = 0;
        double worst = 0.0;
        for (int field = 0; field < fields; field++)
        {
            var run = new int[nodes * (3 + (int)(rng.NextDouble() * 8))];
            for (int p = 0; p < run.Length; p += nodes)
            {
                int blocks = 2 + (int)(rng.NextDouble() * 4);
                for (int i = 0; i < nodes; i++) run[p + i] = (int)(rng.NextDouble() * blocks);
            }
            DenseCoMembership dense = DenseCoMembership.FromPartitions(nodes, [run]);
            double theta = Thetas[field % Thetas.Length];

            var weights = new double[nodes, nodes];
            for (int i = 0; i < nodes; i++)
                for (int j = i + 1; j < nodes; j++)
                    weights[i, j] = dense.At(i, j) - theta;
            double best = double.NegativeInfinity;
            foreach (int[] p in partitions)
            {
                double score = 0.0;
                for (int i = 0; i < nodes; i++)
                    for (int j = i + 1; j < nodes; j++)
                        if (p[i] == p[j]) score += weights[i, j];
                best = Math.Max(best, score);
            }

            var start = new int[nodes];    // no graph here, so the closure's start is the single cluster
            PairLossSearch.Found found = PairLossSearch.Best(dense, theta, start, 4, field);
            Assert.True(found.Score <= best + Exact);
            if (found.Score < best - Exact)
            {
                misses++;
                if (best > 0.0) worst = Math.Max(worst, (best - found.Score) / best);
            }
        }
        Assert.True(misses <= 0.1 * fields, $"{misses} of {fields} fields short of the best partition");
        Assert.True(worst <= 0.05, $"short by {worst} of the best score");
    }

    [Fact]
    public void Search_PartsPlantedBlocks_WhereClosureChains()
    {
        // Blocks A and B of 20 and one bridge: each pair across reads 0.3, the bridge's ends 0.6 and A's other pairs at
        // least 0.7. At θ = ½ the closure takes the bridge and joins all. The rule's score bounds any partition at 178.7,
        // and a partition holding the bridge's ends together gives up 3.8 for 0.1, so the best is exactly A | B, at 178.6.
        const int size = 20, nodes = 2 * size;
        var run = new List<int>();
        for (int copy = 0; copy < 3; copy++) run.AddRange(new int[nodes]);
        for (int copy = 0; copy < 3; copy++)
            for (int i = 0; i < nodes; i++)
                run.Add(i == 0 || i == size ? 0 : i < size ? 1 : 2);
        for (int copy = 0; copy < 4; copy++)
            for (int i = 0; i < nodes; i++)
                run.Add(i < size ? 0 : 1);
        DenseCoMembership dense = DenseCoMembership.FromPartitions(nodes, [run.ToArray()]);

        var edges = new List<(int, int, double)> { (0, size, 1.0) };
        for (int i = 1; i < size; i++)
        {
            edges.Add((i - 1, i, 1.0));
            edges.Add((size + i - 1, size + i, 1.0));
        }
        CsrGraph graph = Graph(nodes, edges.ToArray());
        int[] closure = AffinityThreshold.Connect(graph, dense.EdgeColumn(graph), 0.5).Labels;
        Assert.Single(closure.Distinct());

        PairLossSearch.Found found = PairLossSearch.Best(dense, 0.5, closure, 4, 1);
        Assert.Equal(178.6, found.Score, 1e-9);
        for (int i = 0; i < nodes; i++)
            Assert.Equal(i < size, found.Labels[i] == found.Labels[0]);
    }

    [Fact]
    public void Polish_LeavesANodeAlone_WhenEveryClusterCostsIt()
    {
        // Node 0 reads 0.1 with {1, 2} and 0.2 with {3, 4}, so at θ = ½ its own cluster costs 0.8 and the other 0.6. Alone
        // costs nothing, and moving to the cheaper cluster would still be a loss.
        var run = new List<int>();
        run.AddRange([0, 0, 0, 3, 3]);
        for (int copy = 0; copy < 2; copy++) run.AddRange([0, 1, 1, 0, 0]);
        for (int copy = 0; copy < 7; copy++) run.AddRange([0, 1, 1, 3, 3]);
        DenseCoMembership dense = DenseCoMembership.FromPartitions(5, [run.ToArray()]);

        PairLossSearch.Found found = PairLossSearch.Polished(dense, 0.5, [0, 0, 0, 3, 3]);

        Assert.Equal(1.0, found.Score, Exact);    // {1, 2} and {3, 4}, at ½ each
        Assert.Equal(3, found.Labels.Distinct().Count());
        Assert.DoesNotContain(found.Labels[0], found.Labels.Skip(1));
    }

    [Fact]
    public void Search_IsRepeatable_AndNoWorseThanClosure_OnASampledField()
    {
        CsrGraph graph = SwAccumulatorProbeTests.Random16();
        Accumulator acc = Run(graph, 0.6, 20, new AccumulationSpec { PartitionEvery = 1 }, 2000, 11, 0);
        DenseCoMembership dense = DenseCoMembership.FromPartitions(16, [acc.KeptPartitions!]);
        int[] closure = AffinityThreshold.Connect(graph, dense.EdgeColumn(graph), 0.474).Labels;

        PairLossSearch.Found a = PairLossSearch.Best(dense, 0.474, closure, 4, 7), b = PairLossSearch.Best(dense, 0.474, closure, 4, 7);
        Assert.Equal(a.Labels, b.Labels);
        Assert.Equal(a.Score, b.Score);
        Assert.True(a.Score >= PairLossSearch.Score(dense, closure, 0.474) - Exact);
    }

    // ── The sampler against the enumeration: the measure the checks assume is the one SW samples ──

    [Theory]
    [InlineData(0.17)]
    [InlineData(0.35)]
    public void Sampled_CoMembership_MatchesTheExactMeasure(double temperature)
    {
        CsrGraph graph = Fixture("strong-weak");
        double[,] exact = ExactFk(graph, 1.0 / temperature, 20);
        Accumulator acc = Run(graph, temperature, 20, new AccumulationSpec { PartitionEvery = 1 }, 40000, 29, 0);
        DenseCoMembership dense = DenseCoMembership.FromPartitions(6, [acc.KeptPartitions!]);

        foreach ((int i, int j) in AllPairs(6))
            Assert.InRange(dense.At(i, j) - exact[i, j], -0.02, 0.02);
    }

    // ── The sweep view ──────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Curve_ScoresTheTStacksOwnCut()
    {
        CsrGraph graph = SwAccumulatorProbeTests.Random16();
        List<Accumulator> frames = Sweep(graph, new AccumulationSpec { CoMembership = true, PartitionEvery = 1 }, replicas: 4);

        ClosureGapCurve edges = SweepClosureGap.OverEdges(graph, frames);
        ClosureGapCurve pairs = SweepClosureGap.OverPairs(graph, frames);
        PartitionHierarchy stack = SweepDenseTStack.Build(graph, frames);

        Assert.Equal((20 * 0.5 - 1) / 19, edges.Theta, Exact);
        Assert.Equal(stack.Levels.Select(l => l.Temperature), edges.Temperatures);
        for (int t = 0; t < edges.Temperatures.Length; t++)
        {
            Assert.Equal(stack.Levels[t].Partition.Count, edges.Results[t].ClusterCount);
            Assert.Equal(stack.Levels[t].Partition.Count, pairs.Results[t].ClusterCount);    // every draw kept
            Assert.Equal(34, edges.Results[t].Pairs);
            Assert.Equal(120, pairs.Results[t].Pairs);
        }
    }

    [Fact]
    public void Curve_CrossFit_SplitsReplicasAlternately_AndEachWayRound()
    {
        CsrGraph graph = SwAccumulatorProbeTests.Random16();
        List<Accumulator> frames = Sweep(graph, new AccumulationSpec { CoMembership = true, PartitionEvery = 1 }, replicas: 4);
        ClosureGapCurve edges = SweepClosureGap.OverEdges(graph, frames);
        ClosureGapCurve pairs = SweepClosureGap.OverPairs(graph, frames);

        for (int t = 0; t < edges.Temperatures.Length; t++)
        {
            Accumulator[] at = frames.Where(f => f.Temperature == edges.Temperatures[t]).OrderBy(f => f.ReplicaIndex).ToArray();
            Accumulator[] a = [at[0], at[2]], b = [at[1], at[3]];

            double[] first = Pooled(a), second = Pooled(b);
            double byEdges = 0.5 * (ClosureGap.OverEdges(graph, first, edges.Theta, second).Gap
                                  + ClosureGap.OverEdges(graph, second, edges.Theta, first).Gap);
            Assert.Equal(byEdges, edges.CrossFitGap[t], Exact);

            DenseCoMembership x = DenseCoMembership.FromPartitions(16, a.Select(f => f.KeptPartitions!).ToArray());
            DenseCoMembership y = DenseCoMembership.FromPartitions(16, b.Select(f => f.KeptPartitions!).ToArray());
            double byPairs = 0.5 * (ClosureGap.OverPairs(graph, x, pairs.Theta, y).Gap
                                  + ClosureGap.OverPairs(graph, y, pairs.Theta, x).Gap);
            Assert.Equal(byPairs, pairs.CrossFitGap[t], Exact);
        }
    }

    [Fact]
    public void Curve_DoesNotDependOnTheFramesOrder()
    {
        CsrGraph graph = SwAccumulatorProbeTests.Random16();
        List<Accumulator> frames = Sweep(graph, new AccumulationSpec { CoMembership = true }, replicas: 3);
        var reversed = new List<Accumulator>(frames);
        reversed.Reverse();

        ClosureGapCurve a = SweepClosureGap.OverEdges(graph, frames), b = SweepClosureGap.OverEdges(graph, reversed);

        Assert.Equal(a.Temperatures, b.Temperatures);
        Assert.Equal(a.Results, b.Results);
        Assert.Equal(a.CrossFitGap, b.CrossFitGap);
        Assert.Equal(a.Integrated, b.Integrated);
    }

    [Fact]
    public void Curve_OneReplica_HasNoCrossFit()
    {
        CsrGraph graph = SwAccumulatorProbeTests.Random16();
        List<Accumulator> frames = Sweep(graph, new AccumulationSpec { CoMembership = true, PartitionEvery = 1 }, replicas: 1);

        Assert.All(SweepClosureGap.OverEdges(graph, frames).CrossFitGap, g => Assert.True(double.IsNaN(g)));
        Assert.All(SweepClosureGap.OverPairs(graph, frames).CrossFitGap, g => Assert.True(double.IsNaN(g)));
    }

    [Fact]
    public void Curve_NamesTheOptionItNeeds()
    {
        CsrGraph graph = SwAccumulatorProbeTests.Random16();
        List<Accumulator> frames = Sweep(graph, AccumulationSpec.None, replicas: 1);

        Assert.Contains("CoMembership", Assert.Throws<InvalidOperationException>(() => SweepClosureGap.OverEdges(graph, frames)).Message);
        Assert.Contains("PartitionEvery", Assert.Throws<InvalidOperationException>(() => SweepClosureGap.OverPairs(graph, frames)).Message);
    }

    [Theory]
    [InlineData(0.04, 20)]
    [InlineData(1.01, 20)]
    [InlineData(0.5, 1)]
    public void CoMembershipTheta_RefusesWhatHasNoCostRatio(double theta, int q)
        => Assert.Throws<ArgumentOutOfRangeException>(() => SweepClosureGap.CoMembershipTheta(theta, q));

    // ── Ground truth, and helpers ───────────────────────────────────────────────────────────────

    /// <summary>Every bond configuration's partition and its probability, from the FK weight Π_open (e^{βJ} − 1) · q^k.</summary>
    private static List<(int[] Labels, double Weight)> FkLaw(CsrGraph graph, double beta, int q)
    {
        int n = graph.NodeCount;
        var edges = new List<(int I, int J, double W)>();
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
            edges.Add((edge.Source, edge.Target, Math.Exp(beta * graph.Weights[edge.Slot]) - 1.0));

        var law = new List<(int[] Labels, double Weight)>(1 << edges.Count);
        double total = 0.0;
        var parent = new int[n];
        for (int mask = 0; mask < 1 << edges.Count; mask++)
        {
            for (int i = 0; i < n; i++) parent[i] = i;
            double weight = 1.0;
            for (int e = 0; e < edges.Count; e++)
                if ((mask >> e & 1) == 1)
                {
                    weight *= edges[e].W;
                    parent[Find(parent, edges[e].I)] = Find(parent, edges[e].J);
                }
            var labels = new int[n];
            int k = 0;
            for (int i = 0; i < n; i++)
            {
                labels[i] = Find(parent, i);
                if (labels[i] == i) k++;
            }
            weight *= Math.Pow(q, k);
            total += weight;
            law.Add((labels, weight));
        }
        return law.Select(entry => (entry.Labels, entry.Weight / total)).ToList();
    }

    /// <summary>Pair co-membership under the FK measure.</summary>
    private static double[,] ExactFk(CsrGraph graph, double beta, int q)
    {
        int n = graph.NodeCount;
        var tau = new double[n, n];
        foreach ((int[] labels, double weight) in FkLaw(graph, beta, q))
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    if (labels[i] == labels[j]) tau[i, j] += weight;
        return tau;
    }

    /// <summary>The FK law as kept partitions: each configuration's partition, repeated in proportion to its weight.</summary>
    private static DenseCoMembership Dense(CsrGraph graph, double beta, int q, int scale)
    {
        int n = graph.NodeCount;
        var run = new List<int>();
        foreach ((int[] labels, double weight) in FkLaw(graph, beta, q))
            for (int copy = (int)Math.Round(weight * scale); copy > 0; copy--)
                run.AddRange(labels);
        return DenseCoMembership.FromPartitions(n, [run.ToArray()]);
    }

    private static int Find(int[] parent, int i)
    {
        while (parent[i] != i) i = parent[i];
        return i;
    }

    /// <summary>The level at which closure first joins each pair: the max-min of the field over graph paths.</summary>
    private static double[,] WidestPaths(CsrGraph graph, double[,] tau)
    {
        int n = graph.NodeCount;
        var m = new double[n, n];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
            m[edge.Source, edge.Target] = m[edge.Target, edge.Source] = tau[edge.Source, edge.Target];
        for (int k = 0; k < n; k++)
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    m[i, j] = Math.Max(m[i, j], Math.Min(m[i, k], m[k, j]));
        return m;
    }

    private static IEnumerable<(double[,] Tau, double Theta)> Cases(CsrGraph graph)
    {
        foreach (int q in Qs)
            foreach (double beta in Betas)
            {
                double[,] tau = ExactFk(graph, beta, q);
                foreach (double theta in Thetas)
                    yield return (tau, theta);
            }
    }

    private static double Loss(int[] partition, List<(int I, int J)> pairs, double[,] tau, double theta)
        => pairs.Sum(e => partition[e.I] == partition[e.J] ? theta * (1 - tau[e.I, e.J]) : (1 - theta) * tau[e.I, e.J]);

    private static List<int[]> Partitions(int n)
    {
        var all = new List<int[]>();
        var labels = new int[n];
        void Fill(int i, int blocks)
        {
            if (i == n)
            {
                all.Add((int[])labels.Clone());
                return;
            }
            for (int c = 0; c <= blocks; c++)
            {
                labels[i] = c;
                Fill(i + 1, Math.Max(blocks, c + 1));
            }
        }
        Fill(0, 0);
        return all;
    }

    private static double[,] Tau(DenseCoMembership dense)
    {
        var tau = new double[dense.NodeCount, dense.NodeCount];
        for (int i = 0; i < dense.NodeCount; i++)
            for (int j = 0; j < dense.NodeCount; j++)
                tau[i, j] = dense.At(i, j);
        return tau;
    }

    private static double Midpoint(Func<double, double> f)
    {
        const int steps = 10000;
        double sum = 0.0;
        for (int i = 0; i < steps; i++) sum += f((i + 0.5) / steps);
        return sum / steps;
    }

    private static double Square(double x) => x * x;

    private static List<(int I, int J)> AllPairs(int n)
    {
        var pairs = new List<(int, int)>();
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                pairs.Add((i, j));
        return pairs;
    }

    private static List<(int I, int J)> Edges(CsrGraph graph)
    {
        var edges = new List<(int, int)>();
        foreach (UndirectedEdge edge in graph.UndirectedEdges()) edges.Add((edge.Source, edge.Target));
        return edges;
    }

    private static double[] Column(CsrGraph graph, double[,] tau)
    {
        var column = new double[graph.Targets.Length];
        foreach (UndirectedEdge edge in graph.UndirectedEdges()) column[edge.Slot] = tau[edge.Source, edge.Target];
        return column;
    }

    private static double[] Slots(CsrGraph graph, params ((int I, int J) Edge, double Value)[] values)
    {
        var column = new double[graph.Targets.Length];
        foreach (((int i, int j), double value) in values) column[graph.FindSlot(i, j)] = value;
        return column;
    }

    private static double[] Pooled(Accumulator[] frames)
    {
        var column = new double[frames[0].CoMembershipCount!.Length];
        long draws = frames.Sum(f => (long)f.DrawCount);
        foreach (Accumulator frame in frames)
            for (int e = 0; e < column.Length; e++)
                column[e] += (double)frame.CoMembershipCount![e] / draws;
        return column;
    }

    private static CsrGraph Graph(int nodes, params (int I, int J, double W)[] edges)
        => CsrGraph.FromEdges(edges.Select(e => new Edge(e.I, e.J, e.W)).ToArray(), nodes);

    private static Accumulator Run(CsrGraph graph, double temperature, int q, AccumulationSpec spec, int cycles, int seed, int replica)
        => SwRunner.Run(new SwRunSpec
        {
            Graph = graph, Temperature = temperature, Q = q, Accumulation = spec, Seed = seed, ReplicaIndex = replica,
            Budget = new RunBudget(200, cycles),
        }).Accumulator;

    private static List<Accumulator> Sweep(CsrGraph graph, AccumulationSpec spec, int replicas)
    {
        var frames = new List<Accumulator>();
        int seed = 1000;
        foreach (double t in new[] { 0.1, 0.3, 0.6, 1.2 })
            for (int r = 0; r < replicas; r++)
                frames.Add(Run(graph, t, 20, spec, 300, seed++, r));
        return frames;
    }
}
