using System;
using System.Collections.Generic;
using Clustering.Dendrograms;
using Clustering.Primitives;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Core.Solver;
using Graphs.Primitives;
using Xunit;

using Clustering.Graphical.SPC.Runtime.Core;

using Graphs;

using Clustering.Graphical.SPC.Demos;

namespace Clustering.Graphical.SPC.Tests.Sampler;

/// <summary>
/// Empirical adjudication of the Wang 2020 reduction: the paper's M-draw Monte
/// Carlo estimates a quantity with a closed form, and the clustering it induces
/// is thermal single-linkage. Both reference computations here — the M-draw
/// inverse-transform MC and the Kruskal single-linkage — are reimplemented
/// independently of the PKWang kernel, so agreement is genuine cross-validation,
/// not a tautology (validation independence).
/// </summary>
public sealed class PKWangMeanFieldTests
{
    // Two well-separated triplets {0,1,2} and {3,4,5} joined by one weak bridge
    // (2,3). DISTINCT couplings → unambiguous edge ranking (no ties). Coupling J
    // is higher for stronger/closer pairs, per CsrGraph weight semantics.
    private static CsrGraph BuildGraph()
    {
        var edges = new[]
        {
            new Edge(0, 1, 10.0),
            new Edge(1, 2, 9.0),
            new Edge(0, 2, 5.0),   // redundant within-cluster edge (triangle)
            new Edge(3, 4, 8.0),
            new Edge(4, 5, 7.0),
            new Edge(2, 3, 1.0),   // weak bridge between the two triplets
        };
        return CsrGraph.FromEdges(edges, nodeCount: 6);
    }

    // Independent ascending cumulative ladder, computed straight from the graph
    // (the test's own reimplementation — does not call MeanField.BuildHcum).
    private static List<(int I, int J, double Hcum)> IndependentLadder(CsrGraph g)
    {
        var uniq = new List<(int I, int J, double Coupling)>();
        for (int i = 0; i < g.NodeCount; i++)
        {
            int rowEnd = g.RowPointers[i + 1];
            for (int e = g.RowPointers[i]; e < rowEnd; e++)
                if (g.Targets[e] > i)
                    uniq.Add((i, g.Targets[e], g.Weights[e]));
        }
        uniq.Sort((a, b) => a.Coupling.CompareTo(b.Coupling));

        var ladder = new List<(int, int, double)>(uniq.Count);
        double cum = 0.0;
        foreach (var (i, j, coupling) in uniq)
        {
            cum += coupling;
            ladder.Add((i, j, cum));
        }
        return ladder;
    }

    private static int SlotOf(CsrGraph g, int i, int j)
    {
        int rowEnd = g.RowPointers[i + 1];
        for (int e = g.RowPointers[i]; e < rowEnd; e++)
            if (g.Targets[e] == j) return e;
        throw new InvalidOperationException($"No CSR slot for edge ({i},{j}).");
    }

    // Two label arrays describe the same partition iff they induce the same
    // same-cluster equivalence relation (relabeling-invariant).
    private static bool SamePartition(int[] a, int[] b)
    {
        if (a.Length != b.Length) return false;
        for (int i = 0; i < a.Length; i++)
            for (int j = i + 1; j < a.Length; j++)
                if ((a[i] == a[j]) != (b[i] == b[j]))
                    return false;
        return true;
    }

    [Fact]
    public void ClosedForm_MatchesMonteCarlo()
    {
        CsrGraph g = BuildGraph();
        const double T = 6.0;
        PKWangContext ctx = PKWang.Prepare(g, EdgeWeightKind.Coupling, Field.Mean);

        Affinities closed = PKWang.Solve(ctx, T);
        Affinities mc = WangMonteCarloDemo.Sample(ctx, T, draws: 200_000, seed: 20260607);

        // Method equivalence (draws vs CDF) at every edge — the empirical
        // statement that the paper's Monte Carlo was estimating a closed form —
        // plus a cross-check of the closed form against the analytic
        // 1 - exp(-Hcum/T) using an independently recomputed ladder (pins
        // BuildHcum, not merely the kernel).
        foreach (var (i, j, hcum) in IndependentLadder(g))
        {
            int slot = SlotOf(g, i, j);
            Assert.Equal(mc.G[slot], closed.G[slot], precision: 2);
            Assert.Equal(1.0 - Math.Exp(-hcum / T), closed.G[slot], precision: 12);
        }
    }

    [Fact]
    public void HardCut_IsHcumExceedsTLn2()
    {
        CsrGraph g = BuildGraph();
        const double T = 9.0;
        double cut = T * Math.Log(2.0);

        PKWangContext ctx = PKWang.Prepare(g, EdgeWeightKind.Coupling, Field.Mean);
        Affinities corr = PKWang.Solve(ctx, T);

        foreach (var (i, j, hcum) in IndependentLadder(g))
        {
            bool active = corr.G[SlotOf(g, i, j)] > 0.5;
            Assert.Equal(hcum > cut, active);
        }
    }

    [Theory]
    [InlineData(0.7)]    // T·ln2 ≈ 0.49 — every edge active, one cluster
    [InlineData(4.5)]    // bridge (Hcum=1) cut → {0,1,2} | {3,4,5}
    [InlineData(14.0)]   // also cuts (0,2),(4,5)
    [InlineData(36.0)]   // only the two strongest survive
    [InlineData(72.0)]   // all edges cut → six singletons
    public void MeanField_EqualsSingleLinkage(double T)
    {
        CsrGraph g = BuildGraph();
        Assignment pk = PKWang.Cluster(PKWang.Prepare(g, EdgeWeightKind.Coupling, Field.Mean), T);

        // Independent cross-subsystem oracle: the canonical MST → dendrogram
        // single-linkage, cut to the same cluster count PKWang produced. The
        // theorem says PKWang's family IS the single-linkage family, so the two
        // k-partitions must coincide (distinct couplings ⇒ a unique k-cut).
        int[] singleLinkage = SingleLinkage.FromCouplingGraph(g).CutToK(pk.Count);

        Assert.True(SamePartition(pk.Labels, singleLinkage),
            $"MeanField partition diverged from single-linkage at T={T}.");
    }

    [Fact]
    public void ClusterCount_IsMonotoneInTemperature()
    {
        CsrGraph g = BuildGraph();
        PKWangContext ctx = PKWang.Prepare(g, EdgeWeightKind.Coupling, Field.Mean);

        int prev = 0;
        foreach (double T in new[] { 0.7, 4.5, 14.0, 36.0, 72.0 })
        {
            int count = PKWang.Cluster(ctx, T).Count;
            Assert.True(count >= prev, $"Cluster count fell from {prev} to {count} at T={T}.");
            prev = count;
        }
    }
}
