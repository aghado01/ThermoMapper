using System;
using System.Collections.Generic;
using Clustering.Primitives;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Core.Solver;
using Graphs.Primitives;
using Xunit;

using Clustering.Graphical.SPC.Runtime.Core;

using Graphs;

namespace Clustering.Graphical.SPC.Tests.Sampler;

/// <summary>
/// LocalField (per-site energy ladder) tests: implementation correctness against
/// an independent per-site reimplementation, the <see cref="SymmetrizationRule"/>
/// mode semantics, and a demonstration that LocalField genuinely diverges from
/// MeanField (i.e. it is not a relabeled global sort).
/// </summary>
public sealed class PKWangLocalFieldTests
{
    // Triangle {0,1,2} + bridge to a pair {3,4}; distinct couplings.
    private static CsrGraph TripletGraph()
    {
        var edges = new[]
        {
            new Edge(0, 1, 10.0),
            new Edge(1, 2, 9.0),
            new Edge(0, 2, 5.0),
            new Edge(3, 4, 8.0),
            new Edge(2, 3, 1.0),
        };
        return CsrGraph.FromEdges(edges, nodeCount: 5);
    }

    // Path 0-1-2-3-4 with strictly increasing couplings — local and global
    // edge rankings diverge here, which is what separates LocalField from
    // MeanField.
    private static CsrGraph IncreasingPath()
    {
        var edges = new[]
        {
            new Edge(0, 1, 1.0),
            new Edge(1, 2, 2.0),
            new Edge(2, 3, 3.0),
            new Edge(3, 4, 4.0),
        };
        return CsrGraph.FromEdges(edges, nodeCount: 5);
    }

    private static int SlotOf(CsrGraph g, int i, int j)
    {
        int rowEnd = g.RowPointers[i + 1];
        for (int e = g.RowPointers[i]; e < rowEnd; e++)
            if (g.Targets[e] == j) return e;
        throw new InvalidOperationException($"No CSR slot for edge ({i},{j}).");
    }

    // Independent per-site cumulative ladder (the test's own reimplementation;
    // does not call LocalField.BuildHcum). Returns slot-indexed directed Hcum.
    private static double[] IndependentPerSiteHcum(CsrGraph g)
    {
        var hcum = new double[g.Targets.Length];
        for (int i = 0; i < g.NodeCount; i++)
        {
            int start = g.RowPointers[i];
            int deg = g.RowPointers[i + 1] - start;
            var local = new List<(int Slot, double Coupling)>(deg);
            for (int t = 0; t < deg; t++)
                local.Add((start + t, g.Weights[start + t]));
            local.Sort((a, b) => a.Coupling.CompareTo(b.Coupling));
            double cum = 0.0;
            foreach (var (slot, coupling) in local)
            {
                cum += coupling;
                hcum[slot] = cum;
            }
        }
        return hcum;
    }

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
    public void MutualG_MatchesIndependentPerSiteLadder()
    {
        CsrGraph g = TripletGraph();
        const double T = 6.0;
        double[] hcum = IndependentPerSiteHcum(g);

        PKWangContext ctx = PKWang.Prepare(g, EdgeWeightKind.Coupling, Field.Local, SymmetrizationRule.Mutual);
        Affinities corr = PKWang.Solve(ctx, T);

        for (int i = 0; i < g.NodeCount; i++)
        {
            int rowEnd = g.RowPointers[i + 1];
            for (int e = g.RowPointers[i]; e < rowEnd; e++)
            {
                int j = g.Targets[e];
                if (j <= i) continue;

                double gij = 1.0 - Math.Exp(-hcum[e] / T);                 // i → j
                double gji = 1.0 - Math.Exp(-hcum[SlotOf(g, j, i)] / T);   // j → i
                double expected = Math.Min(gij, gji);                       // Mutual

                Assert.Equal(expected, corr.G[e], precision: 12);
            }
        }
    }

    [Fact]
    public void Modes_OrderClusterCounts()
    {
        CsrGraph g = TripletGraph();
        const double T = 6.0;

        int inclusive = PKWang.Cluster(PKWang.Prepare(g, EdgeWeightKind.Coupling, Field.Local, SymmetrizationRule.Inclusive), T).Count;
        int mean = PKWang.Cluster(PKWang.Prepare(g, EdgeWeightKind.Coupling, Field.Local, SymmetrizationRule.Mean), T).Count;
        int mutual = PKWang.Cluster(PKWang.Prepare(g, EdgeWeightKind.Coupling, Field.Local, SymmetrizationRule.Mutual), T).Count;

        // max G ≥ mean G ≥ min G ⇒ Inclusive keeps the most edges (fewest
        // clusters), Mutual the fewest (most clusters).
        Assert.True(inclusive <= mean, $"Inclusive ({inclusive}) should not exceed Mean ({mean}).");
        Assert.True(mean <= mutual, $"Mean ({mean}) should not exceed Mutual ({mutual}).");
    }

    [Fact]
    public void ClusterCount_IsMonotoneInTemperature()
    {
        CsrGraph g = TripletGraph();
        PKWangContext ctx = PKWang.Prepare(g, EdgeWeightKind.Coupling, Field.Local, SymmetrizationRule.Mutual);

        int prev = 0;
        foreach (double T in new[] { 0.7, 3.0, 8.0, 20.0, 60.0 })
        {
            int count = PKWang.Cluster(ctx, T).Count;
            Assert.True(count >= prev, $"Cluster count fell from {prev} to {count} at T={T}.");
            prev = count;
        }
    }

    [Fact]
    public void LocalField_DivergesFromMeanField()
    {
        CsrGraph g = IncreasingPath();
        const double T = 3.6; // T·ln2 ≈ 2.50

        Assignment mean = PKWang.Cluster(PKWang.Prepare(g, EdgeWeightKind.Coupling, Field.Mean), T);
        Assignment local = PKWang.Cluster(PKWang.Prepare(g, EdgeWeightKind.Coupling, Field.Local, SymmetrizationRule.Mutual), T);

        // MeanField cumulates globally → keeps (1,2),(2,3),(3,4): {0}|{1,2,3,4}.
        // LocalField cumulates per-site → also cuts (1,2): {0}|{1}|{2,3,4}.
        Assert.False(SamePartition(mean.Labels, local.Labels),
            "LocalField should not reduce to MeanField on a graph where local and global rankings differ.");
        Assert.True(local.Count > mean.Count,
            "Per-site cumulation cuts at least as aggressively as the global pool here.");
    }
}
