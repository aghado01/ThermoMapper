using System;
using System.Collections.Generic;
using Clustering.Evaluation.External;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// Rand index and link overlap against independent references: brute-force
/// agreement counts over all node pairs (Rand) and over the edge list a graph was
/// built from (link overlap), plus the identity that ties them — link overlap on
/// the complete graph is the Rand index.
/// </summary>
public sealed class RandAndLinkOverlapTests
{
    // Agreement over all unordered pairs, no contingency table.
    private static double BruteRand(int[] a, int[] b)
    {
        int agree = 0, total = 0;
        for (int i = 0; i < a.Length; i++)
            for (int j = i + 1; j < a.Length; j++)
            {
                total++;
                if ((a[i] == a[j]) == (b[i] == b[j])) agree++;
            }
        return (double)agree / total;
    }

    // Agreement over an explicit edge list, weighted or not, independent of the CSR walk.
    private static double BruteLinkOverlap(Edge[] edges, int[] a, int[] b, bool weighted)
    {
        double agree = 0.0, total = 0.0;
        foreach (Edge e in edges)
        {
            double w = weighted ? e.J : 1.0;
            total += w;
            if ((a[e.Source] == a[e.Target]) == (b[e.Source] == b[e.Target])) agree += w;
        }
        return agree / total;
    }

    private static CsrGraph Complete(int n)
    {
        var edges = new List<Edge>();
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                edges.Add(new Edge(i, j, 1.0));
        return CsrGraph.FromEdges(edges.ToArray(), nodeCount: n);
    }

    // Path 0-1-2-3-4-5 with weak couplings on (1,2) and (2,3).
    private static Edge[] PathEdges() => new[]
    {
        new Edge(0, 1, 1.0),
        new Edge(1, 2, 0.25),
        new Edge(2, 3, 0.25),
        new Edge(3, 4, 1.0),
        new Edge(4, 5, 1.0),
    };

    public static IEnumerable<object[]> LabelPairs => new[]
    {
        new object[] { new[] { 0, 0, 1, 1, 2, 2 }, new[] { 0, 0, 0, 1, 1, 1 } },
        new object[] { new[] { 5, 5, 9, 9, 9, 3, 3 }, new[] { 1, 2, 1, 2, 1, 2, 1 } },        // sparse label ids
        new object[] { new[] { 0, 1, 2, 3, 4 }, new[] { 0, 0, 0, 0, 0 } },                    // singletons vs one cluster
        new object[] { new[] { 2, 2, 2, 7, 7, 1, 1, 1 }, new[] { 2, 2, 2, 7, 7, 1, 1, 1 } },  // identical
    };

    [Theory]
    [MemberData(nameof(LabelPairs))]
    public void RandIndex_MatchesBruteForcePairCount(int[] a, int[] b)
        => Assert.Equal(BruteRand(a, b), new RandIndex().Evaluate(a, b), 12);

    [Theory]
    [MemberData(nameof(LabelPairs))]
    public void RandIndex_IsSymmetric(int[] a, int[] b)
        => Assert.Equal(new RandIndex().Evaluate(a, b), new RandIndex().Evaluate(b, a), 12);

    [Fact]
    public void RandIndex_IgnoresLabelNames()
    {
        int[] a = { 0, 0, 1, 1, 2, 2, 2 };
        int[] relabeled = { 7, 7, 3, 3, 0, 0, 0 };
        Assert.Equal(1.0, new RandIndex().Evaluate(a, relabeled), 12);
    }

    [Theory]
    [MemberData(nameof(LabelPairs))]
    public void LinkOverlap_OnCompleteGraph_EqualsRandIndex(int[] a, int[] b)
        => Assert.Equal(
            new RandIndex().Evaluate(a, b),
            new LinkOverlap().Evaluate(Complete(a.Length), a, b), 12);

    [Fact]
    public void LinkOverlap_ReadsBoundaries_NotVolume()
    {
        // A cuts the path at (2,3), B at (1,2): they disagree on exactly those two
        // edges and agree on the other three, Q = 3/5. The all-pairs Rand index of
        // the same labels is 10/15 — it counts volume, not boundary.
        Edge[] edges = PathEdges();
        var graph = CsrGraph.FromEdges(edges, nodeCount: 6);
        int[] a = { 0, 0, 0, 1, 1, 1 };
        int[] b = { 0, 0, 1, 1, 1, 1 };

        Assert.Equal(0.6, new LinkOverlap().Evaluate(graph, a, b), 12);
        Assert.Equal(10.0 / 15.0, new RandIndex().Evaluate(a, b), 12);
    }

    [Fact]
    public void LinkOverlap_Weighted_CountsEdgesByCoupling()
    {
        // The two disagreeing edges are the weak ones: Q_w = 3 / 3.5.
        Edge[] edges = PathEdges();
        var graph = CsrGraph.FromEdges(edges, nodeCount: 6);
        int[] a = { 0, 0, 0, 1, 1, 1 };
        int[] b = { 0, 0, 1, 1, 1, 1 };

        Assert.Equal(3.0 / 3.5, new LinkOverlap(weighted: true).Evaluate(graph, a, b), 12);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void LinkOverlap_MatchesBruteForceEdgeCount(bool weighted)
    {
        // Irregular graph, sparse label ids.
        Edge[] edges =
        {
            new Edge(0, 1, 0.9), new Edge(0, 2, 0.4), new Edge(1, 2, 0.7),
            new Edge(2, 3, 0.2), new Edge(3, 4, 0.8), new Edge(3, 5, 0.6),
            new Edge(4, 5, 0.5), new Edge(5, 6, 0.3), new Edge(1, 6, 0.1),
        };
        var graph = CsrGraph.FromEdges(edges, nodeCount: 7);
        int[] a = { 4, 4, 4, 8, 8, 8, 2 };
        int[] b = { 1, 1, 3, 3, 3, 1, 1 };

        Assert.Equal(
            BruteLinkOverlap(edges, a, b, weighted),
            new LinkOverlap(weighted).Evaluate(graph, a, b), 12);
    }

    [Fact]
    public void LinkOverlap_RejectsLabelLengthMismatch()
    {
        var graph = CsrGraph.FromEdges(PathEdges(), nodeCount: 6);
        Assert.Throws<ArgumentException>(
            () => new LinkOverlap().Evaluate(graph, new[] { 0, 0, 0 }, new[] { 0, 0, 0, 1, 1, 1 }));
    }
}
