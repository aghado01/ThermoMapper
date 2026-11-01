using System.Linq;
using Graphs;
using Graphs.Primitives;
using Graphs.Proximity;
using Xunit;

namespace VizCore.Tests;

/// <summary>Local Mutual Proximity's two rank rules on a path 0-1-2 (graph-compiler-brief, decision 1b).</summary>
public sealed class LmpRankRuleTests
{
    private static CsrGraph Path() =>
        CsrGraph.FromEdges(new[] { new Edge(0, 1, 0.5), new Edge(1, 2, 0.3) }, 3);

    private static double WeightOf(CsrGraph graph, int from, int to)
    {
        for (int slot = graph.RowPointers[from]; slot < graph.RowPointers[from + 1]; slot++)
            if (graph.Targets[slot] == to) return graph.Weights[slot];
        throw new System.InvalidOperationException($"No edge {from}-{to}.");
    }

    [Fact]
    public void StrictRank_WritesZeroOnEveryEdgeOfAPath_BecauseEachHasALeafEnd()
    {
        CsrGraph scaled = LocalMutualProximity.ApplyLocalScaling(Path(), weightsAreCouplings: true);

        Assert.All(scaled.Weights, w => Assert.Equal(0.0, w));
    }

    [Fact]
    public void InclusiveRank_IsTheEmpiricalCdf_SoNoEdgeIsZero()
    {
        CsrGraph scaled = LocalMutualProximity.ApplyLocalScaling(
            Path(), weightsAreCouplings: true, rankRule: LmpRankRule.Inclusive);

        // Edge 0-1 (0.5): rank 1 at the leaf 0, 2/2 at node 1. Edge 1-2 (0.3): 1/2 at node 1, rank 1 at the leaf 2.
        Assert.Equal(1.0, WeightOf(scaled, 0, 1));
        Assert.Equal(0.5, WeightOf(scaled, 1, 2));
        Assert.Equal(WeightOf(scaled, 0, 1), WeightOf(scaled, 1, 0));
    }

    [Fact]
    public void InclusiveRank_OnDistances_CountsTheEdgesAtLeastAsFar()
    {
        CsrGraph path = CsrGraph.FromEdges(new[] { new Edge(0, 1, 1.0), new Edge(1, 2, 3.0) }, 3);

        CsrGraph scaled = LocalMutualProximity.ApplyLocalScaling(
            path, weightsAreCouplings: false, rankRule: LmpRankRule.Inclusive);

        Assert.Equal(1.0, WeightOf(scaled, 0, 1));    // d = 1: both of node 1's edges are at least as far
        Assert.Equal(0.5, WeightOf(scaled, 1, 2));    // d = 3: only itself at node 1
    }

    [Fact]
    public void ProtectedEdges_KeepTheirWeightUnderEitherRule()
    {
        var protectedPairs = new System.Collections.Generic.HashSet<(int Lo, int Hi)> { (0, 1) };

        foreach (LmpRankRule rule in new[] { LmpRankRule.Strict, LmpRankRule.Inclusive })
        {
            CsrGraph scaled = LocalMutualProximity.ApplyLocalScaling(
                Path(), weightsAreCouplings: true, protectedEdges: protectedPairs, rankRule: rule);
            Assert.Equal(0.5, WeightOf(scaled, 0, 1));
        }
    }
}
