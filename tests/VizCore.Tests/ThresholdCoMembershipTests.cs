using System;
using Clustering.Graphical.SPC.Partitions.Strategies;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// A2 parity: ThresholdCoMembership — BWD1996 eq-4 discriminant
/// δ̄_ij = ((q−1)⟨n_ij⟩+1)/q threshold-and-CC partition strategy.
/// </summary>
public sealed class ThresholdCoMembershipTests
{
    // 2-node graph: one undirected edge (0→1 upper-tri slot=0, 1→0 lower-tri slot=1).
    private static CsrGraph TwoNodeGraph()
        => CsrGraph.FromEdges(new Edge[] { new(0, 1, 1.0) }, 2);

    private static Affinities DummyAffinities(int slotCount)
        => new() { Temperature = 1.0, G = new double[slotCount] };

    // ── Error paths ──────────────────────────────────────────────────────────

    [Fact]
    public void Apply_Throws_WhenCoMembershipIsNull()
    {
        CsrGraph graph = TwoNodeGraph();
        var strategy = new ThresholdCoMembership { Theta = 0.5 };

        Assert.Throws<InvalidOperationException>(() =>
            strategy.Apply(graph, DummyAffinities(graph.Targets.Length), alignments: null, coMembership: null));
    }

    // ── Connectivity ─────────────────────────────────────────────────────────

    [Fact]
    public void Apply_ConnectsNodes_WhenDeltaAboveTheta()
    {
        // Q=20, ⟨n_ij⟩=0.9 → δ̄ = (19·0.9+1)/20 = 18.1/20 = 0.905 > 0.5
        CsrGraph graph = TwoNodeGraph();
        var cm = new CoMembership
        {
            Temperature = 1.0,
            Q           = 20,
            G           = new[] { 0.9, 0.0 },  // slot 0 upper-tri, slot 1 lower-tri (mirror, ignored)
        };
        var result = new ThresholdCoMembership { Theta = 0.5 }
            .Apply(graph, DummyAffinities(graph.Targets.Length), alignments: null, coMembership: cm);

        Assert.Equal(result.Labels[0], result.Labels[1]); // same cluster
    }

    [Fact]
    public void Apply_DisconnectsNodes_WhenDeltaBelowTheta()
    {
        // Q=20, ⟨n_ij⟩=0.1 → δ̄ = (19·0.1+1)/20 = 2.9/20 = 0.145 < 0.5
        CsrGraph graph = TwoNodeGraph();
        var cm = new CoMembership
        {
            Temperature = 1.0,
            Q           = 20,
            G           = new[] { 0.1, 0.0 },
        };
        var result = new ThresholdCoMembership { Theta = 0.5 }
            .Apply(graph, DummyAffinities(graph.Targets.Length), alignments: null, coMembership: cm);

        Assert.NotEqual(result.Labels[0], result.Labels[1]); // different clusters
    }

    // ── Eq-4 formula accuracy ─────────────────────────────────────────────────

    [Theory]
    [InlineData(20,  9.0 / 19.0, 0.5)]             // exact boundary: (19·(9/19)+1)/20 = 10/20 = 0.5
    [InlineData(20,  1.0,        1.0)]              // fully co-clustered: (19·1+1)/20 = 1.0
    [InlineData(20,  0.0,        0.05)]             // pure noise: (19·0+1)/20 = 1/20 = 0.05
    [InlineData(4,   0.75,       0.8125)]           // Q=4, ⟨n⟩=0.75: (3·0.75+1)/4 = 3.25/4 = 0.8125
    [InlineData(2,   0.5,        0.75)]             // Q=2, ⟨n⟩=0.5: (1·0.5+1)/2 = 0.75
    public void Eq4Transform_IsCorrect(int q, double nij, double expectedDelta)
    {
        // Force exactly one edge above/below: use a threshold that gates on
        // whether the transform output matches the formula.
        // Connect iff delta > theta; set theta = expectedDelta ± epsilon to probe.
        CsrGraph graph = TwoNodeGraph();
        var cm = new CoMembership { Temperature = 1.0, Q = q, G = new[] { nij, 0.0 } };

        // Just above expected delta → should connect
        bool connectsAbove = Connects(graph, cm, theta: expectedDelta - 1e-10);
        // Exactly at expected delta → should NOT connect (strict >)
        bool connectsAt    = Connects(graph, cm, theta: expectedDelta);

        Assert.True(connectsAbove,  $"Q={q}, n={nij}: expected delta={expectedDelta} should be > (theta-eps)");
        Assert.False(connectsAt,    $"Q={q}, n={nij}: expected delta={expectedDelta} should not be strictly > delta");
    }

    private static bool Connects(CsrGraph graph, CoMembership cm, double theta)
    {
        var result = new ThresholdCoMembership { Theta = theta }
            .Apply(graph, DummyAffinities(graph.Targets.Length), alignments: null, coMembership: cm);
        return result.Labels[0] == result.Labels[1];
    }
}
