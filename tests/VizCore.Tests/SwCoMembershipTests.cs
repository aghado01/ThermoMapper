using System;
using System.IO;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// A1 parity: co-membership accumulator (⟨n_ij⟩) — the improved estimator that counts
/// draws where two nodes fall in the same bond cluster, regardless of whether the direct
/// bond froze. Key invariant: co-membership ≥ bond-frequency for every edge because
/// transitive co-clustering via multi-hop paths is captured.
/// </summary>
public sealed class SwCoMembershipTests
{
    // Triangle: 0─1─2─0, weight=1.0. All edges have equal affinity.
    // At T=0.3, bond probability = 1 - exp(-1/0.3) ≈ 0.964, so nearly always all bonds form
    // and co-membership ≈ bond-frequency ≈ 1.0. Still useful for gating/shape tests.
    // At T=1.5, bond probability = 1 - exp(-1/1.5) ≈ 0.487 — nodes often co-cluster via paths.
    private static CsrGraph Triangle()
    {
        var edges = new Edge[]
        {
            new(0, 1, 1.0),
            new(1, 2, 1.0),
            new(0, 2, 1.0),
        };
        return CsrGraph.FromEdges(edges, 3);
    }

    private static SwRunSpec MakeSpec(CsrGraph graph, AccumulationSpec accumulation, int draws = 2000)
        => new()
        {
            Graph        = graph,
            Temperature  = 1.5,
            Q            = 4,
            Accumulation = accumulation,
            Seed         = 42,
            Budget       = new RunBudget(BurnIn: 200, Cycles: draws),
        };

    // ── Gating ───────────────────────────────────────────────────────────────

    [Fact]
    public void CoMembership_IsNull_WhenNotRequested()
    {
        Accumulator acc = SwRunner.Run(MakeSpec(Triangle(),
            new AccumulationSpec { Affinities = true })).Accumulator;

        Assert.Null(acc.CoMembershipCount);
    }

    [Fact]
    public void CoMembership_IsNonNull_WhenRequested()
    {
        Accumulator acc = SwRunner.Run(MakeSpec(Triangle(),
            new AccumulationSpec { Affinities = true, CoMembership = true })).Accumulator;

        Assert.NotNull(acc.CoMembershipCount);
        // CSR for the 3-node triangle has 6 directed slots (2 per undirected edge).
        Assert.Equal(acc.BondFormedCount!.Length, acc.CoMembershipCount!.Length);
    }

    // ── Core invariant: co-membership ≥ bond-frequency for every upper-triangle slot ──

    [Fact]
    public void CoMembership_GeqBondFrequency_OnEveryEdge()
    {
        // Proven invariant: if the direct bond froze, Find(i)==Find(j) is guaranteed,
        // so every BondFormedCount[e] increment also increments CoMembershipCount[e].
        Accumulator acc = SwRunner.Run(MakeSpec(Triangle(),
            new AccumulationSpec { Affinities = true, CoMembership = true })).Accumulator;

        CsrGraph graph = Triangle();
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
        {
            int e = edge.Slot;
            Assert.True(acc.CoMembershipCount![e] >= acc.BondFormedCount![e],
                $"Slot {e} (nodes {edge.Source},{edge.Target}): " +
                $"co-membership {acc.CoMembershipCount[e]} < bond-frequency {acc.BondFormedCount[e]}");
        }
    }

    // ── Mint ─────────────────────────────────────────────────────────────────

    [Fact]
    public void ToCoMembership_MintsCorrectly_FromHandCraftedAccumulator()
    {
        // 2-node graph, one undirected edge (0→1 upper-triangle slot = index 0).
        // 3 out of 4 draws they co-clustered → rate = 0.75.
        var accumulator = new Accumulator
        {
            Temperature              = 1.5,
            Q                        = 4,
            DrawCount                = 4,
            Spins                    = new int[2],
            ClusterSizeHistogram     = new int[2],
            RngState0 = 1, RngState1 = 2, RngState2 = 3, RngState3 = 4,
            RunningSumSqClusterSizes = 0, RunningSumSqClusterSizesExcl = 0,
            RunningSumEnergy = 0, RunningSumEnergySq = 0,
            RunningSumMag = 0, RunningSumMagSq = 0,
            CoMembershipCount = new[] { 3, 0 },   // slot 0 = (0→1) upper-tri: 3/4 → 0.75; slot 1 = lower-tri: 0
        };

        CoMembership cm = SwCurrencies.ToCoMembership(accumulator);

        Assert.Equal(1.5,  cm.Temperature);
        Assert.Equal(0.75, cm.G[0], precision: 12);
        Assert.Equal(0.0,  cm.G[1]);
    }

    [Fact]
    public void ToCoMembership_Throws_WhenCountArrayIsNull()
    {
        var accumulator = new Accumulator
        {
            Temperature = 1.0, Q = 4, DrawCount = 10,
            Spins = new int[2], ClusterSizeHistogram = new int[2],
            RngState0 = 1, RngState1 = 2, RngState2 = 3, RngState3 = 4,
            RunningSumSqClusterSizes = 0, RunningSumSqClusterSizesExcl = 0,
            RunningSumEnergy = 0, RunningSumEnergySq = 0,
            RunningSumMag = 0, RunningSumMagSq = 0,
            CoMembershipCount = null,
        };

        Assert.Throws<InvalidOperationException>(() => SwCurrencies.ToCoMembership(accumulator));
    }

    // ── Serialization round-trip ─────────────────────────────────────────────

    [Fact]
    public void CoMembership_SurvivesSerializationRoundTrip()
    {
        Accumulator acc = SwRunner.Run(MakeSpec(Triangle(),
            new AccumulationSpec { Affinities = true, CoMembership = true })).Accumulator;

        using var ms = new MemoryStream();
        AccumulatorSerializer.Instance.WriteTo(acc, ms);
        ms.Position = 0;
        Accumulator round = AccumulatorSerializer.Instance.ReadFrom(ms);

        Assert.Equal(acc.CoMembershipCount, round.CoMembershipCount);
    }

    [Fact]
    public void Serialization_PreservesNullCoMembership_WhenUntracked()
    {
        Accumulator acc = SwRunner.Run(MakeSpec(Triangle(),
            new AccumulationSpec { Affinities = true, CoMembership = false })).Accumulator;

        using var ms = new MemoryStream();
        AccumulatorSerializer.Instance.WriteTo(acc, ms);
        ms.Position = 0;
        Accumulator round = AccumulatorSerializer.Instance.ReadFrom(ms);

        Assert.Null(round.CoMembershipCount);
    }
}
