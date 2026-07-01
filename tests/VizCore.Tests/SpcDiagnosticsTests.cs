using System;
using System.Collections.Generic;
using Clustering.Graphical.SPC;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs.Primitives;
using Graphs.Observables;
using Maths.Information;
using Xunit;

namespace VizCore.Tests;

public sealed class SpcDiagnosticsTests
{
    [Fact]
    public void Shannon_EntropyNats_Counts_DegenerateAndUniformCasesBehaveAsExpected()
    {
        Assert.Equal(0.0, Shannon.EntropyNats(new int[] { 0, 0, 0 }));
        Assert.Equal(0.0, Shannon.EntropyNats(new int[] { 5, 0, 0 }));
        Assert.InRange(Math.Abs(Shannon.EntropyNats(new int[] { 1, 1, 1, 1 }) - Math.Log(4.0)), 0.0, 1e-12);
    }

    [Fact]
    public void Shannon_EntropyNats_Counts_MatchesHandComputedSkewedDistribution()
    {
        double expected = -(0.99 * Math.Log(0.99) + 0.01 * Math.Log(0.01));

        Assert.InRange(Math.Abs(Shannon.EntropyNats(new int[] { 99, 1 }) - expected), 0.0, 1e-12);
    }

    [Fact]
    public void AffinityEntropy_EntropyNats_ComputesOverNormalizedAffinityField()
    {
        var edges = CreateFrame(temperature: 1.0, cycleCount: 4) with
        {
            BondFormedCount = new[] { 2, 0, 2, 0 },
            SpinAgreementCount = new[] { 0, 0, 0, 0 },
        };

        // Mint the currency (G = BondFormedCount/draws), then reduce at the graph tier.
        double entropy = AffinityEntropy.EntropyNats(SwCurrencies.ToAffinities(edges));

        Assert.InRange(Math.Abs(entropy - Math.Log(2.0)), 0.0, 1e-12);
    }

    [Fact]
    public void ByTemperature_Susceptibility_GroupsReplicasAndSortsTemperatures()
    {
        var frames = new List<Accumulator>
        {
            CreateFrame(temperature: 2.0, cycleCount: 5, nodeCount: 3, sumSqClusterSizes: 30.0),
            CreateFrame(temperature: 1.0, cycleCount: 4, nodeCount: 3, sumSqClusterSizes: 12.0),
            CreateFrame(temperature: 1.0, cycleCount: 4, nodeCount: 3, sumSqClusterSizes: 24.0),
        };

        var (temperatures, chi) = SweepCurves.ByTemperature(
            frames, f => Graphs.Models.Potts.Observables.Susceptibility.Reduce(
                f.RunningSumSqClusterSizes, f.DrawCount, f.Spins.Length));

        Assert.Equal(new[] { 1.0, 2.0 }, temperatures);
        Assert.InRange(Math.Abs(chi[0] - 1.5), 0.0, 1e-12);
        Assert.InRange(Math.Abs(chi[1] - 2.0), 0.0, 1e-12);
    }

    [Fact]
    public void LabelEntropy_From_UsesClusterSizeHistogramPerTemperature()
    {
        var frames = new List<Accumulator>
        {
            CreateFrame(temperature: 2.0, clusterSizeHistogram: new[] { 1, 0, 0, 0 }),
            CreateFrame(temperature: 1.0, clusterSizeHistogram: new[] { 2, 0, 0, 0 }),
            CreateFrame(temperature: 1.0, clusterSizeHistogram: new[] { 0, 2, 0, 0 }),
        };

        LabelEntropyCurve curve = LabelEntropyCurve.From(frames);

        Assert.Equal(new[] { 1.0, 2.0 }, curve.Temperatures);
        Assert.InRange(Math.Abs(curve.Entropy[0] - Math.Log(2.0)), 0.0, 1e-12);
        Assert.Equal(0.0, curve.Entropy[1]);
    }

    [Fact]
    public void ByTemperature_SpecificHeat_ComputesBetaSquaredEnergyVariance()
    {
        var frames = new List<Accumulator>
        {
            CreateFrame(temperature: 2.0, cycleCount: 3, runningSumEnergy: 6.0, runningSumEnergySq: 15.0),
            CreateFrame(temperature: 2.0, cycleCount: 3, runningSumEnergy: 6.0, runningSumEnergySq: 15.0),
        };

        var (temperatures, cv) = SweepCurves.ByTemperature(
            frames, f => Graphs.Models.Potts.Observables.SpecificHeat.Reduce(
                f.RunningSumEnergy, f.RunningSumEnergySq, f.DrawCount, f.Temperature));

        Assert.Equal(new[] { 2.0 }, temperatures);
        Assert.InRange(Math.Abs(cv[0] - 0.25), 0.0, 1e-12);
    }

    [Fact]
    public void SweepProfile_From_PopulatesAdditionalChannels()
    {
        var runs = new List<SpcRunResult>
        {
            new SpcRunResult
            {
                Graph = BuildGraph(2, (0, 1, 1.0)),
                Accumulator = new Accumulator
                {
                    Temperature = 1.0,
                    Q = 3,
                    DrawCount = 4,
                    Spins = new[] { 0, 0 },
                    ClusterSizeHistogram = new[] { 2 },
                    RngState0 = 1UL,
                    RngState1 = 2UL,
                    RngState2 = 3UL,
                    RngState3 = 4UL,
                    RunningSumSqClusterSizes = 8.0,
                    RunningSumSqClusterSizesExcl = 0.0,
                    RunningSumEnergy = -4.0,
                    RunningSumEnergySq = 4.0,
                    RunningSumMag = 2.0,
                    RunningSumMagSq = 2.0,
                }
            },
            new SpcRunResult
            {
                Graph = BuildGraph(2, (0, 1, 1.0)),
                Accumulator = new Accumulator
                {
                    Temperature = 2.0,
                    Q = 3,
                    DrawCount = 2,
                    Spins = new[] { 0, 1 },
                    ClusterSizeHistogram = new[] { 1, 1 },
                    RngState0 = 5UL,
                    RngState1 = 6UL,
                    RngState2 = 7UL,
                    RngState3 = 8UL,
                    RunningSumSqClusterSizes = 4.0,
                    RunningSumSqClusterSizesExcl = 0.0,
                    RunningSumEnergy = 0.0,
                    RunningSumEnergySq = 0.0,
                    RunningSumMag = 0.0,
                    RunningSumMagSq = 0.0,
                }
            }
        };

        SweepProfile profile = SweepProfile.From(runs);

        Assert.Equal(2, profile.Count);
        Assert.True(profile.AdditionalChannels.ContainsKey("MeanEnergy"));
        Assert.Equal(new[] { -1.0, 0.0 }, profile.AdditionalChannels["MeanEnergy"]);
        Assert.Equal(new[] { 0.5, 0.0 }, profile.AdditionalChannels["MeanMagnetization"]);
        Assert.Equal(new[] { 0.25, 0.0 }, profile.AdditionalChannels["MagnetizationVariance"]);
    }

    [Fact]
    public void SpcProfileAnalysis_DirectMethods_UseSweepProfileValues()
    {
        var profile = new SweepProfile(
            Temperatures: new[] { 1.0, 2.0, 3.0 },
            Susceptibility: new[] { 0.5, 2.0, 1.0 },
            SpecificHeat: new[] { 0.0, 0.0, 0.0 },
            LabelEntropy: new[] { 0.0, 0.0, 0.0 },
            BondEntropy: null,
            AdditionalChannels: new Dictionary<string, IReadOnlyList<double>>());

        Assert.Equal(2.0, SpcProfileAnalysis.PickPeakTemperature(profile));
        Assert.Equal((2.0, 3.0), SpcProfileAnalysis.ComputeHalfMaximumBand(profile));
        Assert.InRange(SpcProfileAnalysis.ComputeStability(profile), 0.0, 1.0);
    }

    [Fact]
    public void CriticalTemperatureEstimator_Estimate_MatchesUniformWeightSanityCheck()
    {
        CsrGraph graph = BuildGraph(3,
            (0, 1, 1.0),
            (1, 2, 1.0),
            (0, 2, 1.0));

        double estimate = SpcScheduleHelpers.Estimate(graph, q: 4);

        Assert.InRange(Math.Abs(estimate - (2.0 / Math.Log(4.0))), 0.0, 1e-12);
    }

    [Fact]
    public void SwCurrencies_Mint_FeedsRelocatedGraphSignals()
    {
        // 2-node graph, one undirected edge (0,1). CSR slot 0 = the upper-triangle (0→1) slot.
        CsrGraph graph = BuildGraph(2, (0, 1, 1.0));

        Accumulator accumulator = CreateFrame(temperature: 1.5, cycleCount: 4, nodeCount: 2) with
        {
            BondFormedCount    = new[] { 2, 0 },   // bonded on 2 of 4 draws  → affinity 0.5
            SpinAgreementCount = new[] { 4, 0 },   // agreed on every draw    → alignment 1.0
        };

        // The SW realization of the currencies (PKWang mints Affinities directly; no Alignments).
        Affinities affinities = SwCurrencies.ToAffinities(accumulator);
        Alignments alignments = SwCurrencies.ToAlignments(accumulator);

        Assert.Equal(1.5, affinities.Temperature);
        Assert.Equal(new[] { 0.5, 0.0 }, affinities.G);
        Assert.Equal(new[] { 1.0, 0.0 }, alignments.G);

        // The relocated graph-tier signals now have a producer — mint → signal end-to-end.
        double[] degree = new AffinityDegree().Compute(affinities, graph);
        Assert.Equal(new[] { 0.5, 0.5 }, degree);  // 0.5 into both endpoints

        double[] entropyBits = new AffinityBinaryEntropySum().Compute(affinities, graph);
        Assert.InRange(Math.Abs(entropyBits[0] - 1.0), 0.0, 1e-12);  // H₂(0.5) = 1 bit
        Assert.InRange(Math.Abs(entropyBits[1] - 1.0), 0.0, 1e-12);

        double[] centrality = new AlignmentEigenCentrality().Compute(alignments, graph);
        Assert.Equal(2, centrality.Length);
        Assert.True(centrality[0] > 0.0 && centrality[1] > 0.0);
        Assert.InRange(Math.Abs(centrality[0] - centrality[1]), 0.0, 1e-9);  // symmetric → equal entries
    }

    [Fact]
    public void MagnetizationSusceptibility_Reduce_ScalesPooledVarianceByNOverT()
    {
        // ⟨m⟩=0.5, ⟨m²⟩=0.5 → var=0.25; (N/T)·var = (2/1)·0.25 = 0.5
        Assert.Equal(0.5, Graphs.Models.Potts.Observables.MagnetizationSusceptibility.Reduce(0.5, 0.5, 2, 1.0));
        Assert.Equal(0.0, Graphs.Models.Potts.Observables.MagnetizationSusceptibility.Reduce(0.5, 0.5, 0, 1.0));  // N=0
        Assert.Equal(0.0, Graphs.Models.Potts.Observables.MagnetizationSusceptibility.Reduce(0.5, 0.5, 2, 0.0));  // T=0
    }

    [Fact]
    public void SweepProfile_From_EmitsThreeSusceptibilityChannels_KindSelectsPrimary()
    {
        var runs = TwoTemperatureRuns();

        // FK cluster χ = Σ|c|²/(draws·N): T=1 → 8/(4·2)=1.0; T=2 → 4/(2·2)=1.0.
        // FK reduced (giant-excluded sum = 0): 0, 0.
        // Magnetization (N/T)·var: T=1 → (2/1)·0.25=0.5; T=2 → 0.
        SweepProfile fk = SweepProfile.From(runs);  // default kind = FkCluster
        Assert.Equal(new[] { 1.0, 1.0 }, fk.AdditionalChannels["SusceptibilityFkCluster"]);
        Assert.Equal(new[] { 0.0, 0.0 }, fk.AdditionalChannels["SusceptibilityFkReduced"]);
        Assert.Equal(new[] { 0.5, 0.0 }, fk.AdditionalChannels["SusceptibilityMagnetization"]);
        Assert.Equal(fk.AdditionalChannels["SusceptibilityFkCluster"], fk.Susceptibility);

        SweepProfile mag = SweepProfile.From(runs, SusceptibilityKind.Magnetization);
        Assert.Equal(mag.AdditionalChannels["SusceptibilityMagnetization"], mag.Susceptibility);

        SweepProfile reduced = SweepProfile.From(runs, SusceptibilityKind.FkReduced);
        Assert.Equal(reduced.AdditionalChannels["SusceptibilityFkReduced"], reduced.Susceptibility);
    }

    [Fact]
    public void SweepProfile_From_SpecificHeat_PoolsEnergyMomentsAcrossReplicas()
    {
        // Two replicas at T=1, each with ZERO within-replica energy variance but different
        // means (⟨E⟩=1 and ⟨E⟩=3). Reduce-then-average → 0; pool-then-reduce → Var=1, Cv=1.
        // Asserts the pooled value, proving Cv is not averaged per-frame (the #3 follow-up fix).
        var graph = BuildGraph(4, (0, 1, 1.0));
        var runs = new List<SpcRunResult>
        {
            new() { Graph = graph, Accumulator = CreateFrame(temperature: 1.0, cycleCount: 2, runningSumEnergy: 2.0, runningSumEnergySq: 2.0) },
            new() { Graph = graph, Accumulator = CreateFrame(temperature: 1.0, cycleCount: 2, runningSumEnergy: 6.0, runningSumEnergySq: 18.0) },
        };

        SweepProfile profile = SweepProfile.From(runs);

        Assert.Single(profile.SpecificHeat);
        Assert.InRange(Math.Abs(profile.SpecificHeat[0] - 1.0), 0.0, 1e-12);
    }

    [Fact]
    public void SweepProfile_From_BondEntropy_PoolsRatesAcrossReplicas()
    {
        // The channel is DISPERSION entropy (AffinityEntropy = Shannon over
        // the NORMALIZED field; sibling of, and distinct from, the per-edge
        // binary entropy). Replica A concentrates all bond activity on edge
        // 0, replica B on edge 1: per-replica dispersion is 0 for both, so
        // reduce-then-average gives 0 — the POOLED rates [0.5, 0.5] give
        // ln 2. Asserts pool-then-reduce (the Jensen discipline).
        var graph = BuildGraph(4, (0, 1, 1.0), (2, 3, 1.0));
        var runs = new List<SpcRunResult>
        {
            new() { Graph = graph, Accumulator = CreateFrame(temperature: 1.0, cycleCount: 10) with { BondFormedCount = new[] { 10, 0 } } },
            new() { Graph = graph, Accumulator = CreateFrame(temperature: 1.0, cycleCount: 10) with { BondFormedCount = new[] { 0, 10 } } },
        };

        SweepProfile profile = SweepProfile.From(runs);

        Assert.NotNull(profile.BondEntropy);
        Assert.Single(profile.BondEntropy!);

        Assert.InRange(Math.Abs(profile.BondEntropy![0] - Math.Log(2.0)), 0.0, 1e-12);
    }

    private static List<SpcRunResult> TwoTemperatureRuns() => new()
    {
        new SpcRunResult
        {
            Graph = BuildGraph(2, (0, 1, 1.0)),
            Accumulator = new Accumulator
            {
                Temperature = 1.0, Q = 3, DrawCount = 4,
                Spins = new[] { 0, 0 }, ClusterSizeHistogram = new[] { 2 },
                RngState0 = 1UL, RngState1 = 2UL, RngState2 = 3UL, RngState3 = 4UL,
                RunningSumSqClusterSizes = 8.0, RunningSumSqClusterSizesExcl = 0.0,
                RunningSumEnergy = -4.0, RunningSumEnergySq = 4.0,
                RunningSumMag = 2.0, RunningSumMagSq = 2.0,
            }
        },
        new SpcRunResult
        {
            Graph = BuildGraph(2, (0, 1, 1.0)),
            Accumulator = new Accumulator
            {
                Temperature = 2.0, Q = 3, DrawCount = 2,
                Spins = new[] { 0, 1 }, ClusterSizeHistogram = new[] { 1, 1 },
                RngState0 = 5UL, RngState1 = 6UL, RngState2 = 7UL, RngState3 = 8UL,
                RunningSumSqClusterSizes = 4.0, RunningSumSqClusterSizesExcl = 0.0,
                RunningSumEnergy = 0.0, RunningSumEnergySq = 0.0,
                RunningSumMag = 0.0, RunningSumMagSq = 0.0,
            }
        },
    };

    private static Accumulator CreateFrame(
        double temperature,
        int cycleCount = 4,
        int nodeCount = 4,
        double sumSqClusterSizes = 0.0,
        int[]? clusterSizeHistogram = null,
        double runningSumEnergy = 0.0,
        double runningSumEnergySq = 0.0)
    {
        return new Accumulator
        {
            Temperature = temperature,
            Q = 4,
            DrawCount = cycleCount,
            Spins = new int[nodeCount],
            ClusterSizeHistogram = clusterSizeHistogram ?? new int[nodeCount],
            RngState0 = 1,
            RngState1 = 2,
            RngState2 = 3,
            RngState3 = 4,
            RunningSumSqClusterSizes = sumSqClusterSizes,
            RunningSumSqClusterSizesExcl = 0.0,
            RunningSumEnergy = runningSumEnergy,
            RunningSumEnergySq = runningSumEnergySq,
            RunningSumMag = 0.0,
            RunningSumMagSq = 0.0,
        };
    }

    private static CsrGraph BuildGraph(int nodeCount, params (int Source, int Target, double Weight)[] edges)
    {
        var graphEdges = new Edge[edges.Length];
        for (int i = 0; i < edges.Length; i++)
            graphEdges[i] = new Edge(edges[i].Source, edges[i].Target, edges[i].Weight);

        return CsrGraph.FromEdges(graphEdges, nodeCount);
    }
}
