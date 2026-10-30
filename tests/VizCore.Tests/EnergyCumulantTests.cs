using System;
using System.Linq;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Execution;
using Graphs.Models.Potts.Observables;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The energy moments ⟨E⟩–⟨E⁴⟩ over the per-draw energy record, and the energy cumulant V = 1 − ⟨E⁴⟩/(3⟨E²⟩²) built
/// from them and assembled per temperature from the replicas' pooled moments.
/// </summary>
public sealed class EnergyCumulantTests
{
    [Fact]
    public void Moments_ByHand()
    {
        EnergyMoments m = EnergyCumulant.Moments([-1.0, -2.0, -3.0]);

        Assert.Equal(-2.0, m.Mean, 14);
        Assert.Equal(14.0 / 3.0, m.SecondMoment, 14);
        Assert.Equal(-12.0, m.ThirdMoment, 14);
        Assert.Equal(98.0 / 3.0, m.FourthMoment, 14);
    }

    [Fact]
    public void Reduce_IsTwoThirdsForASharpEnergy_AndBelowForTwoPeaks()
    {
        Assert.Equal(2.0 / 3.0, EnergyCumulant.Reduce([-5.0, -5.0, -5.0, -5.0]), 14);

        // Two equal peaks at −1 and −3: ⟨E²⟩ = 5, ⟨E⁴⟩ = 41, so V = 1 − 41/75.
        Assert.Equal(1.0 - 41.0 / 75.0, EnergyCumulant.Reduce([-1.0, -3.0, -1.0, -3.0]), 14);
    }

    [Fact]
    public void Reduce_IsNaN_WithoutAPositiveSecondMoment()
    {
        Assert.True(double.IsNaN(EnergyCumulant.Reduce(ReadOnlySpan<double>.Empty)));
        Assert.True(double.IsNaN(EnergyCumulant.Reduce([0.0, 0.0])));
    }

    [Fact]
    public void Moments_MatchTheRunningSums_BitForBit()
    {
        Accumulator acc = SwRunner.Run(new SwRunSpec
        {
            Graph = BinderCumulantTests.PeriodicLattice(4), Temperature = 1.2, Q = 2, Seed = 7,
            Budget = new RunBudget(100, 500),
        }).Accumulator;

        EnergyMoments m = EnergyCumulant.Moments(acc.EnergyPerDraw);

        Assert.Equal(acc.RunningSumEnergy / acc.DrawCount, m.Mean);
        Assert.Equal(acc.RunningSumEnergySq / acc.DrawCount, m.SecondMoment);
    }

    [Fact]
    public void SweepProfile_ReducesThePooledMoments_Once()
    {
        double[][] series = [[-1.0, -3.0, -2.0, -2.0], [-4.0, -4.0, -3.0, -5.0]];
        SpcRunResult[] runs = [Run(1.0, series[0]), Run(1.0, series[1]), Run(2.0, [-1.0, -1.0, -1.0, -1.0])];

        SweepProfile profile = SweepProfile.From(runs);
        var v = profile.AdditionalChannels["EnergyCumulant"];

        double secondMoment = series.Average(s => s.Average(e => e * e));
        double fourthMoment = series.Average(s => s.Average(e => e * e * e * e));
        Assert.Equal(EnergyCumulant.Reduce(secondMoment, fourthMoment), v[0], 12);
        Assert.Equal(2.0 / 3.0, v[1], 14);
    }

    [Fact]
    public void SweepProfile_IsNaN_WhereAFrameHasNoRecord()
    {
        SpcRunResult bare = Run(1.0, [-1.0, -2.0]);
        bare = bare with { Accumulator = bare.Accumulator with { EnergyPerDraw = [] } };

        SweepProfile profile = SweepProfile.From([bare]);

        Assert.True(double.IsNaN(profile.AdditionalChannels["EnergyCumulant"][0]));
    }

    internal static SpcRunResult Run(double temperature, double[] energyPerDraw) => new()
    {
        Graph = CsrGraph.FromEdges([new Edge(0, 1, 1.0)], 2),
        Accumulator = new Accumulator
        {
            Temperature = temperature, Q = 2, DrawCount = energyPerDraw.Length,
            Spins = new int[2], ClusterSizeHistogram = new long[] { 0, (long)energyPerDraw.Length },
            RngState0 = 1, RngState1 = 2, RngState2 = 3, RngState3 = 4,
            RunningSumSqClusterSizes = 4.0 * energyPerDraw.Length, RunningSumSqClusterSizesExcl = 0,
            RunningSumEnergy = energyPerDraw.Sum(), RunningSumEnergySq = energyPerDraw.Sum(e => e * e),
            RunningSumMag = energyPerDraw.Length, RunningSumMagSq = energyPerDraw.Length,
            EnergyPerDraw = energyPerDraw,
            MagPerDraw = Enumerable.Repeat(1.0, energyPerDraw.Length).ToArray(),
            SumSqClusterSizesPerDraw = Enumerable.Repeat(4L, energyPerDraw.Length).ToArray(),
            MaxClusterSizePerDraw = Enumerable.Repeat(2, energyPerDraw.Length).ToArray(),
        },
    };
}
