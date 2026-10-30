using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Xunit;
using PottsObs = Graphs.Models.Potts.Observables;

namespace VizCore.Tests;

/// <summary>
/// The replica spread (S6): <see cref="SweepCurves.ByTemperatureWithSpread"/> keeps, beside each temperature's replica
/// mean, the replicas' sample standard deviation (NaN for one replica) and their count. Its means are
/// <see cref="SweepCurves.ByTemperature"/>'s bit for bit, which <see cref="SweepProfile"/> still assembles from.
/// </summary>
public sealed class SweepCurvesSpreadTests
{
    [Fact]
    public void Means_AreByTemperaturesBitForBit()
    {
        // Three replicas at each of the triangle's probe temperatures, seeded per (T, replica).
        var frames = (from t in new[] { 0.3, 1.0, 3.0 }
                      from r in Enumerable.Range(0, 3)
                      select SwRunner.Run(new SwRunSpec
                      {
                          Graph = SwAccumulatorProbeTests.Triangle(), Temperature = t, Q = 4, ReplicaIndex = r,
                          Seed = SpcSeedHelper.Derive(SwAccumulatorProbeTests.Seed, t, r),
                          Budget = SwAccumulatorProbeTests.Budget,
                      }).Accumulator).ToList();

        foreach (Func<Accumulator, double> perFrame in new Func<Accumulator, double>[]
                 {
                     f => PottsObs.Susceptibility.Reduce(f.RunningSumSqClusterSizes, f.DrawCount, f.Spins.Length),
                     f => PottsObs.MeanEnergy.Reduce(f.RunningSumEnergy, f.DrawCount),
                     f => PottsObs.Magnetization.Reduce(f.RunningSumMag, f.RunningSumMagSq, f.DrawCount).SecondMoment,
                 })
        {
            var (temperatures, values) = SweepCurves.ByTemperature(frames, perFrame);
            SweepCurve curve = SweepCurves.ByTemperatureWithSpread(frames, perFrame);

            Assert.Equal(temperatures, curve.Temperatures);
            Assert.Equal(values.Select(BitConverter.DoubleToInt64Bits), curve.Means.Select(BitConverter.DoubleToInt64Bits));
            Assert.Equal(new[] { 3, 3, 3 }, curve.Replicas);
            Assert.All(curve.StandardDeviations, sd => Assert.True(sd >= 0.0));
        }
    }

    [Fact]
    public void Spread_IsTheSampleStandardDeviation()
    {
        SweepCurve curve = SweepCurves.ByTemperatureWithSpread(
            [Frame(1.0, 1.0), Frame(2.0, 5.0), Frame(1.0, 2.0), Frame(1.0, 4.0)], f => f.RunningSumEnergy);

        Assert.Equal(new[] { 1.0, 2.0 }, curve.Temperatures);
        Assert.Equal(new[] { 3, 1 }, curve.Replicas);
        Assert.InRange(Math.Abs(curve.Means[0] - 7.0 / 3.0), 0.0, 1e-15);
        // Σ(v − 7/3)² = 42/9 over n − 1 = 2.
        Assert.InRange(Math.Abs(curve.StandardDeviations[0] - Math.Sqrt(7.0 / 3.0)), 0.0, 1e-15);
        Assert.Equal(5.0, curve.Means[1]);
        Assert.True(double.IsNaN(curve.StandardDeviations[1]));
    }

    [Fact]
    public void NoFrames_GiveAnEmptyCurve()
    {
        SweepCurve curve = SweepCurves.ByTemperatureWithSpread(new List<Accumulator>(), f => f.RunningSumEnergy);

        Assert.Empty(curve.Temperatures);
        Assert.Empty(curve.Means);
        Assert.Empty(curve.StandardDeviations);
        Assert.Empty(curve.Replicas);
    }

    private static Accumulator Frame(double temperature, double value) => new()
    {
        Temperature = temperature, Q = 2, DrawCount = 1,
        Spins = new int[2], ClusterSizeHistogram = new long[2],
        RngState0 = 1, RngState1 = 2, RngState2 = 3, RngState3 = 4,
        RunningSumSqClusterSizes = 0, RunningSumSqClusterSizesExcl = 0,
        RunningSumEnergy = value, RunningSumEnergySq = 0,
        RunningSumMag = 0, RunningSumMagSq = 0,
    };
}
