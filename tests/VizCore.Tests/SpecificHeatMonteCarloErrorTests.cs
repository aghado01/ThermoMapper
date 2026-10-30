using System;
using System.Linq;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Execution;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace VizCore.Tests;

/// <summary>
/// C_v's Monte Carlo error, sd(y)/√ESS(y) with y = (E − Ē)²/T², against references from outside the estimator: Gaussian
/// theory for an iid and an AR(1) energy series, and the spread of C_v across independent SW replicas. The tolerances
/// were fixed before the first run.
/// </summary>
public sealed class SpecificHeatMonteCarloErrorTests(ITestOutputHelper output)
{
    [Fact]
    public void IidGaussianEnergy_MatchesTheSampleVariancesError()
    {
        // E ~ N(−50, 2²) iid: Var((E − μ)²) = 2σ⁴, so C_v's error at T = 1 is √(2σ⁴/n). Tolerance 10%.
        const int n = 40_000;
        const double sigma = 2.0;
        var rng = new Xoshiro256PlusPlus(11);
        double[] e = Enumerable.Range(0, n).Select(_ => -50.0 + sigma * Gaussian(rng)).ToArray();

        double error = SpecificHeatMonteCarloError.Of([e], 1.0);

        double expected = Math.Sqrt(2.0 * Math.Pow(sigma, 4) / n);
        output.WriteLine($"iid: error / theory = {error / expected:F4}");
        Assert.InRange(error / expected, 0.9, 1.1);
    }

    [Fact]
    public void Ar1Energy_MatchesTheSquaredSeriesAutocorrelation()
    {
        // E_t = φE_{t−1} + √(1 − φ²)·z_t, stationary N(0, 1). E² has autocorrelation φ^{2k}, so
        // τ_int = (1 + φ²)/(1 − φ²) and C_v's error at T = 1 is √(2τ_int/n). Tolerance 15%.
        const int n = 100_000;
        const double phi = 0.9;
        var rng = new Xoshiro256PlusPlus(12);
        var e = new double[n];
        e[0] = Gaussian(rng);
        double innovation = Math.Sqrt(1.0 - phi * phi);
        for (int t = 1; t < n; t++)
            e[t] = phi * e[t - 1] + innovation * Gaussian(rng);

        double error = SpecificHeatMonteCarloError.Of([e], 1.0);

        double tau = (1.0 + phi * phi) / (1.0 - phi * phi);
        double expected = Math.Sqrt(2.0 * tau / n);
        output.WriteLine($"AR(1): error / theory = {error / expected:F4}");
        Assert.InRange(error / expected, 0.85, 1.15);
    }

    [Fact]
    public void SwReplicas_ErrorMatchesTheSpreadAcrossReplicas()
    {
        // 32 independent replicas of the q = 2 model on an 8 × 8 torus at T = 1.2, just above T_c. Each replica's own
        // error estimate, root-mean-squared, against the standard deviation of their C_v; the pooled error against
        // that deviation over √32. Each ratio within [0.7, 1.4].
        const int replicas = 32;
        const double temperature = 1.2;
        var graph = BinderCumulantTests.PeriodicLattice(8);
        Accumulator[] frames = Enumerable.Range(0, replicas).Select(r => SwRunner.Run(new SwRunSpec
        {
            Graph = graph, Temperature = temperature, Q = 2, ReplicaIndex = r, Seed = 2_000 + r,
            Budget = new RunBudget(500, 4_000),
        }).Accumulator).ToArray();

        double[] cv = frames.Select(f => Variance(f.EnergyPerDraw) / (temperature * temperature)).ToArray();
        double[] own = frames.Select(f => SpecificHeatMonteCarloError.Of([f.EnergyPerDraw], temperature)).ToArray();
        double spread = Math.Sqrt(cv.Sum(c => (c - cv.Average()) * (c - cv.Average())) / (replicas - 1));
        double rmsOwn = Math.Sqrt(own.Average(x => x * x));
        double pooled = SpecificHeatMonteCarloError.Of(frames.Select(f => f.EnergyPerDraw).ToArray(), temperature);

        output.WriteLine($"SW: rms own error / replica spread = {rmsOwn / spread:F4}; pooled / (spread/√R) = {pooled / (spread / Math.Sqrt(replicas)):F4}");
        Assert.InRange(rmsOwn / spread, 0.7, 1.4);
        Assert.InRange(pooled / (spread / Math.Sqrt(replicas)), 0.7, 1.4);
    }

    [Fact]
    public void SweepProfile_CarriesTheError_ByTemperature()
    {
        var graph = BinderCumulantTests.PeriodicLattice(4);
        SpcRunResult[] runs = [.. new[] { 1.5, 1.0 }.SelectMany(t => Enumerable.Range(0, 2).Select(r => new SpcRunResult
        {
            Graph = graph,
            Accumulator = SwRunner.Run(new SwRunSpec
            {
                Graph = graph, Temperature = t, Q = 2, ReplicaIndex = r, Seed = 30 + r, Budget = new RunBudget(50, 400),
            }).Accumulator,
        }))];

        SweepProfile profile = SweepProfile.From(runs);
        var (temperatures, values) = SpecificHeatMonteCarloError.ByTemperature(runs.Select(x => x.Accumulator).ToArray());

        Assert.Equal(new[] { 1.0, 1.5 }, temperatures);
        Assert.Equal(values, profile.AdditionalChannels["SpecificHeatMonteCarloError"]);
        Assert.All(values, v => Assert.True(v > 0.0));
    }

    [Fact]
    public void ConstantEnergy_HasNoError()
        => Assert.Equal(0.0, SpecificHeatMonteCarloError.Of([[-3.0, -3.0, -3.0, -3.0, -3.0]], 1.0));

    [Fact]
    public void UndefinedCases_AreNaN()
    {
        double[] e = [-1.0, -2.0, -3.0, -2.0, -1.0, -3.0];

        Assert.True(double.IsNaN(SpecificHeatMonteCarloError.Of([e], 0.0)));
        Assert.True(double.IsNaN(SpecificHeatMonteCarloError.Of([e, e[..5]], 1.0)));
        Assert.True(double.IsNaN(SpecificHeatMonteCarloError.Of([e[..3]], 1.0)));
        Assert.True(double.IsNaN(SpecificHeatMonteCarloError.Of(Array.Empty<double[]>(), 1.0)));
        // Two replicas frozen at energies symmetric about their mean: y is constant, but E never mixed.
        Assert.True(double.IsNaN(SpecificHeatMonteCarloError.Of([[-1.0, -1.0, -1.0, -1.0], [-3.0, -3.0, -3.0, -3.0]], 1.0)));
        // Three frozen replicas, asymmetric: y is frozen too.
        Assert.True(double.IsNaN(SpecificHeatMonteCarloError.Of(
            [[-1.0, -1.0, -1.0, -1.0], [-1.0, -1.0, -1.0, -1.0], [-4.0, -4.0, -4.0, -4.0]], 1.0)));
    }

    [Fact]
    public void ByTemperature_IsNaN_WhereAFrameHasNoRecord()
    {
        SpcRunResult recorded = EnergyCumulantTests.Run(1.0, [-1.0, -2.0, -3.0, -2.0, -1.0]);
        SpcRunResult bare = EnergyCumulantTests.Run(2.0, [-1.0, -2.0, -3.0, -2.0, -1.0]);
        bare = bare with { Accumulator = bare.Accumulator with { EnergyPerDraw = [] } };

        var (_, values) = SpecificHeatMonteCarloError.ByTemperature([recorded.Accumulator, bare.Accumulator]);

        Assert.False(double.IsNaN(values[0]));
        Assert.True(double.IsNaN(values[1]));
    }

    private static double Variance(double[] x)
    {
        double mean = x.Average();
        return x.Average(v => (v - mean) * (v - mean));
    }

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }
}
