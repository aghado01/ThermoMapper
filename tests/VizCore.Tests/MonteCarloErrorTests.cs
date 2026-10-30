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
/// The sweep curves' Monte Carlo errors, sd/√ESS of a per-draw series (a mean) or of y = (x − x̄)² (a variance), against
/// references from outside the estimator: Gaussian theory for iid and AR(1) series, and the spread of each curve across
/// independent SW replicas. The tolerances were fixed before the first run.
/// </summary>
public sealed class MonteCarloErrorTests(ITestOutputHelper output)
{
    private const int Replicas = 32;
    private const double ReplicaTemperature = 1.2;

    // ── The variance form ────────────────────────────────────────────────────

    [Fact]
    public void Variance_IidGaussian_MatchesTheSampleVariancesError()
    {
        // x ~ N(−50, 2²) iid: Var((x − μ)²) = 2σ⁴, so the variance's error is √(2σ⁴/n). Tolerance 10%.
        const int n = 40_000;
        const double sigma = 2.0;
        double[] x = Iid(n, -50.0, sigma, seed: 11);

        double error = MonteCarloError.OfVariance([x]);

        double expected = Math.Sqrt(2.0 * Math.Pow(sigma, 4) / n);
        output.WriteLine($"variance, iid: error / theory = {error / expected:F4}");
        Assert.InRange(error / expected, 0.9, 1.1);
    }

    [Fact]
    public void Variance_Ar1_MatchesTheSquaredSeriesAutocorrelation()
    {
        // x_t = φx_{t−1} + √(1 − φ²)·z_t, stationary N(0, 1). x² has autocorrelation φ^{2k}, so
        // τ_int = (1 + φ²)/(1 − φ²) and the variance's error is √(2τ_int/n). Tolerance 15%.
        const int n = 100_000;
        const double phi = 0.9;
        double[] x = Ar1(n, phi, seed: 12);

        double error = MonteCarloError.OfVariance([x]);

        double tau = (1.0 + phi * phi) / (1.0 - phi * phi);
        double expected = Math.Sqrt(2.0 * tau / n);
        output.WriteLine($"variance, AR(1): error / theory = {error / expected:F4}");
        Assert.InRange(error / expected, 0.85, 1.15);
    }

    // ── The mean form ────────────────────────────────────────────────────────

    [Fact]
    public void Mean_IidGaussian_MatchesSigmaOverRootN()
    {
        // x ~ N(3, 2²) iid: the mean's error is σ/√n. Tolerance 10%.
        const int n = 40_000;
        const double sigma = 2.0;
        double[] x = Iid(n, 3.0, sigma, seed: 13);

        double error = MonteCarloError.OfMean([x]);

        double expected = sigma / Math.Sqrt(n);
        output.WriteLine($"mean, iid: error / theory = {error / expected:F4}");
        Assert.InRange(error / expected, 0.9, 1.1);
    }

    [Fact]
    public void Mean_Ar1_MatchesTheSeriesAutocorrelation()
    {
        // The AR(1) above: x has autocorrelation φ^k, so τ_int = (1 + φ)/(1 − φ) and the mean's error is √(τ_int/n).
        // Tolerance 15%.
        const int n = 100_000;
        const double phi = 0.9;
        double[] x = Ar1(n, phi, seed: 14);

        double error = MonteCarloError.OfMean([x]);

        double tau = (1.0 + phi) / (1.0 - phi);
        double expected = Math.Sqrt(tau / n);
        output.WriteLine($"mean, AR(1): error / theory = {error / expected:F4}");
        Assert.InRange(error / expected, 0.85, 1.15);
    }

    // ── Against the spread across SW replicas ────────────────────────────────

    [Fact]
    public void SwReplicas_EachCurvesError_MatchesTheSpreadAcrossReplicas()
    {
        // 32 independent replicas of the q = 2 model on an 8 × 8 torus at T = 1.2, just above T_c. For C_v, Var(m) and
        // χ_FK: each replica's own error, root-mean-squared, against the standard deviation of the replicas' values;
        // the pooled error against that deviation over √32. Each ratio within [0.7, 1.4].
        Accumulator[] frames = ReplicaFrames();
        double tSq = ReplicaTemperature * ReplicaTemperature;
        int n = frames[0].Spins.Length;

        Check("C_v",
            frames.Select(f => f.EnergyPerDraw).ToArray(),
            e => Variance(e) / tSq,
            replicas => MonteCarloError.OfVariance(replicas) / tSq);
        Check("Var(m)",
            frames.Select(f => f.MagPerDraw).ToArray(),
            Variance,
            MonteCarloError.OfVariance);
        Check("χ_FK",
            frames.Select(f => f.SumSqClusterSizesPerDraw.Select(s => (double)s / n).ToArray()).ToArray(),
            x => x.Average(),
            MonteCarloError.OfMean);

        void Check(string curve, double[][] series, Func<double[], double> value, Func<double[][], double> error)
        {
            double[] values = series.Select(value).ToArray();
            double spread = Math.Sqrt(values.Sum(v => (v - values.Average()) * (v - values.Average())) / (Replicas - 1));
            double rmsOwn = Math.Sqrt(series.Average(s => Math.Pow(error([s]), 2)));
            double pooled = error(series);

            output.WriteLine($"{curve}: rms own error / replica spread = {rmsOwn / spread:F4}; " +
                             $"pooled / (spread/√R) = {pooled / (spread / Math.Sqrt(Replicas)):F4}");
            Assert.InRange(rmsOwn / spread, 0.7, 1.4);
            Assert.InRange(pooled / (spread / Math.Sqrt(Replicas)), 0.7, 1.4);
        }
    }

    // ── The curves and SweepProfile ──────────────────────────────────────────

    [Fact]
    public void SweepProfile_CarriesEachError_ByTemperature()
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
        Accumulator[] frames = runs.Select(x => x.Accumulator).ToArray();

        SweepProfile profile = SweepProfile.From(runs);

        foreach ((string channel, (double[] Temperatures, double[] Values) curve) in new[]
        {
            ("SpecificHeatMonteCarloError", MonteCarloError.SpecificHeat(frames)),
            ("MagnetizationVarianceMonteCarloError", MonteCarloError.MagnetizationVariance(frames)),
            ("SusceptibilityFkClusterMonteCarloError", MonteCarloError.SusceptibilityFkCluster(frames)),
        })
        {
            Assert.Equal(new[] { 1.0, 1.5 }, curve.Temperatures);
            Assert.Equal(curve.Values, profile.AdditionalChannels[channel]);
            Assert.All(curve.Values, v => Assert.True(v > 0.0, $"{channel}: {v}"));
        }
    }

    [Fact]
    public void ConstantSeries_HaveNoError()
    {
        double[] constant = [-3.0, -3.0, -3.0, -3.0, -3.0];

        Assert.Equal(0.0, MonteCarloError.OfVariance([constant]));
        Assert.Equal(0.0, MonteCarloError.OfMean([constant, constant]));
    }

    [Fact]
    public void UndefinedCases_AreNaN()
    {
        double[] x = [-1.0, -2.0, -3.0, -2.0, -1.0, -3.0];
        double[][] frozenSymmetric = [[-1.0, -1.0, -1.0, -1.0], [-3.0, -3.0, -3.0, -3.0]];
        double[][] frozenAsymmetric = [[-1.0, -1.0, -1.0, -1.0], [-1.0, -1.0, -1.0, -1.0], [-4.0, -4.0, -4.0, -4.0]];

        foreach (Func<double[][], double> error in new Func<double[][], double>[] { MonteCarloError.OfVariance, MonteCarloError.OfMean })
        {
            Assert.True(double.IsNaN(error([x, x[..5]])));
            Assert.True(double.IsNaN(error([x[..3]])));
            Assert.True(double.IsNaN(error([])));
            // Frozen symmetric about the mean: a variance's y would be constant, but x never mixed.
            Assert.True(double.IsNaN(error(frozenSymmetric)));
            Assert.True(double.IsNaN(error(frozenAsymmetric)));
        }

        var frame = EnergyCumulantTests.Run(0.0, x).Accumulator;
        Assert.True(double.IsNaN(MonteCarloError.SpecificHeat([frame]).Values[0]));
    }

    [Fact]
    public void Curves_AreNaN_WhereAFrameHasNoRecord()
    {
        Accumulator recorded = EnergyCumulantTests.Run(1.0, [-1.0, -2.0, -3.0, -2.0, -1.0]).Accumulator;
        Accumulator bare = EnergyCumulantTests.Run(2.0, [-1.0, -2.0, -3.0, -2.0, -1.0]).Accumulator
            with { EnergyPerDraw = [], SumSqClusterSizesPerDraw = [] };

        var (_, cv) = MonteCarloError.SpecificHeat([recorded, bare]);
        var (_, chi) = MonteCarloError.SusceptibilityFkCluster([recorded, bare]);

        Assert.False(double.IsNaN(cv[0]));
        Assert.True(double.IsNaN(cv[1]));
        Assert.Equal(0.0, chi[0]);   // the helper's Σ|c|² is constant
        Assert.True(double.IsNaN(chi[1]));
    }

    // ── Helpers ──────────────────────────────────────────────────────────────

    private static Accumulator[] ReplicaFrames()
    {
        var graph = BinderCumulantTests.PeriodicLattice(8);
        return Enumerable.Range(0, Replicas).Select(r => SwRunner.Run(new SwRunSpec
        {
            Graph = graph, Temperature = ReplicaTemperature, Q = 2, ReplicaIndex = r, Seed = 2_000 + r,
            Budget = new RunBudget(500, 4_000),
        }).Accumulator).ToArray();
    }

    private static double Variance(double[] x)
    {
        double mean = x.Average();
        return x.Average(v => (v - mean) * (v - mean));
    }

    private static double[] Iid(int n, double mean, double sigma, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        return Enumerable.Range(0, n).Select(_ => mean + sigma * Gaussian(rng)).ToArray();
    }

    private static double[] Ar1(int n, double phi, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[n];
        x[0] = Gaussian(rng);
        double innovation = Math.Sqrt(1.0 - phi * phi);
        for (int t = 1; t < n; t++)
            x[t] = phi * x[t - 1] + innovation * Gaussian(rng);
        return x;
    }

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }
}
