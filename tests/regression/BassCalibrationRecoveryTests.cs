using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Recovery of a known θ* (stage 7, C3; D09): a synthetic simulator η(t; θ) = θ₁·sin(2π(t + θ₂/4)) + θ₂·t/2 on 20
/// locations, run at 36 inputs on a 6 × 6 grid of θ ∈ [0.5, 2] × [0, 1]; BASS fitted once on a PCA basis and once directly
/// on the functional response; field data simulated at θ* = (1.3, 0.4) with Gaussian noise, and calibrated by the nested
/// reference sampler under each emulator, then with a discrepancy in the field data and declared, and with two variance
/// groups of different noise. Each case runs over five field-noise seeds and counts the seeds whose central 95% interval
/// holds each θ*ⱼ.
/// </summary>
public sealed class BassCalibrationRecoveryTests
{
    private readonly ITestOutputHelper _out;
    public BassCalibrationRecoveryTests(ITestOutputHelper output) => _out = output;

    private const int M = 20;
    private static readonly double[] T = [.. Enumerable.Range(0, M).Select(l => l / (M - 1.0))];
    private static readonly double[] ThetaStar = [1.3, 0.4];
    private static readonly ContinuousParameter Theta1 = new("theta1", 0.5, 2.0), Theta2 = new("theta2", 0.0, 1.0);

    private static double Eta(double t, double theta1, double theta2) => theta1 * Math.Sin(2 * Math.PI * (t + 0.25 * theta2)) + 0.5 * theta2 * t;

    private static (double[,] Y, BassPredictor[] Predictors) Runs()
    {
        var y = new double[36, M];
        var a = new double[36];
        var b = new double[36];
        for (int i = 0; i < 36; i++)
        {
            a[i] = Theta1.ToPhysical((i / 6) / 5.0);
            b[i] = Theta2.ToPhysical((i % 6) / 5.0);
            for (int l = 0; l < M; l++) y[i, l] = Eta(T[l], a[i], b[i]);
        }
        return (y, [new ContinuousPredictor("theta1", a), new ContinuousPredictor("theta2", b)]);
    }

    private static BassConfig FitConfig(int seed) => BassRunTests.Config(chains: 2, burn: 1500, samples: 100, seed: seed, maxBasis: 12);

    private static readonly Lazy<BassBasisFit> BasisFit = new(() =>
    {
        (double[,] y, BassPredictor[] predictors) = Runs();
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 0.9999 });
        return BassComponents.Fit(d, predictors, FitConfig(4));
    });

    private static readonly Lazy<(BassData Data, BassFit Fit)> DirectFit = new(() =>
    {
        (double[,] y, BassPredictor[] predictors) = Runs();
        var data = new BassData(predictors, new FunctionalResponse(y, [new ContinuousPredictor("t", T)]));
        return (data, Bass.Fit(FitConfig(6), data));
    });

    private static double[] Field(int seed, double[] sd, Func<double, double>? delta = null)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        return [.. Enumerable.Range(0, M).Select(l => Eta(T[l], ThetaStar[0], ThetaStar[1]) + (delta?.Invoke(T[l]) ?? 0.0) + sd[l] * Normal.SampleStandard(rng))];
    }

    private static double WeightedQuantile(double[] values, double[] weights, double p)
    {
        int[] order = [.. Enumerable.Range(0, values.Length).OrderBy(i => values[i])];
        double total = weights.Sum(), cumulative = 0.0;
        foreach (int i in order)
            if ((cumulative += weights[i] / total) >= p) return values[i];
        return values[order[^1]];
    }

    private (int Theta1, int Theta2) Coverage(string label, ICalibrationEmulator emulator, Func<int, double[]> field, GaussianMeasurement groups,
                                              CovarianceFormulation covariance, int seeds = 5)
    {
        int[] covered = new int[2];
        for (int seed = 1; seed <= seeds; seed++)
        {
            var config = new CalibrationConfig
            {
                Parameters = [Theta1, Theta2],
                Measurement = groups,
                Covariance = covariance,
                Sampling = new CalibrationSampling { MasterSeed = seed, BurnIn = 600, DecorrelationEvery = 50 },
                Reference = new ReferenceSampling { OuterDraws = 6, Chains = 2, MaxSamples = 1500, BatchSize = 300, RHatTarget = 1.05, MinEssPerChain = 50 },
            };
            NestedCalibration result = NestedCalibrationSampler.Run(new CalibrationModel(config, field(seed), emulator));
            (IReadOnlyList<CalibrationDraw> draws, double[] weights) = result.Pooled();
            Assert.NotEmpty(draws);
            var line = new List<string>();
            for (int j = 0; j < 2; j++)
            {
                double[] values = [.. draws.Select(d => d.Theta[j])];
                double lo = WeightedQuantile(values, weights, 0.025), hi = WeightedQuantile(values, weights, 0.975);
                double mean = values.Zip(weights, (v, w) => v * w).Sum() / weights.Sum();
                if (ThetaStar[j] >= lo && ThetaStar[j] <= hi) covered[j]++;
                line.Add($"θ{j + 1} {mean:F3} [{lo:F3}, {hi:F3}]");
            }
            _out.WriteLine($"{label}, seed {seed}: {string.Join("; ", line)}; {result.Runs.Count(r => r.Converged)}/{result.Runs.Count} runs converged " +
                           $"({string.Join(",", result.Runs.Select(r => r.SamplesPerChain))} draws per chain)");
        }
        _out.WriteLine($"{label}: θ* covered {covered[0]}/{seeds} and {covered[1]}/{seeds}");
        return (covered[0], covered[1]);
    }

    private static GaussianMeasurement OneGroup() => new([new VarianceGroup([.. Enumerable.Range(0, M)], new InverseGammaVariance(4, 0.05))]);

    [Fact]
    public void BasisEmulator_RecoversThetaStar()
    {
        var emulator = new BassBasisEmulator(BasisFit.Value, [Theta1, Theta2]);
        _out.WriteLine($"basis: rank {BasisFit.Value.Decomposition.Rank}, {emulator.DrawCount} joint draws, residual RMS {BasisFit.Value.Decomposition.Truncation.Rms:E2}");
        double[] sd = [.. Enumerable.Repeat(0.05, M)];
        (int a, int b) = Coverage("basis", emulator, seed => Field(seed, sd), OneGroup(), new FunctionalFormulation());
        Assert.True(a >= 4 && b >= 4);
    }

    [Fact]
    public void DirectEmulator_RecoversThetaStar()
    {
        (BassData data, BassFit fit) = DirectFit.Value;
        var emulator = new BassFitEmulator(fit, data, [Theta1, Theta2]);
        Assert.Equal(M, emulator.Outputs);
        double[] sd = [.. Enumerable.Repeat(0.05, M)];
        (int a, int b) = Coverage("direct", emulator, seed => Field(seed, sd), OneGroup(), new FunctionalFormulation());
        Assert.True(a >= 4 && b >= 4);
    }

    [Fact]
    public void WithADeclaredDiscrepancy_RecoversThetaStar()
    {
        var emulator = new BassBasisEmulator(BasisFit.Value, [Theta1, Theta2]);
        double[] sd = [.. Enumerable.Repeat(0.05, M)];
        double[][] basis = [.. T.Select(t => new[] { 1.0, t - 0.5 })];
        var covariance = new FunctionalFormulation { Discrepancy = new FunctionalDiscrepancy(basis, 0.2) };
        (int a, int b) = Coverage("discrepancy", emulator, seed => Field(seed, sd, t => 0.1 + 0.1 * (t - 0.5)), OneGroup(), covariance);
        Assert.True(a >= 4 && b >= 4);
    }

    [Fact]
    public void WithTwoVarianceGroups_RecoversThetaStar()
    {
        var emulator = new BassBasisEmulator(BasisFit.Value, [Theta1, Theta2]);
        double[] sd = [.. Enumerable.Range(0, M).Select(l => l < M / 2 ? 0.03 : 0.08)];
        var groups = new GaussianMeasurement([new VarianceGroup([.. Enumerable.Range(0, M / 2)], new InverseGammaVariance(4, 0.05)),
                                              new VarianceGroup([.. Enumerable.Range(M / 2, M / 2)], new HalfCauchyDeviation(0.05))]);
        (int a, int b) = Coverage("two groups", emulator, seed => Field(seed, sd), groups, new FunctionalFormulation());
        Assert.True(a >= 4 && b >= 4);
    }
}
