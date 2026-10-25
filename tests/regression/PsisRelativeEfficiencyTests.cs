using System;
using System.Linq;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// PSIS-LOO's relative efficiency and Monte Carlo error (the loo package's r_eff and mcse_elpd_loo), against truths from
/// outside the estimator. The model is y_i ~ N(θ, 1) under a flat prior, so θ | y ~ N(ȳ, 1/n). Its draws are an AR(1)
/// chain started at stationarity, θ_d − ȳ = ρ(θ_{d−1} − ȳ) + √((1 − ρ²)/n)·z, so every draw is exact and only their
/// autocorrelation is chosen. Each unit's r_eff is then exact: with g(θ) = exp(−(y_i − θ)²/2), s² = 1/n and
/// δ = y_i − ȳ, two draws correlated r have E[g g] = ((1 + s²)² − s⁴r²)^(−1/2)·exp(−δ²/(1 + s²(1 + r))), so g's lag-k
/// correlation follows at r = ρᵏ. A unit far from ȳ, where g is nearly linear in θ, has about the AR(1) value
/// (1 − ρ)/(1 + ρ); one at ȳ, where g is quadratic, about (1 − ρ²)/(1 + ρ²). The spread of the estimates across
/// independent runs is the Monte Carlo error itself.
/// </summary>
public sealed class PsisRelativeEfficiencyTests
{
    private readonly ITestOutputHelper _out;
    public PsisRelativeEfficiencyTests(ITestOutputHelper output) => _out = output;

    private const int N = 20;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double[] Data()
    {
        var rng = new Xoshiro256PlusPlus(11);
        return Enumerable.Range(0, N).Select(_ => 0.5 + Gaussian(rng)).ToArray();
    }

    // Chains of exact posterior draws with lag-one autocorrelation rho, scored as log N(y_i; θ_d, 1).
    private static PointwiseLogDensity Density(double[] y, double rho, int chains, int length, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double mean = y.Average(), sd = 1.0 / Math.Sqrt(N), step = sd * Math.Sqrt(1.0 - rho * rho);
        var ll = new double[chains * length, N];
        int row = 0;
        for (int c = 0; c < chains; c++)
        {
            double theta = mean + sd * Gaussian(rng);
            for (int d = 0; d < length; d++, row++)
            {
                theta = mean + rho * (theta - mean) + step * Gaussian(rng);
                for (int i = 0; i < N; i++) ll[row, i] = -0.5 * Math.Log(2.0 * Math.PI) - 0.5 * (y[i] - theta) * (y[i] - theta);
            }
        }
        var report = new[] { new LatentReport("θ", LatentScope.Global, LatentTreatment.Conditioned) };
        return new PointwiseLogDensity(ll, UnitPartition.EachObservation(N), report, Enumerable.Repeat(length, chains).ToArray());
    }

    // The exact r_eff of g(θ) = exp(−(y_i − θ)²/2) under the stationary chain: 1/τ, τ = 1 + 2Σₖ corr(g(θ_t), g(θ_{t+k})).
    private static double ExactRelativeEfficiency(double delta, double rho)
    {
        double s2 = 1.0 / N;
        double Product(double r) =>
            Math.Exp(-delta * delta / (1.0 + s2 * (1.0 + r))) / Math.Sqrt((1.0 + s2) * (1.0 + s2) - s2 * s2 * r * r);
        double mean = Math.Exp(-delta * delta / (2.0 * (1.0 + s2))) / Math.Sqrt(1.0 + s2);
        double variance = Product(1.0) - mean * mean, tau = 1.0, lag = 1.0;
        for (int k = 1; k < 2000; k++)
        {
            lag *= rho;
            tau += 2.0 * (Product(lag) - mean * mean) / variance;
        }
        return 1.0 / tau;
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(0.8)]
    public void RelativeEfficiency_ReadsTheChains(double rho)
    {
        double[] y = Data();
        PointwiseLogDensity density = Density(y, rho, chains: 4, length: 2000, seed: 3);
        double[] rEff = PredictiveChecks.PsisLoo(density).RelativeEfficiency!;
        double[] exact = y.Select(v => ExactRelativeEfficiency(v - y.Average(), rho)).ToArray();
        double[] ratio = rEff.Zip(exact, (r, e) => r / e).ToArray();
        _out.WriteLine($"ρ = {rho}: exact r_eff {exact.Min():F3}–{exact.Max():F3}; estimated/exact {ratio.Min():F3}–{ratio.Max():F3}");
        Assert.All(ratio, q => Assert.InRange(q, 0.85, 1.15));
        // Without chains the draws count as independent.
        Assert.All(PredictiveChecks.PsisLoo(density with { ChainLengths = null }).RelativeEfficiency!, r => Assert.Equal(1.0, r));
    }

    /// <summary>
    /// Across 200 independent runs of four chains of 250 autocorrelated draws, each unit's reported Monte Carlo error
    /// matches the spread of its estimates. Read as independent draws, the same terms understate it by √r_eff, the exact one. The
    /// total's error, √Σ of the units', assumes the units' errors independent, as loo's does; they share the draws, so the
    /// test reports how the total fares without asserting it.
    /// </summary>
    [Fact]
    public void MonteCarloError_MatchesTheSpreadAcrossRuns()
    {
        double[] y = Data();
        const int runs = 200;
        var estimates = new double[runs, N];
        var withChains = new double[runs, N];
        var asIndependent = new double[runs, N];
        var totals = new double[runs];
        var totalErrors = new double[runs];
        for (int r = 0; r < runs; r++)
        {
            PointwiseLogDensity density = Density(y, rho: 0.8, chains: 4, length: 250, seed: 1000 + r);
            PredictiveCheck loo = PredictiveChecks.PsisLoo(density);
            PredictiveCheck independent = PredictiveChecks.PsisLoo(density with { ChainLengths = null });
            for (int i = 0; i < N; i++)
            {
                estimates[r, i] = loo.PerUnit[i];
                withChains[r, i] = loo.MonteCarloErrors![i];
                asIndependent[r, i] = independent.MonteCarloErrors![i];
            }
            totals[r] = loo.Elpd;
            totalErrors[r] = loo.MonteCarloError;
        }

        static double Sd(double[] v)
        {
            double m = v.Average();
            return Math.Sqrt(v.Sum(x => (x - m) * (x - m)) / (v.Length - 1));
        }
        static double Rms(double[] v) => Math.Sqrt(v.Average(x => x * x));
        double[] Column(double[,] a, int i) => Enumerable.Range(0, runs).Select(r => a[r, i]).ToArray();

        var ratio = new double[N];
        var naive = new double[N];
        for (int i = 0; i < N; i++)
        {
            double spread = Sd(Column(estimates, i));
            ratio[i] = Rms(Column(withChains, i)) / spread;
            naive[i] = Rms(Column(asIndependent, i)) / spread / Math.Sqrt(ExactRelativeEfficiency(y[i] - y.Average(), 0.8));
        }
        _out.WriteLine($"reported/actual, with r_eff: {ratio.Min():F2}–{ratio.Max():F2}; as independent, over √r_eff: {naive.Min():F2}–{naive.Max():F2}");
        _out.WriteLine($"total: reported {totalErrors.Average():F4} against an actual spread of {Sd(totals):F4}");
        Assert.All(ratio, q => Assert.InRange(q, 0.8, 1.25));
        Assert.All(naive, q => Assert.InRange(q, 0.8, 1.25));
    }

    [Fact]
    public void MonteCarloError_IsWithheld_AboveTheThreshold()
    {
        Assert.Equal(0.7, PredictiveChecks.ParetoKThreshold(4000));
        Assert.Equal(1.0 - 1.0 / Math.Log10(400), PredictiveChecks.ParetoKThreshold(400), 1e-15);
        // One unit with a heavy-tailed ratio: its log terms are Cauchy-like across the draws.
        var rng = new Xoshiro256PlusPlus(5);
        var ll = new double[400, 2];
        for (int d = 0; d < 400; d++)
        {
            ll[d, 0] = -1.0 + 0.01 * Gaussian(rng);
            ll[d, 1] = -2.0 * Math.Log(1.0 + Math.Abs(Math.Tan(Math.PI * (rng.NextDouble() - 0.5))));
        }
        var report = new[] { new LatentReport("θ", LatentScope.Global, LatentTreatment.Conditioned) };
        PredictiveCheck loo = PredictiveChecks.PsisLoo(new PointwiseLogDensity(ll, UnitPartition.EachObservation(2), report));
        Assert.True(loo.ParetoK![1] > PredictiveChecks.ParetoKThreshold(400), $"k̂ {loo.ParetoK[1]:F2}");
        Assert.True(double.IsNaN(loo.MonteCarloError));
        Assert.All(loo.MonteCarloErrors!, e => Assert.True(double.IsFinite(e)));
        Assert.True(double.IsNaN(PredictiveChecks.Waic(new PointwiseLogDensity(ll, UnitPartition.EachObservation(2), report)).MonteCarloError));
    }
}
