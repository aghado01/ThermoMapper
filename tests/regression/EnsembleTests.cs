using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The chain-ensemble layer: the R̂ and ESS formulae against hand-computed / synthetic chains, a converged
/// multi-chain run (R̂ near 1, consensus map sized, pooled fit denoises), the R̂ adaptive stop, and the facade.
/// </summary>
public sealed class EnsembleTests
{
    private readonly ITestOutputHelper _out;
    public EnsembleTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double Mse(double[] a, double[] b)
    {
        double s = 0.0;
        for (int i = 0; i < a.Length; i++) { double d = a[i] - b[i]; s += d * d; }
        return s / a.Length;
    }

    [Fact]
    public void RHat_IdenticalChains_BelowOne_ShiftedChains_AboveOne()
    {
        // Four chains of {1,2,3,4,5}: each mean 3, each sample variance 2.5.
        Assert.Equal(Math.Sqrt(0.8), ChainDiagnostics.RHat([3.0, 3.0, 3.0, 3.0], [2.5, 2.5, 2.5, 2.5], 5).Value, 12);
        // The same chains shifted by 0, 1, 2 and 3.
        Assert.True(ChainDiagnostics.RHat([3.0, 4.0, 5.0, 6.0], [2.5, 2.5, 2.5, 2.5], 5).Value > 1.1);
    }

    [Fact]
    public void Ensemble_Converges_RHatNearOne_ConsensusSized_AndDenoises()
    {
        var rng = new Xoshiro256PlusPlus(seed: 31);
        int n = 100;
        var x = new double[n];
        var f = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            f[i] = Math.Sin(2.0 * Math.PI * x[i]);
            y[i] = f[i] + 0.15 * Gaussian(rng);
        }

        var ensemble = new BarsEnsemble(
            new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(4.0), new LocalBetaKernel(50.0));
        BarsResult result = ensemble.Run(x, y, grid: x, chains: 4, masterSeed: 7, burn: 1500, samples: 1500);

        double mseFit = Mse(result.Fit, f);
        _out.WriteLine($"R̂(fit)={result.RHatFitMax:F3} R̂(k)={result.Knots.RHat:F3} ESS(k)={result.Knots.BulkEss:F0} " +
                       $"meanK={result.MeanKnots:F2} used={result.SamplesUsed} mseFit={mseFit:F5}");

        Assert.Equal(DiagnosticState.Estimated, result.RHatFitMax.State);
        Assert.True(result.RHatFitMax.MeetsTarget(1.3));
        Assert.True(result.Knots.RHat.MeetsTarget(1.3));
        Assert.Equal(x.Length, result.RHatFit.Length);          // consensus map sized to the grid
        Assert.True(result.Knots.BulkEss.Value > 0.0 && result.Knots.BulkEss.Value <= 4 * 1500);
        Assert.Equal(1500, result.SamplesUsed);
        Assert.True(mseFit < 0.5 * Mse(y, f));
        Assert.Equal(4, result.ChainSeeds.Length);
    }

    [Fact]
    public void AdaptiveStop_HaltsBeforeMax_OnConvergence()
    {
        var rng = new Xoshiro256PlusPlus(seed: 12);
        int n = 100;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.15 * Gaussian(rng);
        }

        var ensemble = new BarsEnsemble(
            new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(4.0), new LocalBetaKernel(50.0));
        // The knot count's autocorrelation time here is about 450 steps (batch means agree), so the default gate of
        // 100 effective draws per chain would need about 45,000 draws per chain; a gate of 10 keeps the test short.
        BarsResult result = ensemble.Run(x, y, grid: x, chains: 4, masterSeed: 7, burn: 1000, samples: 40_000,
                                         rHatTarget: 1.5, batchSize: 500, minEssPerChain: 10.0);

        _out.WriteLine($"used={result.SamplesUsed}/40000 R̂(fit)={result.RHatFitMax:F3} R̂(k)={result.Knots.RHat:F3} " +
                       $"ESS(k) bulk={result.Knots.BulkEss:F0} tail={result.Knots.TailEss:F0}");
        Assert.True(result.SamplesUsed < 40_000, $"adaptive stop used all {result.SamplesUsed} samples");
        // The stop waits for the knot count's bulk and tail ESS, not for a draw count.
        Assert.True(result.Knots.Ess.MeetsMinimum(4 * 10.0), $"stopped at ESS {result.Knots.Ess:F0}");
    }

    [Fact]
    public void Facade_BarsRun_Denoises()
    {
        var rng = new Xoshiro256PlusPlus(seed: 88);
        int n = 100;
        var x = new double[n];
        var f = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            f[i] = Math.Sin(2.0 * Math.PI * x[i]);
            y[i] = f[i] + 0.15 * Gaussian(rng);
        }

        var config = new BarsConfig { Prior = new PoissonPrior(4.0), Chains = 3, MasterSeed = 5, BurnIn = 1000, MaxSamples = 1500 };
        BarsResult result = Bars.Run(config, x, y, x);

        Assert.True(Mse(result.Fit, f) < 0.5 * Mse(y, f));
        Assert.Equal(3, result.ChainSeeds.Length);
    }
}
