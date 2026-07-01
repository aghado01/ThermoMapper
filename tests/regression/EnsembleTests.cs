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
        var sums = new[] { 15.0, 15.0, 15.0, 15.0 };      // four chains of {1,2,3,4,5}
        var sqs = new[] { 55.0, 55.0, 55.0, 55.0 };
        Assert.Equal(Math.Sqrt(0.8), ChainDiagnostics.RHat(sums, sqs, 5), 6);

        var sumsShift = new[] { 15.0, 20.0, 25.0, 30.0 }; // shifted by a constant each
        var sqsShift = new[] { 55.0, 90.0, 135.0, 190.0 };
        Assert.True(ChainDiagnostics.RHat(sumsShift, sqsShift, 5) > 1.1);
    }

    [Fact]
    public void Ess_IidChains_NearTotal_AutocorrelatedChains_Reduced()
    {
        var rng = new Xoshiro256PlusPlus(seed: 1);
        int c = 4, n = 600;
        var iid = new double[c][];
        for (int j = 0; j < c; j++)
        {
            var a = new double[n];
            for (int i = 0; i < n; i++) a[i] = Gaussian(rng);
            iid[j] = a;
        }
        var ar = new double[c][];
        for (int j = 0; j < c; j++)
        {
            var a = new double[n];
            double prev = 0.0;
            for (int i = 0; i < n; i++) { prev = 0.8 * prev + Gaussian(rng); a[i] = prev; }
            ar[j] = a;
        }

        double total = c * (double)n;
        double essIid = ChainDiagnostics.Ess(iid);
        double essAr = ChainDiagnostics.Ess(ar);
        _out.WriteLine($"ESS iid={essIid:F0}/{total} ar1={essAr:F0}/{total}");

        Assert.True(essIid > 0.5 * total, $"iid ESS {essIid} should be near total {total}");
        Assert.True(essAr < 0.5 * total, $"AR(1) ESS {essAr} should be well below total {total}");
        Assert.True(essAr < essIid);
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
        _out.WriteLine($"R̂(fit)={result.RHatFitMax:F3} R̂(k)={result.RHatKnots:F3} ESS(k)={result.EssKnots:F0} " +
                       $"meanK={result.MeanKnots:F2} used={result.SamplesUsed} mseFit={mseFit:F5}");

        Assert.True(result.RHatFitMax < 1.3);
        Assert.True(result.RHatKnots < 1.3);
        Assert.Equal(x.Length, result.RHatFit.Length);          // consensus map sized to the grid
        Assert.True(result.EssKnots > 0.0 && result.EssKnots <= 4 * 1500);
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
        BarsResult result = ensemble.Run(x, y, grid: x, chains: 4, masterSeed: 7, burn: 1000, samples: 6000,
                                         rHatTarget: 1.5, batchSize: 500);

        _out.WriteLine($"used={result.SamplesUsed}/6000 R̂(fit)={result.RHatFitMax:F3} R̂(k)={result.RHatKnots:F3}");
        Assert.True(result.SamplesUsed < 6000, $"adaptive stop used all {result.SamplesUsed} samples");
        Assert.True(result.SamplesUsed >= 500);
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
