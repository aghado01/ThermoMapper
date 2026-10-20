using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// R05: mutable proposal state belongs to one chain. Adaptive τ is requested declaratively and built per chain;
/// a tunable kernel instance cannot be shared across chains or tempering replicas.
/// </summary>
public sealed class KernelOwnershipTests
{
    [Fact]
    public void SuppliedAdaptiveKernel_IsRejected()
    {
        Assert.Throws<ArgumentException>(() =>
            new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(4.0),
                             new AdaptiveLocalBetaKernel(50.0)));
    }

    [Fact]
    public void SuppliedAdaptiveKernel_WithTempering_IsRejected()
    {
        Assert.Throws<ArgumentException>(() =>
            new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(4.0),
                             new AdaptiveLocalBetaKernel(50.0), temperLevels: 3));
    }

    [Fact]
    public void AdaptiveTau_WithTempering_IsRejected()
    {
        Assert.Throws<ArgumentException>(() =>
            new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(4.0),
                             new LocalBetaKernel(50.0), adaptKernelTau: 50.0, temperLevels: 3));
    }

    [Fact]
    public void ConcurrentChains_WithAdaptiveTau_AreDeterministic()
    {
        (double[] x, double[] y) = Data();
        BarsResult first = Run(x, y);
        for (int i = 0; i < 3; i++)
        {
            BarsResult again = Run(x, y);
            Assert.Equal(first.MeanKnots, again.MeanKnots);
            Assert.Equal(first.Fit, again.Fit);
            Assert.Equal(first.AcceptanceRate, again.AcceptanceRate);
        }
    }

    private static BarsResult Run(double[] x, double[] y) =>
        new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(4.0),
                         new LocalBetaKernel(50.0), adaptKernelTau: 50.0)
            .Run(x, y, grid: x, chains: 4, masterSeed: 9, burn: 1000, samples: 300);

    private static (double[] X, double[] Y) Data()
    {
        var rng = new Xoshiro256PlusPlus(seed: 23);
        const int n = 60;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            double u1 = 1.0 - rng.NextDouble(), u2 = 1.0 - rng.NextDouble();
            y[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.15 * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }
        return (x, y);
    }
}
