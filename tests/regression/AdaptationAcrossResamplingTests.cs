using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// R07: adaptive τ must see every relocate attempt, whatever the resampling period. A resampler that keeps the
/// weights leaves the target unchanged, so a run with it must reproduce a run without it exactly — for resampling
/// periods that do and don't divide the 200-step adaptation window, and for burn-ins that end mid-window.
/// </summary>
public sealed class AdaptationAcrossResamplingTests
{
    [Theory]
    [InlineData(25, 1000)]
    [InlineData(30, 1000)]
    [InlineData(25, 950)]
    [InlineData(30, 950)]
    public void AdaptiveTau_IsUnaffectedByAResamplerThatKeepsTheWeights(int weightEvery, int burn)
    {
        (double[] x, double[] y) = Data();

        BarsResult plain = Ensemble(resampler: null)
            .Run(x, y, grid: x, chains: 2, masterSeed: 5, burn: burn, samples: 200);
        BarsResult resampled = Ensemble(new KeepWeights())
            .Run(x, y, grid: x, chains: 2, masterSeed: 5, burn: burn, samples: 200, weightEvery: weightEvery);

        Assert.Equal(plain.MeanKnots, resampled.MeanKnots);
        Assert.Equal(plain.Fit, resampled.Fit);
        Assert.Equal(plain.AcceptanceRate, resampled.AcceptanceRate);
    }

    private static BarsEnsemble Ensemble(IWeightResampler? resampler) =>
        new(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(4.0), new LocalBetaKernel(50.0),
            resampler: resampler, adaptKernelTau: 50.0);

    private static (double[] X, double[] Y) Data()
    {
        var rng = new Xoshiro256PlusPlus(seed: 17);
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

    // Unit weights every time, drawing nothing from the RNG: the target after a resample equals the one before.
    private sealed class KeepWeights : IWeightResampler
    {
        public double[] Resample(KnotConfig config, double[]? currentWeights, double[] x, double[] y,
                                 IBasis basis, IObservationModel model, Xoshiro256PlusPlus rng)
        {
            var w = new double[y.Length];
            Array.Fill(w, 1.0);
            return w;
        }
    }
}
