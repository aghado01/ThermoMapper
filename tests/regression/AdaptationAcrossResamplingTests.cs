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

        // The knot trajectory must be identical. The fit agrees to rounding: unit weights take the weighted marginal's
        // path, which factors Zᵀ(W + I/m)Z rather than rescaling ZᵀZ.
        Assert.Equal(plain.MeanKnots, resampled.MeanKnots);
        Assert.Equal(plain.AcceptanceRate, resampled.AcceptanceRate);
        for (int i = 0; i < plain.Fit.Length; i++) Assert.Equal(plain.Fit[i], resampled.Fit[i], 1e-12);
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

    // Unit weights every time: the target after a resample equals the one before. The weight block draws from its
    // own stream, so the chain's trajectory is untouched.
    private sealed class KeepWeights : IWeightResampler
    {
        public double[] Resample(ReadOnlySpan<double> residuals, double sigma2, Xoshiro256PlusPlus rng)
        {
            var w = new double[residuals.Length];
            Array.Fill(w, 1.0);
            return w;
        }

        // Unit weights are the Normal error model itself.
        public double MarginalLogDensity(double residual, double sigma2)
            => -0.5 * Math.Log(2.0 * Math.PI * sigma2) - 0.5 * residual * residual / sigma2;
    }
}
