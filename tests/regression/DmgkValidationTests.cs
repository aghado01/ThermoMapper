using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The DMGK miniboss: fit the full free-knot sampler (LocalBeta locality + weighted-Normal marginal +
/// Poisson prior) to noisy data from three function shapes — smooth, a sharp localized peak, and a jump —
/// and check the posterior-mean fit denoises (closer to truth than the raw data). The smooth/peak cases get
/// strong-denoise thresholds; the jump is softer because v1 has no coincident-knot mechanism yet, so it
/// approximates the discontinuity with a steep continuous spline. Each case runs under a Conway–Maxwell–Poisson count
/// prior at ν = 1 (Poisson) and ν = 2 (the law BARS ran under before the R01 measure fix), at the same thresholds.
/// </summary>
public sealed class DmgkValidationTests
{
    private readonly ITestOutputHelper _out;
    public DmgkValidationTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    /// <summary>Run the sampler and return the chain-averaged posterior-mean fit evaluated at <paramref name="x"/>,
    /// and the mean knot count.</summary>
    private static (double[] Fit, double MeanKnots) PosteriorMeanFit(
        double[] x, double[] y, IComplexityPrior prior, double tau, int seed, int burn, int samples)
    {
        var basis = new SplineBasis(degree: 3);
        var model = new WeightedNormalModel();
        var target = new SplineTarget(basis, model, prior, x, y);
        var kernel = new LocalBetaKernel(tau);
        var moves = new IRjMove<KnotConfig>[]
        {
            new KnotBirthMove(kernel), new KnotDeathMove(kernel), new KnotRelocateMove(kernel),
        };
        var chain = new ReversibleJumpChain<KnotConfig>(
            moves, target, new KnotConfig(Array.Empty<double>()), new Xoshiro256PlusPlus(seed));

        for (int i = 0; i < burn; i++) chain.Step();

        var sum = new double[x.Length];
        double knots = 0.0;
        for (int s = 0; s < samples; s++)
        {
            KnotConfig cfg = chain.Step();
            knots += cfg.Count;
            double[,] zTrain = basis.Design(cfg, x);
            double[] c = model.PosteriorMeanCoefficients(zTrain, y, null);
            for (int g = 0; g < x.Length; g++)
            {
                double f = 0.0;
                for (int j = 0; j < c.Length; j++) f += zTrain[g, j] * c[j];
                sum[g] += f;
            }
        }
        for (int g = 0; g < sum.Length; g++) sum[g] /= samples;
        return (sum, knots / samples);
    }

    private static double Mse(double[] a, double[] b)
    {
        double s = 0.0;
        for (int i = 0; i < a.Length; i++) { double d = a[i] - b[i]; s += d * d; }
        return s / a.Length;
    }

    private void RunCase(string name, Func<double, double> truth, double sigma, int n,
                         double lambda, double nu, double maxFitToDataRatio)
    {
        (double mseFit, double mseData, double meanKnots) =
            Fit(truth, sigma, n, new ConwayMaxwellPoissonPrior(lambda, nu), dataSeed: 2026, chainSeed: 99);
        _out.WriteLine($"{name} ν={nu}: mseFit={mseFit:F5} mseData={mseData:F5} ratio={mseFit / mseData:F3} k̄={meanKnots:F2}");

        Assert.True(mseFit < maxFitToDataRatio * mseData,
            $"{name} ν={nu}: fit MSE {mseFit:F5} not below {maxFitToDataRatio}× data MSE {mseData:F5}");
    }

    private static (double MseFit, double MseData, double MeanKnots) Fit(
        Func<double, double> truth, double sigma, int n, IComplexityPrior prior, int dataSeed, int chainSeed)
    {
        var rng = new Xoshiro256PlusPlus(dataSeed);
        var x = new double[n];
        var f = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            f[i] = truth(x[i]);
            y[i] = f[i] + sigma * Gaussian(rng);
        }

        (double[] fit, double meanKnots) = PosteriorMeanFit(x, y, prior, tau: 50.0, seed: chainSeed, burn: 1500, samples: 3000);
        return (Mse(fit, f), Mse(y, f), meanKnots);
    }

    private static readonly Func<double, double> SmoothTruth = t => Math.Sin(2.0 * Math.PI * t);
    private static readonly Func<double, double> PeakTruth = t => Math.Exp(-150.0 * (t - 0.4) * (t - 0.4));
    private static readonly Func<double, double> JumpTruth = t => (t < 0.5 ? 0.2 : 0.8) + 0.3 * t;

    [Theory]
    [InlineData(1.0)]
    [InlineData(2.0)]
    public void Smooth_StrongDenoise(double nu)
        => RunCase("smooth", SmoothTruth, sigma: 0.15, n: 120, lambda: 4.0, nu, maxFitToDataRatio: 0.5);

    [Theory]
    [InlineData(1.0)]
    [InlineData(2.0)]
    public void SharpPeak_Denoise(double nu)
        => RunCase("peak", PeakTruth, sigma: 0.10, n: 150, lambda: 5.0, nu, maxFitToDataRatio: 0.6);

    [Theory]
    [InlineData(1.0)]
    [InlineData(2.0)]
    public void Discontinuity_StillFits(double nu)
        => RunCase("jump", JumpTruth, sigma: 0.08, n: 150, lambda: 6.0, nu, maxFitToDataRatio: 1.0);
}
