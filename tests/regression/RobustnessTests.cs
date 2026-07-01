using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Scale-mixture robustness: with gross outliers injected, the Student-t weight-Gibbs fit beats the plain
/// least-squares (unweighted) fit on MSE-to-truth — the weights downweight the outliers the LS fit chases.
/// </summary>
public sealed class RobustnessTests
{
    private readonly ITestOutputHelper _out;
    public RobustnessTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void Robust_DownweightsOutliers_BeatsLeastSquares()
    {
        var rng = new Xoshiro256PlusPlus(seed: 44);
        int n = 100;
        var x = new double[n];
        var f = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            f[i] = Math.Sin(2.0 * Math.PI * x[i]);
            y[i] = f[i] + 0.08 * Gaussian(rng);
        }
        foreach (int b in new[] { 8, 21, 37, 50, 63, 79, 88, 95 })
            y[b] = f[b] + 3.5;   // gross outliers

        var basis = new SplineBasis(degree: 3);
        var model = new WeightedNormalModel();
        var prior = new PoissonPrior(4.0);
        var kernel = new LocalBetaKernel(50.0);

        BarsResult ls = new BarsEnsemble(basis, model, prior, kernel)
            .Run(x, y, grid: x, chains: 2, masterSeed: 9, burn: 1000, samples: 1500);
        BarsResult robust = new BarsEnsemble(basis, model, prior, kernel, new StudentTWeights(4.0))
            .Run(x, y, grid: x, chains: 2, masterSeed: 9, burn: 1000, samples: 1500);

        double mseLs = Mse(ls.Fit, f);
        double mseRobust = Mse(robust.Fit, f);
        _out.WriteLine($"mseLS={mseLs:F5} mseRobust={mseRobust:F5} ratio={mseRobust / mseLs:F3}");

        Assert.True(mseRobust < mseLs, $"robust MSE {mseRobust:F5} should beat least-squares {mseLs:F5}");
    }

    [Fact]
    public void ContaminatedNormal_RejectsOutliers_BeatsLeastSquares()
    {
        var rng = new Xoshiro256PlusPlus(seed: 44);
        int n = 100;
        var x = new double[n];
        var f = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            f[i] = Math.Sin(2.0 * Math.PI * x[i]);
            y[i] = f[i] + 0.08 * Gaussian(rng);
        }
        foreach (int b in new[] { 8, 21, 37, 50, 63, 79, 88, 95 })
            y[b] = f[b] + 3.5;

        var basis = new SplineBasis(3);
        var model = new WeightedNormalModel();
        var prior = new PoissonPrior(4.0);
        var kernel = new LocalBetaKernel(50.0);

        BarsResult ls = new BarsEnsemble(basis, model, prior, kernel)
            .Run(x, y, grid: x, chains: 2, masterSeed: 9, burn: 1000, samples: 1500);
        BarsResult robust = new BarsEnsemble(basis, model, prior, kernel, new ContaminatedNormalWeights(0.1, 12.0))
            .Run(x, y, grid: x, chains: 2, masterSeed: 9, burn: 1000, samples: 1500);

        double mseLs = Mse(ls.Fit, f);
        double mseRobust = Mse(robust.Fit, f);
        _out.WriteLine($"[contam] mseLS={mseLs:F5} mseRobust={mseRobust:F5} ratio={mseRobust / mseLs:F3}");

        Assert.True(mseRobust < mseLs, $"contaminated-normal MSE {mseRobust:F5} should beat least-squares {mseLs:F5}");
    }

    private static double Mse(double[] a, double[] b)
    {
        double s = 0.0;
        for (int i = 0; i < a.Length; i++) { double d = a[i] - b[i]; s += d * d; }
        return s / a.Length;
    }
}
