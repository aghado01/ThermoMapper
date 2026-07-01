using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The step-function carrier: <see cref="StepBasis"/> produces segment-indicator designs, and the same engine
/// + moves + marginal + ensemble recover a true step (denoising) — the carrier for level/integer observables
/// over a swept parameter (b₁(T)) and true discontinuities.
/// </summary>
public sealed class StepBasisTests
{
    private readonly ITestOutputHelper _out;
    public StepBasisTests(ITestOutputHelper output) => _out = output;

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
    public void Design_And_Evaluate_AreSegmentIndicators()
    {
        var basis = new StepBasis();
        var config = new KnotConfig(new[] { 0.5 });   // one changepoint → two segments
        Assert.Equal(2, basis.Dimension(1));

        double[,] z = basis.Design(config, new[] { 0.2, 0.8 });
        Assert.Equal(1.0, z[0, 0]);
        Assert.Equal(0.0, z[0, 1]);
        Assert.Equal(0.0, z[1, 0]);
        Assert.Equal(1.0, z[1, 1]);

        Assert.Equal(3.0, basis.Evaluate(config, new[] { 3.0, 7.0 }, 0.2));
        Assert.Equal(7.0, basis.Evaluate(config, new[] { 3.0, 7.0 }, 0.8));
    }

    [Fact]
    public void StepBasis_RecoversStepFunction()
    {
        var rng = new Xoshiro256PlusPlus(seed: 21);
        int n = 120;
        var x = new double[n];
        var f = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            f[i] = x[i] < 0.5 ? 0.0 : 1.0;
            y[i] = f[i] + 0.1 * Gaussian(rng);
        }

        var config = new BarsConfig
        {
            Basis = new StepBasis(),
            Prior = new PoissonPrior(2.0),
            Chains = 3,
            MasterSeed = 5,
            BurnIn = 1000,
            MaxSamples = 1500,
        };
        BarsResult result = Bars.Run(config, x, y, x);

        double mseFit = Mse(result.Fit, f);
        double mseData = Mse(y, f);
        _out.WriteLine($"step mseFit={mseFit:F5} mseData={mseData:F5} meanK={result.MeanKnots:F2}");

        Assert.True(mseFit < 0.5 * mseData, $"step fit MSE {mseFit:F5} not below half data MSE {mseData:F5}");
        Assert.True(result.MeanKnots >= 0.5, $"should find ~1 changepoint, got {result.MeanKnots:F2}");
    }
}
