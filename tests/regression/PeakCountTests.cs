using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The significant-peak count functional: a two-bump curve yields a peak-count posterior near 2 and a single
/// bump near 1, with prominence filtering the noise wiggles — the "how many transitions" readout for a
/// multimodal SPC response.
/// </summary>
public sealed class PeakCountTests
{
    private readonly ITestOutputHelper _out;
    public PeakCountTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double PeakCount(Func<double, double> truth, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        int n = 150;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = truth(x[i]) + 0.05 * Gaussian(rng);
        }
        var config = new BarsConfig
        {
            Prior = new PoissonPrior(5.0),
            Chains = 4,
            MasterSeed = 3,
            BurnIn = 1500,
            MaxSamples = 2000,
            PeakProminence = 0.15,
        };
        return Bars.Run(config, x, y, x).PeakCountMean;
    }

    [Fact]
    public void TwoBumps_PeakCountNearTwo()
    {
        double count = PeakCount(t => Math.Exp(-150.0 * (t - 0.3) * (t - 0.3))
                                    + Math.Exp(-150.0 * (t - 0.7) * (t - 0.7)), seed: 11);
        _out.WriteLine($"two-bump peakCountMean={count:F2}");
        Assert.InRange(count, 1.5, 2.5);
    }

    [Fact]
    public void SingleBump_PeakCountNearOne()
    {
        double count = PeakCount(t => Math.Exp(-150.0 * (t - 0.5) * (t - 0.5)), seed: 12);
        _out.WriteLine($"one-bump peakCountMean={count:F2}");
        Assert.InRange(count, 0.7, 1.3);
    }
}
