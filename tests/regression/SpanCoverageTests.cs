using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The span-coverage readout π(T): a two-bump curve yields coverage plateaus around each transition (the
/// structural FWHM-esque spans) with a dip in the stable regime between them — the matching-free interval
/// readout that would aim an SPC scheduler's second pass.
/// </summary>
public sealed class SpanCoverageTests
{
    private readonly ITestOutputHelper _out;
    public SpanCoverageTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static int Nearest(double[] grid, double t)
    {
        int best = 0;
        double bd = double.PositiveInfinity;
        for (int i = 0; i < grid.Length; i++) { double d = Math.Abs(grid[i] - t); if (d < bd) { bd = d; best = i; } }
        return best;
    }

    [Fact]
    public void TwoBumps_CoverageCoversTransitionsAndDipsInRegime()
    {
        var rng = new Xoshiro256PlusPlus(11);
        int n = 150;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Exp(-150.0 * (x[i] - 0.3) * (x[i] - 0.3))
                 + Math.Exp(-150.0 * (x[i] - 0.7) * (x[i] - 0.7))
                 + 0.05 * Gaussian(rng);
        }
        var config = new BarsConfig
        {
            Prior = new PoissonPrior(5.0),
            Chains = 4,
            MasterSeed = 3,
            BurnIn = 1500,
            MaxSamples = 2000,
            PeakProminence = 0.15,
            SpanDropFraction = 0.5,
        };

        BarsResult r = Bars.Run(config, x, y, x);

        double left = r.SpanCoverage[Nearest(x, 0.3)];
        double right = r.SpanCoverage[Nearest(x, 0.7)];
        double valley = r.SpanCoverage[Nearest(x, 0.5)];
        _out.WriteLine($"π(0.3)={left:F2} π(0.5)={valley:F2} π(0.7)={right:F2}");

        // Both transition spans are covered with high posterior probability...
        Assert.True(left > 0.5, $"left transition coverage {left:F2} should be high");
        Assert.True(right > 0.5, $"right transition coverage {right:F2} should be high");
        // ...and the stable regime between them is a coverage dip (the spans don't reach the valley centre).
        Assert.True(valley < left && valley < right, $"valley coverage {valley:F2} should dip below the peaks");

        // π is a coverage probability in [0,1] for well-separated peaks (no within-draw span overlap).
        foreach (double pi in r.SpanCoverage) Assert.InRange(pi, 0.0, 1.0 + 1e-9);
    }
}
