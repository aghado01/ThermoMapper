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

        // π is a coverage probability in [0,1].
        foreach (double pi in r.SpanCoverage) Assert.InRange(pi, 0.0, 1.0 + 1e-9);
    }

    [Fact]
    public void SpanCoverage_CountsTheUnionOfADrawsSpans()
    {
        // Nested [0.35, 0.45] ⊂ [0.30, 0.60], plus a disjoint [0.80, 0.90]: each covered grid point gets exactly 1.
        PeakSpan[] spans =
        [
            new(0.50, 1.0, 0.30, 0.60, false, false),
            new(0.40, 0.9, 0.35, 0.45, false, false),
            new(0.85, 0.5, 0.80, 0.90, false, false),
        ];
        var grid = new double[21];
        for (int i = 0; i < grid.Length; i++) grid[i] = i * 0.05;
        var sum = new double[grid.Length];

        BarsEnsemble.AddSpanCoverage(spans, grid, sum);

        for (int i = 0; i < grid.Length; i++)
        {
            double t = grid[i];
            bool covered = (t >= 0.30 && t <= 0.60) || (t >= 0.80 && t <= 0.90);
            Assert.Equal(covered ? 1.0 : 0.0, sum[i]);
        }
    }

    [Fact]
    public void NestedPeaks_CoverageStaysAProbability()
    {
        // R14: two peaks joined by a high col. The higher peak (0.54) has its half-prominence level below the col, so
        // its span runs over the col and contains the lower peak's span: a draw covers the region around 0.4 twice.
        // π(T) is P(T lies in some significant span), so each draw counts a grid point at most once.
        var rng = new Xoshiro256PlusPlus(12);
        int n = 200;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Exp(-0.5 * Math.Pow((x[i] - 0.40) / 0.08, 2))
                 + 0.8 * Math.Exp(-0.5 * Math.Pow((x[i] - 0.54) / 0.03, 2))
                 + 0.02 * Gaussian(rng);
        }
        var config = new BarsConfig
        {
            Prior = new PoissonPrior(6.0),
            Chains = 2,
            MasterSeed = 4,
            BurnIn = 1500,
            MaxSamples = 1500,
            PeakProminence = 0.1,
            SpanDropFraction = 0.5,
        };

        BarsResult r = Bars.Run(config, x, y, x);

        double inner = r.SpanCoverage[Nearest(x, 0.40)];
        double outer = r.SpanCoverage[Nearest(x, 0.54)];
        double max = 0.0;
        foreach (double pi in r.SpanCoverage) max = Math.Max(max, pi);
        _out.WriteLine($"π(0.40)={inner:F3} π(0.54)={outer:F3} max π={max:F3} peaks/draw={r.PeakCountMean:F2}");
        foreach (PeakMode m in r.PeakModes)
            _out.WriteLine($"mode {m.Location:F3} [{m.LocationLo:F3}, {m.LocationHi:F3}] span [{m.SpanLeft:F3}, {m.SpanRight:F3}] mass {m.Mass:F2}");

        Assert.True(r.PeakCountMean > 1.5, $"both peaks should be significant in most draws ({r.PeakCountMean:F2})");
        Assert.True(inner > 0.9 && outer > 0.9, $"both peaks should be covered: π(0.40)={inner:F3}, π(0.54)={outer:F3}");
        foreach (double pi in r.SpanCoverage) Assert.InRange(pi, 0.0, 1.0 + 1e-9);
    }
}
