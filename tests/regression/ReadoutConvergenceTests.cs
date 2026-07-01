using System;
using System.Linq;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The per-mode readout convergence: <see cref="BarsResult.PeakModes"/> resolves each transition from the pooled
/// λ/π fields, and on a single-peak curve its dominant mode reconciles with <see cref="BarsResult.Peak"/> — the
/// global <see cref="PeakPosterior"/> as the k = 1 special case, not a second source of truth. On a two-peak
/// curve it generalizes to both transitions (where the single global peak cannot).
/// </summary>
public sealed class ReadoutConvergenceTests
{
    private readonly ITestOutputHelper _out;
    public ReadoutConvergenceTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static BarsResult Fit(Func<double, double> truth, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        int n = 150;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++) { x[i] = (i + 0.5) / n; y[i] = truth(x[i]) + 0.05 * Gaussian(rng); }
        var config = new BarsConfig
        {
            Prior = new PoissonPrior(5.0),
            Chains = 4, MasterSeed = 3, BurnIn = 1500, MaxSamples = 2000, PeakProminence = 0.15,
        };
        return Bars.Run(config, x, y, x);
    }

    [Fact]
    public void SingleBump_DominantModeReconcilesWithPeakPosterior()
    {
        BarsResult r = Fit(t => Math.Exp(-150.0 * (t - 0.5) * (t - 0.5)), seed: 12);
        Assert.NotEmpty(r.PeakModes);
        PeakMode m = r.PeakModes[0];   // dominant
        _out.WriteLine($"mode loc={m.Location:F3} [{m.LocationLo:F3},{m.LocationHi:F3}] span=[{m.SpanLeft:F3},{m.SpanRight:F3}] R̂={m.LocationRHat:F3} | peak={r.Peak.LocationMean:F3}");

        // k = 1: the dominant mode IS the global peak — same location, and its credible interval brackets it.
        Assert.Equal(r.Peak.LocationMean, m.Location, 1);
        Assert.InRange(r.Peak.LocationMean, m.LocationLo - 1e-9, m.LocationHi + 1e-9);
        // The mode carries a real structural span around its location.
        Assert.True(m.SpanLeft < m.Location && m.Location < m.SpanRight, "dominant mode should have a span");
        // ...and a cross-chain convergence diagnostic (matching-free): ≈1 for a well-mixed single peak.
        Assert.True(!double.IsNaN(m.LocationRHat) && m.LocationRHat < 1.5, $"dominant mode R̂={m.LocationRHat:F3} should be ≈1");
    }

    [Fact]
    public void TwoBumps_ResolvesBothModesTheGlobalPeakCannot()
    {
        BarsResult r = Fit(t => Math.Exp(-150.0 * (t - 0.3) * (t - 0.3))
                              + Math.Exp(-150.0 * (t - 0.7) * (t - 0.7)), seed: 11);
        _out.WriteLine($"modes: {string.Join(", ", r.PeakModes.Select(m => $"{m.Location:F2}(m={m.Mass:F2})"))}");

        Assert.True(r.PeakModes.Count >= 2, $"expected ≥2 modes, got {r.PeakModes.Count}");
        Assert.Contains(r.PeakModes, m => Math.Abs(m.Location - 0.3) < 0.06);
        Assert.Contains(r.PeakModes, m => Math.Abs(m.Location - 0.7) < 0.06);
    }
}
