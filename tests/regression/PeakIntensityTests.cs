using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The peak-intensity readout λ(T): a two-bump curve yields a grid-aligned intensity with mass at BOTH
/// transitions (not just the global one); its total equals the peak-count mean by construction, and each
/// bump carries ≈ one transition's worth of mass — the matching-free multi-peak posterior (MP-4).
/// </summary>
public sealed class PeakIntensityTests
{
    private readonly ITestOutputHelper _out;
    public PeakIntensityTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void TwoBumps_IntensityResolvesBothTransitions()
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
        };

        BarsResult r = Bars.Run(config, x, y, x);

        double total = 0.0, left = 0.0, right = 0.0;
        for (int i = 0; i < n; i++)
        {
            total += r.PeakIntensity[i];
            if (x[i] >= 0.2 && x[i] <= 0.4) left += r.PeakIntensity[i];
            if (x[i] >= 0.6 && x[i] <= 0.8) right += r.PeakIntensity[i];
        }
        _out.WriteLine($"λ total={total:F3} (count={r.PeakCountMean:F3}) leftMass={left:F2} rightMass={right:F2}");

        // Invariant: ∑ λ == PeakCountMean (the same significant-peak set feeds both reductions).
        Assert.Equal(r.PeakCountMean, total, 6);
        // Both transitions are resolved, not just the global one: each window carries ≈ one peak's worth.
        Assert.True(left > 0.5, $"left transition mass {left:F2} should be ≳ 1");
        Assert.True(right > 0.5, $"right transition mass {right:F2} should be ≳ 1");
    }
}
