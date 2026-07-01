using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The exact-peak readout: the closed-form argmax matches a fine grid scan (no optimizer slop), and the
/// chain-pooled peak posterior recovers a known peak location with a credible interval that brackets it —
/// the T_c-as-distribution payoff.
/// </summary>
public sealed class PeakTests
{
    private readonly ITestOutputHelper _out;
    public PeakTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void Extrema_ClosedForm_MatchesFineGridScan()
    {
        var basis = new SplineBasis(3);
        var config = new KnotConfig(new[] { 0.3, 0.55, 0.8 });
        var rng = new Xoshiro256PlusPlus(seed: 2);
        int nu = basis.Dimension(config.Count);
        var coef = new double[nu];
        for (int j = 0; j < nu; j++) coef[j] = Gaussian(rng);

        (double loc, double height) = SplineExtrema.Argmax(config, coef, basis);

        double bx = 0.0, bf = double.NegativeInfinity;
        const int n = 40001;
        for (int i = 0; i < n; i++)
        {
            double x = (double)i / (n - 1);
            double f = basis.Evaluate(config, coef, x);
            if (f > bf) { bf = f; bx = x; }
        }

        _out.WriteLine($"closed-form=({loc:F6},{height:F6}) grid=({bx:F6},{bf:F6})");
        Assert.Equal(bf, height, 4);
        Assert.True(Math.Abs(loc - bx) < 1e-3, $"closed-form loc {loc} vs grid {bx}");
    }

    [Theory]
    [InlineData(4)]   // quartic — beyond the 4-eval cubic fit; exercises the general reconstruction
    [InlineData(5)]   // quintic
    public void Extrema_MatchesFineGridScan_ForHigherDegreeSplines(int degree)
    {
        var basis = new SplineBasis(degree);
        var config = new KnotConfig(new[] { 0.3, 0.55, 0.8 });
        var rng = new Xoshiro256PlusPlus(seed: 5 + degree);
        int nu = basis.Dimension(config.Count);
        var coef = new double[nu];
        for (int j = 0; j < nu; j++) coef[j] = Gaussian(rng);

        (double loc, double height) = SplineExtrema.Argmax(config, coef, basis);

        double bx = 0.0, bf = double.NegativeInfinity;
        const int n = 40001;
        for (int i = 0; i < n; i++)
        {
            double x = (double)i / (n - 1);
            double f = basis.Evaluate(config, coef, x);
            if (f > bf) { bf = f; bx = x; }
        }

        _out.WriteLine($"deg={degree} closed-form=({loc:F6},{height:F6}) grid=({bx:F6},{bf:F6})");
        Assert.Equal(bf, height, 4);
        Assert.True(Math.Abs(loc - bx) < 1e-3, $"deg {degree}: closed-form loc {loc} vs grid {bx}");
    }

    [Fact]
    public void PeakPosterior_RecoversKnownPeak_WithBracketingInterval()
    {
        var rng = new Xoshiro256PlusPlus(seed: 7);
        int n = 150;
        const double truePeak = 0.4;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Exp(-120.0 * (x[i] - truePeak) * (x[i] - truePeak)) + 0.05 * Gaussian(rng);
        }

        var config = new BarsConfig { Prior = new PoissonPrior(5.0), Chains = 4, MasterSeed = 3, BurnIn = 1500, MaxSamples = 2000 };
        BarsResult result = Bars.Run(config, x, y, x);
        PeakPosterior p = result.Peak;

        _out.WriteLine($"peak={p.LocationMean:F3} [{p.LocationLo:F3},{p.LocationHi:F3}] " +
                       $"height={p.HeightMean:F3} R̂={p.LocationRHat:F3} ESS={p.LocationEss:F0}");

        Assert.True(Math.Abs(p.LocationMean - truePeak) < 0.06, $"peak mean {p.LocationMean:F3} off {truePeak}");
        Assert.True(p.LocationLo <= truePeak && truePeak <= p.LocationHi, "95% interval should bracket the true peak");
        Assert.True(p.LocationLo < p.LocationHi);
    }
}
