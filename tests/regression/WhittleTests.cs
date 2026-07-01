using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The Whittle spectral model: fit an exponential periodogram from a known log-spectrum with a peak and
/// recover the dominant frequency — the spectral-density-estimation case on the same engine.
/// </summary>
public sealed class WhittleTests
{
    private readonly ITestOutputHelper _out;
    public WhittleTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void WhittleModel_RecoversSpectralPeak()
    {
        var rng = new Xoshiro256PlusPlus(seed: 55);
        int m = 200;
        const double peakFreq = 0.3;
        var omega = new double[m];
        var periodogram = new double[m];
        for (int i = 0; i < m; i++)
        {
            omega[i] = (i + 0.5) / m;
            double logSpectrum = -0.5 + 1.6 * Math.Exp(-80.0 * (omega[i] - peakFreq) * (omega[i] - peakFreq));
            double f = Math.Exp(logSpectrum);
            double u = 1.0 - rng.NextDouble();
            periodogram[i] = -f * Math.Log(u);   // Exponential(mean f)
        }

        var config = new BarsConfig
        {
            Model = new WhittleModel(),
            Prior = new PoissonPrior(5.0),
            Chains = 3,
            MasterSeed = 3,
            BurnIn = 1000,
            MaxSamples = 1500,
        };
        BarsResult result = Bars.Run(config, omega, periodogram, omega);

        _out.WriteLine($"whittle peak={result.Peak.LocationMean:F3} " +
                       $"[{result.Peak.LocationLo:F3},{result.Peak.LocationHi:F3}] meanK={result.MeanKnots:F2}");

        Assert.InRange(result.Peak.LocationMean, peakFreq - 0.1, peakFreq + 0.1);
    }
}
