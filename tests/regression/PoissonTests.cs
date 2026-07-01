using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The Poisson observation model: fit count data from a known log-intensity bump and recover the intensity
/// peak — the spike-train / density case, the same engine with a Poisson likelihood instead of Normal.
/// </summary>
public sealed class PoissonTests
{
    private readonly ITestOutputHelper _out;
    public PoissonTests(ITestOutputHelper output) => _out = output;

    // Knuth's Poisson sampler (adequate for the moderate λ here).
    private static int PoissonSample(Xoshiro256PlusPlus rng, double lambda)
    {
        double l = Math.Exp(-lambda), p = 1.0;
        int k = 0;
        do { k++; p *= rng.NextDouble(); } while (p > l);
        return k - 1;
    }

    [Fact]
    public void PoissonModel_RecoversIntensityPeak()
    {
        var rng = new Xoshiro256PlusPlus(seed: 33);
        int n = 150;
        const double truePeak = 0.4;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            double logIntensity = 1.2 + 1.3 * Math.Exp(-60.0 * (x[i] - truePeak) * (x[i] - truePeak));
            y[i] = PoissonSample(rng, Math.Exp(logIntensity));
        }

        var config = new BarsConfig
        {
            Model = new PoissonModel(),
            Prior = new PoissonPrior(5.0),
            Chains = 3,
            MasterSeed = 3,
            BurnIn = 1000,
            MaxSamples = 1500,
        };
        BarsResult result = Bars.Run(config, x, y, x);

        _out.WriteLine($"poisson peak={result.Peak.LocationMean:F3} " +
                       $"[{result.Peak.LocationLo:F3},{result.Peak.LocationHi:F3}] meanK={result.MeanKnots:F2}");

        Assert.InRange(result.Peak.LocationMean, truePeak - 0.1, truePeak + 0.1);
    }
}
