using System;
using System.Linq;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// A BARS run surfaces per-move acceptance telemetry (birth / death / relocate) in <c>BarsResult.MoveStats</c>,
/// banked across robust-path chain rebuilds — the "which move is limiting mixing" view the DMGK / adaptive-τ work
/// makes actionable.
/// </summary>
public sealed class MoveTelemetryTests
{
    private readonly ITestOutputHelper _out;
    public MoveTelemetryTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void Ensemble_SurfacesPerMoveAcceptance()
    {
        var noise = new Xoshiro256PlusPlus(seed: 19);
        int n = 70;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++) { x[i] = (i + 0.5) / n; y[i] = Math.Sin(3.0 * Math.PI * x[i]) + 0.1 * Gaussian(noise); }

        BarsResult r = new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(4.0), new LocalBetaKernel(50.0))
            .Run(x, y, grid: x, chains: 3, masterSeed: 2, burn: 500, samples: 1000);

        Assert.NotEmpty(r.MoveStats);
        foreach (MoveStat m in r.MoveStats)
        {
            Assert.True(m.Attempts >= 0 && m.Accepted >= 0 && m.Accepted <= m.Attempts);
            double rate = m.Attempts > 0 ? (double)m.Accepted / m.Attempts : 0.0;
            _out.WriteLine($"[move] {m.Key,-9} attempts={m.Attempts,7} acc={rate:F3}");
            Assert.InRange(rate, 0.0, 1.0);
        }

        // The standard palette is present and exercised.
        foreach (string key in new[] { "birth", "death", "relocate" })
        {
            MoveStat m = r.MoveStats.First(s => s.Key == key);
            Assert.True(m.Attempts > 0, $"move '{key}' should have been attempted");
        }
    }
}
