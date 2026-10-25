using System;
using System.Linq;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// BARS's chains start on the posterior's support. An overdispersed start can hold more knots than a few observations
/// admit; a chain started there leaves only by a proposal onto the support, and where every neighbour is off it too it
/// never moves, so its "draws" were fits at configurations the posterior gives no mass, NaN where the design was rank
/// deficient. Such a start now sheds knots until it is admissible. When even the zero-knot design is refused, the run
/// fails with the reason instead of returning a NaN curve with peak modes.
/// </summary>
public sealed class BarsStartTests
{
    private static (double[] X, double[] Y) Curve(int n)
    {
        double[] x = Enumerable.Range(0, n).Select(i => i / (n - 1.0)).ToArray();
        return (x, x.Select(v => Math.Sin(3.0 * v) + 0.1 * Math.Cos(17.0 * v)).ToArray());
    }

    private static BarsEnsemble Bars() =>
        new(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(4.0), new LocalBetaKernel(50.0));

    private static readonly double[] Grid = Enumerable.Range(0, 11).Select(i => i / 10.0).ToArray();

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    public void NoAdmissibleConfiguration_IsRefused_WithTheReason(int n)
    {
        (double[] x, double[] y) = Curve(n);
        var e = Assert.Throws<InvalidOperationException>(() => Bars().Run(x, y, Grid, chains: 4, masterSeed: 1, burn: 50, samples: 50));
        Assert.Contains($"{n} observations", e.Message);
        Assert.Contains(nameof(DesignRejection.ThinSupport), e.Message);
    }

    [Theory]
    [InlineData(6)]
    [InlineData(9)]
    public void StartsOffTheSupport_MoveOntoIt(int n)
    {
        (double[] x, double[] y) = Curve(n);
        // Four chains with up to four starting knots: on so few points the larger starts are off the support.
        BarsResult r = Bars().Run(x, y, Grid, chains: 4, masterSeed: 1, burn: 200, samples: 200, startDispersion: 4);
        Assert.All(r.Draws.SelectMany(c => c), d =>
        {
            Assert.True(double.IsFinite(d.Sigma2), $"σ² {d.Sigma2} at {d.Knots.Count} knots");
            Assert.All(d.Coefficients, b => Assert.True(double.IsFinite(b)));
        });
        Assert.All(r.Fit, f => Assert.True(double.IsFinite(f)));
        Assert.All(r.FitSd, s => Assert.True(double.IsFinite(s)));
        Assert.True(r.AcceptanceRate > 0.1, $"acceptance {r.AcceptanceRate:F3}");
    }
}
