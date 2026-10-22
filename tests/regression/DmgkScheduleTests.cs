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
/// The DMGK birth/death schedule (prior-aware, state-dependent move-selection weights) recovers the curve through the
/// full ensemble. That it samples the same posterior as constant weights is checked exactly by the count-law oracle
/// (<see cref="CountLawTests"/>). Mixing is logged, not asserted. At equal c the schedule never raises a birth/death
/// flow (<see cref="DmgkSchedule"/>), and over 12 master seeds on this problem, before the R01 measure fix, DMGK had
/// the higher ESS(k) in 5, with single-seed ESS(k) from 6 to 223 (preparation brief, BARS mixing).
/// </summary>
public sealed class DmgkScheduleTests
{
    private readonly ITestOutputHelper _out;
    public DmgkScheduleTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double Mse(double[] a, double[] b)
    {
        double s = 0.0;
        for (int i = 0; i < a.Length; i++) { double d = a[i] - b[i]; s += d * d; }
        return s / a.Length;
    }

    [Fact]
    public void Dmgk_RecoversTheCurve()
    {
        var rng = new Xoshiro256PlusPlus(seed: 314);
        int n = 90;
        var x = new double[n];
        var f = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            f[i] = Math.Sin(3.0 * Math.PI * x[i]);   // ~1.5 cycles → a few interior knots
            y[i] = f[i] + 0.10 * Gaussian(rng);
        }

        var basis = new SplineBasis(3);
        var prior = new PoissonPrior(4.0);
        var kernel = new LocalBetaKernel(50.0);

        BarsResult equal = new BarsEnsemble(basis, new WeightedNormalModel(), prior, kernel)
            .Run(x, y, grid: x, chains: 4, masterSeed: 3, burn: 1500, samples: 3000);
        BarsResult dmgk = new BarsEnsemble(basis, new WeightedNormalModel(), prior, kernel, dmgkConstant: 0.4)
            .Run(x, y, grid: x, chains: 4, masterSeed: 3, burn: 1500, samples: 3000);

        double mseEqual = Mse(equal.Fit, f);
        double mseDmgk = Mse(dmgk.Fit, f);
        _out.WriteLine($"[dmgk] equal: mse={mseEqual:F5} k̄={equal.MeanKnots:F2} ESS(k)={equal.Knots.BulkEss:F0} acc={equal.AcceptanceRate:F3}");
        _out.WriteLine($"[dmgk] dmgk : mse={mseDmgk:F5} k̄={dmgk.MeanKnots:F2} ESS(k)={dmgk.Knots.BulkEss:F0} acc={dmgk.AcceptanceRate:F3}");

        Assert.True(mseDmgk < 0.01, $"DMGK fit MSE {mseDmgk:F5} should recover the curve");
    }

    [Fact]
    public void GreenConstant_IsTheLargestCWithinTheBudget()
    {
        // Poisson(1): the worst k is 1, with min(1, 1/2) + min(1, 1/1) = 1.5, so c = 0.9/1.5 = 0.6, above DMGK's old cap.
        Assert.Equal(0.6, DmgkSchedule.GreenConstant(new PoissonPrior(1.0), maxKnots: 40), 12);
        // Poisson(4): the worst k is 4, with min(1, 4/5) + min(1, 4/4) = 1.8, so c = 0.5.
        Assert.Equal(0.5, DmgkSchedule.GreenConstant(new PoissonPrior(4.0), maxKnots: 40), 12);
        // A cap at 5: interior k have 1 + 1 = 2, so c = 0.45; the cap itself has no birth.
        Assert.Equal(0.45, DmgkSchedule.GreenConstant(new UniformPrior(5), maxKnots: 40), 12);

        // At Green's c every reachable k keeps b_k + d_k within the budget, so the weights sum to one.
        var kernel = new LocalBetaKernel(50.0);
        var prior = new PoissonPrior(1.0);
        IRjMove<KnotConfig>[] moves = DmgkSchedule.Wrap(new KnotBirthMove(kernel), new KnotDeathMove(kernel), new KnotRelocateMove(kernel),
                                                        prior, DmgkSchedule.GreenConstant(prior, maxKnots: 40));
        for (int k = 0; k <= 40; k++)
        {
            var state = new KnotConfig(Enumerable.Range(0, k).Select(i => (i + 0.5) / (k + 1)).ToArray());
            double b = moves[0].Weight(state), d = moves[1].Weight(state), r = moves[2].Weight(state);
            Assert.True(b + d <= 0.9 + 1e-12, $"k = {k}: b + d = {b + d}");
            Assert.Equal(1.0, b + d + r, 12);
        }
    }
}
