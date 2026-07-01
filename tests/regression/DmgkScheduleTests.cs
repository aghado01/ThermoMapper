using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The DMGK birth/death schedule (prior-aware, state-dependent move-selection weights) must sample the same
/// posterior as equal weights — that it does is the proof the Green-general selection ratio is unbiased — while
/// mixing the dimension at least as well, since folding the prior into the proposal lets it cancel in the
/// acceptance so dimension moves turn on likelihood evidence.
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
    public void Dmgk_SamplesSamePosterior_AndMixesDimension()
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
        _out.WriteLine($"[dmgk] equal: mse={mseEqual:F5} k̄={equal.MeanKnots:F2} ESS(k)={equal.EssKnots:F0} acc={equal.AcceptanceRate:F3}");
        _out.WriteLine($"[dmgk] dmgk : mse={mseDmgk:F5} k̄={dmgk.MeanKnots:F2} ESS(k)={dmgk.EssKnots:F0} acc={dmgk.AcceptanceRate:F3}");

        // Correctness: both recover the curve, and DMGK does not bias the dimension posterior.
        Assert.True(mseDmgk < 0.01, $"DMGK fit MSE {mseDmgk:F5} should recover the curve");
        Assert.True(Math.Abs(dmgk.MeanKnots - equal.MeanKnots) < 0.6,
            $"DMGK mean knots {dmgk.MeanKnots:F2} should match equal-weight {equal.MeanKnots:F2} (unbiased)");

        // Mixing: DMGK mixes the knot count better — folding the prior into the proposal lets it cancel in the
        // acceptance, so dimension moves turn on likelihood evidence (here ~8× the ESS of equal weights).
        Assert.True(dmgk.EssKnots > equal.EssKnots,
            $"DMGK ESS(k) {dmgk.EssKnots:F0} should exceed equal-weight {equal.EssKnots:F0}");
    }
}
