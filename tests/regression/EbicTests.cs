using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The EBIC complexity prior: the γ knob is flat at 0 and penalizes growing model space at 1 (pinned against
/// the closed form), and a γ=1 fit selects fewer knots than γ=0 on noisy over-knotting-prone data.
/// </summary>
public sealed class EbicTests
{
    private readonly ITestOutputHelper _out;
    public EbicTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void GammaZero_IsFlat_GammaOne_PenalizesGrowingK()
    {
        var flat = new EbicPrior(gamma: 0.0, candidateCount: 200);
        Assert.Equal(0.0, flat.LogPrior(0), 12);
        Assert.Equal(0.0, flat.LogPrior(5), 12);
        Assert.Equal(0.0, flat.LogPrior(20), 12);

        var penalised = new EbicPrior(gamma: 1.0, candidateCount: 200);
        Assert.True(penalised.LogPrior(2) > penalised.LogPrior(5));    // −log C(n,k) strictly decreasing for k « n
        Assert.True(penalised.LogPrior(5) > penalised.LogPrior(10));
        Assert.Equal(double.NegativeInfinity, penalised.LogPrior(201));
    }

    [Fact]
    public void HigherGamma_SelectsFewerKnots()
    {
        (BarsResult flat, BarsResult strict) = FlatAndStrict(noiseSeed: 17, masterSeed: 3);

        _out.WriteLine($"meanK γ=0: {flat.MeanKnots:F2}, γ=1: {strict.MeanKnots:F2}");
        Assert.True(strict.MeanKnots < flat.MeanKnots,
            $"γ=1 mean knots {strict.MeanKnots:F2} should be below γ=0 {flat.MeanKnots:F2}");
    }

    // The test's computation at any seeds, which the seed sweep reruns (ComparativeSweep).
    internal static (BarsResult Flat, BarsResult Strict) FlatAndStrict(int noiseSeed, int masterSeed)
    {
        var rng = new Xoshiro256PlusPlus(seed: noiseSeed);
        int n = 90;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.25 * Gaussian(rng);   // noisy → over-knotting tempting
        }

        var basis = new SplineBasis(degree: 3);
        var model = new WeightedNormalModel();
        var kernel = new LocalBetaKernel(50.0);

        // Burn-in 2,000. Births land near existing knots, so a chain started with knots far apart can hold one for
        // about 1,500 steps although deaths would raise the posterior; at 800, γ = 1 sometimes kept more knots than
        // γ = 0 (4 of 100 runs in the step 3.5 seed sweep).
        const int burn = 2000;
        BarsResult flat = new BarsEnsemble(basis, model, new EbicPrior(0.0, 200), kernel)
            .Run(x, y, grid: x, chains: 2, masterSeed: masterSeed, burn: burn, samples: 1200);
        BarsResult strict = new BarsEnsemble(basis, model, new EbicPrior(1.0, 200), kernel)
            .Run(x, y, grid: x, chains: 2, masterSeed: masterSeed, burn: burn, samples: 1200);
        return (flat, strict);
    }
}
