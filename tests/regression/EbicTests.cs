using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The EBIC count law in He, Yang and Kang's convention (step 3.8): log p(k) = (1 − γ)·log C(C, k), pinned against the
/// closed form. γ = 0 is BIC's uniform prior over models, γ = 1 is flat, γ = 2 penalizes; and a γ = 2 fit selects fewer
/// knots than a γ = 1 fit on noisy over-knotting-prone data.
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
    public void GammaOne_IsFlat_GammaTwo_PenalizesGrowingK()
    {
        var flat = new EbicPrior(gamma: 1.0, candidateCount: 200);
        Assert.Equal(0.0, flat.LogPrior(0), 12);
        Assert.Equal(0.0, flat.LogPrior(5), 12);
        Assert.Equal(0.0, flat.LogPrior(20), 12);

        var penalised = new EbicPrior(gamma: 2.0, candidateCount: 200);
        Assert.True(penalised.LogPrior(2) > penalised.LogPrior(5));    // −log C(C,k) strictly decreasing for k « C
        Assert.True(penalised.LogPrior(5) > penalised.LogPrior(10));
        Assert.Equal(double.NegativeInfinity, penalised.LogPrior(201));
        Assert.Equal(double.NegativeInfinity, penalised.LogPrior(-1));
    }

    [Fact]
    public void Law_IsOneMinusGamma_TimesTheLogModelCount()
    {
        // log C(10, 3) = log 120.
        Assert.Equal(Math.Log(120.0), new EbicPrior(0.0, 10).LogPrior(3), 12);           // BIC's uniform prior over models
        Assert.Equal(0.5 * Math.Log(120.0), new EbicPrior(0.5, 10).LogPrior(3), 12);
        Assert.Equal(-Math.Log(120.0), new EbicPrior(2.0, 10).LogPrior(3), 12);
        // Below γ = 1 the law peaks at C/2; above it, the law is symmetric and smallest there.
        var favoring = new EbicPrior(0.5, 10);
        Assert.True(favoring.LogPrior(5) > favoring.LogPrior(4) && favoring.LogPrior(5) > favoring.LogPrior(6));
        var penalizing = new EbicPrior(2.0, 10);
        Assert.Equal(penalizing.LogPrior(2), penalizing.LogPrior(8), 12);
        Assert.Equal(0.0, penalizing.LogPrior(10), 12);   // it rises again to the far end, which C ≥ 2·k_max keeps away
        Assert.Throws<ArgumentOutOfRangeException>(() => new EbicPrior(-0.1, 10));
        Assert.Throws<ArgumentOutOfRangeException>(() => new EbicPrior(double.NaN, 10));
        _ = new EbicPrior(2.5, 10);   // past the authors' recommended [0.5, 2], still a valid law
    }

    [Fact]
    public void HigherGamma_SelectsFewerKnots()
    {
        (BarsResult flat, BarsResult strict) = FlatAndStrict(noiseSeed: 17, masterSeed: 3);

        _out.WriteLine($"meanK γ=1: {flat.MeanKnots:F2}, γ=2: {strict.MeanKnots:F2}");
        Assert.True(strict.MeanKnots < flat.MeanKnots,
            $"γ=2 mean knots {strict.MeanKnots:F2} should be below γ=1 {flat.MeanKnots:F2}");
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

        const int burn = 800;
        BarsResult flat = new BarsEnsemble(basis, model, new EbicPrior(1.0, 200), kernel)
            .Run(x, y, grid: x, chains: 2, masterSeed: masterSeed, burn: burn, samples: 1200);
        BarsResult strict = new BarsEnsemble(basis, model, new EbicPrior(2.0, 200), kernel)
            .Run(x, y, grid: x, chains: 2, masterSeed: masterSeed, burn: burn, samples: 1200);
        return (flat, strict);
    }
}
