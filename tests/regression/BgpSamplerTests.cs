using System;
using Maths.Regression.Bgp;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The full Bayesian-GP instrument (Tang, Wu, Cheng &amp; Dunson 2025) end-to-end on the Swiss Roll — a 2-D manifold
/// embedded in R³. The empirical-Bayes bandwidth prior + Metropolis sampler recover a smooth function on the
/// manifold while never being told the intrinsic dimension d=2; that the estimator lands at a length scale giving
/// low out-of-sample error is the dimension-adaptivity payoff. The manifold dimension is an independent ground
/// truth — fixed by construction, never seen by the prior.
/// </summary>
public sealed class BgpSamplerTests
{
    private readonly ITestOutputHelper _out;
    public BgpSamplerTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    // Swiss roll: a 2-D manifold (unroll coord × height) embedded in R³. Returns raw ambient coords + the
    // intrinsic-smooth target f*.
    private static (double[,] X, double[] F) SwissRoll(int n, Xoshiro256PlusPlus rng)
    {
        var x = new double[n, 3];
        var f = new double[n];
        for (int i = 0; i < n; i++)
        {
            double u = rng.NextDouble();             // unroll coordinate
            double v = rng.NextDouble();             // height coordinate
            double phi = 1.5 * Math.PI * (1.0 + 2.0 * u);
            x[i, 0] = phi * Math.Cos(phi);
            x[i, 1] = phi * Math.Sin(phi);
            x[i, 2] = 12.0 * v;
            f[i] = Math.Sin(2.0 * Math.PI * u) + 0.5 * Math.Cos(2.0 * Math.PI * v);   // smooth on the manifold
        }
        return (x, f);
    }

    // Min-max normalize each ambient coordinate into [0,1] using the supplied bounds (paper assumes X ⊂ [0,1]^D).
    private static void Normalize(double[,] x, double[] lo, double[] span)
    {
        int n = x.GetLength(0), d = x.GetLength(1);
        for (int i = 0; i < n; i++)
            for (int k = 0; k < d; k++) x[i, k] = (x[i, k] - lo[k]) / span[k];
    }

    [Fact]
    public void Bgp_RecoversManifoldFunction_WithoutKnowingDimension()
    {
        var rng = new Xoshiro256PlusPlus(seed: 42);
        const int nTrain = 200, nTest = 300;
        (double[,] xTr, double[] fTr) = SwissRoll(nTrain, rng);
        (double[,] xTe, double[] fTe) = SwissRoll(nTest, rng);

        // Shared [0,1]³ frame from the training bounds.
        var lo = new double[3]; var span = new double[3];
        for (int k = 0; k < 3; k++)
        {
            double mn = double.MaxValue, mx = double.MinValue;
            for (int i = 0; i < nTrain; i++) { mn = Math.Min(mn, xTr[i, k]); mx = Math.Max(mx, xTr[i, k]); }
            lo[k] = mn; span[k] = mx - mn;
        }
        Normalize(xTr, lo, span);
        Normalize(xTe, lo, span);

        const double sigma = 0.1;
        var y = new double[nTrain];
        for (int i = 0; i < nTrain; i++) y[i] = fTr[i] + sigma * Gaussian(rng);

        var gp = new GpRegression(xTr, y, sigma * sigma, new SquaredExponentialKernel());
        var prior = new EmpiricalBayesBandwidthPrior(xTr, gp);          // d=2 never supplied
        var sampler = new BgpSampler(gp, prior);
        BgpResult r = sampler.Run(xTe, draws: 300, burn: 200, proposalSd: 0.5, seed: 1);

        double mse = 0.0;
        for (int i = 0; i < nTest; i++) { double e = r.PosteriorMean[i] - fTe[i]; mse += e * e; }
        mse /= nTest;

        double varF = 0.0, meanF = 0.0;
        for (int i = 0; i < nTest; i++) meanF += fTe[i];
        meanF /= nTest;
        for (int i = 0; i < nTest; i++) varF += (fTe[i] - meanF) * (fTe[i] - meanF);
        varF /= nTest;

        _out.WriteLine($"[bgp-swiss] t̂={r.BandwidthMean:F3} [{r.BandwidthLo:F3},{r.BandwidthHi:F3}] " +
                       $"accept={r.AcceptanceRate:F3} test MSE={mse:F4} (Var f*={varF:F3})");

        // Recovers the manifold function: out-of-sample error far below the function's own variance.
        Assert.True(mse < 0.1 * varF, $"test MSE {mse:F4} should be ≪ Var f* {varF:F3}");
        Assert.InRange(r.AcceptanceRate, 0.1, 0.85);                    // healthy RW mixing
        Assert.True(r.BandwidthMean > prior.LowerBound && r.BandwidthMean <= 1.0, "t̂ inside the prior band");
        Assert.True(r.BandwidthLo < r.BandwidthHi, "non-degenerate bandwidth posterior");
    }
}
