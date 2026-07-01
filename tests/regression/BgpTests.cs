using System;
using Maths.Regression.Bgp;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The Bayesian-GP core (Tang, Wu, Cheng &amp; Dunson 2025): the kernel-affinity statistic v̂_n(t) scales as
/// t^{d/2} in the intrinsic dimension d (the mechanism the empirical-Bayes bandwidth prior exploits — it reads d
/// off this slope without estimating it), and the dense conjugate fit recovers a smooth function. The affinity
/// test is an independent mechanism check: the manifold dimension is fixed by construction and the statistic
/// never sees it.
/// </summary>
public sealed class BgpTests
{
    private readonly ITestOutputHelper _out;
    public BgpTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    // Least-squares slope of log v̂_n(t) vs log t over a geometric t-grid — the empirical exponent.
    private static double AffinityLogLogSlope(GpRegression gp, double tLo, double tHi, int points)
    {
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (int i = 0; i < points; i++)
        {
            double t = tLo * Math.Pow(tHi / tLo, i / (double)(points - 1));
            double x = Math.Log(t), y = Math.Log(gp.KernelAffinity(t));
            sx += x; sy += y; sxx += x * x; sxy += x * y;
        }
        return (points * sxy - sx * sy) / (points * sxx - sx * sx);
    }

    [Theory]
    [InlineData(1, 0.5)]   // points on a line embedded in R² ⇒ v̂_n ∼ t^{1/2}
    [InlineData(2, 1.0)]   // points on a square embedded in R³ ⇒ v̂_n ∼ t^{1}
    public void KernelAffinity_ScalesAsIntrinsicDimension(int d, double expectedSlope)
    {
        var rng = new Xoshiro256PlusPlus(seed: 20 + d);
        const int n = 500, ambient = 3;
        var x = new double[n, ambient];
        for (int i = 0; i < n; i++)
            for (int k = 0; k < d; k++) x[i, k] = rng.NextDouble();   // remaining ambient coords stay 0

        var gp = new GpRegression(x, new double[n], sigma2: 1.0, new SquaredExponentialKernel());
        double slope = AffinityLogLogSlope(gp, tLo: 0.004, tHi: 0.04, points: 12);

        _out.WriteLine($"[bgp-affinity] d={d}: v̂_n ∼ t^slope, slope={slope:F3} (expected {expectedSlope:F2} = d/2)");
        Assert.True(Math.Abs(slope - expectedSlope) < 0.2, $"affinity slope {slope:F3} should track d/2 = {expectedSlope:F2}");
    }

    [Fact]
    public void Fit_RecoversSmoothFunction_AndPeaksMarginalAtSensibleBandwidth()
    {
        var rng = new Xoshiro256PlusPlus(seed: 7);
        const int n = 150;
        var x = new double[n, 2];
        var f = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = rng.NextDouble();
            x[i, 1] = rng.NextDouble();
            f[i] = Math.Exp(-5.0 * ((x[i, 0] - 0.5) * (x[i, 0] - 0.5) + (x[i, 1] - 0.5) * (x[i, 1] - 0.5)));
            y[i] = f[i] + 0.05 * Gaussian(rng);
        }

        var gp = new GpRegression(x, y, sigma2: 0.05 * 0.05, new SquaredExponentialKernel());

        // Pick the bandwidth by maximizing the marginal evidence over a grid (what the sampler will do continuously).
        GpFit best = gp.Fit(0.05);
        for (int i = 0; i < 16; i++)
        {
            double t = 0.01 * Math.Pow(50.0, i / 15.0);   // 0.01 … 0.5
            GpFit fit = gp.Fit(t);
            if (fit.LogMarginal > best.LogMarginal) best = fit;
        }

        double[] fhat = gp.PredictMean(best, x);
        double mse = 0.0;
        for (int i = 0; i < n; i++) { double e = fhat[i] - f[i]; mse += e * e; }
        mse /= n;

        _out.WriteLine($"[bgp-fit] argmax-marginal t̂={best.Bandwidth:F3} logZ={best.LogMarginal:F1} in-sample MSE={mse:F5}");

        Assert.True(double.IsFinite(best.LogMarginal), "marginal evidence must be finite");
        Assert.True(mse < 0.005, $"GP posterior mean should recover the smooth function (MSE {mse:F5})");
        Assert.InRange(best.Bandwidth, 0.01, 0.5);
    }

    private static double Tn(int n, int d, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[n, 3];
        for (int i = 0; i < n; i++)
            for (int k = 0; k < d; k++) x[i, k] = rng.NextDouble();
        var gp = new GpRegression(x, new double[n], sigma2: 1.0, new SquaredExponentialKernel());
        return new EmpiricalBayesBandwidthPrior(x, gp, k: 5).Tn;
    }

    [Fact]
    public void Tn_ScalesAsInverseIntrinsicDimension()
    {
        // The averaged k-NN distance T_n ∼ n^{−1/d} (fixed k) — the other half of the dimension-implicit mechanism.
        const int d = 2;
        int[] ns = { 200, 400, 800, 1600 };
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        foreach (int n in ns)
        {
            double x = Math.Log(n), y = Math.Log(Tn(n, d, seed: 31 + n));
            sx += x; sy += y; sxx += x * x; sxy += x * y;
        }
        double slope = (ns.Length * sxy - sx * sy) / (ns.Length * sxx - sx * sx);

        _out.WriteLine($"[bgp-Tn] d={d}: T_n ∼ n^slope, slope={slope:F3} (expected {-1.0 / d:F2} = −1/d)");
        Assert.True(Math.Abs(slope - (-1.0 / d)) < 0.15, $"T_n slope {slope:F3} should track −1/d = {-1.0 / d:F2}");
    }

    [Fact]
    public void EbPrior_SupportedOnAdaptiveBand_AndPeaksInside()
    {
        var rng = new Xoshiro256PlusPlus(seed: 5);
        const int n = 300;
        var x = new double[n, 3];
        for (int i = 0; i < n; i++) { x[i, 0] = rng.NextDouble(); x[i, 1] = rng.NextDouble(); }  // d=2 in R³

        var gp = new GpRegression(x, new double[n], sigma2: 1.0, new SquaredExponentialKernel());
        var prior = new EmpiricalBayesBandwidthPrior(x, gp);

        _out.WriteLine($"[bgp-prior] k={prior.K} T_n={prior.Tn:F4} support=({prior.LowerBound:F5}, 1]");

        // Zero (−∞ log) outside the support, finite within.
        Assert.True(double.IsNegativeInfinity(prior.LogDensity(prior.LowerBound)), "0 at the support floor");
        Assert.True(double.IsNegativeInfinity(prior.LogDensity(2.0)), "0 above t=1");
        Assert.True(double.IsFinite(prior.LogDensity(0.5 * (prior.LowerBound + 1.0))), "finite inside the band");

        // The density peaks strictly inside the band (the adaptive scale), not at an edge.
        double argmax = double.NaN, best = double.NegativeInfinity;
        for (int i = 0; i <= 200; i++)
        {
            double t = prior.LowerBound + (1.0 - prior.LowerBound) * i / 200.0;
            double lp = prior.LogDensity(t);
            if (lp > best) { best = lp; argmax = t; }
        }
        _out.WriteLine($"[bgp-prior] argmax t̂={argmax:F4}");
        Assert.True(argmax > prior.LowerBound && argmax < 1.0, $"prior should peak inside the band (argmax {argmax:F4})");
    }
}
