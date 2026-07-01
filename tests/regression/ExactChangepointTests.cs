using System;
using Maths.Regression.Changepoint;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The exact (dynamic-programming) piecewise-constant inference (Hutter 2005) recovers the change-point structure
/// without MCMC: on a 3-segment signal it MAPs the right number of segments, the break probabilities spike at the
/// true boundaries, and the exact identity Σ P(break after i) = E[#breaks] holds.
/// </summary>
public sealed class ExactChangepointTests
{
    private readonly ITestOutputHelper _out;
    public ExactChangepointTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void Infer_RecoversSegmentCount_AndBreaks()
    {
        var rng = new Xoshiro256PlusPlus(seed: 71);
        const int n = 60;
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            double level = i < 20 ? 0.0 : i < 40 ? 2.0 : -1.0;   // breaks after positions 20 and 40 (1-based)
            y[i] = level + 0.30 * Gaussian(rng);
        }

        double mean = 0.0; foreach (double v in y) mean += v; mean /= n;
        ChangepointResult r = ExactChangepoint.Infer(y, sigma2: 0.1, nu: mean, rho2: 2.0, maxSegments: 8);

        _out.WriteLine($"[changepoint] k̂={r.MapSegments} P(k|y)=[{string.Join(", ", Array.ConvertAll(r.SegmentNumberPosterior, p => p.ToString("F3")))}]");
        _out.WriteLine($"[changepoint] breaks @20={r.BreakProbability[19]:F3} @40={r.BreakProbability[39]:F3} logEvidence={r.LogEvidence:F2}");

        // 3 segments recovered.
        Assert.Equal(3, r.MapSegments);

        // Break probability spikes at the two true boundaries (after positions 20 and 40), flat elsewhere.
        Assert.True(r.BreakProbability[19] > 0.7, $"break after 20 should be near-certain (got {r.BreakProbability[19]:F3})");
        Assert.True(r.BreakProbability[39] > 0.7, $"break after 40 should be near-certain (got {r.BreakProbability[39]:F3})");
        double maxSpurious = 0.0;
        for (int i = 0; i < r.BreakProbability.Length; i++)
            if (Math.Abs(i - 19) > 1 && Math.Abs(i - 39) > 1)
                maxSpurious = Math.Max(maxSpurious, r.BreakProbability[i]);
        Assert.True(maxSpurious < 0.3, $"no spurious breaks away from the true boundaries (max {maxSpurious:F3})");

        // Exact identity: Σ_i P(break after i) = E[#breaks] = Σ_k (k−1) P(k|y).
        double sumBreaks = 0.0;
        foreach (double p in r.BreakProbability) sumBreaks += p;
        double expectedBreaks = 0.0;
        for (int s = 1; s <= r.SegmentNumberPosterior.Length; s++) expectedBreaks += (s - 1) * r.SegmentNumberPosterior[s - 1];
        Assert.Equal(expectedBreaks, sumBreaks, 6);

        Assert.True(double.IsFinite(r.LogEvidence));
    }

    [Fact]
    public void RegressionCurve_RecoversLevels_AndJumpsAtBreaks()
    {
        var rng = new Xoshiro256PlusPlus(seed: 88);
        const int n = 60;
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            double level = i < 20 ? 0.0 : i < 40 ? 2.0 : -1.0;
            y[i] = level + 0.30 * Gaussian(rng);
        }
        double mean = 0.0; foreach (double v in y) mean += v; mean /= n;
        ChangepointResult r = ExactChangepoint.Infer(y, sigma2: 0.1, nu: mean, rho2: 2.0, maxSegments: 8);

        _out.WriteLine($"[curve] seg1≈{r.RegressionCurve[5]:F2} seg2≈{r.RegressionCurve[30]:F2} seg3≈{r.RegressionCurve[50]:F2}; " +
                       $"break20: {r.RegressionCurve[19]:F2}→{r.RegressionCurve[20]:F2}; std mid={r.RegressionStd[30]:F3}");

        // Mid-segment levels recovered (posterior-mean curve ≈ true levels).
        Assert.InRange(r.RegressionCurve[5], -0.3, 0.3);    // segment 1 ≈ 0
        Assert.InRange(r.RegressionCurve[30], 1.7, 2.3);    // segment 2 ≈ 2
        Assert.InRange(r.RegressionCurve[50], -1.3, -0.7);  // segment 3 ≈ −1

        // The curve jumps at the boundary (between points 20 and 21) rather than blurring across it.
        double jump = Math.Abs(r.RegressionCurve[20] - r.RegressionCurve[19]);
        Assert.True(jump > 1.5, $"curve should jump at the break, not blur it (got {jump:F3})");

        // Error band: finite, and modest within a clear segment.
        Assert.True(r.RegressionStd[30] > 0.0 && r.RegressionStd[30] < 0.5, $"mid-segment std {r.RegressionStd[30]:F3}");
    }

    [Fact]
    public void Cauchy_IsRobust_WhereGaussianIsFooled()
    {
        var rng = new Xoshiro256PlusPlus(seed: 123);
        const int n = 60;
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            double level = i < 30 ? 0.0 : 3.0;       // one true break after position 30
            y[i] = level + 0.30 * Gaussian(rng);
        }
        foreach (int o in new[] { 10, 20, 45 }) y[o] += 10.0;   // gross outliers inside the segments

        ChangepointResult gauss = ExactChangepoint.Infer(y, sigma2: 0.1, nu: 1.5, rho2: 4.0, maxSegments: 10);
        ChangepointResult cauchy = ExactChangepoint.Infer(
            y, new CauchySegmentModel(y, scale: 0.5, nu: 1.5, rho2: 4.0), maxSegments: 10);

        _out.WriteLine($"[robust] gaussian k̂={gauss.MapSegments}; cauchy k̂={cauchy.MapSegments} break@30={cauchy.BreakProbability[29]:F3}");

        // Cauchy recovers the true 2-segment structure despite the outliers.
        Assert.Equal(2, cauchy.MapSegments);
        Assert.True(cauchy.BreakProbability[29] > 0.5, $"true break after 30 (got {cauchy.BreakProbability[29]:F3})");
        double spurious = 0.0;
        for (int i = 0; i < cauchy.BreakProbability.Length; i++)
            if (Math.Abs(i - 29) > 1) spurious = Math.Max(spurious, cauchy.BreakProbability[i]);
        Assert.True(spurious < 0.3, $"no spurious breaks at the outliers (max {spurious:F3})");

        // The Gaussian model is fooled — it isolates the outliers into extra segments.
        Assert.True(gauss.MapSegments > cauchy.MapSegments,
            $"gaussian should over-segment around outliers (k̂={gauss.MapSegments}) vs cauchy ({cauchy.MapSegments})");
    }
}
