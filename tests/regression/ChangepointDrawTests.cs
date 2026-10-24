using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Changepoint;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The retained changepoint posterior (<see cref="ChangepointResult.Posterior"/>) against brute force: on a series short
/// enough to enumerate every segmentation, with each segment's Gaussian evidence from an explicit covariance, the
/// segment-count posterior, break marginals, segment probabilities, curve and evidence agree with the dynamic program,
/// and the exact draws reproduce the posterior over segmentations and the curve's moments.
/// </summary>
public sealed class ChangepointDrawTests
{
    private const double Sigma2 = 0.3, Nu = 0.5, Rho2 = 1.5;
    private static readonly double[] Short = { 0.1, -0.3, 0.4, 1.6, 1.2, 2.1, 0.9 };

    // log N(y_{a..b}; ν1, σ²I + ρ²11ᵀ) (0-based, inclusive) through an explicit Cholesky factor.
    private static double LogSegmentEvidence(double[] y, int first, int last)
    {
        int len = last - first + 1;
        var l = new double[len, len];
        for (int i = 0; i < len; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = (i == j ? Sigma2 : 0.0) + Rho2;
                for (int k = 0; k < j; k++) s -= l[i, k] * l[j, k];
                l[i, j] = i == j ? Math.Sqrt(s) : s / l[j, j];
            }
        var z = new double[len];
        double logDet = 0.0, quad = 0.0;
        for (int i = 0; i < len; i++)
        {
            double s = y[first + i] - Nu;
            for (int k = 0; k < i; k++) s -= l[i, k] * z[k];
            z[i] = s / l[i, i];
            logDet += 2.0 * Math.Log(l[i, i]);
            quad += z[i] * z[i];
        }
        return -0.5 * (len * Math.Log(2.0 * Math.PI) + logDet + quad);
    }

    private static double LogFactorial(int m)
    {
        double s = 0.0;
        for (int i = 2; i <= m; i++) s += Math.Log(i);
        return s;
    }

    // Every segmentation of the short series, as its break mask (bit i: a break after index i), with its unnormalized log
    // posterior: −log K − log C(n−1, k−1) + Σ log D(segment), and k uniform on 1..n.
    private static (int Mask, double LogWeight)[] Enumerate(double[] y)
    {
        int n = y.Length, gaps = n - 1;
        var all = new (int, double)[1 << gaps];
        for (int mask = 0; mask < 1 << gaps; mask++)
        {
            int k = 1 + System.Numerics.BitOperations.PopCount((uint)mask);
            double w = -Math.Log(n) - (LogFactorial(gaps) - LogFactorial(k - 1) - LogFactorial(gaps - k + 1));
            int first = 0;
            for (int i = 0; i < n; i++)
                if (i == n - 1 || (mask >> i & 1) == 1) { w += LogSegmentEvidence(y, first, i); first = i + 1; }
            all[mask] = (mask, w);
        }
        return all;
    }

    private static double[] Normalize((int Mask, double LogWeight)[] all, out double logTotal)
    {
        double top = all.Max(e => e.LogWeight);
        double sum = all.Sum(e => Math.Exp(e.LogWeight - top));
        logTotal = top + Math.Log(sum);
        double lt = logTotal;
        return all.Select(e => Math.Exp(e.LogWeight - lt)).ToArray();
    }

    private static int MaskOf(ChangepointDraw d)
    {
        int mask = 0;
        for (int s = 0; s + 1 < d.SegmentEnds.Length; s++) mask |= 1 << d.SegmentEnds[s];
        return mask;
    }

    private static IEnumerable<(int First, int Last)> Segments(int mask, int n)
    {
        int first = 0;
        for (int i = 0; i < n; i++)
            if (i == n - 1 || (mask >> i & 1) == 1) { yield return (first, i); first = i + 1; }
    }

    [Fact]
    public void DynamicProgram_MatchesBruteForce()
    {
        double[] y = Short;
        int n = y.Length;
        ChangepointResult r = ExactChangepoint.Infer(y, Sigma2, Nu, Rho2, maxSegments: n);
        (int Mask, double LogWeight)[] all = Enumerate(y);
        double[] p = Normalize(all, out double logTotal);

        Assert.Equal(logTotal, r.LogEvidence, 10);
        Assert.Equal(logTotal, r.Posterior.LogEvidence, 10);
        for (int k = 1; k <= n; k++)
            Assert.Equal(all.Where(e => 1 + System.Numerics.BitOperations.PopCount((uint)e.Mask) == k).Sum(e => p[e.Mask]),
                         r.SegmentNumberPosterior[k - 1], 12);
        for (int i = 0; i < n - 1; i++)
            Assert.Equal(all.Where(e => (e.Mask >> i & 1) == 1).Sum(e => p[e.Mask]), r.BreakProbability[i], 12);

        var curve = new double[n];
        for (int first = 0; first < n; first++)
            for (int last = first; last < n; last++)
            {
                Assert.Equal(LogSegmentEvidence(y, first, last), r.Posterior.SegmentLogEvidence(first, last), 10);
                double ps = all.Where(e => Segments(e.Mask, n).Contains((first, last))).Sum(e => p[e.Mask]);
                Assert.Equal(ps, r.Posterior.SegmentProbability(first, last), 12);
                double sum = 0.0;
                for (int i = first; i <= last; i++) sum += y[i];
                double level = (Nu / Rho2 + sum / Sigma2) / (1.0 / Rho2 + (last - first + 1) / Sigma2);
                for (int i = first; i <= last; i++) curve[i] += ps * level;
            }
        for (int i = 0; i < n; i++) Assert.Equal(curve[i], r.RegressionCurve[i], 12);
    }

    [Fact]
    public void Draws_FollowTheExactPosterior_OverSegmentationsAndLevels()
    {
        double[] y = Short;
        int n = y.Length;
        ChangepointResult r = ExactChangepoint.Infer(y, Sigma2, Nu, Rho2, maxSegments: n);
        double[] p = Normalize(Enumerate(y), out _);

        const int draws = 100_000;
        IReadOnlyList<ChangepointDraw> sample = r.Posterior.Sample(draws, seed: 17);
        var counts = new int[p.Length];
        var sum = new double[n];
        var sumSq = new double[n];
        foreach (ChangepointDraw d in sample)
        {
            Assert.Equal(n - 1, d.SegmentEnds[^1]);
            Assert.Equal(d.SegmentEnds.Length, d.Levels.Length);
            counts[MaskOf(d)]++;
            double[] c = d.Curve();
            for (int i = 0; i < n; i++) { sum[i] += c[i]; sumSq[i] += c[i] * c[i]; }
        }
        for (int mask = 0; mask < p.Length; mask++)
        {
            double f = counts[mask] / (double)draws, se = Math.Sqrt(p[mask] * (1.0 - p[mask]) / draws);
            Assert.True(Math.Abs(f - p[mask]) <= 5.0 * se + 1e-9, $"segmentation {mask}: {f} vs {p[mask]}");
        }
        for (int i = 0; i < n; i++)
        {
            double mean = sum[i] / draws, variance = sumSq[i] / draws - mean * mean, sd = r.RegressionStd[i];
            Assert.True(Math.Abs(mean - r.RegressionCurve[i]) <= 5.0 * sd / Math.Sqrt(draws), $"point {i}: {mean} vs {r.RegressionCurve[i]}");
            Assert.InRange(variance / (sd * sd), 0.97, 1.03);   // a sample variance's relative SE is about √(2/N) ≈ 0.45%
        }
    }

    [Fact]
    public void CauchyDraws_MatchTheQuadraturesMoments()
    {
        var rng = new Xoshiro256PlusPlus(seed: 5);
        int n = 30;
        var y = new double[n];
        for (int i = 0; i < n; i++) y[i] = (i < 12 ? 0.0 : i < 22 ? 2.0 : -1.0) + 0.3 * Math.Tan(Math.PI * (rng.NextDouble() - 0.5));
        ChangepointResult r = ExactChangepoint.Infer(y, new CauchySegmentModel(y, 0.3, 0.3, 2.0), maxSegments: 8);

        const int draws = 20_000;
        var sum = new double[n];
        foreach (ChangepointDraw d in r.Posterior.Sample(draws, seed: 3))
        {
            double[] c = d.Curve();
            for (int i = 0; i < n; i++) sum[i] += c[i];
        }
        for (int i = 0; i < n; i++)
            Assert.True(Math.Abs(sum[i] / draws - r.RegressionCurve[i]) <= 5.0 * r.RegressionStd[i] / Math.Sqrt(draws),
                        $"point {i}: {sum[i] / draws} vs {r.RegressionCurve[i]}");
    }

    [Fact]
    public void EachPoint_LiesInExactlyOneSegment()
    {
        var rng = new Xoshiro256PlusPlus(seed: 9);
        int n = 40;
        var y = new double[n];
        for (int i = 0; i < n; i++) y[i] = (i < 15 ? 0.0 : 1.5) + 0.4 * Math.Sqrt(-2.0 * Math.Log(1.0 - rng.NextDouble())) * Math.Cos(2.0 * Math.PI * rng.NextDouble());
        (double s2, double nu, double rho2) = ExactChangepoint.DefaultHyperparameters(y);
        foreach (ChangepointResult r in new[] { ExactChangepoint.Infer(y, s2, nu, rho2, 8),
                                                ExactChangepoint.Infer(y, new CauchySegmentModel(y, 0.4, nu, rho2), 8) })
            for (int i = 0; i < n; i++)
            {
                double covered = 0.0;
                for (int first = 0; first <= i; first++)
                    for (int last = i; last < n; last++) covered += r.Posterior.SegmentProbability(first, last);
                Assert.Equal(1.0, covered, 12);
            }
    }

    [Fact]
    public void Draws_AreReproducibleFromTheirSeed()
    {
        ChangepointResult r = ExactChangepoint.Infer(Short, Sigma2, Nu, Rho2, maxSegments: 4);
        IReadOnlyList<ChangepointDraw> a = r.Posterior.Sample(50, seed: 8), b = r.Posterior.Sample(50, seed: 8);
        for (int d = 0; d < 50; d++)
        {
            Assert.Equal(a[d].SegmentEnds, b[d].SegmentEnds);
            Assert.Equal(a[d].Levels, b[d].Levels);
            Assert.True(a[d].SegmentEnds.Length <= 4);
        }
        Assert.Empty(r.Posterior.Sample(0, seed: 1));
    }
}
