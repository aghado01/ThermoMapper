using System;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// The Weibull distribution: its density integrates to its CDF and reduces to the exponential at shape 1; at shape ½ it
/// is the law of τ² for an exponential τ; its tails keep their digits in log space; its quantile inverts its CDF; and its
/// draws have its moments and quantiles.
/// </summary>
public sealed class WeibullTests
{
    [Theory]
    [InlineData(0.5, 2.0)]
    [InlineData(1.0, 0.7)]
    [InlineData(3.0, 1.5)]
    public void Density_IntegratesToTheCdf(double shape, double scale)
    {
        Weibull w = Weibull.ShapeScale(shape, scale);
        foreach (double x in new[] { 0.3, 1.0, 4.0 })
        {
            // ∫₀ˣ f by the substitution t = x·s², which removes the k < 1 singularity at 0: f(t) dt = f(xs²)·2xs ds.
            const int m = 200_000;
            double sum = 0.0;
            for (int i = 0; i < m; i++)
            {
                double s = (i + 0.5) / m;
                sum += Math.Exp(w.LogPdf(x * s * s)) * 2.0 * x * s / m;
            }
            Assert.Equal(w.Cdf(x), sum, 6);
        }
    }

    [Fact]
    public void Shapes_OneAndAHalf_AreTheExponentialAndTheSquaredExponential()
    {
        Weibull exponential = Weibull.ShapeScale(1.0, 2.5);
        foreach (double x in new[] { 0.1, 1.0, 7.0 })
            Assert.Equal(-Math.Log(2.5) - x / 2.5, exponential.LogPdf(x), 12);

        // τ ~ Exp(rate θ) gives P(τ² ≤ x) = 1 − e^(−θ√x), which is W(½, 1/θ²).
        const double theta = 0.8;
        Weibull squared = Weibull.ShapeScale(0.5, 1.0 / (theta * theta));
        foreach (double x in new[] { 0.01, 0.5, 3.0, 20.0 })
            Assert.Equal(1.0 - Math.Exp(-theta * Math.Sqrt(x)), squared.Cdf(x), 12);
        Assert.Equal(2.0 / (theta * theta), squared.Mean, 10);
        Assert.Equal(20.0 / Math.Pow(theta, 4), squared.Variance, 8);
    }

    [Fact]
    public void Tails_KeepTheirDigits_AndTheQuantileInverts()
    {
        Weibull w = Weibull.ShapeScale(0.5, 3.0);
        Assert.Equal(0.5 * Math.Log(1e-20 / 3.0), w.LogCdf(1e-20), 9);   // F(x) ≈ (x/λ)^k near 0
        Assert.Equal(-Math.Sqrt(1e6 / 3.0), w.LogSurvival(1e6), 9);
        Assert.Equal(w.Survival(2.0) - w.Survival(5.0), w.IntervalProbability(2.0, 5.0), 14);
        foreach (double p in new[] { 1e-12, 0.01, 0.5, 0.99 })
            Assert.Equal(p, w.Cdf(w.Quantile(p)), 12);
        Assert.Equal(w.Quantile(0.5), w.Median, 12);
        Assert.Equal(0.0, w.Mode);
        Assert.Equal(1.5 * Math.Pow(2.0 / 3.0, 1.0 / 3.0), Weibull.ShapeScale(3.0, 1.5).Mode, 12);
    }

    [Fact]
    public void Draws_HaveTheMomentsAndQuantiles()
    {
        Weibull w = Weibull.ShapeScale(2.0, 1.3);
        var rng = new Xoshiro256PlusPlus(seed: 4);
        double[] x = Enumerable.Range(0, 400_000).Select(_ => w.Sample(rng)).OrderBy(v => v).ToArray();
        double mean = x.Average(), variance = x.Average(v => (v - mean) * (v - mean));
        Assert.InRange(mean, w.Mean * 0.995, w.Mean * 1.005);
        Assert.InRange(variance, w.Variance * 0.99, w.Variance * 1.01);
        foreach (double p in new[] { 0.05, 0.5, 0.95 })
            Assert.InRange(x[(int)(p * x.Length)], w.Quantile(p - 0.002), w.Quantile(p + 0.002));
    }

    [Fact]
    public void InvalidParameters_AreRefused()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Weibull.ShapeScale(0.0, 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Weibull.ShapeScale(1.0, double.PositiveInfinity));
        Assert.Throws<ArgumentOutOfRangeException>(() => Weibull.ShapeScale(1.0, 1.0).Quantile(1.5));
    }
}
