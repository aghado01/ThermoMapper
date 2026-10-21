using System;
using Maths.Distributions;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// The normal quantile Φ⁻¹ (AS 241), checked against tabulated values and, independently, against the log-space CDF
/// <see cref="Normal.StandardLogCdf"/>, which shares no code with it: <c>ln Φ(Φ⁻¹(p)) = ln p</c> from the centre to
/// p = 10⁻³⁰⁰, in both tails.
/// </summary>
public sealed class NormalQuantileTests
{
    [Theory]
    [InlineData(0.5, 0.0, 0.0)]
    [InlineData(0.975, 1.959963984540054, 1e-15)]
    [InlineData(0.025, -1.959963984540054, 1e-15)]
    [InlineData(0.95, 1.6448536269514722, 1e-15)]
    [InlineData(0.999, 3.090232306167813, 1e-14)]
    [InlineData(1e-10, -6.361340902404056, 1e-9)]
    public void Quantile_MatchesTabulatedValues(double p, double expected, double tolerance)
    {
        Assert.Equal(expected, Normal.StandardQuantile(p), tolerance);
    }

    [Fact]
    public void Quantile_InvertsTheLogCdf_IntoBothTails()
    {
        for (int e = 1; e <= 300; e++)
            foreach (double mantissa in new[] { 1.0, 3.7, 7.1 })
            {
                double p = mantissa * Math.Pow(10.0, -e);
                if (p >= 0.5) continue;
                double z = Normal.StandardQuantile(p);
                Assert.Equal(Math.Log(p), Normal.StandardLogCdf(z), 1e-11 * Math.Max(1.0, -Math.Log(p)));
            }
        for (double p = 0.01; p < 1.0; p += 0.01)
        {
            Assert.Equal(Math.Log(p), Normal.StandardLogCdf(Normal.StandardQuantile(p)), 1e-14);
            Assert.Equal(-Normal.StandardQuantile(p), Normal.StandardQuantile(1.0 - p), 1e-13);   // 1 − p is exact enough here
        }
    }

    [Fact]
    public void Quantile_AtAndBeyondTheEnds()
    {
        Assert.Equal(double.NegativeInfinity, Normal.StandardQuantile(0.0));
        Assert.Equal(double.PositiveInfinity, Normal.StandardQuantile(1.0));
        Assert.True(double.IsNaN(Normal.StandardQuantile(-0.1)));
        Assert.True(double.IsNaN(Normal.StandardQuantile(1.1)));
        Assert.True(double.IsNaN(Normal.StandardQuantile(double.NaN)));
        Assert.Equal(2.0 + 3.0 * 1.959963984540054, Normal.MeanSd(2.0, 3.0).Quantile(0.975), 1e-14);
    }
}
