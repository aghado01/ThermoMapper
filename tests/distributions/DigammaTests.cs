using System;
using Maths.Distributions;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// The digamma function: closed values at rational points (Gauss's digamma theorem), the recurrence and the reflection,
/// the small-x Laurent series, the large-x asymptote, and the derivative of ln Γ by central differences.
/// </summary>
public sealed class DigammaTests
{
    private const double EulerGamma = 0.57721566490153286061;

    [Theory]
    [InlineData(1.0, -EulerGamma)]
    [InlineData(2.0, 1.0 - EulerGamma)]
    [InlineData(0.5, -EulerGamma - 1.3862943611198906)]                                          // −γ − 2 ln 2
    [InlineData(0.25, -EulerGamma - 1.5707963267948966 - 2.0794415416798357)]                     // −γ − π/2 − 3 ln 2
    [InlineData(0.75, -EulerGamma + 1.5707963267948966 - 2.0794415416798357)]                     // −γ + π/2 − 3 ln 2
    [InlineData(10.0, 2.2517525890667211)]
    [InlineData(100.0, 4.6001618527380874)]
    public void ClosedValues(double x, double expected)
        => Assert.True(Math.Abs(SpecialFunctions.Digamma(x) - expected) <= 2e-15 * Math.Max(1.0, Math.Abs(expected)),
                       $"ψ({x}) = {SpecialFunctions.Digamma(x):R} against {expected:R}");

    [Fact]
    public void Recurrence_AndReflection()
    {
        foreach (double x in new[] { 1e-6, 0.003, 0.37, 1.4616321449683623, 2.5, 7.9, 14.99, 15.0, 31.4, 1e4 })
        {
            double step = SpecialFunctions.Digamma(x + 1.0) - SpecialFunctions.Digamma(x) - 1.0 / x;
            Assert.True(Math.Abs(step) <= 1e-14 * Math.Max(1.0, 1.0 / x), $"x = {x}: recurrence off by {step:E2}");
        }
        foreach (double x in new[] { 0.1, 0.3, 0.7, 2.6, -0.4, -3.75 })
        {
            double reflected = SpecialFunctions.Digamma(1.0 - x) - SpecialFunctions.Digamma(x) - Math.PI / Math.Tan(Math.PI * x);
            Assert.True(Math.Abs(reflected) <= 1e-13, $"x = {x}: reflection off by {reflected:E2}");
        }
        Assert.True(double.IsNaN(SpecialFunctions.Digamma(0.0)));
        Assert.True(double.IsNaN(SpecialFunctions.Digamma(-2.0)));
        Assert.Equal(double.PositiveInfinity, SpecialFunctions.Digamma(double.PositiveInfinity));
    }

    [Fact]
    public void SmallAndLargeArguments()
    {
        // ψ(x) = −1/x − γ + ζ(2)x − ζ(3)x² + ζ(4)x³ − …
        double x = 1e-3;
        double series = -1.0 / x - EulerGamma + 1.6449340668482264 * x - 1.2020569031595943 * x * x + 1.0823232337111382 * x * x * x;
        Assert.True(Math.Abs(SpecialFunctions.Digamma(x) - series) <= 1e-12 * Math.Abs(series));

        // ψ(x) ≈ ln x − 1/(2x) − 1/(12x²) for large x.
        double big = 1e6;
        Assert.True(Math.Abs(SpecialFunctions.Digamma(big) - (Math.Log(big) - 0.5 / big - 1.0 / (12 * big * big))) <= 1e-15 * Math.Log(big));
    }

    [Fact]
    public void IsTheDerivativeOfLogGamma()
    {
        foreach (double x in new[] { 0.2, 0.9, 3.3, 12.0, 40.0 })
        {
            double h = 1e-5 * x;
            double difference = (SpecialFunctions.LogGamma(x + h) - SpecialFunctions.LogGamma(x - h)) / (2 * h);
            Assert.True(Math.Abs(difference - SpecialFunctions.Digamma(x)) <= 1e-7 * Math.Max(1.0, Math.Abs(difference)),
                        $"x = {x}: {SpecialFunctions.Digamma(x):R} against {difference:R}");
        }
    }
}
