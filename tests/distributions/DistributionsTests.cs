using System;
using Maths.Distributions;
using Maths.Rng;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// Pins the distribution primitives: the log-gamma normalizer against known factorials, the Beta/Gamma
/// log-densities against closed forms, and the samplers against their analytic means. These back the BARS
/// locality proposal and the scale-mixture robustness weights.
/// </summary>
public sealed class DistributionsTests
{
    [Fact]
    public void LogGamma_AtIntegers_MatchesLogFactorial()
    {
        Assert.Equal(0.0, SpecialFunctions.LogGamma(1.0), 9);            // Γ(1) = 0! = 1
        Assert.Equal(0.0, SpecialFunctions.LogGamma(2.0), 9);            // Γ(2) = 1! = 1
        Assert.Equal(Math.Log(24.0), SpecialFunctions.LogGamma(5.0), 9); // Γ(5) = 4! = 24
        Assert.Equal(Math.Log(120.0), SpecialFunctions.LogGamma(6.0), 9);// Γ(6) = 5! = 120
    }

    [Fact]
    public void Beta_LogPdf_MatchesClosedForm()
    {
        // Beta(2,2) pdf at 0.5: x(1-x)/B(2,2) = 0.25 / (1/6) = 1.5
        Assert.Equal(Math.Log(1.5), Beta.LogPdf(0.5, 2.0, 2.0), 9);
        Assert.Equal(double.NegativeInfinity, Beta.LogPdf(0.0, 2.0, 2.0));
        Assert.Equal(double.NegativeInfinity, Beta.LogPdf(1.0, 2.0, 2.0));
    }

    [Fact]
    public void Gamma_LogPdf_MatchesClosedForm()
    {
        // Gamma(1, scale) is Exponential(rate 1/scale): pdf(x) = (1/scale) e^(-x/scale).
        double scale = 2.0, x = 1.5;
        double expected = -Math.Log(scale) - x / scale;
        Assert.Equal(expected, Gamma.LogPdf(x, 1.0, scale), 9);
        Assert.Equal(double.NegativeInfinity, Gamma.LogPdf(0.0, 3.0));
    }

    [Fact]
    public void Beta_SampleMean_ApproachesAnalyticMean()
    {
        var rng = new Xoshiro256PlusPlus(seed: 20260612);
        double a = 2.0, b = 5.0;
        int n = 200_000;
        double sum = 0.0;
        for (int i = 0; i < n; i++) sum += Beta.Sample(rng, a, b);
        double mean = sum / n;
        Assert.Equal(a / (a + b), mean, 2);   // analytic mean 2/7 ≈ 0.2857, 2-decimal tolerance
    }

    [Fact]
    public void Gamma_SampleMean_ApproachesAnalyticMean()
    {
        var rng = new Xoshiro256PlusPlus(seed: 777);
        double shape = 0.5, scale = 3.0;   // shape < 1 exercises the Stuart boost
        int n = 200_000;
        double sum = 0.0;
        for (int i = 0; i < n; i++) sum += Gamma.Sample(rng, shape, scale);
        double mean = sum / n;
        Assert.Equal(shape * scale, mean, 1);   // analytic mean 1.5, 1-decimal tolerance
    }
}
