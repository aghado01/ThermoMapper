using System;
using Maths.Distributions;
using Maths.Regression;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// GBASS's <c>nw_est_mom</c> and <c>nw_gamma_prior</c> corrected (E1's slice 3, slice 0's defects): the NIG by its moments recovers an NIG from
/// that NIG's moments taken by quadrature of its density, and from a large sample; the prior on γ puts its two masses where asked, read by
/// quadrature of the normal's part on γ &gt; 0, and GBASS's default targets are refused with the most a normal reaches.
/// </summary>
public sealed class NormalWaldCalibrationTests
{
    [Theory]
    [InlineData(0.3, 1.0, 1.0, 2.0)]
    [InlineData(-1.0, 0.5, -2.0, 3.0)]
    [InlineData(2.0, 4.0, 0.3, 1.5)]
    public void Moments_RecoverTheNig_FromItsDensitysMoments(double mu, double delta, double beta, double alpha)
    {
        double gamma = Math.Sqrt(alpha * alpha - beta * beta);
        // X − μ = β·V + √V·Z with V ~ IG(δ/γ, δ²) = GIG(−½, γ², δ²).
        var law = GeneralizedHyperbolic.Of(-0.5, gamma * gamma, delta * delta, beta, 1.0);
        double Moment(int k, double centre) => Quadrature.LogIntegral(x => law.LogPdf(x) + k * Math.Log(Math.Abs(x - centre)), 0.0, lower: centre)
                                               is var right && Quadrature.LogIntegral(x => law.LogPdf(x) + k * Math.Log(Math.Abs(x - centre)), 0.0, upper: centre) is var left
                                               ? Math.Exp(right) + (k % 2 == 0 ? 1.0 : -1.0) * Math.Exp(left) : 0.0;
        double mean = Moment(1, 0.0), variance = Moment(2, mean), third = Moment(3, mean), fourth = Moment(4, mean);
        NormalInverseGaussian fit = NormalWaldCalibration.Moments(mu + mean, variance, third / Math.Pow(variance, 1.5), fourth / (variance * variance));
        Assert.Equal(mu, fit.Mu, 7);
        Assert.Equal(delta, fit.Delta, 7);
        Assert.Equal(beta, fit.Beta, 7);
        Assert.Equal(alpha, fit.Alpha, 7);
        // In the core's terms: σ = δ, β′ = βδ, γ′ = γδ, and the steepness (1 + δγ)^(−½).
        Assert.Equal(delta * delta, fit.Sigma2, 6);
        Assert.Equal(gamma * delta, fit.CoreGamma, 6);
        Assert.Equal(1.0 / Math.Sqrt(1.0 + delta * gamma), fit.Steepness, 7);
    }

    [Fact]
    public void Moments_FromASample_AndRefusedOutsideTheNigs()
    {
        var law = GeneralizedHyperbolic.Of(-0.5, 4.0 - 0.36, 1.0, 0.6, 1.0);   // α = 2, β = 0.6, δ = 1
        var rng = new Xoshiro256PlusPlus(20261114);
        var sample = new double[400_000];
        for (int i = 0; i < sample.Length; i++) sample[i] = law.Sample(rng);
        NormalInverseGaussian fit = NormalWaldCalibration.Moments(sample);
        Assert.InRange(fit.Alpha, 1.8, 2.2);
        Assert.InRange(fit.Beta, 0.5, 0.7);
        Assert.InRange(fit.Delta, 0.9, 1.1);
        Assert.Throws<ArgumentException>(() => NormalWaldCalibration.Moments(0.0, 1.0, 1.0, 3.5));   // 3·0.5 < 4·1
    }

    [Theory]
    [InlineData(0.5, 0.9, 0.5, 0.03)]
    [InlineData(0.25, 0.75, 0.5, 0.02)]
    [InlineData(0.3, 0.6, 0.2, 0.15)]
    public void GammaPrior_PutsItsMassesWhereAsked(double q1, double q2, double p1, double p2)
    {
        (double m, double s) = NormalWaldCalibration.GammaPrior(q1, q2, p1, p2);
        double g1 = 1.0 / (q1 * q1) - 1.0, g2 = 1.0 / (q2 * q2) - 1.0;
        double Density(double g) => -0.5 * Math.Pow((g - m) / s, 2.0);
        double total = Quadrature.LogIntegral(Density, Math.Max(m, 0.0), lower: 0.0);
        double above = Math.Exp(Quadrature.LogIntegral(Density, Math.Max(m, g1), lower: g1) - total);
        double below = Math.Exp(Quadrature.LogIntegral(Density, Math.Clamp(m, 0.0, g2), lower: 0.0, upper: g2) - total);
        Assert.Equal(p1, above, 8);
        Assert.Equal(p2, below, 8);
    }

    [Fact]
    public void GammaPrior_RefusesGbasssDefaults_WithTheMostANormalReaches()
    {
        var refusal = Assert.Throws<ArgumentException>(() => NormalWaldCalibration.GammaPrior(0.1, 0.9, 0.5, 0.05));
        Assert.Contains("the most it puts below", refusal.Message);
        // At half the mass above g₁ the bound is the exponential limit 1 − 2^(−g₂/g₁): 0.0016 at GBASS's levels.
        double g1 = 99.0, g2 = 1.0 / 0.81 - 1.0;
        Assert.Throws<ArgumentException>(() => NormalWaldCalibration.GammaPrior(0.1, 0.9, 0.5, 1.05 * (1.0 - Math.Pow(2.0, -g2 / g1))));
        (double m, double s) = NormalWaldCalibration.GammaPrior(0.1, 0.9, 0.5, 0.9 * (1.0 - Math.Pow(2.0, -g2 / g1)));
        Assert.True(m < 0.0 && s > 0.0, $"near the bound the mean runs negative: ({m}, {s})");
    }
}
