using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// Validates the proposal seam: the mixture-density log-sum-exp, the Uniform baseline, and an end-to-end
/// LocalBeta chain. The Hastings correctness of birth/death is structural (the death ratio is the negative
/// of the reverse birth's by construction); these pin the density math and that the locality path runs.
/// </summary>
public sealed class ProposalTests
{
    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void MixtureDensity_SingleCenter_EqualsKernelDensity()
    {
        var kernel = new LocalBetaKernel(tau: 50.0);
        double x = 0.37, center = 0.5;
        Assert.Equal(kernel.LogDensity(x, center),
                     ProposalMath.LogMixtureDensity(x, new[] { center }, kernel), 12);
    }

    [Fact]
    public void MixtureDensity_NoCenters_IsUniform()
    {
        var kernel = new LocalBetaKernel();
        Assert.Equal(0.0, ProposalMath.LogMixtureDensity(0.5, Array.Empty<double>(), kernel), 12);
    }

    [Fact]
    public void UniformKernel_GivesZeroLogDensity()
    {
        var kernel = new UniformKernel();
        Assert.Equal(0.0, kernel.LogDensity(0.3, 0.7));
        Assert.Equal(0.0, ProposalMath.LogMixtureDensity(0.3, new[] { 0.2, 0.6 }, kernel), 12);
    }

    [Fact]
    public void Chain_LocalBeta_RunsEndToEnd()
    {
        var basis = new SplineBasis(degree: 3);
        var noise = new Xoshiro256PlusPlus(seed: 5);
        int m = 60;
        var x = new double[m];
        for (int i = 0; i < m; i++) x[i] = (i + 0.5) / m;

        var trueConfig = new KnotConfig(new[] { 0.4 });
        double[,] zt = basis.Design(trueConfig, x);
        int nut = zt.GetLength(1);
        var y = new double[m];
        for (int i = 0; i < m; i++)
        {
            double s = 0.0;
            for (int j = 0; j < nut; j++) s += zt[i, j] * (1.0 + 0.3 * j);
            y[i] = s + 0.05 * Gaussian(noise);
        }

        var target = new SplineTarget(basis, new WeightedNormalModel(), new PoissonPrior(3.0), x, y);
        var kernel = new LocalBetaKernel(tau: 50.0);
        var moves = new IRjMove<KnotConfig>[]
        {
            new KnotBirthMove(kernel), new KnotDeathMove(kernel), new KnotRelocateMove(kernel),
        };
        var chain = new ReversibleJumpChain<KnotConfig>(
            moves, target, new KnotConfig(Array.Empty<double>()), new Xoshiro256PlusPlus(seed: 101));

        for (int s = 0; s < 500; s++) chain.Step();

        Assert.True(chain.Attempts > 0);
        Assert.InRange((double)chain.Accepted / chain.Attempts, 0.0, 1.0);
    }
}
