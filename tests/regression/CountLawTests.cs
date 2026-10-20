using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// R01 / Q3: BARS's knot measure. Knot locations are the order statistics of k uniforms (density k! on the ordered
/// simplex), so with the likelihood held constant the knot count must follow the count prior p(k) exactly. Checked
/// under both move schedules, constant weights and DMGK's, which settles that the two sample the same posterior;
/// for a hard cap; and for the Conway–Maxwell–Poisson tails. The kernels mix within the batches: under a flat
/// target the default LocalBeta(50) wanders slowly enough that 40 batches understate its error.
/// </summary>
public sealed class CountLawTests
{
    private readonly ITestOutputHelper _out;
    public CountLawTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData("poisson", false, "uniform")]
    [InlineData("poisson", true, "uniform")]
    [InlineData("poisson", false, "beta")]
    [InlineData("poisson", true, "beta")]
    [InlineData("uniform", false, "beta")]
    [InlineData("uniform", true, "beta")]
    [InlineData("cmp-nu2", false, "beta")]
    [InlineData("cmp-nu2", true, "beta")]
    [InlineData("cmp-nu0.5", false, "beta")]
    public void KnotCount_FollowsTheCountPrior(string priorName, bool dmgk, string kernelName)
    {
        (IComplexityPrior prior, int bins) = priorName switch
        {
            "poisson" => ((IComplexityPrior)new PoissonPrior(4.0), 12),
            "uniform" => (new UniformPrior(6), 9),
            "cmp-nu2" => (new ConwayMaxwellPoissonPrior(4.0, 2.0), 7),
            "cmp-nu0.5" => (new ConwayMaxwellPoissonPrior(1.5, 0.5), 12),
            _ => throw new ArgumentOutOfRangeException(nameof(priorName)),
        };
        IKnotKernel kernel = kernelName == "uniform" ? new UniformKernel() : new LocalBetaKernel(5.0);

        IRjMove<KnotConfig>[] moves = dmgk
            ? DmgkSchedule.Wrap(new KnotBirthMove(kernel), new KnotDeathMove(kernel), new KnotRelocateMove(kernel), prior)
            : [new KnotBirthMove(kernel), new KnotDeathMove(kernel), new KnotRelocateMove(kernel)];
        var target = new SplineTarget(new SplineBasis(3), new FlatModel(), prior, x: [0.5], y: [0.0]);

        CountLaw.Result r = CountLaw.Sample(moves, target, new KnotConfig([]), s => s.Count, prior.LogPrior,
                                            bins, burn: 10_000, steps: 400_000, batches: 40, seed: 31);
        string label = $"{priorName} dmgk={dmgk} kernel={kernelName}";
        _out.WriteLine($"[count-law] {label}:{r.Describe()}");

        Assert.True(r.Agrees(z: 5.0, floor: 0.002), $"{label}:{r.Describe()}");
    }

    [Fact]
    public void ConwayMaxwellPoisson_AtNuOne_IsPoisson()
    {
        var cmp = new ConwayMaxwellPoissonPrior(4.0, 1.0);
        var poisson = new PoissonPrior(4.0);
        double offset = cmp.LogPrior(0) - poisson.LogPrior(0);
        for (int k = 0; k <= 30; k++) Assert.Equal(poisson.LogPrior(k) + offset, cmp.LogPrior(k), 10);
        Assert.Equal(double.NegativeInfinity, cmp.LogPrior(-1));
    }

    [Fact]
    public void ConwayMaxwellPoisson_RejectsAnImproperGeometric()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new ConwayMaxwellPoissonPrior(1.0, 0.0));
        _ = new ConwayMaxwellPoissonPrior(0.5, 0.0);
    }

    [Fact]
    public void Target_IsMinusInfinityOffTheOpenOrderedSimplex()
    {
        var target = new SplineTarget(new SplineBasis(3), new FlatModel(), new PoissonPrior(4.0), x: [0.5], y: [0.0]);
        Assert.Equal(double.NegativeInfinity, target.LogPosterior(new KnotConfig([0.0, 0.5])));
        Assert.Equal(double.NegativeInfinity, target.LogPosterior(new KnotConfig([0.5, 1.0])));
        Assert.Equal(double.NegativeInfinity, target.LogPosterior(new KnotConfig([0.3, 0.3])));
        Assert.Equal(double.NegativeInfinity, target.LogPosterior(new KnotConfig([0.6, 0.3])));
        Assert.True(double.IsFinite(target.LogPosterior(new KnotConfig([0.3, 0.6]))));
    }

    // A constant likelihood: the target reduces to the count prior and the location measure.
    private sealed class FlatModel : IObservationModel
    {
        public double LogMarginalLikelihood(double[,] design, double[] y, double[]? weights) => 0.0;
        public int EffectiveDimension(double[,] design, double[]? weights) => design.GetLength(1);
        public double[] PosteriorMeanCoefficients(double[,] design, double[] y, double[]? weights)
            => new double[design.GetLength(1)];
    }
}
