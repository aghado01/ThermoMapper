using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// R01 / Q3: BARS's knot measure. Knot locations are the order statistics of k uniforms (density k! on the ordered
/// simplex), so with the likelihood held constant the knot count must follow the count prior p(k) exactly. Checked
/// under both move schedules, constant weights and DMGK's, which settles that the two sample the same posterior;
/// for a hard cap; for the Conway–Maxwell–Poisson tails; and for the EBIC law on both sides of γ = 1, including its
/// rise past C/2 when γ &gt; 1, which its support [0, C] keeps. The kernels mix within the batches: under a flat
/// target the default LocalBeta(50) wanders slowly enough that 40 batches understate its error. These targets put one
/// data point under the basis, so the support gate is off: they check the bare knot prior. With the gate on, the law
/// is the prior restricted to supported configurations, checked against direct sampling.
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
    [InlineData("ebic-0.5", false, "beta")]
    [InlineData("ebic-2", false, "beta")]
    [InlineData("ebic-2", true, "beta")]
    public void KnotCount_FollowsTheCountPrior(string priorName, bool dmgk, string kernelName)
    {
        (IComplexityPrior prior, int bins) = priorName switch
        {
            "poisson" => ((IComplexityPrior)new PoissonPrior(4.0), 12),
            "uniform" => (new UniformPrior(6), 9),
            "cmp-nu2" => (new ConwayMaxwellPoissonPrior(4.0, 2.0), 7),
            "cmp-nu0.5" => (new ConwayMaxwellPoissonPrior(1.5, 0.5), 12),
            "ebic-0.5" => (new EbicPrior(0.5, 12), 14),   // peaks at C/2 = 6; zero past C = 12
            "ebic-2" => (new EbicPrior(2.0, 6), 8),       // U-shaped: falls to C/2 = 3, rises to C = 6, zero past it
            _ => throw new ArgumentOutOfRangeException(nameof(priorName)),
        };
        IKnotKernel kernel = kernelName == "uniform" ? new UniformKernel() : new LocalBetaKernel(5.0);

        IRjMove<KnotConfig>[] moves = dmgk
            ? DmgkSchedule.Wrap(new KnotBirthMove(kernel), new KnotDeathMove(kernel), new KnotRelocateMove(kernel), prior)
            : [new KnotBirthMove(kernel), new KnotDeathMove(kernel), new KnotRelocateMove(kernel)];
        var target = new SplineTarget(new SplineBasis(3), new FlatModel(), prior, x: [0.5], y: [0.0], minColumnMass: 0.0);

        CountLaw.Result r = CountLaw.Sample(moves, target, new KnotConfig([]), s => s.Count, prior.LogPrior,
                                            bins, burn: 10_000, steps: 400_000, batches: 40, seed: 31);
        string label = $"{priorName} dmgk={dmgk} kernel={kernelName}";
        _out.WriteLine($"[count-law] {label}:{r.Describe()}");

        Assert.True(r.Agrees(z: 5.0, floor: 0.002), $"{label}:{r.Describe()}");
    }

    [Theory]
    [InlineData(0.1, "uniform")]
    [InlineData(0.1, "design")]
    [InlineData(1.0, "design")]
    public void MixtureBirths_KeepTheCountLaw(double exploration, string baseName)
    {
        // Births with an exploratory share ε from a base measure (step 3.9): any proposal density leaves the target alone
        // when the death carries it in reverse. The design base is a histogram of 40 points packed into [0.6, 0.9], with
        // a uniform floor: deliberately uneven. The local kernel is the default τ = 50.
        var prior = new PoissonPrior(4.0);
        double[] packed = new double[40];
        for (int i = 0; i < packed.Length; i++) packed[i] = 0.6 + 0.3 * (i + 0.5) / packed.Length;
        IKnotBase knotBase = baseName == "design" ? new DesignHistogramBase(packed) : UniformKnotBase.Instance;
        var density = new KnotBirthDensity(new LocalBetaKernel(50.0), exploration, knotBase);
        var target = new SplineTarget(new SplineBasis(3), new FlatModel(), prior, x: [0.5], y: [0.0], minColumnMass: 0.0);

        CountLaw.Result r = CountLaw.Sample<KnotConfig>(
            [new KnotBirthMove(density), new KnotDeathMove(density), new KnotRelocateMove(density.Kernel)], target,
            new KnotConfig([]), s => s.Count, prior.LogPrior, 12, burn: 10_000, steps: 400_000, batches: 40, seed: 31);
        string label = $"ε={exploration} base={baseName}";
        _out.WriteLine($"[count-law-mixture] {label}:{r.Describe()}");
        Assert.True(r.Agrees(z: 5.0, floor: 0.002), $"{label}:{r.Describe()}");
    }

    [Fact]
    public void SupportGate_RestrictsTheKnotPrior_AndNothingElse()
    {
        // Under a constant likelihood with the gate on, the target is the knot prior restricted to configurations whose
        // every basis function has column mass ≥ 1. The knot count's law is then p(k)·s(k)/Z, where s(k) is the chance
        // that k sorted uniforms pass the gate: a prior-support restriction leaves the knot measure alone. s(k) is
        // estimated by direct sampling, apart from the chain.
        var prior = new PoissonPrior(4.0);
        var basis = new SplineBasis(3);
        double[] x = new double[30];
        for (int i = 0; i < x.Length; i++) x[i] = (i + 0.5) / x.Length;
        var target = new SplineTarget(basis, new FlatModel(), prior, x, new double[x.Length]);
        Assert.Equal(SplineTarget.DefaultMinColumnMass, target.MinColumnMass);

        const int bins = 12, draws = 50_000;
        var rng = new Xoshiro256PlusPlus(seed: 23);
        var logPass = new double[40];
        for (int k = 0; k < logPass.Length; k++)
        {
            int pass = 0;
            for (int d = 0; d < draws; d++)
            {
                var knots = new double[k];
                for (int i = 0; i < k; i++) knots[i] = rng.NextDouble();
                Array.Sort(knots);
                double[,] z = basis.Design(new KnotConfig(knots), x);
                bool supported = true;
                for (int j = 0; j < z.GetLength(1) && supported; j++)
                {
                    double mass = 0.0;
                    for (int i = 0; i < x.Length; i++) mass += z[i, j];
                    supported = mass >= 1.0;
                }
                if (supported) pass++;
            }
            logPass[k] = Math.Log((double)pass / draws);
        }
        _out.WriteLine($"[count-law-gate] s(k) for k = 0..8: {string.Join(", ", Array.ConvertAll(logPass[..9], v => Math.Exp(v).ToString("F3")))}");

        IKnotKernel kernel = new LocalBetaKernel(5.0);
        CountLaw.Result r = CountLaw.Sample<KnotConfig>(
            [new KnotBirthMove(kernel), new KnotDeathMove(kernel), new KnotRelocateMove(kernel)], target, new KnotConfig([]),
            s => s.Count, k => prior.LogPrior(k) + (k < logPass.Length ? logPass[k] : double.NegativeInfinity),
            bins, burn: 10_000, steps: 400_000, batches: 40, seed: 31);
        _out.WriteLine($"[count-law-gate]{r.Describe()}");
        Assert.True(r.Agrees(z: 5.0, floor: 0.003), r.Describe());
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
        var target = new SplineTarget(new SplineBasis(3), new FlatModel(), new PoissonPrior(4.0), x: [0.5], y: [0.0], minColumnMass: 0.0);
        Assert.Equal(double.NegativeInfinity, target.LogPosterior(new KnotConfig([0.0, 0.5])));
        Assert.Equal(double.NegativeInfinity, target.LogPosterior(new KnotConfig([0.5, 1.0])));
        Assert.Equal(double.NegativeInfinity, target.LogPosterior(new KnotConfig([0.3, 0.3])));
        Assert.Equal(double.NegativeInfinity, target.LogPosterior(new KnotConfig([0.6, 0.3])));
        Assert.True(double.IsFinite(target.LogPosterior(new KnotConfig([0.3, 0.6]))));
    }

    [Theory]
    [InlineData(TemperingPath.Likelihood)]
    [InlineData(TemperingPath.Powered)]
    public void TemperedRungs_FollowTheirPathsCountLaw(TemperingPath path)
    {
        // Under a constant likelihood the likelihood path's every rung samples the knot prior itself. The powered path's
        // rung β samples (p(k)·k!)^β on the ordered simplex, whose volume is 1/k!, so for a Poisson(λ) count prior its
        // knot count is Poisson(λ^β): a powered density on a space of varying dimension depends on the reference measure
        // of each subspace, and its hot rungs leave the prior (Q5).
        const double lambda = 4.0;
        var prior = new PoissonPrior(lambda);
        IKnotKernel kernel = new LocalBetaKernel(5.0);
        double[] betas = [1.0, 0.5, 0.25];
        var target = new SplineTarget(new SplineBasis(3), new FlatModel(), prior, x: [0.5], y: [0.0], minColumnMass: 0.0);
        Func<int, double> LawAt(int rung) => path == TemperingPath.Likelihood
            ? prior.LogPrior
            : new PoissonPrior(Math.Pow(lambda, betas[rung])).LogPrior;

        CountLaw.Result[] rungs = CountLaw.SampleLadder(
            _ => [new KnotBirthMove(kernel), new KnotDeathMove(kernel), new KnotRelocateMove(kernel)], target,
            new KnotConfig([]), s => s.Count, betas, path, LawAt, bins: 10, burn: 10_000, steps: 300_000, batches: 30,
            seed: 41);

        for (int r = 0; r < rungs.Length; r++)
        {
            string label = $"{path} β={betas[r]}";
            _out.WriteLine($"[count-law-ladder] {label}:{rungs[r].Describe()}");
            Assert.True(rungs[r].Agrees(z: 5.0, floor: 0.002), $"{label}:{rungs[r].Describe()}");
        }
    }

    [Theory]
    [InlineData(TemperingPath.Likelihood)]
    [InlineData(TemperingPath.Powered)]
    public void ReferenceRung_FollowsItsPathsCountLaw(TemperingPath path)
    {
        // The exact reference rung at β = 0 (OQ2): on the likelihood path the knot prior, Poisson(λ); on the powered path
        // the normalized reference measure, Poisson(1) = Poisson(λ⁰), with uniform locations on both. Every rung above it
        // keeps its path's law.
        const double lambda = 4.0;
        var prior = new PoissonPrior(lambda);
        IKnotKernel kernel = new LocalBetaKernel(5.0);
        double[] betas = [1.0, 0.5, 0.25, 0.0];
        var target = new SplineTarget(new SplineBasis(3), new FlatModel(), prior, x: [0.5], y: [0.0], minColumnMass: 0.0);
        Func<int, double> LawAt(int rung) => path == TemperingPath.Likelihood
            ? prior.LogPrior
            : new PoissonPrior(Math.Pow(lambda, betas[rung])).LogPrior;

        CountLaw.Result[] rungs = CountLaw.SampleLadder(
            _ => [new KnotBirthMove(kernel), new KnotDeathMove(kernel), new KnotRelocateMove(kernel)], target,
            new KnotConfig([]), s => s.Count, betas, path, LawAt, bins: 10, burn: 10_000, steps: 300_000, batches: 30,
            seed: 43, reference: target.Reference(path));

        for (int r = 0; r < rungs.Length; r++)
        {
            string label = $"{path} β={betas[r]}";
            _out.WriteLine($"[count-law-reference] {label}:{rungs[r].Describe()}");
            Assert.True(rungs[r].Agrees(z: 5.0, floor: 0.002), $"{label}:{rungs[r].Describe()}");
        }
    }

    [Fact]
    public void Reference_DrawsOnlyAdmittedConfigurations()
    {
        // With the support gate on, the reference is the b = 0 law restricted by rejection: every draw is admitted.
        double[] x = new double[30];
        for (int i = 0; i < x.Length; i++) x[i] = (i + 0.5) / x.Length;
        var target = new SplineTarget(new SplineBasis(3), new FlatModel(), new PoissonPrior(4.0), x, new double[x.Length]);
        var rng = new Xoshiro256PlusPlus(seed: 29);
        foreach (TemperingPath path in new[] { TemperingPath.Likelihood, TemperingPath.Powered })
        {
            IReferenceSampler<KnotConfig> reference = target.Reference(path);
            Assert.Equal(path, reference.Path);
            for (int d = 0; d < 2_000; d++)
            {
                (KnotConfig config, ITargetEvaluation e) = reference.Draw(rng);
                Assert.True(double.IsFinite(new Temperature(0.0, path).LogDensity(e)));
                Assert.True(double.IsFinite(target.LogPosterior(config)));
            }
        }
        var normal = new SplineTarget(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(4.0), x, x);
        Assert.Throws<InvalidOperationException>(() => normal.Reference(TemperingPath.Likelihood));
        _ = normal.Reference(TemperingPath.Powered);
    }

    [Fact]
    public void Target_OffTheCountPriorsSupport_BuildsNoDesign()
    {
        // A zero count prior decides the evaluation before any design is built or factored.
        var target = new SplineTarget(new SplineBasis(3), new ThrowingModel(), new UniformPrior(1), x: [0.5], y: [0.0], minColumnMass: 0.0);
        ITargetEvaluation e = target.Evaluate(new KnotConfig([0.3, 0.6]));
        Assert.Equal(double.NegativeInfinity, e.LogPrior);
        Assert.Equal(double.NegativeInfinity, e.LogLikelihood(1.0));
        Assert.Throws<InvalidOperationException>(() => target.Evaluate(new KnotConfig([0.3])));
    }

    private sealed class ThrowingModel : IObservationModel
    {
        public double LogMarginalLikelihood(double[,] design, double[] y, double[]? weights)
            => throw new InvalidOperationException("The model was consulted.");
        public int EffectiveDimension(double[,] design, double[]? weights) => design.GetLength(1);
        public double[] PosteriorMeanCoefficients(double[,] design, double[] y, double[]? weights)
            => new double[design.GetLength(1)];
    }

    // A constant likelihood at every temperature: the target reduces to the count prior and the location measure.
    private sealed class FlatModel : IObservationModel
    {
        public double LogMarginalLikelihood(double[,] design, double[] y, double[]? weights) => 0.0;
        public bool CanTemperLikelihood => true;
        public bool CanTemperToPrior => true;
        public ITemperedMarginal TemperedMarginal(double[,] design, double[] y, double[]? weights) => Constant.Zero;
        public int EffectiveDimension(double[,] design, double[]? weights) => design.GetLength(1);
        public double[] PosteriorMeanCoefficients(double[,] design, double[] y, double[]? weights)
            => new double[design.GetLength(1)];
    }

    private sealed class Constant : ITemperedMarginal
    {
        public static readonly Constant Zero = new();
        public double LogMarginalLikelihood(double inverseTemperature) => 0.0;
    }
}
