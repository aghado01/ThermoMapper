using System;
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Tempering inside BARS: each "chain" of the ensemble a parallel-tempering ladder whose cold replica the readouts
/// read unchanged, and its incompatibility with weight resampling. The ladder's own behaviour is checked in
/// Maths.Samplers.Tests.
/// </summary>
public sealed class TemperingTests
{
    private readonly ITestOutputHelper _out;
    public TemperingTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void TemperedEnsemble_RecoversPeak_ViaColdReplica()
    {
        var rng = new Xoshiro256PlusPlus(seed: 11);
        int n = 150;
        const double truePeak = 0.4;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Exp(-120.0 * (x[i] - truePeak) * (x[i] - truePeak)) + 0.05 * Gaussian(rng);
        }

        // Each "chain" is a 5-replica tempering ladder; the readout reads its cold (β = 1) replica unchanged.
        BarsResult result = new BarsEnsemble(
                new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0),
                temperLevels: 5, temperBetaMin: 0.05)
            .Run(x, y, grid: x, chains: 4, masterSeed: 3, burn: 1500, samples: 2000);
        PeakPosterior p = result.Peak;

        _out.WriteLine($"[tempered-ens] peak={p.LocationMean:F3} [{p.LocationLo:F3},{p.LocationHi:F3}] " +
                       $"R̂={p.LocationConvergence.RHat:F3} ESS={p.LocationConvergence.BulkEss:F0} coldAcc={result.AcceptanceRate:F3}");

        Assert.True(Math.Abs(p.LocationMean - truePeak) < 0.06, $"tempered peak {p.LocationMean:F3} off {truePeak}");
        Assert.True(p.LocationLo <= truePeak && truePeak <= p.LocationHi, "95% interval should bracket the peak");
        Assert.InRange(result.AcceptanceRate, 1e-6, 1.0);   // cold replica genuinely mixed
    }

    [Fact]
    public void TemperedRun_ReportsItsLadders()
    {
        var rng = new Xoshiro256PlusPlus(seed: 21);
        int n = 100;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.15 * Gaussian(rng);
        }
        BarsEnsemble Ensemble(bool adapt) => new(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0),
                                                  new LocalBetaKernel(50.0), temperLevels: 4, temperAdaptLadder: adapt);

        BarsResult adapted = Ensemble(adapt: true).Run(x, y, grid: x, chains: 2, masterSeed: 4, burn: 600, samples: 600);
        TemperingReadout t = adapted.Tempering!;
        _out.WriteLine($"[ladder] trips/scan {t.RoundTripsPerScan:F4} (ideal {t.IdealRoundTripRate:F4}), Λ̂ {t.Barrier:F3}, " +
                       $"rejection [{string.Join(", ", Array.ConvertAll(t.PairRejection, r => r.ToString("F3")))}]");
        Assert.Equal(TemperingPath.Likelihood, t.Path);
        Assert.Equal(SwapSchedule.DeterministicEvenOdd, t.Schedule);
        Assert.Equal(3, t.PairRejection.Length);
        Assert.All(t.PairRejection, r => Assert.InRange(r, 0.0, 1.0));
        Assert.Equal(t.PairRejection[0] + t.PairRejection[1] + t.PairRejection[2], t.Barrier, 12);
        Assert.InRange(t.RoundTripsPerScan, 0.0, t.IdealRoundTripRate);   // real local moves fall short of exact draws
        Assert.Equal(2, t.Ladders.Length);
        double[] geometric = ParallelTempering<KnotConfig>.GeometricLadder(4, 0.05);
        foreach (double[] ladder in t.Ladders)
        {
            Assert.Equal(1.0, ladder[0]);
            Assert.Equal(0.05, ladder[^1], 15);
            Assert.NotEqual(geometric, ladder);   // the interior rungs moved in burn-in
        }

        TemperingReadout fixedLadder = Ensemble(adapt: false).Run(x, y, grid: x, chains: 2, masterSeed: 4, burn: 300, samples: 300).Tempering!;
        Assert.All(fixedLadder.Ladders, ladder => Assert.Equal(geometric, ladder));

        Assert.Null(new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0))
            .Run(x, y, grid: x, chains: 2, masterSeed: 4, burn: 100, samples: 100).Tempering);
    }

    [Theory]
    [InlineData(TemperingPath.Powered, TemperedForm.Integrated, false)]
    [InlineData(TemperingPath.Likelihood, TemperedForm.Linear, false)]
    [InlineData(TemperingPath.Likelihood, TemperedForm.Integrated, true)]
    public void ReferenceRungLadder_RecoversThePeak(TemperingPath path, TemperedForm form, bool proper)
    {
        // Each ladder ends at an exact reference rung (OQ2): the powered path's Poisson(1) reference for the default
        // Normal model, or the knot prior on the likelihood path, for the Normal model's linear form, or for its integrated
        // form under the proper scaled σ² prior.
        var rng = new Xoshiro256PlusPlus(seed: 11);
        int n = 150;
        const double truePeak = 0.4;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Exp(-120.0 * (x[i] - truePeak) * (x[i] - truePeak)) + 0.05 * Gaussian(rng);
        }
        BarsResult result = new BarsEnsemble(
                new SplineBasis(3), new WeightedNormalModel(form: form, sigma2Prior: proper ? Sigma2Prior.ScaledToResponse() : null),
                new PoissonPrior(5.0), new LocalBetaKernel(50.0),
                temperLevels: 5, temperBetaMin: 0.1, temperPath: path, temperReferenceRung: true)
            .Run(x, y, grid: x, chains: 4, masterSeed: 3, burn: 1500, samples: 2000);
        PeakPosterior p = result.Peak;
        TemperingReadout t = result.Tempering!;
        _out.WriteLine($"[reference-rung] {path}/{form}{(proper ? "/proper" : "")}: peak={p.LocationMean:F3} [{p.LocationLo:F3},{p.LocationHi:F3}] " +
                       $"knots={result.MeanKnots:F2} trips/scan {t.RoundTripsPerScan:F4} (ideal {t.IdealRoundTripRate:F4}), Λ̂ {t.Barrier:F3}");

        Assert.True(Math.Abs(p.LocationMean - truePeak) < 0.06, $"peak {p.LocationMean:F3} off {truePeak}");
        Assert.All(t.Ladders, ladder => Assert.Equal(0.0, ladder[^1]));
        Assert.True(t.RoundTripsPerScan > 0.0);
    }

    [Fact]
    public void ReferenceRung_IsRejected_WhereItCannotBeDrawn()
    {
        IBasis basis = new SplineBasis(3);
        var prior = new PoissonPrior(5.0);
        var kernel = new LocalBetaKernel(50.0);
        Assert.Throws<ArgumentException>(() => new BarsEnsemble(basis, new WeightedNormalModel(), prior, kernel,
            temperLevels: 4, temperReferenceRung: true));                       // the integrated form refuses b = 0
        Assert.Throws<ArgumentOutOfRangeException>(() => new BarsEnsemble(basis, new WeightedNormalModel(), prior, kernel,
            temperLevels: 2, temperPath: TemperingPath.Powered, temperReferenceRung: true));
        Assert.Throws<ArgumentException>(() => new BarsEnsemble(basis, new WeightedNormalModel(), prior, kernel, temperReferenceRung: true));
        _ = new BarsEnsemble(basis, new WeightedNormalModel(form: TemperedForm.Linear), prior, kernel, temperLevels: 4, temperReferenceRung: true);
        _ = new BarsEnsemble(basis, new WeightedNormalModel(sigma2Prior: Sigma2Prior.ScaledToResponse()), prior, kernel, temperLevels: 4,
                             temperReferenceRung: true);                       // the integrated form under a proper σ² prior
        _ = new BarsEnsemble(basis, new PoissonModel(), prior, kernel, temperLevels: 4, temperReferenceRung: true);
    }

    [Fact]
    public void BimodalFixture_IsMirrorSymmetric()
    {
        // The fixture's ground truth, each mode at exactly half the mass, rests on the target's invariance under
        // ξ ↦ 1 − ξ: a configuration and its mirror evaluate alike in both tempered forms, and the mirror's curve peaks
        // at the mirrored location, so the peak's side of ½ has a symmetric law.
        (double[] x, double[] y) = BimodalFixture.Data(200);
        var basis = new SplineBasis(3);
        var rng = new Xoshiro256PlusPlus(seed: 5);
        foreach (TemperedForm form in new[] { TemperedForm.Integrated, TemperedForm.Linear })
        {
            var model = new WeightedNormalModel(form: form);
            var target = new SplineTarget(basis, model, BimodalFixture.Prior(), x, y);
            int admitted = 0;
            for (int trial = 0; trial < 300; trial++)
            {
                var knots = new double[1 + trial % BimodalFixture.Cap];
                for (int j = 0; j < knots.Length; j++) knots[j] = rng.NextDouble();
                Array.Sort(knots);
                var config = new KnotConfig(knots);
                KnotConfig mirror = BimodalFixture.Mirror(config);
                ITargetEvaluation e = target.Evaluate(config), m = target.Evaluate(mirror);
                Assert.Equal(e.LogPrior, m.LogPrior);
                if (double.IsNegativeInfinity(e.LogLikelihood(1.0)))
                {
                    Assert.True(double.IsNegativeInfinity(m.LogLikelihood(1.0)));
                    continue;
                }
                admitted++;
                foreach (double b in new[] { 1.0, 0.3 })
                    Assert.Equal(e.LogLikelihood(b), m.LogLikelihood(b), 1e-9 * Math.Abs(e.LogLikelihood(b)));

                double[] coef = model.PosteriorMeanCoefficients(basis.Design(config, x), y, null);
                double[] mirrorCoef = model.PosteriorMeanCoefficients(basis.Design(mirror, x), y, null);
                (double location, double height) = SplineExtrema.Argmax(config, coef, basis);
                double[,] at = basis.Design(mirror, [1.0 - location]);
                double mirrored = 0.0;
                for (int j = 0; j < mirrorCoef.Length; j++) mirrored += at[0, j] * mirrorCoef[j];
                Assert.Equal(height, SplineExtrema.Argmax(mirror, mirrorCoef, basis).Height, 1e-9);
                Assert.Equal(height, mirrored, 1e-9);
            }
            Assert.True(admitted > 150, $"only {admitted} configurations admitted");
        }
    }

    [Fact]
    public void LikelihoodPath_IsRejected_ForAModelWithoutATemperedMarginal()
    {
        Assert.Throws<ArgumentException>(() => new BarsEnsemble(
            new SplineBasis(3), new ConstantModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0), temperLevels: 3));
        _ = new BarsEnsemble(new SplineBasis(3), new ConstantModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0),
                             temperLevels: 3, temperPath: TemperingPath.Powered);
    }

    [Fact]
    public void PoissonLadder_OnTheLikelihoodPath_RecoversThePeak()
    {
        // Counts at the noiseless intensity exp(1 + 1.5·exp(−80(x − 0.4)²)); the Laplace marginal tempered linearly.
        int n = 120;
        const double truePeak = 0.4;
        var x = new double[n];
        var counts = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            counts[i] = Math.Round(Math.Exp(1.0 + 1.5 * Math.Exp(-80.0 * (x[i] - truePeak) * (x[i] - truePeak))));
        }

        BarsResult result = new BarsEnsemble(
                new SplineBasis(3), new PoissonModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0), temperLevels: 3)
            .Run(x, counts, grid: x, chains: 2, masterSeed: 5, burn: 300, samples: 500);
        PeakPosterior p = result.Peak;

        _out.WriteLine($"[tempered-poisson] peak={p.LocationMean:F3} [{p.LocationLo:F3},{p.LocationHi:F3}] knots={result.MeanKnots:F2}");
        Assert.True(Math.Abs(p.LocationMean - truePeak) < 0.05, $"tempered Poisson peak {p.LocationMean:F3} off {truePeak}");
    }

    // A model that implements only the required members, so it has no tempered marginal.
    private sealed class ConstantModel : IObservationModel
    {
        public double LogMarginalLikelihood(double[,] design, double[] y, double[]? weights) => 0.0;
        public int EffectiveDimension(double[,] design, double[]? weights) => design.GetLength(1);
        public double[] PosteriorMeanCoefficients(double[,] design, double[] y, double[]? weights)
            => new double[design.GetLength(1)];
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(1.0)]
    public void HotEnd_OutsideTheOpenUnitInterval_IsRejected(double betaMin)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new BarsEnsemble(
            new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0),
            temperLevels: 3, temperBetaMin: betaMin));
    }

    [Fact]
    public void Tempering_IsRejected_AlongsideResampler()
    {
        Assert.Throws<ArgumentException>(() => new BarsEnsemble(
            new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0),
            resampler: new StudentTWeights(4.0), temperLevels: 5));
    }
}
