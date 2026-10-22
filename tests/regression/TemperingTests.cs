using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers;
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
