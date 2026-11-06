using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The generalized beta prime on the BASS core (E1's slice 4), the horseshoe and the Strawderman–Berger: their declaration; the latent block,
/// whose slice step moves each v from its current value, under σ² integrated and under a drawn skew with σ² on the record; σ² under the
/// half-Cauchy, <c>hbass</c>'s prior; retention and resume bit for bit; and the readers, the pointwise and candidate densities refused until
/// opted into and then the law's integral, the one the other's column at the training cases, the predictive's quantiles from the mixture's
/// table, and the moments the horseshoe lacks.
/// </summary>
public sealed class BassGeneralizedBetaPrimeTests(ITestOutputHelper output)
{
    private static readonly BassNoise Skewed = BassNoise.Horseshoe(new BassNoisePrior(0.0, 0.5));

    [Fact]
    public void Config_DeclaresTheKind()
    {
        BassConfig horseshoe = BassStudentTTests.Config(BassNoise.Horseshoe());
        Assert.Contains("\"Noise\":{\"Kind\":\"GeneralizedBetaPrime\",\"VarianceScale\":1,\"Reference\":\"Moments\",\"BetaPrimePrior\":{\"P\":1,\"A\":0.5,\"B\":0.5}}",
                        horseshoe.ToJson());
        Assert.Equal(horseshoe, BassConfig.FromJson(horseshoe.ToJson()));
        BassConfig sb = BassStudentTTests.Config(BassNoise.StrawdermanBerger(new BassNoisePrior(0.2, 0.5)));
        Assert.Contains("\"BetaPrimePrior\":{\"P\":0.5,\"A\":1,\"B\":0.5}", sb.ToJson());
        Assert.Equal(sb, BassConfig.FromJson(sb.ToJson()));
        Assert.IsType<GbpMixing>(BassNoise.Horseshoe().Mixing);
        Assert.False(BassNoise.Horseshoe().Skewed);
        Assert.True(Skewed.Skewed && Skewed.SkewDrawn);
        // Another kind's config does not name the GBP.
        Assert.DoesNotContain("BetaPrimePrior", BassStudentTTests.Config(BassNoise.NormalWald()).ToJson());

        Assert.Throws<ArgumentException>(() => new BassNoise(BassNoiseKind.GeneralizedBetaPrime));
        Assert.Throws<ArgumentException>(() => new BassNoise(BassNoiseKind.NormalWald, betaPrimePrior: GbpParameters.Horseshoe));
        Assert.Throws<ArgumentOutOfRangeException>(() => BassNoise.GeneralizedBetaPrime(new GbpParameters(0.0, 0.5, 0.5)));
        Assert.Throws<ArgumentOutOfRangeException>(() => BassNoise.GeneralizedBetaPrime(new GbpParameters(1.0, 0.5, double.PositiveInfinity)));
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void ASweep_StepsEachV_FromItsCurrentValue(bool skewDrawn)
    {
        BassData data = BassQuantileTests.Data();
        BassNoise noise = skewDrawn ? Skewed : BassNoise.Horseshoe();
        BassModel model = BassModel.Create(BassStudentTTests.Config(noise) with { Sigma2 = Sigma2Prior.HalfCauchy() }, data);
        Assert.NotNull(model.Target.Sigma2Rate);
        var start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(1)).State;
        var chain = new BassChain(model, start.Terms, start.LogTau, new Xoshiro256PlusPlus(2), new Xoshiro256PlusPlus(3));
        string[] blocks = skewDrawn ? ["moves", "tau", "coefficients", "latent scales", "scale", "noise shape", "sigma2 rate"]
                                    : ["moves", "tau", "coefficients", "latent scales", "sigma2 rate"];
        Assert.Equal(blocks, chain.Scan.BlockStats().Select(s => s.Key));
        double[] previous = chain.Current.Fit.Working!.Latents!.ValueArray;
        int moved = 0, total = 0;
        var betas = new List<double>();
        for (int k = 0; k < 150; k++)
        {
            chain.Step();
            ObservationLatents latents = chain.Current.Fit.Working!.Latents!;
            for (int i = 0; i < latents.Values.Length; i++)
            {
                Assert.True(latents.Values[i] > 0.0 && double.IsFinite(latents.Values[i]));
                moved += latents.Values[i] != previous[i] ? 1 : 0;
                total++;
            }
            previous = latents.ValueArray;
            betas.Add(latents.Skew);
            if (skewDrawn) Assert.NotNull(latents.Sigma2);
            else Assert.Null(latents.Sigma2);
        }
        Assert.True(moved > 0.95 * total, $"v moved {moved} of {total}");
        if (skewDrawn) Assert.True(betas.Distinct().Count() > 100, "β moves");
        foreach (MoveStat s in chain.Scan.BlockStats().Skip(2)) output.WriteLine($"{s.Key}: {s.Accepted} of {s.Attempts}");
    }

    public static IEnumerable<object[]> Noises() =>
    [
        [BassNoise.Horseshoe()],
        [Skewed],
        [BassNoise.StrawdermanBerger()],
    ];

    [Theory]
    [MemberData(nameof(Noises))]
    public void ARetainedDraw_RebuildsBitForBit(BassNoise noise)
    {
        BassData data = BassQuantileTests.Data();
        BassModel model = BassModel.Create(BassStudentTTests.Config(noise), data);
        var start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(1)).State;
        var chain = new BassChain(model, start.Terms, start.LogTau, new Xoshiro256PlusPlus(2), new Xoshiro256PlusPlus(3));
        var terms = new BassTermTable();
        var compose = new Xoshiro256PlusPlus(4);
        for (int k = 0; k < 40; k++)
        {
            chain.Step();
            BassDraw draw = chain.Compose(compose, terms);
            WorkingResponse held = chain.Current.Fit.Working!;
            Assert.Same(held.Latents, draw.Latents);
            WorkingResponse rebuilt = draw.WorkingOn(model);
            Assert.True(rebuilt.Values.SequenceEqual(held.Values));
            Assert.True(rebuilt.Weights.SequenceEqual(held.Weights));
            Assert.Equal(held.Latents!.Skew, rebuilt.Latents!.Skew);
        }
    }

    [Fact]
    public void AHorseshoeFit_ResumesBitForBit()
    {
        BassData data = BassQuantileTests.Data();
        BassConfig whole = BassStudentTTests.Config(Skewed, samples: 60) with { Sigma2 = Sigma2Prior.HalfCauchy() };
        BassConfig half = whole with { Sampling = whole.Sampling with { Samples = 30 } };
        BassFit a = Bass.Fit(whole, data);
        BassFit first = Bass.Fit(half, data, new BassRunOptions { Checkpoint = true });
        BassFit b = Bass.Resume(first.Checkpoint!, half, data);
        BassRunTests.AssertSameFit(a, b);
        for (int c = 0; c < a.Draws.Count; c++)
            for (int d = 0; d < a.Draws[c].Count; d++)
            {
                ObservationLatents la = a.Draws[c][d].Latents!, lb = b.Draws[c][d].Latents!;
                Assert.True(la.Values.SequenceEqual(lb.Values));
                Assert.Equal(la.Sigma2, lb.Sigma2);
                Assert.Equal(la.Skew, lb.Skew);
                Assert.Equal(a.Draws[c][d].Sigma2Rate, b.Draws[c][d].Sigma2Rate);
            }
    }

    [Fact]
    public void TheReaders_TakeTheHorseshoesLaw_ThePointwiseOptedInto()
    {
        BassData data = BassQuantileTests.Data();
        BassNoise noise = BassNoise.Horseshoe();
        BassFit fit = Bass.Fit(BassStudentTTests.Config(noise) with { Sigma2 = Sigma2Prior.HalfCauchy() }, data);
        var refused = Assert.Throws<NotSupportedException>(() => fit.PointwiseLogDensity(data));
        Assert.Contains("quadrature: true", refused.Message);
        PointwiseLogDensity pointwise = fit.PointwiseLogDensity(data, quadrature: true);
        BassPrediction p = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Seed = 2 });
        ReadOnlySpan<double> y = data.Response.Values;
        List<BassDraw> draws = [.. fit.Draws.SelectMany(chain => chain)];
        IBassNoiseLaw law = new VarianceMeanMixture(noise);
        Assert.True(law.DensityByQuadrature);
        Assert.False(new VarianceMeanMixture(BassNoise.NormalWald()).Law.DensityByQuadrature);
        var mixing = new GbpMixing(GbpParameters.Horseshoe);
        for (int d = 0; d < p.DrawCount; d += 11)
        {
            double sigma2 = draws[d].Sigma2;
            Assert.NotNull(draws[d].Sigma2Rate);
            Assert.Equal(sigma2, p.NoiseVariances[d]);
            Assert.Equal(0.0, p.MeanShifts[d]);
            Assert.Equal(double.PositiveInfinity, law.Variance(sigma2, draws[d].Latents));
            for (int i = 0; i < y.Length; i += 9)
                Assert.Equal(mixing.MarginalLogDensity(y[i] - p.Latent[d, i, 0], sigma2, 0.0), pointwise.Values[d, i], 10);
        }
        // The candidate density takes the same opt-in, and at the training cases it is the pointwise column.
        double[][] grids = [.. y.ToArray().Select(v => new[] { v })];
        Assert.Throws<NotSupportedException>(() => fit.CandidateLogDensity(data, PredictionCases.Of(data), grids));
        CandidateLogDensity candidates = fit.CandidateLogDensity(data, PredictionCases.Of(data), grids, quadrature: true);
        for (int s = 0; s < pointwise.Draws; s += 13)
            for (int i = 0; i < y.Length; i += 7)
                Assert.Equal(pointwise.Values[s, i], candidates.Values[s, i], 12);
        double[,] lo = p.Quantile(0.05), hi = p.Quantile(0.95);
        for (int i = 0; i < y.Length; i++) Assert.True(double.IsFinite(lo[i, 0]) && double.IsFinite(hi[i, 0]) && lo[i, 0] < hi[i, 0]);
        // The horseshoe has no variance, so the training fit states no noise sd.
        Assert.Null(fit.TrainingFit(data).NoiseSd);
    }
}
