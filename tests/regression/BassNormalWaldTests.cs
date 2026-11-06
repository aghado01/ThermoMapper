using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The Normal-Wald, the skewed t and the generalized hyperbolic on the BASS core (E1's slice 3): their declaration; the sweep's noise-shape
/// block, which draws the skew β and the Normal-Wald's γ and holds them on the record with σ²; σ² integrated again where the skew is held at
/// 0; retention and resume bit for bit with the shape; and the readers, the noise law's density, its mean and variance, and its CDF, the
/// mixture over v's nodes, against the density's own integral.
/// </summary>
public sealed class BassNormalWaldTests(ITestOutputHelper output)
{
    private static readonly BassNoise Shaped = BassNoise.NormalWald(new BassNoisePrior(2.0, 1.0), new BassNoisePrior(0.0, 0.7));

    [Fact]
    public void Config_DeclaresTheKinds()
    {
        BassConfig nw = BassStudentTTests.Config(Shaped);
        Assert.Contains("\"Noise\":{\"Kind\":\"NormalWald\",\"VarianceScale\":1,\"Reference\":\"Moments\",\"SkewPrior\":{\"Mean\":0,\"Sd\":0.7},\"GammaPrior\":{\"Mean\":2,\"Sd\":1}}", nw.ToJson());
        Assert.Equal(nw, BassConfig.FromJson(nw.ToJson()));
        Assert.True(Shaped.Skewed && Shaped.SkewDrawn);
        BassNoise plain = BassNoise.NormalWald();
        Assert.Equal(BassNoisePrior.Held(1.0), plain.GammaPrior);
        Assert.False(plain.Skewed);
        BassConfig skewedT = BassStudentTTests.Config(BassNoise.SkewedT(5, new BassNoisePrior(0.2, 0.5)));
        Assert.Equal(skewedT, BassConfig.FromJson(skewedT.ToJson()));
        BassConfig gh = BassStudentTTests.Config(BassNoise.GeneralizedHyperbolic(new GigParameters(0.7, 1.2, 0.8), new BassNoisePrior(0.3, 0.5)));
        Assert.Contains("\"MixingPrior\":{\"P\":0.7,\"A\":1.2,\"B\":0.8}", gh.ToJson());
        Assert.Equal(gh, BassConfig.FromJson(gh.ToJson()));
        // The start of γ under a prior on γ > 0 is that part's mean.
        double m = 2.0, s = 1.0, mills = Math.Exp(-0.5 * 4.0 - 0.5 * Math.Log(2.0 * Math.PI) - Normal.StandardLogCdf(m / s));
        Assert.Equal(m + s * mills, Shaped.GammaStart!.Value, 14);

        Assert.Throws<ArgumentException>(() => BassNoise.SkewedT(5, BassNoisePrior.Held(0.0)));
        Assert.Throws<ArgumentException>(() => new BassNoise(BassNoiseKind.AsymmetricLaplace, varianceScale: 2 / 0.25, quantile: 0.5, skewPrior: new BassNoisePrior(0, 1)));
        Assert.Throws<ArgumentException>(() => new BassNoise(BassNoiseKind.StudentT, 3, gammaPrior: BassNoisePrior.Held(1.0)));
        Assert.Throws<ArgumentException>(() => new BassNoise(BassNoiseKind.GeneralizedHyperbolic));
        Assert.Throws<ArgumentOutOfRangeException>(() => BassNoise.NormalWald(BassNoisePrior.Held(0.0)));
        Assert.Throws<ArgumentException>(() => BassNoise.GeneralizedHyperbolic(new GigParameters(0.5, 0.0, 1.0)));
        Assert.Throws<ArgumentOutOfRangeException>(() => new BassNoisePrior(0.0, -1.0));
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void ASweep_DrawsTheShape_AndHoldsItOnTheRecord(bool skewDrawn)
    {
        BassData data = BassQuantileTests.Data();
        BassNoise noise = skewDrawn ? Shaped : BassNoise.NormalWald(new BassNoisePrior(2.0, 1.0));
        BassModel model = BassModel.Create(BassStudentTTests.Config(noise), data);
        var start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(1)).State;
        var chain = new BassChain(model, start.Terms, start.LogTau, new Xoshiro256PlusPlus(2), new Xoshiro256PlusPlus(3));
        string[] blocks = skewDrawn ? ["moves", "tau", "coefficients", "latent scales", "scale", "noise shape"]
                                    : ["moves", "tau", "coefficients", "latent scales", "noise shape"];
        Assert.Equal(blocks, chain.Scan.BlockStats().Select(s => s.Key));
        var betas = new List<double>();
        var gammas = new List<double>();
        for (int k = 0; k < 150; k++)
        {
            chain.Step();
            WorkingResponse working = chain.Current.Fit.Working!;
            ObservationLatents latents = working.Latents!;
            gammas.Add(latents.Gamma!.Value);
            betas.Add(latents.Skew);
            Assert.True(latents.Gamma > 0.0);
            if (!skewDrawn)
            {
                // A skew held at 0 is a pure scale mixture: σ² is integrated out of the structure step, so the record holds none.
                Assert.Null(latents.Sigma2);
                Assert.Equal(0.0, latents.Skew);
                continue;
            }
            double sigma2 = latents.Sigma2!.Value, drift = Math.Sqrt(sigma2) * latents.Skew;
            for (int i = 0; i < 5; i++) Assert.Equal(drift * latents.Values[i], working.Offset[i]);
            BassEvaluation e = model.Target.Evaluate(chain.Current);
            Assert.Equal(e.Marginal().LogMarginalGiven(sigma2), e.LogLikelihood(1.0), 12);
        }
        Assert.True(gammas.Distinct().Count() > 100, "γ moves");
        if (skewDrawn) Assert.True(betas.Distinct().Count() > 100, "β moves");
        foreach (MoveStat s in chain.Scan.BlockStats().Skip(2)) output.WriteLine($"{s.Key}: {s.Accepted} of {s.Attempts}");
        output.WriteLine($"γ mean {gammas.Average():F3}, β mean {betas.Average():F3}");
    }

    public static IEnumerable<object[]> Noises() =>
    [
        [Shaped],
        [BassNoise.SkewedT(5, new BassNoisePrior(0.0, 0.5), 0.6)],
        [BassNoise.GeneralizedHyperbolic(new GigParameters(0.7, 1.2, 0.8), new BassNoisePrior(0.3, 0.5))],
    ];

    [Theory]
    [MemberData(nameof(Noises))]
    public void ARetainedDraw_RebuildsBitForBit_WithItsShape(BassNoise noise)
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
            Assert.Null(draw.Working);
            Assert.Same(held.Latents, draw.Latents);
            WorkingResponse rebuilt = draw.WorkingOn(model);
            Assert.True(rebuilt.Values.SequenceEqual(held.Values));
            Assert.True(rebuilt.Offset.SequenceEqual(held.Offset));
            Assert.True(rebuilt.Weights.SequenceEqual(held.Weights));
            Assert.Equal(held.Weighting, rebuilt.Weighting);
            Assert.Equal(held.Latents!.Gamma, rebuilt.Latents!.Gamma);
            Assert.Equal(held.Latents.Skew, rebuilt.Latents.Skew);
        }
    }

    [Fact]
    public void ANormalWaldFit_ResumesBitForBit()
    {
        BassData data = BassQuantileTests.Data();
        BassConfig whole = BassStudentTTests.Config(Shaped, samples: 60), half = whole with { Sampling = whole.Sampling with { Samples = 30 } };
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
                Assert.Equal(la.Gamma, lb.Gamma);
            }
        foreach (BassChainAccount account in a.Diagnostics.Chains)
            output.WriteLine(string.Join("; ", account.Blocks!.Skip(2).Select(s => $"{s.Key} {s.Accepted} of {s.Attempts}")));
    }

    [Fact]
    public void TheReaders_TakeTheNormalWaldsLaw()
    {
        BassData data = BassQuantileTests.Data();
        BassFit fit = Bass.Fit(BassStudentTTests.Config(Shaped), data);
        BassPrediction p = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Seed = 2 });
        PointwiseLogDensity pointwise = fit.PointwiseLogDensity(data);
        ReadOnlySpan<double> y = data.Response.Values;
        List<BassDraw> draws = [.. fit.Draws.SelectMany(chain => chain)];
        IBassNoiseLaw law = new VarianceMeanMixture(Shaped);
        Assert.Equal(draws.Count, p.DrawCount);
        for (int d = 0; d < p.DrawCount; d += 7)
        {
            ObservationLatents latents = draws[d].Latents!;
            double sigma2 = latents.Sigma2!.Value, gamma = latents.Gamma!.Value, beta = latents.Skew;
            Assert.Equal(sigma2, p.NoiseVariances[d]);
            // E[v] = 1/γ and Var[v] = 1/γ³ under the inverse Gaussian.
            Assert.Equal(Math.Sqrt(sigma2) * beta / gamma, p.MeanShifts[d], 12);
            Assert.Equal(sigma2 * (beta * beta / (gamma * gamma * gamma) + 1.0 / gamma), law.Variance(sigma2, latents), 12);
            var mixing = GigMixing.NormalWald(gamma);
            for (int i = 0; i < y.Length; i += 9)
                Assert.Equal(mixing.MarginalLogDensity(y[i] - p.Latent[d, i, 0], sigma2, Math.Sqrt(sigma2) * beta), pointwise.Values[d, i], 10);
        }
        Assert.NotNull(fit.TrainingFit(data).NoiseSd);
    }

    // The noise at a draw is the generalized hyperbolic. Its tails by the integral over v against the mass its closed density puts below
    // its quantiles from 0.001 to 0.999, read from the side whose tail is lighter; the table a reader takes, against that integral out to
    // the 10⁻⁹ tails and at the levels.
    [Theory]
    [InlineData(-0.5, 4.0, 1.0, 0.0, 0.8)]
    [InlineData(-0.5, 0.09, 1.0, 1.5, 0.3)]
    [InlineData(-2.5, 0.0, 5.0, -0.7, 1.0)]
    [InlineData(-1.5, 0.0, 3.0, 0.4, 2.0)]
    [InlineData(0.7, 1.2, 0.8, 0.3, 0.5)]
    [InlineData(2.0, 3.0, 0.0, -0.4, 1.0)]
    public void TheNoisesTails_AreItsDensitysIntegral(double p, double a, double b, double beta, double sigma2)
    {
        var mixing = new GigMixing(new GigParameters(p, a, b));
        double m = Math.Sqrt(sigma2) * beta;
        GeneralizedHyperbolic noise = mixing.Noise(sigma2, m);
        double integral = 0.0, level = 0.0, table = 0.0;
        foreach (double probability in new[] { 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999 })
        {
            double x = noise.Quantile(probability);
            // From the lighter side: under a skewed t the drift's side has a power tail, which the reference quadrature truncates.
            bool fromBelow = m > 0.0 || (m == 0.0 && probability <= 0.5);
            double reference = fromBelow
                ? Math.Exp(Quadrature.LogIntegral(r => mixing.MarginalLogDensity(r, sigma2, m), x, upper: x))
                : 1.0 - Math.Exp(Quadrature.LogIntegral(r => mixing.MarginalLogDensity(r, sigma2, m), x, lower: x));
            double exact = fromBelow ? Math.Exp(noise.IntegratedLogCdf(x)) : 1.0 - Math.Exp(noise.IntegratedLogSurvival(x));
            integral = Math.Max(integral, Math.Abs(exact - reference));
            level = Math.Max(level, Math.Abs(reference - probability));
        }
        double low = noise.Quantile(1e-9), high = noise.Quantile(1.0 - 1e-9);
        for (int k = 0; k <= 40; k++)
        {
            double x = low + (high - low) * k / 40.0, below = noise.IntegratedLogCdf(x), above = noise.IntegratedLogSurvival(x);
            table = Math.Max(table, Math.Abs(noise.LogCdf(x) - below) / Math.Max(1.0, Math.Abs(below)));
            table = Math.Max(table, Math.Abs(noise.LogSurvival(x) - above) / Math.Max(1.0, Math.Abs(above)));
        }
        output.WriteLine($"GIG({p}, {a}, {b}), β {beta}, σ² {sigma2}: integral {integral:G3}, table {table:G3}, levels {level:G3}");
        Assert.True(integral <= 1e-8, $"the integral against the density's mass, worst {integral:G3}");
        Assert.True(table <= 5e-4, $"the table against the integral, worst relative gap in the log tails {table:G3}");
        Assert.True(level <= 1e-4, $"the table's quantiles against the density's mass, worst {level:G3}");
    }
}
