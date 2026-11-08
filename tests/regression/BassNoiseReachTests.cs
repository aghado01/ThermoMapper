using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Maths.Samplers.Smc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The engine's reach under a noise (E1's slice 5): ladders and SMC on the likelihood path in the integrated form, whose rungs temper the
/// likelihood given the latents, and the forms that couple them refused with the reason; the prior's exact draw with the latents (a
/// reference rung's and SMC's start); a rung at b = 0 that draws every block from its prior, read through SMC's own kernel; tempered and
/// SMC fits with their evidence; and the direct functional response and categorical predictors under a latent per observation.
/// </summary>
public sealed class BassNoiseReachTests(ITestOutputHelper output)
{
    private const double KolmogorovCut = 2.6;

    private static readonly BassNoise Skewed = BassNoise.NormalWald(new BassNoisePrior(1.5, 0.5), new BassNoisePrior(0.0, 0.7));

    private static BassConfig Config(BassNoise noise) => BassStudentTTests.Config(noise) with
    {
        Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 6 },
        Sigma2 = Sigma2Prior.InverseGamma(3, 0.5),
    };

    private static string Refused(BassConfig config) => Assert.Throws<ArgumentException>(config.Validate).Message;

    [Fact]
    public void Interweaving_IsOptIn_SerializedOnlyWhenOn_AndRefusedWithAHeldSigma2()
    {
        // E1's slice 7: the switch stays out of the serialized form when off, so every existing fit's fingerprint stands.
        BassConfig plain = Config(Skewed), woven = Config(Skewed with { Interweave = true });
        Assert.DoesNotContain("Interweave", plain.ToJson());
        Assert.Contains("\"Interweave\":true", woven.ToJson());
        Assert.Equal(woven, BassConfig.FromJson(woven.ToJson()));
        BassData data = BassQuantileTests.Data();
        Assert.NotEqual(Bass.Fingerprint(plain, data), Bass.Fingerprint(woven, data));
        Assert.Contains("Noise.Interweave", Assert.Throws<ArgumentException>(() => (woven with { Sigma2 = Sigma2Prior.Held(0.3) }).Validate()).Message);

        // The block runs last in the scale slot, and moves σ and v: a woven chain's draws are not the plain chain's.
        BassFit a = Bass.Fit(plain, data), b = Bass.Fit(woven, data);
        MoveStat woveBlock = b.Diagnostics.Chains[0].Blocks!.Single(k => k.Key == "interweaving");
        output.WriteLine($"interweaving: {woveBlock.Accepted} of {woveBlock.Attempts} accepted");
        Assert.True(woveBlock.Attempts > 0 && woveBlock.Accepted > woveBlock.Attempts / 2);
        Assert.NotEqual(a.Draws[0][^1].Latents!.Sigma2, b.Draws[0][^1].Latents!.Sigma2);
    }

    [Fact]
    public void Config_TakesLaddersAndSmc_InTheIntegratedForm()
    {
        BassConfig t = Config(BassNoise.StudentT(4));
        (t with { Sampling = t.Sampling with { Tempering = new TemperingSpec { ReferenceRung = true } } }).Validate();
        (t with { Sampling = t.Sampling with { Smc = new SmcSpec() } }).Validate();
        Assert.Contains("Sampling.Tempering.Form", Refused(t with { Sampling = t.Sampling with { Tempering = new TemperingSpec { Form = TemperedForm.Linear } } }));
        Assert.Contains("Sampling.Tempering.Path", Refused(t with { Sampling = t.Sampling with { Tempering = new TemperingSpec { Path = TemperingPath.Powered } } }));
        Assert.Contains("Sampling.Smc.Form", Refused(t with { Sampling = t.Sampling with { Smc = new SmcSpec { Form = TemperedForm.Linear } } }));
        var flat = new BassFlatHistogramSpec { Run = new Maths.Samplers.FlatHistogram.FlatHistogramSpec { LearningIterations = 100, ProductionDraws = 100 } };
        Assert.Contains("FlatHistogram", Refused(t with { Sampling = t.Sampling with { FlatHistogram = flat } }));
        // Under a skew σ² sits on the record, and b = 0 draws it from its prior.
        BassConfig skewed = Config(Skewed) with { Sigma2 = Sigma2Prior.Jeffreys };
        Assert.Contains("Sigma2", Refused(skewed with { Sampling = skewed.Sampling with { Smc = new SmcSpec() } }));
        Assert.Contains("Sigma2", Refused(skewed with { Sampling = skewed.Sampling with { Tempering = new TemperingSpec { ReferenceRung = true } } }));
        (skewed with { Sampling = skewed.Sampling with { Tempering = new TemperingSpec { Levels = 3, HotEnd = 0.2 } } }).Validate();
    }

    [Fact]
    public void ThePriorsDraw_TakesTheLatentsFromTheirPrior()
    {
        BassData data = BassQuantileTests.Data();
        BassModel model = BassModel.Create(Config(Skewed), data);
        var sampler = new BassPriorSampler(model);
        var rng = new Xoshiro256PlusPlus(20261107);
        const int draws = 3_000;
        var beta = new double[draws];
        var gamma = new double[draws];
        var sigma2 = new double[draws];
        var u = new double[draws];   // v₀'s law at its γ, read back through its CDF: uniform
        for (int k = 0; k < draws; k++)
        {
            BassState state = sampler.DrawPrior(rng).State;
            ObservationLatents latents = state.Fit.Working!.Latents!;
            beta[k] = latents.Skew;
            gamma[k] = latents.Gamma!.Value;
            sigma2[k] = latents.Sigma2!.Value;
            u[k] = InverseGaussianCdf(latents.Values[0], 1.0 / gamma[k]);
        }
        var betaLaw = Normal.MeanSd(0.0, 0.7);
        var gammaLaw = Normal.MeanSd(1.5, 0.5);
        double gammaMass = gammaLaw.Survival(0.0);
        InverseGamma sigmaLaw = model.Target.Sigma2!.Value;
        double[] ks =
        [
            Kolmogorov(beta, betaLaw.Cdf), Kolmogorov(gamma, g => (gammaLaw.Cdf(g) - gammaLaw.Cdf(0.0)) / gammaMass),
            Kolmogorov(sigma2, sigmaLaw.Cdf), Kolmogorov(u, x => x),
        ];
        output.WriteLine($"√n·D: β {ks[0]:F2}, γ {ks[1]:F2}, σ² {ks[2]:F2}, v₀ through its CDF {ks[3]:F2}; prior attempts {sampler.Attempts} for {draws}");
        Assert.All(ks, s => Assert.True(s <= KolmogorovCut, $"√n·D {s:F2}"));
        // A chain's start keeps the model's working response: the latents at their start.
        Assert.Same(model.Working.Latents, sampler.DrawState(rng).State.Fit.Working!.Latents);
    }

    // One sweep at b of D04's chain on SMC's particle, the kernel a stage runs.
    private static BassState Sweep(BassModel model, BassSmcKernel kernel, BassParticle particle, double b, Xoshiro256PlusPlus rng)
    {
        var target = new BassParticleTarget(model.Target);
        return kernel.Move(particle, target.Evaluate(particle), new SmcStage<BassParticle>(target, new Temperature(b, TemperingPath.Likelihood)), rng).State.State;
    }

    public static IEnumerable<object[]> Noises() =>
    [
        [BassNoise.StudentT(4)],
        [Skewed],
        [BassNoise.Horseshoe(new BassNoisePrior(0.0, 0.5))],
        // The interweaving step at b = 0 moves σ with v under their prior (E1's slice 7).
        [BassNoise.StudentT(4) with { Interweave = true }],
        [Skewed with { Interweave = true }],
        [BassNoise.Horseshoe(new BassNoisePrior(0.0, 0.5)) with { Interweave = true }],
    ];

    [Theory]
    [MemberData(nameof(Noises))]
    public void ARungAtZero_DrawsEveryBlockFromItsPrior(BassNoise noise)
    {
        // Exact prior draws (SMC's particles), each moved by one sweep at b = 0, against fresh prior draws: a rung at b = 0 leaves the
        // prior in place, the latents and σ² on the record with it.
        BassData data = BassQuantileTests.Data();
        BassModel model = BassModel.Create(Config(noise), data);
        var particles = new BassParticles(model);
        var kernel = new BassSmcKernel(model);
        var rng = new Xoshiro256PlusPlus(20261108);
        const int draws = 2_000;
        double[] Read(BassState s)
        {
            ObservationLatents l = s.Fit.Working!.Latents!;
            return [s.Count, s.LogTau, Math.Log(l.Values[0]), l.Values.ToArray().Average(Math.Log), l.Skew, Math.Log(l.Sigma2 ?? 1.0)];
        }
        double[][] fresh = [.. Enumerable.Range(0, draws).Select(_ => Read(particles.Draw(rng).State))];
        double[][] moved = [.. Enumerable.Range(0, draws).Select(_ => Read(Sweep(model, kernel, particles.Draw(rng), 0.0, rng)))];
        string[] names = ["M", "ln τ", "ln v₀", "mean ln v", "β", "ln σ²"];
        var gaps = new List<string>();
        for (int g = 0; g < names.Length; g++)
        {
            double[] a = [.. fresh.Select(r => r[g])], c = [.. moved.Select(r => r[g])];
            if (a.Distinct().Count() == 1 && c.Distinct().Count() == 1) continue;
            double statistic = TwoSample(a, c);
            output.WriteLine($"{names[g],-10} √(n/2)·D {statistic:F2}");
            if (statistic > KolmogorovCut) gaps.Add($"{names[g]} {statistic:F2}");
        }
        foreach (MoveStat s in kernel.Account.Moves) output.WriteLine($"{s.Key}: {s.Accepted} of {s.Attempts}");
        Assert.Empty(gaps);
    }

    public static IEnumerable<object[]> Fits() =>
    [
        ["t", "tempered"],
        ["skewed", "tempered"],
        ["horseshoe", "tempered"],
        ["t", "smc"],
        ["skewed", "smc"],
    ];

    [Theory]
    [MemberData(nameof(Fits))]
    public void ATemperedOrSmcFit_RunsTheNoise_AndReportsItsEvidence(string law, string sampler)
    {
        BassData data = BassQuantileTests.Data();
        BassNoise noise = law switch { "t" => BassNoise.StudentT(4), "skewed" => Skewed, _ => BassNoise.Horseshoe(new BassNoisePrior(0.0, 0.5)) };
        BassConfig chains = Config(noise) with { Sampling = new BassSampling { Chains = 2, MasterSeed = 23, BurnIn = 400, Samples = 300, Thin = 2 } };
        BassConfig config = sampler == "tempered"
            ? chains with { Sampling = chains.Sampling with { Tempering = new TemperingSpec { Levels = 5, HotEnd = 0.05, ReferenceRung = true }, ReportEvidence = true } }
            : chains with { Sampling = chains.Sampling with { Smc = new SmcSpec { Particles = 300, Replicates = 2, MutationSteps = 3 }, ReportEvidence = true } };
        BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 2 });
        BassFit reference = Bass.Fit(chains, data, new BassRunOptions { MaxDegreeOfParallelism = 2 });
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        Assert.All(draws, d => Assert.Equal(data.Response.Observations, d.Latents!.Values.Length));
        Assert.True(draws.Select(d => d.Latents!.Values[0]).Distinct().Count() > draws.Length / 4, "the draws hold their own latents");
        double Mean(BassFit f, Func<BassDraw, double> read)
        {
            BassDraw[] all = [.. f.Draws.SelectMany(c => c)];
            double[]? w = f.DrawWeights;
            return w is null ? all.Average(read) : all.Select((d, i) => w[i] * read(d)).Sum() / w.Sum();
        }
        double s2 = Mean(fit, d => Math.Log(d.Sigma2)), s2Reference = Mean(reference, d => Math.Log(d.Sigma2));
        double m = Mean(fit, d => d.Structure.TermIds.Length), mReference = Mean(reference, d => d.Structure.TermIds.Length);
        BassEvidence evidence = fit.Evidence ?? throw new InvalidOperationException("No evidence reported.");
        output.WriteLine($"{law}, {sampler}: ln σ² {s2:F3} against the chains' {s2Reference:F3}, M {m:F2} against {mReference:F2}; " +
                         $"ln p(y) {evidence.LogEvidence:F3} ± {evidence.MonteCarloError:F3}");
        Assert.True(double.IsFinite(evidence.LogEvidence) && double.IsFinite(evidence.MonteCarloError));
        Assert.True(Math.Abs(s2 - s2Reference) <= 0.25, $"ln σ² {s2:F3} against {s2Reference:F3}");
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void ATemperedNoiseFit_ResumesBitForBit(bool interweave)
    {
        // Every rung's state carries its latents on its working response, through swaps and the checkpoint; with the interweaving step too.
        BassData data = BassQuantileTests.Data();
        BassConfig whole = Config(Skewed with { Interweave = interweave }) with
        {
            Sampling = new BassSampling { Chains = 2, MasterSeed = 29, BurnIn = 100, Samples = 120, Thin = 2,
                                          Tempering = new TemperingSpec { Levels = 4, HotEnd = 0.1, ReferenceRung = true }, ReportEvidence = true },
        };
        BassConfig half = whole with { Sampling = whole.Sampling with { Samples = 60 } };
        BassFit full = Bass.Fit(whole, data);
        BassFit first = Bass.Fit(half, data, new BassRunOptions { Checkpoint = true });
        BassFit resumed = Bass.Resume(first.Checkpoint!, half, data);
        BassRunTests.AssertSameFit(full, resumed);
        for (int c = 0; c < full.Draws.Count; c++)
            for (int d = 0; d < full.Draws[c].Count; d++)
            {
                ObservationLatents a = full.Draws[c][d].Latents!, b = resumed.Draws[c][d].Latents!;
                Assert.True(a.Values.SequenceEqual(b.Values));
                Assert.Equal(a.Sigma2, b.Sigma2);
                Assert.Equal(a.Skew, b.Skew);
                Assert.Equal(a.Gamma, b.Gamma);
            }
        Assert.Equal(full.Evidence!.LogEvidence, resumed.Evidence!.LogEvidence);
    }

    [Fact]
    public void ADirectFunctionalResponse_TakesALatentPerObservation()
    {
        BassData data = SensitivityFixture.Functional();
        BassConfig config = BassStudentTTests.Config(BassNoise.StudentT(4), burnIn: 300, samples: 100);
        BassFit fit = Bass.Fit(config, data);
        int n = data.Response.Cases, m = data.Response.Locations;
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        Assert.All(draws, d => Assert.Equal(n * m, d.Latents!.Values.Length));
        // The pointwise density's unit is a curve: its observations' law densities summed.
        BassModel model = BassModel.Create(fit.Config, data);
        PointwiseLogDensity pointwise = fit.PointwiseLogDensity(data);
        Assert.Equal(n, pointwise.Units.Units);
        IBassNoiseLaw law = model.Observation.Law!;
        var surface = new BassDrawFit(model, fit.Terms);
        var f = new double[n * m];
        ReadOnlySpan<double> y = data.Response.Values;
        for (int d = 0; d < draws.Length; d += 17)
        {
            surface.Evaluate(draws[d], f);
            for (int i = 0; i < n; i += 7)
            {
                double sum = 0.0;
                for (int l = 0; l < m; l++) sum += law.LogDensity(y[i * m + l] - f[i * m + l], draws[d].Sigma2, draws[d].LatentsOn(model));
                Assert.Equal(sum, pointwise.Values[d, i], 10);
            }
        }
        BassPrediction p = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Seed = 3 });
        double[,] lo = p.Quantile(0.05), hi = p.Quantile(0.95);
        int covered = 0;
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++) covered += y[i * m + l] >= lo[i, l] && y[i * m + l] <= hi[i, l] ? 1 : 0;
        output.WriteLine($"functional t(4): {draws.Length} draws, the 90% predictive covers {covered} of {n * m} training values; " +
                         $"latent block {string.Join("; ", fit.Diagnostics.Chains.SelectMany(c => c.Blocks!).Where(b => b.Key == "latent scales").Select(b => $"{b.Accepted} of {b.Attempts}"))}");
        Assert.InRange(covered / (double)(n * m), 0.8, 1.0);
    }

    [Fact]
    public void CategoricalPredictors_TakeTheNoise()
    {
        BassData data = SensitivityFixture.Scalar();
        BassFit fit = Bass.Fit(BassStudentTTests.Config(Skewed, burnIn: 300, samples: 100), data);
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        Assert.All(draws, d => Assert.NotNull(d.Latents!.Sigma2));
        BassPrediction p = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Predictive = false });
        double[,] mean = p.Mean();
        ReadOnlySpan<double> y = data.Response.Values;
        double squares = 0.0;
        for (int i = 0; i < y.Length; i++) squares += (y[i] - mean[i, 0]) * (y[i] - mean[i, 0]);
        output.WriteLine($"categorical, skewed: residual rms {Math.Sqrt(squares / y.Length):F4} over {y.Length}");
        Assert.True(Math.Sqrt(squares / y.Length) < 0.2);
    }

    // The inverse Gaussian of mean μ and shape 1, GIG(−½, 1/μ², 1): Φ(√(1/x)(x/μ − 1)) + e^(2/μ)·Φ(−√(1/x)(x/μ + 1)).
    private static double InverseGaussianCdf(double x, double mu)
    {
        double root = Math.Sqrt(1.0 / x);
        return Math.Exp(Normal.StandardLogCdf(root * (x / mu - 1.0))) + Math.Exp(2.0 / mu + Normal.StandardLogCdf(-root * (x / mu + 1.0)));
    }

    private static double Kolmogorov(double[] xs, Func<double, double> cdf)
    {
        double[] sorted = [.. xs.Order()];
        double gap = 0.0;
        for (int i = 0; i < sorted.Length; i++)
        {
            double f = cdf(sorted[i]);
            gap = Math.Max(gap, Math.Max(Math.Abs(f - (double)i / sorted.Length), Math.Abs(f - (i + 1.0) / sorted.Length)));
        }
        return Math.Sqrt(sorted.Length) * gap;
    }

    // √(nm/(n + m))·D for two samples.
    private static double TwoSample(double[] a, double[] b)
    {
        double[] x = [.. a.Order()], y = [.. b.Order()];
        int i = 0, j = 0;
        double gap = 0.0;
        while (i < x.Length && j < y.Length)
        {
            double at = Math.Min(x[i], y[j]);
            while (i < x.Length && x[i] <= at) i++;
            while (j < y.Length && y[j] <= at) j++;
            gap = Math.Max(gap, Math.Abs((double)i / x.Length - (double)j / y.Length));
        }
        return Math.Sqrt(x.Length * (double)y.Length / (x.Length + y.Length)) * gap;
    }
}
