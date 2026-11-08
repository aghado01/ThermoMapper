using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Student t on the BASS core (E1's slice 1, robust BASS): v ~ IG(ν/2, ν/2), the core fitting y under weights 1/(c·v).
/// <list type="bullet">
/// <item><b>The mixing law</b> (<see cref="StudentTMixing"/>): its marginal against quadrature of the mixture over v, on chosen cases and a
/// seeded random sweep with far ranges; its conditional's draws against the exact Gamma moments at several powers b; its moments; the
/// log-gamma ratio across its switch.</item>
/// <item><b>v's conditional is <see cref="StudentTWeights"/>' draw,</b> number for number from the same stream.</item>
/// <item><b>The config</b> declares the noise, a Gaussian config's form and fingerprint unchanged, and the slices not built yet are
/// refused with the reason.</item>
/// <item><b>The chain</b> runs the coefficients' block then the latent scales' in their slots, and the latents ride the working
/// response; a fit resumes bit for bit; its account counts the blocks.</item>
/// <item><b>The readers</b> take the law: the pointwise density integrates each vᵢ, prediction's noise and quantiles are t's, and the
/// training fit's noise sd is the law's, null where infinite.</item>
/// </list>
/// The gates that need long runs (Geweke, recovery with coverage, large ν, the cost row) are <see cref="BassStudentTQualificationTests"/>'.
/// </summary>
public sealed class BassStudentTTests(ITestOutputHelper output)
{
    // ---- the mixing law ------------------------------------------------------------------------------------------------------------

    // ln ∫ N(r | 0, s²v)·IG(v; ν/2, ν/2) dv by quadrature over u = ln v. With h = ν/2, the inverse gamma's log density times the
    // Jacobian v is h·ln h − ln Γ(h) − h·u − h·e⁻ᵘ; from h = 10³ on it is taken as ½ln h − ½ln 2π − 1/(12h) + 1/(360h³) − h·(e⁻ᵘ − 1 + u),
    // Stirling's form, since the direct one cancels terms of size h and loses ~10⁻¹⁶·h ln h.
    private static double MixtureLogDensity(double r, double s2, double nu)
    {
        double half = 0.5 * nu, logNorm = half * Math.Log(half) - SpecialFunctions.LogGamma(half);
        double Mixing(double u)
        {
            if (half < 1e3) return logNorm - half * u - half * Math.Exp(-u);
            double g = Math.Abs(u) < 1e-2 ? u * u * (0.5 + u * (-1.0 / 6 + u * (1.0 / 24 + u * (-1.0 / 120 + u / 720)))) : Math.Exp(-u) - 1 + u;
            return 0.5 * Math.Log(half) - 0.5 * Math.Log(2 * Math.PI) - 1 / (12 * half) + 1 / (360 * half * half * half) - half * g;
        }
        double H(double u)
        {
            double v = Math.Exp(u);
            return -0.5 * Math.Log(2 * Math.PI * s2) - 0.5 * u - r * r / (2 * s2 * v)   // N(r | 0, s²v)
                   + Mixing(u);                                                         // IG(v; ν/2, ν/2)·v
        }
        // The peak of the integrand in u: v solving the quadratic of its stationary point.
        double guess = Math.Log((half + r * r / (2 * s2)) / (half + 0.5));
        return Quadrature.LogIntegral(H, guess);
    }

    [Fact]
    public void Marginal_MatchesTheMixture_OnChosenCasesAndASweep()
    {
        (double R, double S2, double Nu)[] chosen =
        [
            (0.0, 1.0, 4.0), (1.3, 1.0, 4.0), (-7.5, 0.3, 3.0), (25.0, 1.0, 1.0), (0.2, 4.0, 0.3), (1e-6, 1e-8, 5.0), (3e4, 2e6, 2.0),
            (1.0, 1.0, 1e6), (2.0, 1.0, 1e9), (0.5, 1e-3, 1e12), (40.0, 1.0, 0.05),
        ];
        var rng = new Xoshiro256PlusPlus(20261105);
        var cases = new List<(double R, double S2, double Nu)>(chosen);
        for (int k = 0; k < 400; k++)
        {
            double nu = Math.Pow(10, -1.3 + 12.3 * rng.NextDouble());        // 0.05 to 10¹¹
            double s2 = Math.Pow(10, -12 + 24 * rng.NextDouble());
            double z = Math.Pow(10, -8 + 12 * rng.NextDouble()) * (rng.NextDouble() < 0.5 ? -1 : 1);
            cases.Add((z * Math.Sqrt(s2), s2, nu));
        }
        double worst = 0, bars = 0, barsNu = 0;
        foreach ((double r, double s2, double nu) in cases)
        {
            double law = new StudentTMixing(nu).MarginalLogDensity(r, s2, 0.0), mixture = MixtureLogDensity(r, s2, nu);
            double gap = Math.Abs(law - mixture) / Math.Max(1.0, Math.Abs(mixture));
            worst = Math.Max(worst, gap);
            Assert.True(gap < 1e-10, $"t({nu}) at r = {r}, s² = {s2}: the law {law:R}, the mixture {mixture:R}");
            // A reading, not a check: BARS's StudentTWeights takes the direct forms, which slice 6 routes through this law.
            double direct = Math.Abs(new StudentTWeights(nu).MarginalLogDensity(r, s2) - mixture) / Math.Max(1.0, Math.Abs(mixture));
            if (direct > bars) (bars, barsNu) = (direct, nu);
        }
        output.WriteLine($"{cases.Count} cases, the worst relative gap {worst:E2}; StudentTWeights' direct forms, worst {bars:E2} at ν = {barsNu:G3}");
    }

    [Fact]
    public void LogGammaHalfRatio_HoldsItsDigits_AcrossItsSwitch()
    {
        foreach (double x in new[] { 999.0, 1000.0, 1001.0, 1500.0, 3000.0 })
        {
            double direct = SpecialFunctions.LogGamma(x + 0.5) - SpecialFunctions.LogGamma(x);
            Assert.Equal(direct, StudentTMixing.LogGammaHalfRatio(x), 11);
        }
        // Far out: ½ln x − 1/(8x) + 1/(192x³), whose next term is 10⁻²⁶ at x = 10⁶.
        foreach (double x in new[] { 1e6, 1e9, 1e14 })
            Assert.Equal(0.5 * Math.Log(x) - 1.0 / (8 * x) + 1.0 / (192 * x * x * x), StudentTMixing.LogGammaHalfRatio(x), 15);
    }

    [Theory]
    [InlineData(4.0, 1.7, 1.0)]
    [InlineData(4.0, 1.7, 0.4)]
    [InlineData(0.6, 0.02, 1.0)]
    [InlineData(30.0, 9.0, 1.0)]
    [InlineData(2.5, 0.0, 0.0)]
    public void Conditional_DrawsTheGammaOfTheReciprocal(double nu, double z2, double b)
    {
        // λ = 1/v ~ Gamma((ν + b)/2, (ν + b·z²)/2): its mean and its log's mean, against 40,000 draws.
        var law = new StudentTMixing(nu);
        var rng = new Xoshiro256PlusPlus(7);
        const int n = 40_000;
        double shape = 0.5 * (nu + b), rate = 0.5 * (nu + b * z2);
        var lambdas = new double[n];
        for (int k = 0; k < n; k++) lambdas[k] = 1.0 / law.SampleConditional(Math.Sqrt(z2) * 1.3, 1.69, 0.0, b, rng);
        double mean = lambdas.Average(), sd = Math.Sqrt(shape) / rate;
        double logMean = lambdas.Select(x => Math.Log(x)).Average(), logSd = Math.Sqrt(Trigamma(shape));
        output.WriteLine($"t({nu}), z² {z2}, b {b}: E[λ] {mean:G6} vs {shape / rate:G6}, E[ln λ] {logMean:G6} vs {SpecialFunctions.Digamma(shape) - Math.Log(rate):G6}");
        Assert.True(Math.Abs(mean - shape / rate) < 5 * sd / Math.Sqrt(n));
        Assert.True(Math.Abs(logMean - (SpecialFunctions.Digamma(shape) - Math.Log(rate))) < 5 * logSd / Math.Sqrt(n));
    }

    private static double Trigamma(double x)
    {
        double sum = 0.0;
        for (; x < 20; x++) sum += 1.0 / (x * x);
        return sum + 1 / x + 1 / (2 * x * x) + 1 / (6 * x * x * x) - 1 / (30 * Math.Pow(x, 5)) + 1 / (42 * Math.Pow(x, 7));
    }

    [Fact]
    public void Moments_AreTheInverseGammas_AndInfiniteWhereItIs()
    {
        Assert.Equal(5.0 / 3.0, new StudentTMixing(5).Mean, 15);
        Assert.Equal(2 * 25.0 / (9 * 1), new StudentTMixing(5).Variance, 13);
        Assert.Equal(double.PositiveInfinity, new StudentTMixing(2).Mean);
        Assert.Equal(double.PositiveInfinity, new StudentTMixing(4).Variance);
        // A drift is the skewed t (E1's slice 3): GIG mixing at t's prior, draw for draw.
        var skewed = new GigMixing(new StudentTMixing(4).Prior);
        Assert.Equal(skewed.SampleConditional(1.0, 1.0, 0.3, 1.0, new Xoshiro256PlusPlus(1)),
                     new StudentTMixing(4).SampleConditional(1.0, 1.0, 0.3, 1.0, new Xoshiro256PlusPlus(1)));
    }

    [Theory]
    [InlineData(0.5)]
    [InlineData(3.0)]
    [InlineData(4.0)]
    [InlineData(30.0)]
    public void VConditional_IsStudentTWeightsDraw_NumberForNumber(double nu)
    {
        var rng = new Xoshiro256PlusPlus(11);
        double[] residuals = [.. Enumerable.Range(0, 500).Select(_ => 3 * Normal.SampleStandard(rng) * Math.Exp(rng.NextDouble()))];
        const double sigma2 = 0.83;
        double[] weights = new StudentTWeights(nu).Resample(residuals, sigma2, new Xoshiro256PlusPlus(5));
        var law = new StudentTMixing(nu);
        var stream = new Xoshiro256PlusPlus(5);
        for (int i = 0; i < residuals.Length; i++)
        {
            Assert.Equal(1.0 / weights[i], law.SampleConditional(residuals[i], sigma2, 0.0, 1.0, stream));
            Assert.Equal(new StudentTWeights(nu).MarginalLogDensity(residuals[i], sigma2), law.MarginalLogDensity(residuals[i], sigma2, 0.0), 11);
        }
    }

    // ---- the config, the chain and the run -----------------------------------------------------------------------------------------

    internal static BassData Data(int n = 80, int seed = 3, double nu = 3.0)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var t = StudentT.Standard(nu);
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            y[i] = Math.Sin(2 * Math.PI * x0[i]) + x1[i] + 0.3 * t.Sample(rng);
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
    }

    internal static BassConfig Config(BassNoise? noise, int burnIn = 200, int samples = 60, int chains = 2) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 12 },
        Noise = noise,
        Sampling = new BassSampling { Chains = chains, MasterSeed = 17, BurnIn = burnIn, Samples = samples, Thin = 2 },
    };

    [Fact]
    public void Config_DeclaresTheNoise_AndAGaussianFormStands()
    {
        BassConfig gaussian = Config(null), t = Config(BassNoise.StudentT(3, 1.0 / 3));
        Assert.DoesNotContain("Noise", gaussian.ToJson());
        Assert.Contains("\"Noise\":{\"Kind\":\"StudentT\",\"DegreesOfFreedom\":3,", t.ToJson());
        Assert.Equal(t, BassConfig.FromJson(t.ToJson()));
        BassData data = Data();
        Assert.NotEqual(Bass.Fingerprint(gaussian, data), Bass.Fingerprint(t, data));
        Assert.IsType<VarianceMeanMixture>(BassModel.Create(t, data).Observation);
        Assert.Same(GaussianObservation.Instance, BassModel.Create(gaussian, data).Observation);

        // A form that couples the latents, and a reference that cannot stand, are refused with the reason (slice 5 takes the rest).
        Assert.Contains("Sampling.Tempering.Form",
                        Assert.Throws<ArgumentException>(() => (t with { Sampling = t.Sampling with { Tempering = new TemperingSpec { Form = Spline.Bars.TemperedForm.Linear } } }).Validate()).Message);
        Assert.Contains("Noise.Reference", Assert.Throws<ArgumentException>(() => Config(BassNoise.StudentT(3, reference: BassNoiseReference.None)).Validate()).Message);
        Assert.IsType<VarianceMeanMixture>(BassModel.Create(t, SensitivityFixture.Functional()).Observation);
        Assert.Throws<ArgumentException>(() => BassModel.Create(t, data, GaussianObservation.Instance));
    }

    [Fact]
    public void ASweep_RunsTheBlocksInTheirSlots_AndTheLatentsRideTheWorkingResponse()
    {
        BassData data = Data();
        var noise = BassNoise.StudentT(4, 0.5);
        BassModel model = BassModel.Create(Config(noise), data);
        Assert.Equal(1.0 / 0.5, model.Working.Weights[0]);
        Assert.All(model.Working.Latents!.Values.ToArray(), v => Assert.Equal(1.0, v));
        var start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(1)).State;
        var chain = new BassChain(model, start.Terms, start.LogTau, new Xoshiro256PlusPlus(2), new Xoshiro256PlusPlus(3));
        Assert.Equal(["moves", "tau", "coefficients", "latent scales"], chain.Scan.BlockStats().Select(s => s.Key));
        for (int k = 0; k < 50; k++) chain.Step();
        WorkingResponse working = chain.Current.Fit.Working!;
        ReadOnlySpan<double> v = working.Latents!.Values, w = working.Weights;
        Assert.Equal(data.Response.Observations, v.Length);
        for (int i = 0; i < v.Length; i++) Assert.Equal(1.0 / (0.5 * v[i]), w[i]);
        Assert.True(working.Offset.IsEmpty);
        Assert.True(working.Values.SequenceEqual(data.Response.Values));
        IReadOnlyList<Maths.Samplers.Rjmcmc.MoveStat> blocks = chain.Scan.BlockStats();
        Assert.Equal(50, blocks[2].Attempts);
        Assert.Equal(50, blocks[3].Attempts);
        output.WriteLine($"latent block: {blocks[3].Accepted} of {blocks[3].Attempts} moved; v from {v.ToArray().Min():G3} to {v.ToArray().Max():G3}");
    }

    [Fact]
    public void AFit_ResumesBitForBit_AndCountsItsBlocks()
    {
        BassData data = Data();
        BassConfig whole = Config(BassNoise.StudentT(4), samples: 60), half = Config(BassNoise.StudentT(4), samples: 30);
        BassFit a = Bass.Fit(whole, data);
        BassFit first = Bass.Fit(half, data, new BassRunOptions { Checkpoint = true });
        BassFit b = Bass.Resume(first.Checkpoint!, half, data);
        BassRunTests.AssertSameFit(a, b);
        for (int c = 0; c < a.Draws.Count; c++)
            for (int d = 0; d < a.Draws[c].Count; d++)
                Assert.True(a.Draws[c][d].Latents!.Values.SequenceEqual(b.Draws[c][d].Latents!.Values));
        foreach (BassChainAccount account in a.Diagnostics.Chains)
        {
            IReadOnlyList<Maths.Samplers.Rjmcmc.MoveStat> blocks = account.Blocks!;
            Assert.Equal(["moves", "tau", "coefficients", "latent scales"], blocks.Select(s => s.Key));
            Assert.All(blocks.Skip(1), s => Assert.Equal(whole.Sampling.Samples * whole.Sampling.Thin, s.Attempts));
            output.WriteLine($"latent scales: {blocks[3].Accepted} of {blocks[3].Attempts} followed, {blocks[3].Attempts - blocks[3].Accepted} refused");
        }
    }

    // ---- the readers ---------------------------------------------------------------------------------------------------------------

    private static readonly Lazy<(BassData Data, BassFit Fit)> TFit = new(() =>
    {
        BassData data = Data();
        return (data, Bass.Fit(Config(BassNoise.StudentT(3, 0.7)), data));
    });

    [Fact]
    public void Pointwise_IsTheLawsDensity_EachLatentIntegrated()
    {
        (BassData data, BassFit fit) = TFit.Value;
        PointwiseLogDensity pointwise = fit.PointwiseLogDensity(data);
        // Each vᵢ is integrated against its prior, which no other observation moves given the globals: an in-sample term, so WAIC reads it
        // and PSIS-LOO has an lppd to subtract from.
        Assert.Contains(pointwise.Latents, l => l.Scope == LatentScope.Observation && l.Treatment == LatentTreatment.Integrated);
        Assert.False(pointwise.LeavesUnitsOut);
        Assert.True(double.IsFinite(PredictiveChecks.Waic(pointwise).Elpd));
        Assert.True(double.IsFinite(PredictiveChecks.PsisLoo(pointwise).EffectiveParameters));
        BassPrediction at =fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Predictive = false });
        ReadOnlySpan<double> y = data.Response.Values;
        for (int d = 0; d < at.DrawCount; d++)
        {
            var t = StudentT.LocationScale(3, 0, Math.Sqrt(at.NoiseVariances[d] * 0.7));
            for (int i = 0; i < y.Length; i++)
                Assert.Equal(t.LogPdf(y[i] - at.Latent[d, i, 0]), pointwise.Values[d, i], 10);
        }
        Assert.Throws<NotSupportedException>(() => fit.PointwiseLogDensity(data, content: PointwiseContent.CollapsedPredictive));
    }

    [Fact]
    public void Prediction_DrawsTheLawsNoise_AndSolvesItsMixture()
    {
        (BassData data, BassFit fit) = TFit.Value;
        BassPrediction p = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Seed = 3 });
        Assert.All(p.MeanShifts, s => Assert.Equal(0.0, s));
        // The standardized noise (ŷ − f)/√(σ²c) carries t₃'s probabilities at its deciles.
        var z = new List<double>();
        for (int d = 0; d < p.DrawCount; d++)
            for (int i = 0; i < p.Cases; i++) z.Add((p.Predictive![d, i, 0] - p.Latent[d, i, 0]) / Math.Sqrt(p.NoiseVariances[d] * 0.7));
        z.Sort();
        var t3 = StudentT.Standard(3);
        for (int k = 1; k <= 9; k++)
        {
            double q = z[(int)(0.1 * k * z.Count)], f = t3.Cdf(q);
            Assert.True(Math.Abs(f - 0.1 * k) < 5 * Math.Sqrt(0.1 * k * (1 - 0.1 * k) / z.Count), $"decile {k}: F = {f}");
        }
        // One draw's mixture is its own t: the quantile is f + √(σ²c)·t₃'s.
        BassPrediction one = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Draws = [5], Predictive = false });
        double[,] q90 = one.Quantile(0.9);
        double scale = Math.Sqrt(one.NoiseVariances[0] * 0.7);
        for (int i = 0; i < one.Cases; i++) Assert.Equal(one.Latent[0, i, 0] + scale * t3.Quantile(0.9), q90[i, 0], 9);
    }

    [Fact]
    public void TrainingFit_StatesTheLawsNoiseSd_AndNoneWhereItIsInfinite()
    {
        (BassData data, BassFit fit) = TFit.Value;
        BassTrainingFit training = fit.TrainingFit(data);
        BassPrediction p = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Predictive = false });
        double expected = Enumerable.Range(0, p.DrawCount).Sum(d => p.Weights[d] * Math.Sqrt(p.NoiseVariances[d] * 0.7 * 3.0));   // c·E[v] = 0.7·3
        Assert.Equal(expected, training.NoiseSd!.Value, 12);
        Assert.Null(training.Offset);

        BassData small = Data(n: 40);
        BassFit cauchy = Bass.Fit(Config(BassNoise.StudentT(2), burnIn: 50, samples: 10), small);
        Assert.Null(cauchy.TrainingFit(small).NoiseSd);
    }
}

/// <summary>ln ∫ exp(h(u)) du over (lower, upper) for a unimodal h, by Gauss–Legendre panels a standard deviation of h's curvature wide,
/// walked out from its peak (or the bound nearest it) until a panel adds under e⁻⁴⁰; each panel summed about its own largest term.</summary>
internal static class Quadrature
{
    public static double LogIntegral(Func<double, double> h, double guess, double lower = double.NegativeInfinity, double upper = double.PositiveInfinity)
    {
        double peak = Peak(h, guess, lower, upper), e = 1e-4 * Math.Max(1.0, Math.Abs(peak));
        double curvature = -(h(peak + e) - 2.0 * h(peak) + h(peak - e)) / (e * e);
        double width = curvature > 4.0 ? 1.0 / Math.Sqrt(curvature) : 0.5;
        var panels = new List<double>();
        foreach (int side in new[] { 1, -1 })
        {
            double best = double.NegativeInfinity, step = width, a = peak;
            for (int k = 0; k < 200000; k++)
            {
                if (side > 0 ? a >= upper : a <= lower) break;
                double b = side > 0 ? Math.Min(a + step, upper) : Math.Max(a - step, lower);
                double panel = Panel(h, Math.Min(a, b), Math.Max(a, b));
                panels.Add(panel);
                best = Math.Max(best, panel);
                if (panel < best - 40.0 && h(b) < h(a)) break;
                a = b;
                if (k > 20) step = Math.Min(0.5, step * 1.05);
            }
        }
        double top = panels.Max();
        return top + Math.Log(panels.Sum(p => Math.Exp(p - top)));
    }

    private static double Peak(Func<double, double> h, double guess, double lower, double upper)
    {
        double Clamp(double u) => Math.Min(Math.Max(u, lower), upper);
        double x = Clamp(guess), step = 0.25;
        int dir = h(Clamp(x + 1e-6)) >= h(x) ? 1 : -1;
        double a = x, b = x;
        for (int k = 0; k < 200; k++)
        {
            double next = Clamp(b + dir * step);
            if (next == b || h(next) < h(b)) break;
            a = b;
            b = next;
            step *= 2.0;
        }
        double lo = Math.Min(a, Clamp(b + dir * step)), hi = Math.Max(a, Clamp(b + dir * step)), g = (Math.Sqrt(5.0) - 1.0) / 2.0;
        for (int k = 0; k < 300 && hi - lo > 1e-13 * Math.Max(1.0, Math.Abs(lo)); k++)
        {
            double c = hi - g * (hi - lo), d = lo + g * (hi - lo);
            if (h(c) >= h(d)) hi = d; else lo = c;
        }
        return 0.5 * (lo + hi);
    }

    private static double Panel(Func<double, double> h, double a, double b)
    {
        GaussLegendre rule = GaussLegendre.Of(40);
        double mid = 0.5 * (a + b), half = 0.5 * (b - a);
        var terms = new double[rule.Count];
        for (int i = 0; i < rule.Count; i++) terms[i] = Math.Log(rule.Weights[i] * half) + h(mid + half * rule.Nodes[i]);
        double top = terms.Max();
        return double.IsNegativeInfinity(top) ? top : top + Math.Log(terms.Sum(t => Math.Exp(t - top)));
    }
}
