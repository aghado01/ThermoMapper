using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The asymmetric Laplace on the BASS core (E1's slice 2, quantile regression), and the primitives it brings.
/// <list type="bullet">
/// <item><b>The GIG and the modified half-normal</b>: their draws carry each law's probabilities at their deciles, the law integrated
/// here, on chosen cases and a seeded random sweep with far ranges, every branch of the modified half-normal, the GIG's limits and its
/// truncated draws.</item>
/// <item><b>The asymmetric Laplace</b>: its density is the exponential mixture's, its tails and quantiles round-trip, its draws on and off
/// an interval carry the law.</item>
/// <item><b>Exponential mixing</b>: its marginal against quadrature of the mixture, any drift; its conditional's draws against the GIG.</item>
/// <item><b>The sweep</b>: σ² rides the record under the skew, the target reads it, and σ²'s block runs in the scale slot.</item>
/// <item><b>Retention</b>: a draw keeps its latents alone, and the working response rebuilt from them equals the one its chain held, bit for
/// bit; a fit resumes bit for bit.</item>
/// </list>
/// The long gates are <see cref="BassQuantileQualificationTests"/>'.
/// </summary>
public sealed class BassQuantileTests(ITestOutputHelper output)
{
    // ---- the laws by quadrature ---------------------------------------------------------------------------------------------------

    // The binomial z of a decile k/10 against a law's CDF there.
    private static double DecileZ(double cdf, int k, int n) => (cdf - 0.1 * k) / Math.Sqrt(0.1 * k * (1 - 0.1 * k) / n);

    // Draws' worst decile z against the law of log kernel `logKernel` in u = ln x (Jacobian included), its CDF by quadrature.
    private static double WorstDecile(double[] draws, Func<double, double> logInU, double guess, double lower = double.NegativeInfinity, double upper = double.PositiveInfinity)
    {
        Array.Sort(draws);
        double logZ = Quadrature.LogIntegral(logInU, guess, lower, upper), worst = 0;
        for (int k = 1; k <= 9; k++)
        {
            double q = draws[(int)(0.1 * k * draws.Length)];
            double cdf = Math.Exp(Quadrature.LogIntegral(logInU, Math.Min(guess, Math.Log(q)), lower, Math.Log(q)) - logZ);
            worst = Math.Max(worst, Math.Abs(DecileZ(cdf, k, draws.Length)));
        }
        return worst;
    }

    private static IEnumerable<(double P, double A, double B)> GigCases()
    {
        (double, double, double)[] chosen =
        [
            (0.5, 2.3, 0.04), (0.5, 2.3, 1e-14), (-1, 0.64, 1.3), (2.5, 1e-3, 1e3), (-8, 0.5, 30), (20, 50, 0.01), (-20, 0.01, 50),
            (0, 1, 1), (0.3, 1e-6, 1e-6), (1e4, 1, 1), (-1e4, 1, 1), (0.5, 1e8, 1e-8), (1, 2, 0), (-2, 0, 3),
        ];
        foreach (var c in chosen) yield return c;
        var rng = new Xoshiro256PlusPlus(20261106);
        for (int k = 0; k < 120; k++)
            yield return (Math.Sign(rng.NextDouble() - 0.5) * Math.Pow(10, -3 + 6 * rng.NextDouble()), Math.Pow(10, -8 + 16 * rng.NextDouble()),
                          Math.Pow(10, -8 + 16 * rng.NextDouble()));
    }

    [Fact]
    public void GigDraws_CarryTheLaw_OnChosenCasesAndASweep()
    {
        var rng = new Xoshiro256PlusPlus(5);
        double worst = 0;
        int n = 0;
        foreach ((double p, double a, double b) in GigCases())
        {
            var gig = GeneralizedInverseGaussian.Of(p, a, b);
            double[] draws;
            try { draws = [.. Enumerable.Range(0, 4000).Select(_ => gig.Sample(rng))]; }
            catch (ArithmeticException e) { throw new InvalidOperationException($"GIG({p:R}, {a:R}, {b:R})", e); }
            Assert.All(draws, x => Assert.True(x > 0 && double.IsFinite(x), $"GIG({p}, {a}, {b}) drew {x}"));
            double H(double u) => gig.LogKernel(Math.Exp(u)) + u;
            double z = WorstDecile(draws, H, Math.Log(gig.Mode > 0 ? gig.Mode : 2 * p / a));
            worst = Math.Max(worst, z);
            Assert.True(z < 5.0, $"GIG({p}, {a}, {b}): decile z {z:F2}");
            n++;
        }
        output.WriteLine($"{n} GIG cases, 4,000 draws each: worst decile z {worst:F2}");
    }

    [Theory]
    [InlineData(0.5, 2.3, 0.04, 0.2, 0.9)]
    [InlineData(-3, 1.0, 4.0, 2.5, 8.0)]
    [InlineData(4, 1.0, 1.0, 1e-3, 0.5)]
    [InlineData(4, 1.0, 1.0, 30, 1e9)]
    public void GigDraws_OnAnInterval_CarryTheTruncatedLaw(double p, double a, double b, double lower, double upper)
    {
        var gig = GeneralizedInverseGaussian.Of(p, a, b);
        var rng = new Xoshiro256PlusPlus(9);
        double[] draws = [.. Enumerable.Range(0, 4000).Select(_ => gig.SampleOn(rng, lower, upper))];
        Assert.All(draws, x => Assert.InRange(x, lower, upper));
        double H(double u) => gig.LogKernel(Math.Exp(u)) + u;
        double z = WorstDecile(draws, H, Math.Log(Math.Clamp(gig.Mode, lower, upper)), Math.Log(lower), Math.Log(upper));
        output.WriteLine($"GIG({p}, {a}, {b}) on [{lower}, {upper}]: decile z {z:F2}");
        Assert.True(z < 5.0);
    }

    private static IEnumerable<(double Alpha, double Beta, double Gamma)> MhnCases()
    {
        (double, double, double)[] chosen =
        [
            (520, 300, -15), (5.5, 2, 3), (40, 1, 30), (305, 12.5, 40), (3, 1e-3, 0.5), (1, 0.7, 4), (1, 2, -3), (1, 2, 0),
            (0.6, 1.3, 2.1), (0.3, 2, 0.5), (0.8, 0.5, -2), (0.2, 1e3, -50), (1e6, 1e4, 1e3), (2e4, 1e-4, -1e-2),
        ];
        foreach (var c in chosen) yield return c;
        var rng = new Xoshiro256PlusPlus(20261107);
        for (int k = 0; k < 120; k++)
        {
            double alpha = Math.Pow(10, -1 + 6 * rng.NextDouble()), beta = Math.Pow(10, -4 + 8 * rng.NextDouble());
            double gamma = Math.Sign(rng.NextDouble() - 0.5) * Math.Pow(10, -3 + 5 * rng.NextDouble());
            if (alpha < 1 && gamma > 0 && gamma * gamma / beta > 1e4) gamma = Math.Sqrt(beta) * 10;   // the series' reach
            yield return (alpha, beta, gamma);
        }
    }

    [Fact]
    public void ModifiedHalfNormalDraws_CarryTheLaw_OnEveryBranch()
    {
        var rng = new Xoshiro256PlusPlus(6);
        var worst = new Dictionary<string, double>();
        foreach ((double alpha, double beta, double gamma) in MhnCases())
        {
            var mhn = ModifiedHalfNormal.Of(alpha, beta, gamma);
            double[] draws;
            try { draws = [.. Enumerable.Range(0, 4000).Select(_ => mhn.Sample(rng))]; }
            catch (ArithmeticException e) { throw new InvalidOperationException($"MHN({alpha:R}, {beta:R}, {gamma:R})", e); }
            Assert.All(draws, x => Assert.True(x > 0 && double.IsFinite(x)));
            double H(double u) => mhn.LogKernel(Math.Exp(u)) + u;
            double guess = Math.Log(Math.Max((Math.Abs(gamma) + Math.Sqrt(gamma * gamma + 8 * beta * alpha)) / (4 * beta), 1e-300));
            double z = WorstDecile(draws, H, guess);
            string branch = alpha >= 1 ? "log-concave in x" : gamma <= 0 ? "log-concave in ln x" : "series";
            worst[branch] = Math.Max(worst.GetValueOrDefault(branch), z);
            Assert.True(z < 5.0, $"MHN({alpha}, {beta}, {gamma}), {branch}: decile z {z:F2}");
        }
        foreach ((string branch, double z) in worst) output.WriteLine($"{branch}: worst decile z {z:F2}");
        Assert.Equal(3, worst.Count);
    }

    [Theory]
    [InlineData(0.9)]
    [InlineData(0.5)]
    [InlineData(0.1)]
    [InlineData(0.01)]
    public void AsymmetricLaplace_IsTheExponentialMixture_AndItsTailsRoundTrip(double q)
    {
        var ald = AsymmetricLaplace.Of(0.3, 1.7, q);
        double theta = (1 - 2 * q) / (q * (1 - q)), tau2 = 2 / (q * (1 - q));
        foreach (double x in new[] { -40.0, -3.1, 0.0, 0.3, 0.31, 2.0, 55.0 })
        {
            // The mixture: x = μ + σ(θv + τ√v·ζ), v ~ Exp(1), so its density is ExponentialMixing's marginal in (x − μ)/σ, over σ.
            double mixture = ExponentialMixing.Instance.MarginalLogDensity((x - 0.3) / 1.7, tau2, theta) - Math.Log(1.7);
            Assert.Equal(mixture, ald.LogPdf(x), 12);
            Assert.Equal(1.0, Math.Exp(ald.LogCdf(x)) + Math.Exp(ald.LogSurvival(x)), 13);
        }
        foreach (double p in new[] { 1e-12, 0.01, q, 0.5, 0.99, 1 - 1e-9 })
            Assert.Equal(p, ald.Cdf(ald.Quantile(p)), 9);
        Assert.Equal(0.3, ald.Quantile(q), 13);
        var rng = new Xoshiro256PlusPlus(3);
        double[] draws = [.. Enumerable.Range(0, 20000).Select(_ => ald.Sample(rng)).Order()];
        for (int k = 1; k <= 9; k++) Assert.True(Math.Abs(DecileZ(ald.Cdf(draws[2000 * k]), k, draws.Length)) < 5.0);
        Assert.True(Math.Abs(draws.Average() - ald.Mean) < 5 * Math.Sqrt(ald.Variance / draws.Length));
        foreach ((double lo, double hi) in new[] { (-1.0, 2.0), (1.0, 3.0), (-9.0, -2.0) })
        {
            double[] on = [.. Enumerable.Range(0, 4000).Select(_ => ald.SampleOn(rng, lo, hi)).Order()];
            Assert.All(on, v => Assert.InRange(v, lo, hi));
            double mass = ald.IntervalProbability(lo, hi);
            for (int k = 1; k <= 9; k++)
                Assert.True(Math.Abs(DecileZ(ald.IntervalProbability(lo, on[400 * k]) / mass, k, on.Length)) < 5.0, $"on [{lo}, {hi}], decile {k}");
        }
    }

    // ln ∫ N(r | m·v, s²·v)e^(−v) dv by quadrature over u = ln v.
    private static double ExponentialMixtureLogDensity(double r, double s2, double m)
    {
        double H(double u)
        {
            double v = Math.Exp(u);
            return -0.5 * Math.Log(2 * Math.PI * s2 * v) - (r - m * v) * (r - m * v) / (2 * s2 * v) - v + u;
        }
        double k = Math.Sqrt(m * m + 2 * s2), guess = r == 0 ? Math.Log(0.5 * s2 / (k * k)) : Math.Log(Math.Abs(r) / k);
        return Quadrature.LogIntegral(H, guess);
    }

    [Fact]
    public void ExponentialMixing_MarginalIsTheMixture_AndItsConditionalTheGig()
    {
        var rng = new Xoshiro256PlusPlus(20261108);
        var cases = new List<(double R, double S2, double M)> { (0.0, 1.0, 0.0), (1.3, 0.7, -2.0), (-4.0, 22.2, 8.9), (1e-6, 1e-4, 0.3), (50, 1, -1) };
        for (int k = 0; k < 300; k++)
        {
            double s2 = Math.Pow(10, -6 + 12 * rng.NextDouble()), m = Math.Sign(rng.NextDouble() - 0.5) * Math.Sqrt(s2) * Math.Pow(10, -3 + 4 * rng.NextDouble());
            cases.Add((Math.Sign(rng.NextDouble() - 0.5) * Math.Sqrt(s2) * Math.Pow(10, -4 + 6 * rng.NextDouble()), s2, m));
        }
        double worst = 0;
        foreach ((double r, double s2, double m) in cases)
        {
            double law = ExponentialMixing.Instance.MarginalLogDensity(r, s2, m), mixture = ExponentialMixtureLogDensity(r, s2, m);
            double gap = Math.Abs(law - mixture) / Math.Max(1.0, Math.Abs(mixture));
            worst = Math.Max(worst, gap);
            Assert.True(gap < 1e-10, $"r {r}, s² {s2}, m {m}: the law {law:R}, the mixture {mixture:R}");
        }
        output.WriteLine($"{cases.Count} cases, the worst relative gap {worst:E2}");
        // The conditional at b = 1 and 0.5 against its GIG, through the GIG's own law by quadrature.
        foreach ((double r, double s2, double m, double b) in new[] { (1.3, 0.7, -2.0, 1.0), (0.0, 1.0, 0.5, 1.0), (-4.0, 22.2, 8.9, 0.5) })
        {
            var gig = GeneralizedInverseGaussian.Of(1 - 0.5 * b, 2 + b * m * m / s2, b * r * r / s2);
            double[] draws = [.. Enumerable.Range(0, 4000).Select(_ => ExponentialMixing.Instance.SampleConditional(r, s2, m, b, rng))];
            double H(double u) => gig.LogKernel(Math.Exp(u)) + u;
            Assert.True(WorstDecile(draws, H, Math.Log(Math.Max(gig.Mode, 0.1))) < 5.0);
        }
    }

    // ---- the config and the sweep --------------------------------------------------------------------------------------------------

    internal static BassData Data(int n = 80, int seed = 4)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            y[i] = Math.Sin(2 * Math.PI * x0[i]) + x1[i] + 0.3 * (Gamma.ShapeRate(2, 1).Sample(rng) - 2);
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
    }

    [Fact]
    public void Config_DeclaresTheQuantileNoise()
    {
        BassConfig config = BassStudentTTests.Config(BassNoise.AsymmetricLaplace(0.9));
        Assert.Contains("\"Noise\":{\"Kind\":\"AsymmetricLaplace\",\"VarianceScale\":", config.ToJson());
        Assert.Equal(config, BassConfig.FromJson(config.ToJson()));
        Assert.Equal(1 / 0.9 - 1 / 0.1, config.Noise!.Skew, 14);
        Assert.Throws<ArgumentException>(() => new BassNoise(BassNoiseKind.AsymmetricLaplace, varianceScale: 1.0, quantile: 0.9));
        Assert.Throws<ArgumentException>(() => BassNoise.StudentT(3, reference: BassNoiseReference.Quantile));
    }

    [Fact]
    public void ASweep_HoldsSigma2OnTheRecord_AndTheTargetReadsIt()
    {
        BassData data = Data();
        BassModel model = BassModel.Create(BassStudentTTests.Config(BassNoise.AsymmetricLaplace(0.25)), data);
        var start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(1)).State;
        var chain = new BassChain(model, start.Terms, start.LogTau, new Xoshiro256PlusPlus(2), new Xoshiro256PlusPlus(3));
        Assert.Equal(["moves", "tau", "coefficients", "latent scales", "scale"], chain.Scan.BlockStats().Select(s => s.Key));
        var w = new List<double>();
        for (int k = 0; k < 100; k++)
        {
            chain.Step();
            WorkingResponse working = chain.Current.Fit.Working!;
            double sigma2 = working.Latents!.Sigma2!.Value;
            w.Add(sigma2);
            // The offset is σβv, and the target's marginal is the one given that σ².
            double drift = Math.Sqrt(sigma2) * working.Latents.Skew;
            for (int i = 0; i < 5; i++) Assert.Equal(drift * working.Latents.Values[i], working.Offset[i]);
            BassEvaluation e = model.Target.Evaluate(chain.Current);
            Assert.Equal(e.Marginal().LogMarginalGiven(sigma2), e.LogLikelihood(1.0), 12);
        }
        Assert.True(w.Distinct().Count() > 50, "σ² moves");
        foreach (MoveStat s in chain.Scan.BlockStats().Skip(2)) output.WriteLine($"{s.Key}: {s.Accepted} of {s.Attempts}");
    }

    // ---- retention -----------------------------------------------------------------------------------------------------------------

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void ARetainedDraw_KeepsItsLatentsAlone_AndItsViewRebuildsBitForBit(bool quantile)
    {
        BassData data = Data();
        BassNoise noise = quantile ? BassNoise.AsymmetricLaplace(0.7) : BassNoise.StudentT(4, 0.5);
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
            Assert.NotSame(held, rebuilt);
            Assert.True(rebuilt.Values.SequenceEqual(held.Values));
            Assert.True(rebuilt.Offset.SequenceEqual(held.Offset));
            Assert.True(rebuilt.Weights.SequenceEqual(held.Weights));
            Assert.Equal(held.ResponseMean, rebuilt.ResponseMean);
            Assert.Equal(held.Weighting, rebuilt.Weighting);
            if (quantile) Assert.Equal(held.Latents!.Sigma2, draw.Sigma2);
        }
    }

    [Fact]
    public void AQuantileFit_ResumesBitForBit()
    {
        BassData data = Data();
        BassConfig whole = BassStudentTTests.Config(BassNoise.AsymmetricLaplace(0.8), samples: 60), half = whole with { Sampling = whole.Sampling with { Samples = 30 } };
        BassFit a = Bass.Fit(whole, data);
        BassFit first = Bass.Fit(half, data, new BassRunOptions { Checkpoint = true });
        BassFit b = Bass.Resume(first.Checkpoint!, half, data);
        BassRunTests.AssertSameFit(a, b);
        for (int c = 0; c < a.Draws.Count; c++)
            for (int d = 0; d < a.Draws[c].Count; d++)
            {
                Assert.True(a.Draws[c][d].Latents!.Values.SequenceEqual(b.Draws[c][d].Latents!.Values));
                Assert.Equal(a.Draws[c][d].Latents!.Sigma2, b.Draws[c][d].Latents!.Sigma2);
            }
        foreach (BassChainAccount account in a.Diagnostics.Chains)
            output.WriteLine(string.Join("; ", account.Blocks!.Skip(2).Select(s => $"{s.Key} {s.Accepted} of {s.Attempts}")));
    }

    [Fact]
    public void TheReaders_TakeTheQuantileLaw()
    {
        BassData data = Data();
        const double q = 0.8;
        BassFit fit = Bass.Fit(BassStudentTTests.Config(BassNoise.AsymmetricLaplace(q)), data);
        BassPrediction p = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Seed = 2 });
        double beta = 1 / q - 1 / (1 - q);
        for (int d = 0; d < p.DrawCount; d++) Assert.Equal(Math.Sqrt(p.NoiseVariances[d]) * beta, p.MeanShifts[d], 13);
        // One draw's predictive quantile at q is its surface: f is the q-quantile.
        BassPrediction one = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Draws = [3], Predictive = false });
        double[,] atQ = one.Quantile(q);
        for (int i = 0; i < one.Cases; i++) Assert.Equal(one.Latent[0, i, 0], atQ[i, 0], 9);
        // The pointwise density is the asymmetric Laplace's at the draw's σ.
        PointwiseLogDensity pointwise = fit.PointwiseLogDensity(data);
        ReadOnlySpan<double> y = data.Response.Values;
        for (int d = 0; d < p.DrawCount; d += 7)
        {
            var ald = AsymmetricLaplace.Of(0, Math.Sqrt(p.NoiseVariances[d]), q);
            for (int i = 0; i < y.Length; i++) Assert.Equal(ald.LogPdf(y[i] - p.Latent[d, i, 0]), pointwise.Values[d, i], 10);
        }
        var training = fit.TrainingFit(data);
        Assert.NotNull(training.Offset);
        double expected = Enumerable.Range(0, p.DrawCount).Sum(d => p.Weights[d] * Math.Sqrt(p.NoiseVariances[d] * (beta * beta + 2 / (q * (1 - q)))));
        Assert.Equal(expected, training.NoiseSd!.Value, 12);
    }
}
