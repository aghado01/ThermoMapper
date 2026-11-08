using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// KHAOS on the engine's reach (slice 5), its long gates, behind <see cref="BassQualification"/>; without it each returns at once. Every
/// fit takes KHAOS's defaults: Legendre terms of degree up to 8 and order up to 2, the modified g-prior (ζ = 1, g₀² ~ IG(0.001, 1,000)),
/// h₁ = 4, h₂ = 4/n, the zero centre and Jeffreys σ², on inputs declared on [0, 1]. The truth is sin(2πx₀) + x₁, outside the model.
/// <list type="bullet">
/// <item><b>Conformal coverage on a Legendre fit:</b> the full conformal Bayes set (CB-1) on held-out cases under the contaminated noise
/// of <see cref="BassCandidateTests"/>, which no law here holds, against its nominal 0.9.</item>
/// <item><b>A Student t KHAOS fit</b> on noise with t₂ tails: its posterior mean surface against the truth, beside the Gaussian fit's
/// on the same data.</item>
/// <item><b>A quantile KHAOS fit</b> at q = 0.8 on heteroscedastic noise: the held-out fraction below its posterior mean surface, and
/// its distance from the true 0.8-quantile surface beside the Gaussian fit's.</item>
/// </list>
/// The noise laws read no reference: the zero centre fixes the level, and Jeffreys reads no scale.
/// </summary>
public sealed class KhaosReachQualificationTests(ITestOutputHelper output)
{
    private static BassData Data(double[] x0, double[] x1, double[] y)
    {
        var domain = new InputDomain(0.0, 1.0);
        return new BassData([new ContinuousPredictor("x0", x0, domain), new ContinuousPredictor("x1", x1, domain)], new ScalarResponse(y));
    }

    internal static BassConfig Config(BassNoise? noise, int cases, int seed, int burnIn = 2_000, int samples = 1_000) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 2, Legendre = new LegendreBasis { MaxDegree = 8 } },
        Count = new BassCountPrior(4, 4.0 / cases),
        Coefficients = BassCoefficientPrior.ModifiedG(TauPrior.Gamma(0.001, 1000.0)),
        Intercept = GPriorCentre.Zero,
        Sigma2 = Sigma2Prior.Jeffreys,
        Noise = noise,
        Sampling = new BassSampling { Chains = 2, MasterSeed = seed, BurnIn = burnIn, Samples = samples },
    };

    private static double Truth(double x0, double x1) => Math.Sin(2 * Math.PI * x0) + x1;

    // ---- conformal coverage on a Legendre fit -----------------------------------------------------------------------------------------

    [Fact]
    public void ConformalCoverage_OnALegendreFit()
    {
        // As BassCandidateTests' held-out coverage: n 99, so m 100 and k 10 at α 0.1, nominal 0.900; twenty fits, fifty held-out cases each.
        if (!BassQualification.Enabled) return;
        const int n = 99, tests = 50, fits = 20;
        int[] seeds = SeedTree.Derive(2026_1108, fits);
        var coverage = new double[fits];
        double width = 0.0;
        int open = 0;
        double[] grid = [.. Enumerable.Range(0, 101).Select(g => -3.0 + 0.08 * g)];
        for (int r = 0; r < fits; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            (double[] x0, double[] x1, double[] y) = BassCandidateTests.Draw(n, rng);
            (double[] t0, double[] t1, double[] yt) = BassCandidateTests.Draw(tests, rng);
            BassData data = Data(x0, x1, y);
            BassFit fit = Bass.Fit(Config(null, n, seeds[r], burnIn: 1_000, samples: 500), data);
            var domain = new InputDomain(0.0, 1.0);
            PredictionCases cases = PredictionCases.Of([new ContinuousPredictor("x0", [.. t0, t0[0]], domain), new ContinuousPredictor("x1", [.. t1, t1[0]], domain)]);
            ConformalPrediction p = ConformalBayes.FullStep(fit.PointwiseLogDensity(data), fit.CandidateLogDensity(data, cases, [.. yt.Select(v => new[] { v }), grid]));
            coverage[r] = Enumerable.Range(0, tests).Count(t => p.Sets[t].Accepts(0)) / (double)tests;
            ConformalSet wide = p.Sets[tests];
            width += wide.Width;
            if (wide.OpenBelow || wide.OpenAbove) open++;
        }
        double mean = coverage.Average(), se = Math.Sqrt(coverage.Sum(c => (c - mean) * (c - mean)) / (fits - 1) / fits);
        output.WriteLine($"KHAOS, Gaussian, contaminated noise: coverage {mean:F4} (nominal 0.9000, se {se:F4} over {fits} fits; binomial alone " +
                         $"{Math.Sqrt(0.09 / (fits * tests)):F4}); mean width {width / fits:F3}; {open} of {fits} sets open at the grid's ends");
        Assert.True(Math.Abs(mean - 0.9) <= 3.5 * Math.Max(se, Math.Sqrt(0.09 / (fits * tests))), $"coverage {mean:F4}");
        Assert.Equal(0, open);
    }

    // ---- a Student t fit and a quantile fit ------------------------------------------------------------------------------------------

    // Each fit's posterior mean surface at the held-out cases.
    private static double[] Surface(BassFit fit, BassData data, double[] t0, double[] t1, int seed)
    {
        var domain = new InputDomain(0.0, 1.0);
        BassPrediction p = fit.Predict(data, PredictionCases.Of([new ContinuousPredictor("x0", t0, domain), new ContinuousPredictor("x1", t1, domain)]),
                                       new BassPredictOptions { Seed = seed, Predictive = false });
        double[,] mean = p.Mean();
        return [.. Enumerable.Range(0, t0.Length).Select(i => mean[i, 0])];
    }

    private static double Rms(double[] a, double[] b) => Math.Sqrt(a.Zip(b, (u, v) => (u - v) * (u - v)).Average());

    [Fact]
    public void StudentT_RecoversTheSurface_UnderHeavyTails()
    {
        // n 150, noise 0.25·t₂; sixteen data sets, each fitted Gaussian and t(3), against the truth at 200 held-out cases. t₂'s tails make the
        // paired gap heavy-tailed (six sets read 1.3 SE with the t fit better on five), so the gate is the sign: better on at least 12 of 16.
        if (!BassQualification.Enabled) return;
        const int n = 150, held = 200, sets = 16;
        int[] seeds = SeedTree.Derive(2026_1111, sets);
        var gaussian = new double[sets];
        var student = new double[sets];
        for (int r = 0; r < sets; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            double T2() => Normal.SampleStandard(rng) / Math.Sqrt(-Math.Log(1.0 - rng.NextDouble()));   // χ²₂/2 = Exp(1)
            double[] x0 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())], x1 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            double[] y = [.. Enumerable.Range(0, n).Select(i => Truth(x0[i], x1[i]) + 0.25 * T2())];
            double[] t0 = [.. Enumerable.Range(0, held).Select(_ => rng.NextDouble())], t1 = [.. Enumerable.Range(0, held).Select(_ => rng.NextDouble())];
            double[] truth = [.. Enumerable.Range(0, held).Select(i => Truth(t0[i], t1[i]))];
            BassData data = Data(x0, x1, y);
            gaussian[r] = Rms(Surface(Bass.Fit(Config(null, n, seeds[r]), data), data, t0, t1, seeds[r]), truth);
            student[r] = Rms(Surface(Bass.Fit(Config(BassNoise.StudentT(3, reference: BassNoiseReference.None), n, seeds[r]), data), data, t0, t1, seeds[r]), truth);
            output.WriteLine($"set {r}: the surface's RMS error against the truth, Gaussian {gaussian[r]:F4}, t(3) {student[r]:F4}");
        }
        double[] gap = [.. gaussian.Zip(student, (g, t) => g - t)];
        double mean = gap.Average(), se = Math.Sqrt(gap.Sum(d => (d - mean) * (d - mean)) / (sets - 1) / sets);
        int wins = gap.Count(d => d > 0.0);
        output.WriteLine($"t(3) against Gaussian KHAOS under 0.25·t₂ noise: better on {wins} of {sets} sets; mean RMS error {student.Average():F4} against " +
                         $"{gaussian.Average():F4}, paired gap {mean:F4} ± {se:F4}");
        Assert.True(wins >= 12, $"the t fit better on {wins} of {sets} sets");
    }

    [Fact]
    public void Quantile_FitsItsQuantileSurface()
    {
        // n 300, noise (0.1 + 0.4x₁)ζ; the 0.8-quantile surface is f + (0.1 + 0.4x₁)z₀.₈. Four data sets, 300 held-out cases each.
        if (!BassQualification.Enabled) return;
        const int n = 300, held = 300, sets = 4;
        const double q = 0.8;
        double z = Normal.StandardQuantile(q);
        int[] seeds = SeedTree.Derive(2026_1110, sets);
        int below = 0;
        double quantileError = 0.0, gaussianError = 0.0;
        for (int r = 0; r < sets; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            double[] x0 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())], x1 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            double[] y = [.. Enumerable.Range(0, n).Select(i => Truth(x0[i], x1[i]) + (0.1 + 0.4 * x1[i]) * Normal.SampleStandard(rng))];
            double[] t0 = [.. Enumerable.Range(0, held).Select(_ => rng.NextDouble())], t1 = [.. Enumerable.Range(0, held).Select(_ => rng.NextDouble())];
            double[] yt = [.. Enumerable.Range(0, held).Select(i => Truth(t0[i], t1[i]) + (0.1 + 0.4 * t1[i]) * Normal.SampleStandard(rng))];
            double[] surface = [.. Enumerable.Range(0, held).Select(i => Truth(t0[i], t1[i]) + (0.1 + 0.4 * t1[i]) * z)];
            BassData data = Data(x0, x1, y);
            double[] fitted = Surface(Bass.Fit(Config(BassNoise.AsymmetricLaplace(q, reference: BassNoiseReference.None), n, seeds[r]), data), data, t0, t1, seeds[r]);
            double[] mean = Surface(Bass.Fit(Config(null, n, seeds[r]), data), data, t0, t1, seeds[r]);
            int under = Enumerable.Range(0, held).Count(i => yt[i] <= fitted[i]);
            below += under;
            quantileError += Rms(fitted, surface) / sets;
            gaussianError += Rms(mean, surface) / sets;
            output.WriteLine($"set {r}: {under} of {held} held-out cases below the fitted 0.8-quantile surface; its RMS distance from the true one " +
                             $"{Rms(fitted, surface):F4}, the Gaussian mean surface's {Rms(mean, surface):F4}");
        }
        int count = sets * held;
        double fraction = below / (double)count, se = Math.Sqrt(q * (1 - q) / count);
        output.WriteLine($"quantile KHAOS at q = {q}: {fraction:F4} of {count} held-out cases below its surface (± {se:F4}); RMS distance from the true " +
                         $"quantile surface {quantileError:F4}, the Gaussian mean surface's {gaussianError:F4}");
        Assert.True(Math.Abs(fraction - q) <= 3 * se + 0.02, $"fraction below {fraction:F4}");
        Assert.True(quantileError < 0.5 * gaussianError, $"the quantile surface's error {quantileError:F4} against the mean's {gaussianError:F4}");
    }
}
