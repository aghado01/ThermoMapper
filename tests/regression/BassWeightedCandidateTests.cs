using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's candidate density on a fit with case weights (E16's seam; <see cref="BassFit.CandidateLogDensity"/>): each new case's weight
/// is read at the fit's meaning as the training rows' are. At the training cases, with their own weights, the candidates are the
/// pointwise columns bit for bit under either meaning, the log CDF included, and the layer takes the pair. Refused: a missing weight, one
/// that is not positive and finite, a count that is not one per case, and a weight given to a fit whose terms read none. With the
/// qualification switch (<c>BASS_QUALIFICATION</c> = 1): held-out coverage of the full step on precision-weighted fits, from a truth
/// outside the model.
/// </summary>
public sealed class BassWeightedCandidateTests(ITestOutputHelper output)
{
    private static double[] WeightsFor(WeightMeaning meaning, int n)
        => meaning == WeightMeaning.Replication ? [.. BassCaseWeightsTests.Counts(n).Select(c => (double)c)] : BassWeightedCoreTests.Weights(n);

    [Theory]
    [InlineData(WeightMeaning.Precision)]
    [InlineData(WeightMeaning.Replication)]
    public void AtTrainingCases_WithTheirWeights_ItIsThePointwiseColumn(WeightMeaning meaning)
    {
        BassData plain = BassStudentTTests.Data();
        double[] w = WeightsFor(meaning, plain.Response.Observations);
        BassData data = BassCaseWeightsTests.Weighted(plain, w, meaning);
        BassFit fit = Bass.Fit(BassStudentTTests.Config(null), data);
        PointwiseLogDensity pointwise = fit.PointwiseLogDensity(data, cdf: true);
        double[] y = data.Response.Values.ToArray();
        CandidateLogDensity candidates = fit.CandidateLogDensity(data, PredictionCases.Of(data), [.. y.Select(v => new[] { v })], cdf: true, caseWeights: w);
        Assert.Equal(pointwise.Draws, candidates.Draws);
        Assert.Equal(pointwise.Latents, candidates.Latents);
        Assert.Equal(pointwise.ChainLengths, candidates.ChainLengths);
        for (int s = 0; s < pointwise.Draws; s++)
            for (int i = 0; i < y.Length; i++)
            {
                Assert.Equal(pointwise.Values[s, i], candidates.Values[s, i]);
                Assert.Equal(pointwise.LogCdf![s, i], candidates.LogCdf![s, i]);
            }
        // The weights enter: the unit-weight reading of the same draws differs.
        BassFit unit = Bass.Fit(BassStudentTTests.Config(null), plain);
        Assert.NotEqual(unit.PointwiseLogDensity(plain).Values[0, 0], pointwise.Values[0, 0]);
        // The layer takes the pair, as the full step's sets and its predictive distribution.
        Assert.Equal(y.Length, ConformalBayes.FullStep(pointwise, candidates).Sets.Count);
        Assert.Equal(y.Length, ConformalBayes.Distribution(pointwise, candidates).Distributions.Count);
        output.WriteLine($"{meaning}: {pointwise.Draws} draws × {y.Length} cases, every term and log CDF equal");
    }

    // A tempered ladder's b = 1 replica and SMC's replicate populations compose their draws on the model's own reading too, so their
    // candidates are the pointwise columns as a chain's are.
    [Theory]
    [InlineData("tempering")]
    [InlineData("smc")]
    public void TheOtherSamplers_DrawOnTheModelsReading(string sampler)
    {
        BassData plain = BassStudentTTests.Data(n: 50);
        double[] w = WeightsFor(WeightMeaning.Precision, 50);
        BassData data = BassCaseWeightsTests.Weighted(plain, w, WeightMeaning.Precision);
        BassConfig config = BassStudentTTests.Config(null, burnIn: 60, samples: 20);
        config = config with
        {
            Sampling = sampler == "smc"
                ? config.Sampling with { Smc = new SmcSpec { Particles = 120, Replicates = 2, MutationSteps = 2, TargetConditionalEss = 0.9 } }
                : config.Sampling with { Tempering = new TemperingSpec { InverseTemperatures = [1.0, 0.5, 0.2] } },
        };
        BassFit fit = Bass.Fit(config, data);
        PointwiseLogDensity pointwise = fit.PointwiseLogDensity(data, cdf: true);
        double[] y = data.Response.Values.ToArray();
        CandidateLogDensity candidates = fit.CandidateLogDensity(data, PredictionCases.Of(data), [.. y.Select(v => new[] { v })], cdf: true, caseWeights: w);
        Assert.Equal(pointwise.Latents, candidates.Latents);
        for (int s = 0; s < pointwise.Draws; s++)
            for (int i = 0; i < y.Length; i++)
            {
                Assert.Equal(pointwise.Values[s, i], candidates.Values[s, i]);
                Assert.Equal(pointwise.LogCdf![s, i], candidates.LogCdf![s, i]);
            }
        Assert.Equal(y.Length, ConformalBayes.FullStep(pointwise, candidates).Sets.Count);
    }

    [Fact]
    public void Refuses_AMissingNonPositiveOrExtraWeight()
    {
        BassData plain = BassStudentTTests.Data(n: 40);
        double[] w = BassWeightedCoreTests.Weights(40);
        BassData data = BassCaseWeightsTests.Weighted(plain, w, WeightMeaning.Precision);
        BassConfig config = BassStudentTTests.Config(null, burnIn: 20, samples: 5, chains: 1);
        BassFit fit = Bass.Fit(config, data);
        PredictionCases cases = PredictionCases.Of(data);
        double[][] grids = [.. Enumerable.Range(0, 40).Select(_ => new[] { 0.0, 1.0 })];
        // Missing.
        Assert.Throws<ArgumentException>(() => fit.CandidateLogDensity(data, cases, grids));
        // Not positive and finite.
        foreach (double bad in new[] { 0.0, -2.0, double.NaN, double.PositiveInfinity })
        {
            double[] given = [.. w];
            given[3] = bad;
            ArgumentException e = Assert.Throws<ArgumentException>(() => fit.CandidateLogDensity(data, cases, grids, caseWeights: given));
            Assert.Contains("Case 3", e.Message);
        }
        // Not one per case.
        Assert.Throws<ArgumentException>(() => fit.CandidateLogDensity(data, cases, grids, caseWeights: w[1..]));
        Assert.Throws<ArgumentException>(() => fit.CandidateLogDensity(data, cases, grids, caseWeights: [.. w, 1.0]));
        // Extra: a fit whose terms read none, Gaussian or under a noise.
        BassFit unweighted = Bass.Fit(config, plain);
        Assert.Throws<ArgumentException>(() => unweighted.CandidateLogDensity(plain, PredictionCases.Of(plain), grids, caseWeights: w));
        BassFit noise = Bass.Fit(BassStudentTTests.Config(BassNoise.StudentT(3, 0.7), burnIn: 20, samples: 5, chains: 1), plain);
        Assert.Throws<ArgumentException>(() => noise.CandidateLogDensity(plain, PredictionCases.Of(plain), grids, caseWeights: w));
        // Given, they are read.
        Assert.Equal(40, fit.CandidateLogDensity(data, cases, grids, caseWeights: w).Cases);
    }

    // ---- with the qualification switch -------------------------------------------------------------------------------------------

    // The truth: f = sin(2πx₀) + x₁, and noise at each case's stated precision w, log-uniform on [¼, 4], whose shape none of the model's
    // Normal terms holds: (0.2ζ plus, one time in seven, 0.8·Exp(1)) / √w.
    private static (double[] X0, double[] X1, double[] W, double[] Y) Draw(int n, Xoshiro256PlusPlus rng)
    {
        var x0 = new double[n];
        var x1 = new double[n];
        var w = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            w[i] = Math.Exp((2.0 * rng.NextDouble() - 1.0) * Math.Log(4.0));
            double e = 0.2 * Normal.SampleStandard(rng) + (rng.NextDouble() < 1.0 / 7 ? -0.8 * Math.Log(1 - rng.NextDouble()) : 0.0);
            y[i] = Math.Sin(2 * Math.PI * x0[i]) + x1[i] + e / Math.Sqrt(w[i]);
        }
        return (x0, x1, w, y);
    }

    [Fact]
    public void HeldOutCoverage_PrecisionWeights()
    {
        // CB-1 over fits and held-out cases, as BassCandidateTests reads the unweighted laws: n 99, so m 100 and k 10 at α 0.1, nominal
        // 0.900. Twenty fits, fifty held-out cases each with its own weight, the set's width read on twenty of them. The cases are
        // exchangeable as (x, w, y) triples, so the guarantee is marginal over the weights; coverage by the weight's half is reported.
        if (!BassQualification.Enabled) return;
        const int n = 99, tests = 50, widths = 20, fits = 20;
        double[] grid = [.. Enumerable.Range(0, 151).Select(g => -4.0 + 0.08 * g)];
        int[] seeds = SeedTree.Derive(2026_1110, fits);
        var coverage = new double[fits];
        double width = 0.0, total = 0.0, slowest = 0.0;
        int open = 0, split = 0, heavy = 0, heavyCovered = 0, light = 0, lightCovered = 0;
        double heavyWidth = 0.0, lightWidth = 0.0;
        int heavyWidths = 0, lightWidths = 0;
        for (int r = 0; r < fits; r++)
        {
            long t0Tick = System.Diagnostics.Stopwatch.GetTimestamp();
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            (double[] x0, double[] x1, double[] w, double[] y) = Draw(n, rng);
            (double[] t0, double[] t1, double[] wt, double[] yt) = Draw(tests, rng);
            BassData data = BassCaseWeightsTests.Weighted(
                new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y)), w, WeightMeaning.Precision);
            BassFit fit = Bass.Fit(BassStudentTTests.Config(null, burnIn: 1000, samples: 500), data);
            PredictionCases cases = PredictionCases.Of([new ContinuousPredictor("x0", [.. t0, .. t0[..widths]]), new ContinuousPredictor("x1", [.. t1, .. t1[..widths]])]);
            double[][] grids = [.. yt.Select(v => new[] { v }), .. Enumerable.Range(0, widths).Select(_ => grid)];
            ConformalPrediction p = ConformalBayes.FullStep(fit.PointwiseLogDensity(data),
                fit.CandidateLogDensity(data, cases, grids, caseWeights: [.. wt, .. wt[..widths]]));
            int covered = 0;
            for (int t = 0; t < tests; t++)
            {
                ConformalSet at = p.Sets[t];
                Assert.Equal(100, at.Units);
                Assert.Equal(10, at.Level);
                bool hit = at.Accepts(0);
                if (hit) covered++;
                if (wt[t] >= 1.0) { light++; if (hit) lightCovered++; }
                else { heavy++; if (hit) heavyCovered++; }
            }
            for (int t = 0; t < widths; t++)
            {
                ConformalSet wide = p.Sets[tests + t];
                width += wide.Width;
                if (wt[t] >= 1.0) { lightWidth += wide.Width; lightWidths++; }
                else { heavyWidth += wide.Width; heavyWidths++; }
                if (wide.OpenBelow || wide.OpenAbove) open++;
                if (wide.Intervals.Count > 1) split++;
            }
            coverage[r] = covered / (double)tests;
            double seconds = (System.Diagnostics.Stopwatch.GetTimestamp() - t0Tick) / (double)System.Diagnostics.Stopwatch.Frequency;
            slowest = Math.Max(slowest, seconds);
            total += seconds;
        }
        double mean = coverage.Average(), se = Math.Sqrt(coverage.Sum(c => (c - mean) * (c - mean)) / (fits - 1) / fits);
        int count = fits * tests;
        output.WriteLine($"precision weights: coverage {mean:F4} (nominal 0.9000, se {se:F4} over {fits} fits; binomial alone {Math.Sqrt(0.09 / count):F4}), "
                         + $"mean width {width / (fits * widths):F3}, {split} of {fits * widths} sets of two or more intervals, {open} open at the grid's ends");
        output.WriteLine($"by the weight's half: w ≥ 1 covers {lightCovered / (double)light:F4} of {light} at mean width {lightWidth / Math.Max(1, lightWidths):F3}; "
                         + $"w < 1 covers {heavyCovered / (double)heavy:F4} of {heavy} at mean width {heavyWidth / Math.Max(1, heavyWidths):F3}; {total:F0} s over the fits, the slowest {slowest:F0} s");
        Assert.True(Math.Abs(mean - 0.9) <= 3.5 * Math.Max(se, Math.Sqrt(0.09 / count)), $"coverage {mean:F4}");
        Assert.Equal(0, open);
    }
}
