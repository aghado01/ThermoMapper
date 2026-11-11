using System;
using System.Collections.Generic;
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
/// BASS's candidate density (E16's seam, the conformal layer's slice 1). For every noise law, at the training cases it is the
/// pointwise density's column, with the same draws, chains and latent report, so the layer takes the pair. It refuses what is
/// not one observation read as a unit is. With the qualification switch (<c>BASS_QUALIFICATION</c> = 1): held-out coverage of
/// Gaussian, t, quantile and Normal-Wald fits on data from outside every one of them, with each set's width beside the model's own
/// interval; and the cost row, per candidate and per grid, at the cost harness's base and at S 4,000.
/// </summary>
public sealed class BassCandidateTests(ITestOutputHelper output)
{
    public static TheoryData<string> Laws => ["gaussian", "t", "quantile", "normal-wald", "skewed t", "generalized hyperbolic"];

    internal static BassNoise? Noise(string law) => law switch
    {
        "gaussian" => null,
        "t" => BassNoise.StudentT(3, 0.7),
        "quantile" => BassNoise.AsymmetricLaplace(0.5),
        "normal-wald" => BassNoise.NormalWald(new BassNoisePrior(2.0, 1.0), new BassNoisePrior(0.0, 0.7)),
        "skewed t" => BassNoise.SkewedT(5, new BassNoisePrior(0.2, 0.5)),
        "generalized hyperbolic" => BassNoise.GeneralizedHyperbolic(new GigParameters(0.7, 1.2, 0.8), new BassNoisePrior(0.3, 0.5)),
        _ => throw new ArgumentOutOfRangeException(nameof(law)),
    };

    [Theory]
    [MemberData(nameof(Laws))]
    public void AtTrainingCases_ItIsThePointwiseColumn(string law)
    {
        BassData data = BassStudentTTests.Data();
        BassFit fit = Bass.Fit(BassStudentTTests.Config(Noise(law)), data);
        PointwiseLogDensity pointwise = fit.PointwiseLogDensity(data);
        double[] y = data.Response.Values.ToArray();
        CandidateLogDensity candidates = fit.CandidateLogDensity(data, PredictionCases.Of(data), [.. y.Select(v => new[] { v })]);
        Assert.Equal(pointwise.Draws, candidates.Draws);
        Assert.Equal(pointwise.Latents, candidates.Latents);
        Assert.Equal(pointwise.ChainLengths, candidates.ChainLengths);
        Assert.Null(candidates.Weights);
        if (law != "gaussian")
            Assert.Contains(candidates.Latents, l => l.Scope == LatentScope.Observation && l.Treatment == LatentTreatment.Integrated);
        double worst = 0.0;
        for (int s = 0; s < pointwise.Draws; s++)
            for (int i = 0; i < y.Length; i++)
            {
                double a = pointwise.Values[s, i], b = candidates.Values[s, i];
                worst = Math.Max(worst, Math.Abs(a - b) / Math.Max(1.0, Math.Abs(a)));
            }
        output.WriteLine($"{law}: {pointwise.Draws} draws × {y.Length} cases, largest relative gap {worst:E2}");
        Assert.True(worst <= 1e-12, $"{law}: {worst:E2}");
        // The layer takes the pair: no refusal, and a set at every case.
        ConformalPrediction prediction = ConformalBayes.FullStep(pointwise, candidates);
        Assert.Equal(y.Length, prediction.Sets.Count);
    }

    [Fact]
    public void AWeightedPopulation_CarriesItsWeights()
    {
        // An SMC fit's replicate populations: the candidates carry the pointwise density's weights and no chains, and the layer
        // multiplies its add-one-in weights by them.
        BassData data = BassStudentTTests.Data();
        BassConfig config = BassStudentTTests.Config(null);
        BassFit fit = Bass.Fit(config with { Sampling = config.Sampling with { Smc = new SmcSpec { Particles = 150, Replicates = 2, MutationSteps = 2, TargetConditionalEss = 0.9 } } }, data);
        PointwiseLogDensity pointwise = fit.PointwiseLogDensity(data);
        double[] y = data.Response.Values.ToArray();
        CandidateLogDensity candidates = fit.CandidateLogDensity(data, PredictionCases.Of(data), [.. y.Select(v => new[] { v })]);
        Assert.NotNull(candidates.Weights);
        Assert.Equal(pointwise.Weights, candidates.Weights);
        Assert.Null(candidates.ChainLengths);
        for (int s = 0; s < pointwise.Draws; s++)
            for (int i = 0; i < y.Length; i++) Assert.Equal(pointwise.Values[s, i], candidates.Values[s, i]);
        Assert.Equal(y.Length, ConformalBayes.FullStep(pointwise, candidates).Sets.Count);
    }

    [Fact]
    public void Refuses_WhatIsNotOneObservationReadAsAUnitIs()
    {
        BassConfig config = BassStudentTTests.Config(null);
        BassData data = BassStudentTTests.Data();
        double[] y = data.Response.Values.ToArray();
        double[][] grids = [.. y.Select(v => new[] { v })];
        // Weights the terms read, with none given for the new cases (BassCaseWeightsTests reads them given).
        BassModel weighted = BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, weights: Enumerable.Repeat(2.0, y.Length).ToArray())));
        Assert.Throws<ArgumentException>(() => BassCandidates.LogDensity(weighted, null!, null!, null, PredictionCases.Of(data), grids, null, UnseenLevel.Refuse));
        // A curve as the unit.
        var t = new double[] { 0.0, 0.5, 1.0 };
        var x = new double[12];
        var curves = new double[12, 3];
        for (int i = 0; i < 12; i++) { x[i] = i / 11.0; for (int r = 0; r < 3; r++) curves[i, r] = x[i] + t[r]; }
        var functional = new BassData([new ContinuousPredictor("x", x)], new FunctionalResponse(curves, [new ContinuousPredictor("t", t)]));
        Assert.Throws<NotSupportedException>(() => BassCandidates.LogDensity(BassModel.Create(config, functional), null!, null!, null,
            PredictionCases.Of(functional), [.. Enumerable.Range(0, 12).Select(_ => new[] { 0.0 })], null, UnseenLevel.Refuse));
        // One grid per case.
        BassFit fit = Bass.Fit(config, data);
        Assert.Throws<ArgumentException>(() => fit.CandidateLogDensity(data, PredictionCases.Of(data), grids[1..]));
    }

    // ---- with the qualification switch -------------------------------------------------------------------------------------------

    // The truth: f = sin(2πx₀) + x₁ under a contaminated noise none of the four laws holds, 0.2ζ plus, one time in seven, 0.8·Exp(1).
    internal static (double[] X0, double[] X1, double[] Y) Draw(int n, Xoshiro256PlusPlus rng)
    {
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            double e = 0.2 * Normal.SampleStandard(rng) + (rng.NextDouble() < 1.0 / 7 ? -0.8 * Math.Log(1 - rng.NextDouble()) : 0.0);
            y[i] = Math.Sin(2 * Math.PI * x0[i]) + x1[i] + e;
        }
        return (x0, x1, y);
    }

    [Fact]
    public void HeldOutCoverage_Gaussian() => HeldOutCoverage("gaussian");

    [Fact]
    public void HeldOutCoverage_StudentT() => HeldOutCoverage("t");

    [Fact]
    public void HeldOutCoverage_Quantile() => HeldOutCoverage("quantile");

    [Fact]
    public void HeldOutCoverage_NormalWald() => HeldOutCoverage("normal-wald");

    // Each law its own fact, so a qualification run can take one at a time.
    private void HeldOutCoverage(string law)
    {
        // CB-1 over fits and held-out cases: n 99, so m 100 and k 10 at α 0.1, nominal 0.900. Twenty fits, fifty held-out cases each,
        // the set's width read on twenty of them; the se is the fits' spread, which holds the training data's share of the variance.
        // The model's own 90% interval is read from its predictive draws, one per retained draw.
        if (!BassQualification.Enabled) return;
        const int n = 99, tests = 50, widths = 20, fits = 20;
        double[] grid = [.. Enumerable.Range(0, 101).Select(g => -3.0 + 0.08 * g)];
        int[] seeds = SeedTree.Derive(2026_1107, fits);
        var coverage = new double[fits];
        double width = 0.0, modelCovered = 0.0, modelWidth = 0.0, lowEss = 0, minEss = double.PositiveInfinity, slowest = 0.0, total = 0.0;
        int open = 0, split = 0;
        for (int r = 0; r < fits; r++)
        {
            long t0Tick = System.Diagnostics.Stopwatch.GetTimestamp();
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            (double[] x0, double[] x1, double[] y) = Draw(n, rng);
            (double[] t0, double[] t1, double[] yt) = Draw(tests, rng);
            var data = new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
            BassFit fit = Bass.Fit(BassStudentTTests.Config(Noise(law), burnIn: 1000, samples: 500), data);
            PredictionCases cases = PredictionCases.Of([new ContinuousPredictor("x0", [.. t0, .. t0[..widths]]), new ContinuousPredictor("x1", [.. t1, .. t1[..widths]])]);
            double[][] grids = [.. yt.Select(v => new[] { v }), .. Enumerable.Range(0, widths).Select(_ => grid)];
            ConformalPrediction p = ConformalBayes.FullStep(fit.PointwiseLogDensity(data), fit.CandidateLogDensity(data, cases, grids));
            BassPrediction own = fit.Predict(data, PredictionCases.Of([new ContinuousPredictor("x0", t0), new ContinuousPredictor("x1", t1)]),
                                             new BassPredictOptions { Seed = seeds[r] });
            int covered = 0;
            for (int t = 0; t < tests; t++)
            {
                ConformalSet at = p.Sets[t];
                Assert.Equal(100, at.Units);
                Assert.Equal(10, at.Level);
                if (at.Accepts(0)) covered++;
                double ess = at.EffectiveSampleSizes![0];
                minEss = Math.Min(minEss, ess);
                if (ess < 50) lowEss++;
                double[] draws = [.. Enumerable.Range(0, own.DrawCount).Select(d => own.Predictive![d, t, 0]).Order()];
                double lo = PosteriorInterval.Quantile(draws, 0.05), hi = PosteriorInterval.Quantile(draws, 0.95);
                if (yt[t] >= lo && yt[t] <= hi) modelCovered++;
                modelWidth += hi - lo;
            }
            for (int t = 0; t < widths; t++)
            {
                ConformalSet wide = p.Sets[tests + t];
                width += wide.Width;
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
        output.WriteLine($"{law}: coverage {mean:F4} (nominal 0.9000, se {se:F4} over {fits} fits; binomial alone {Math.Sqrt(0.09 / count):F4}), mean width {width / (fits * widths):F3}, "
                         + $"{split} of {fits * widths} sets of two or more intervals, {open} open at the grid's ends; ESS at the true value below 50 at {lowEss} of {count}, least {minEss:F1}");
        output.WriteLine($"{law}: the model's own 90% predictive interval covers {modelCovered / count:F4}, mean width {modelWidth / count:F3}; {total:F0} s over the fits, the slowest {slowest:F0} s");
        Assert.True(Math.Abs(mean - 0.9) <= 3.5 * Math.Max(se, Math.Sqrt(0.09 / count)), $"{law}: coverage {mean:F4}");
        Assert.Equal(0, open);
    }

    [Theory]
    [InlineData("gaussian", 500)]
    [InlineData("t", 500)]
    [InlineData("quantile", 500)]
    [InlineData("normal-wald", 500)]
    [InlineData("generalized hyperbolic", 500)]
    [InlineData("gaussian", 2000)]
    [InlineData("generalized hyperbolic", 2000)]
    public void CostRow(string law, int samples)
    {
        // The candidate density for grids of 200 values at two cases (one case reads as a constant predictor), then the layer's step on
        // them, per grid, at the cost harness's base point
        // (rich surface, n 1,000, p 10, M_max 20, K 2) under each law, S = 2 chains × samples. Cycle clock, best of five, indicative.
        if (!BassQualification.Enabled) return;
        BassCost.Point point = BassCost.Base;
        BassData data = BassCost.Data(point);
        BassFit fit = Bass.Fit(BassCost.Config(point, burnIn: 1000, samples: samples) with { Noise = Noise(law) }, data);
        PointwiseLogDensity pointwise = fit.PointwiseLogDensity(data);
        double[] grid = [.. Enumerable.Range(0, 200).Select(g => -3.0 + 0.03 * g)];
        var predictors = new List<BassPredictor>();
        foreach (BassPredictor predictor in data.Predictors)
            if (predictor is ContinuousPredictor c) predictors.Add(new ContinuousPredictor(c.Name, c.Values.Span[..2].ToArray()));
        PredictionCases two = PredictionCases.Of(predictors);
        CandidateLogDensity candidates = null!;
        double best = double.PositiveInfinity, step = double.PositiveInfinity;
        for (int r = 0; r < 6; r++)
        {
            long t0 = BassCost.Tick();
            candidates = fit.CandidateLogDensity(data, two, [grid, grid]);
            long t1 = BassCost.Tick();
            ConformalBayes.FullStep(pointwise, candidates);
            long t2 = BassCost.Tick();
            if (r == 0) continue;
            best = Math.Min(best, (t1 - t0) / BassCost.TicksPerSecond * 1e3 / 2);
            step = Math.Min(step, (t2 - t1) / BassCost.TicksPerSecond * 1e3 / 2);
        }
        output.WriteLine($"{law}, S {candidates.Draws}, n {data.Response.Observations}, G {grid.Length} ({(BassCost.Cycles ? "cycles" : "wall")} clock, indicative): "
                         + $"the candidate density {best:F2} ms per grid, {1e3 * best / grid.Length:F2} µs per candidate; the layer's step {step:F2} ms per grid");
    }
}
