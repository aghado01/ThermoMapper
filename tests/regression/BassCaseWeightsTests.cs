using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bass;
using Maths.Samplers;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Case weights on the data (<see cref="BassData.Weights"/>), the public entry to the weighted core.
/// <list type="bullet">
/// <item><b>Replication.</b> Integer weights fit through <see cref="Bass.Fit"/> as the duplicated rows: two chains' draws state by state,
/// the pointwise terms as the copies' sums, and the tempered marginal and its slope at loaded structures.</item>
/// <item><b>Precision.</b> The data's weights are the core's weighted target bit for bit, under either prior Gram, and a fit's resolved
/// record reads them.</item>
/// <item><b>The fingerprint</b> takes the weights, their meaning and their prior, so a resume on other weights is refused; data that
/// declare none keep theirs.</item>
/// <item><b>Refusals</b>, each with its reason: weights that are not one per case, not positive and finite, or empty; the frozen Gram;
/// and a noise, random effects, a curve or another observation model beside them.</item>
/// </list>
/// </summary>
public sealed class BassCaseWeightsTests(ITestOutputHelper output)
{
    internal static int[] Counts(int n) => [.. Enumerable.Range(0, n).Select(i => 1 + (i * 7 + i / 3) % 3)];

    internal static BassData Weighted(BassData data, double[] weights, WeightMeaning meaning, WeightedPriorGram gram = WeightedPriorGram.Weighted)
        => new(data.Predictors, data.Response, data.Groupings, new CaseWeights(weights, meaning, gram));

    // The rows repeated, each its count's times, in order.
    internal static BassData Duplicated(BassData data, int[] counts)
    {
        IEnumerable<T> Repeat<T>(IEnumerable<T> values) => values.SelectMany((v, i) => Enumerable.Repeat(v, counts[i]));
        BassPredictor[] predictors = [.. data.Predictors.Select(p => (BassPredictor)new ContinuousPredictor(p.Name, [.. Repeat(((ContinuousPredictor)p).Values.ToArray())]))];
        return new BassData(predictors, new ScalarResponse([.. Repeat(data.Response.Values.ToArray())]));
    }

    private static double[] AsDoubles(int[] counts) => [.. counts.Select(c => (double)c)];

    internal static BassConfig FitConfig(int chains = 2, int burnIn = 300, int samples = 50, int thin = 4)
    {
        BassConfig config = BassSweepTests.Config(maxBasis: 6);
        return config with { Sampling = config.Sampling with { Chains = chains, MasterSeed = 9, BurnIn = burnIn, Samples = samples, Thin = thin } };
    }

    // ---- Replication ----

    [Fact]
    public void ReplicationWeights_FitAsDuplicatedRows_ThroughBassFit()
    {
        BassData plain = BassWeightedCoreTests.Scattered(60, 11);
        int[] counts = Counts(plain.Response.Observations);
        BassData weighted = Weighted(plain, AsDoubles(counts), WeightMeaning.Replication), repeated = Duplicated(plain, counts);
        BassConfig config = FitConfig();
        BassFit a = Bass.Fit(config, weighted), b = Bass.Fit(config, repeated);
        Assert.Equal(b.Resolved.SupportMinimum, a.Resolved.SupportMinimum, 1e-12);
        Assert.Equal(b.Resolved.Scale.UnitObservations, a.Resolved.Scale.UnitObservations, 1e-12);
        Assert.Equal(b.Draws.Count, a.Draws.Count);
        int moved = 0;
        for (int c = 0; c < a.Draws.Count; c++)
        {
            Assert.Equal(b.Draws[c].Count, a.Draws[c].Count);
            BassTerm[] last = [];
            for (int d = 0; d < a.Draws[c].Count; d++)
            {
                BassDraw p = a.Draws[c][d], q = b.Draws[c][d];
                BassTerm[] terms = [.. p.Structure.TermIds.Select(id => a.Terms[id])];
                Assert.Equal(q.Structure.TermIds.Select(id => b.Terms[id]), terms);
                if (!terms.SequenceEqual(last)) moved++;
                last = terms;
                BassWeightedCoreTests.Close(q.Tau, p.Tau, 1e-8, $"chain {c}, draw {d}: τ");
                BassWeightedCoreTests.Close(q.Sigma2, p.Sigma2, 1e-8, $"chain {c}, draw {d}: σ²");
                Assert.Equal(q.Coefficients.Length, p.Coefficients.Length);
                for (int j = 0; j < p.Coefficients.Length; j++)
                    BassWeightedCoreTests.Close(q.Coefficients[j], p.Coefficients[j], 1e-7, $"chain {c}, draw {d}: a_{j}");
            }
        }
        Assert.True(moved > 10, $"The retained structures changed {moved} times.");

        // Each weighted row's term is its copies' sum.
        PointwiseLogDensity x = a.PointwiseLogDensity(weighted), y = b.PointwiseLogDensity(repeated);
        int[] first = [.. counts.Select((_, i) => counts[..i].Sum())];
        double worst = 0.0;
        for (int s = 0; s < x.Draws; s++)
            for (int i = 0; i < counts.Length; i++)
            {
                double sum = 0.0;
                for (int k = 0; k < counts[i]; k++) sum += y.Values[s, first[i] + k];
                worst = Math.Max(worst, Math.Abs(sum - x.Values[s, i]) / Math.Max(1.0, Math.Abs(sum)));
            }
        output.WriteLine($"{x.Draws} draws: largest relative gap between a row's term and its copies' {worst:E2}; the structure moved {moved} times");
        Assert.True(worst <= 1e-6, $"{worst:E2}");
    }

    // The tempered marginal and its slope at loaded structures, as the core's model-level check reads them, entered through the data.
    [Theory]
    [InlineData("g", GPriorCentre.Mean)]
    [InlineData("g", GPriorCentre.Flat)]
    [InlineData("ridge", GPriorCentre.Zero)]
    public void ReplicationWeights_OnTheData_EqualDuplicatedRows_ThroughTheMarginal(string form, GPriorCentre intercept)
    {
        BassData data = BassSweepTests.Data(large: true);
        int[] counts = Counts(data.Response.Observations);
        BassData weighted = Weighted(data, AsDoubles(counts), WeightMeaning.Replication), repeated = Duplicated(data, counts);
        CoefficientForm coefficient = form == "g" ? new GForm() : new RidgeForm();
        foreach (Sigma2Prior sigma2 in new[] { Sigma2Prior.ScaledToResponse(0.5), Sigma2Prior.Jeffreys, Sigma2Prior.Held(0.4) })
        {
            BassConfig config = BassWeightedCoreTests.Config(coefficient, intercept, sigma2, 0.2);
            BassModel a = BassModel.Create(config, weighted), b = BassModel.Create(config, repeated);
            Assert.Equal(b.Resolved.CentreLevel, a.Resolved.CentreLevel, 1e-12);
            Assert.Equal(b.Resolved.ResponseMeanSquare, a.Resolved.ResponseMeanSquare, 1e-12);
            foreach (int terms in new[] { 0, 2, 5 })
            {
                BassState p = BassMixing.Loaded(a, BassWeightedCoreTests.Terms[..terms], Math.Log(0.2), out _)!;
                BassState q = BassMixing.Loaded(b, BassWeightedCoreTests.Terms[..terms], Math.Log(0.2), out _)!;
                foreach (double t in new[] { 1.0, 0.5, 0.05 })
                {
                    if (!a.Target.Defines(t)) continue;
                    ITargetEvaluation u = a.Target.Evaluate(p), v = b.Target.Evaluate(q);
                    BassWeightedCoreTests.Close(v.LogLikelihood(t), u.LogLikelihood(t), 1e-10, $"ℓ_{t}, {sigma2}, M = {terms}");
                    BassWeightedCoreTests.Close(((BassEvaluation)v).LogLikelihoodSlope(t), ((BassEvaluation)u).LogLikelihoodSlope(t), 1e-9, $"its slope at {t}, {sigma2}, M = {terms}");
                }
            }
        }
    }

    // ---- Precision ----

    [Theory]
    [InlineData(WeightedPriorGram.Weighted)]
    [InlineData(WeightedPriorGram.Unweighted)]
    public void PrecisionWeights_OnTheData_AreTheCoresWeightedTarget(WeightedPriorGram gram)
    {
        BassData data = BassSweepTests.Data(large: true);
        double[] w = BassWeightedCoreTests.Weights(data.Response.Observations);
        BassData weighted = Weighted(data, w, WeightMeaning.Precision, gram);
        var core = new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, weights: w, meaning: WeightMeaning.Precision, priorGram: gram));
        foreach (CoefficientForm form in new CoefficientForm[] { new GForm(), new RidgeForm() })
            foreach (GPriorCentre intercept in new[] { GPriorCentre.Mean, GPriorCentre.Zero, GPriorCentre.Flat })
            {
                BassConfig config = BassWeightedCoreTests.Config(form, intercept, Sigma2Prior.ScaledToResponse(0.5), 0.2);
                BassModel a = BassModel.Create(config, weighted), b = BassModel.Create(config, data, core);
                Assert.Equal(b.Working.Weights.ToArray(), a.Working.Weights.ToArray());
                Assert.Equal(WeightMeaning.Precision, a.Working.Meaning);
                Assert.Equal(gram, a.Working.PriorGram);
                Assert.Equal(b.Resolved.CentreLevel, a.Resolved.CentreLevel);
                Assert.Equal(b.Resolved.ResponseMeanSquare, a.Resolved.ResponseMeanSquare);
                Assert.Equal(b.Resolved.WeightedMeanSquare, a.Resolved.WeightedMeanSquare);
                Assert.Equal(b.Resolved.Weighting, a.Resolved.Weighting);
                foreach (int terms in new[] { 0, 3, 5 })
                {
                    BassState p = BassMixing.Loaded(a, BassWeightedCoreTests.Terms[..terms], Math.Log(0.2), out _)!;
                    BassState q = BassMixing.Loaded(b, BassWeightedCoreTests.Terms[..terms], Math.Log(0.2), out _)!;
                    foreach (double t in new[] { 1.0, 0.4 })
                    {
                        GPriorMarginal m = a.Target.Marginal(p.Fit, 0.2, t), n = b.Target.Marginal(q.Fit, 0.2, t);
                        Assert.Equal(n.S, m.S);
                        Assert.Equal(n.LogDeterminantTerm, m.LogDeterminantTerm);
                        Assert.Equal(n.LogMarginal(b.Target.Sigma2), m.LogMarginal(a.Target.Sigma2));
                        Assert.Equal(b.Target.Evaluate(q).LogLikelihood(t), a.Target.Evaluate(p).LogLikelihood(t));
                    }
                }
            }
    }

    [Fact]
    public void PrecisionWeights_AFitReadsThem()
    {
        BassData data = BassWeightedCoreTests.Scattered(60, 11);
        double[] w = BassWeightedCoreTests.Weights(data.Response.Observations);
        BassData weighted = Weighted(data, w, WeightMeaning.Precision);
        BassConfig config = FitConfig(burnIn: 100, samples: 20);
        BassFit fit = Bass.Fit(config, weighted), plain = Bass.Fit(config, data);
        BassModel core = BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, weights: w, priorGram: WeightedPriorGram.Weighted)));
        Assert.Equal(core.Resolved.Weighting, fit.Resolved.Weighting);
        Assert.Equal(core.Resolved.CentreLevel, fit.Resolved.CentreLevel);
        Assert.NotEqual(plain.Resolved.CentreLevel, fit.Resolved.CentreLevel);
        // The pointwise terms are N(y | f, σ²/w) at each draw's fit, written out here.
        PointwiseLogDensity terms = fit.PointwiseLogDensity(weighted);
        var at = new BassDrawFit(BassModel.Create(config, weighted), fit.Terms);
        double[] y = data.Response.Values.ToArray(), f = new double[y.Length];
        int s = 0;
        double worst = 0.0;
        foreach (BassDraw draw in fit.Draws.SelectMany(chain => chain))
        {
            at.Evaluate(draw, f);
            for (int i = 0; i < y.Length; i++)
            {
                double e = y[i] - f[i], expected = -0.5 * Math.Log(2.0 * Math.PI * draw.Sigma2 / w[i]) - 0.5 * w[i] * e * e / draw.Sigma2;
                worst = Math.Max(worst, Math.Abs(expected - terms.Values[s, i]) / Math.Max(1.0, Math.Abs(expected)));
            }
            s++;
        }
        Assert.Equal(terms.Draws, s);
        Assert.True(worst <= 1e-12, $"{worst:E2}");
    }

    // ---- The fingerprint ----

    [Fact]
    public void TheFingerprint_TakesTheWeights_SoAResumeCannotChangeThem()
    {
        BassData data = BassWeightedCoreTests.Scattered(40, 5);
        int n = data.Response.Observations;
        double[] w = BassWeightedCoreTests.Weights(n);
        BassConfig config = FitConfig(chains: 2, burnIn: 40, samples: 10, thin: 1);
        string plain = Bass.Fingerprint(config, data);
        Assert.Equal(plain, Bass.Fingerprint(config, new BassData(data.Predictors, data.Response, null, null)));
        string weighted = Bass.Fingerprint(config, Weighted(data, w, WeightMeaning.Precision));
        Assert.NotEqual(plain, weighted);
        Assert.Equal(weighted, Bass.Fingerprint(config, Weighted(data, [.. w], WeightMeaning.Precision)));
        Assert.NotEqual(weighted, Bass.Fingerprint(config, Weighted(data, w, WeightMeaning.Replication)));
        Assert.NotEqual(weighted, Bass.Fingerprint(config, Weighted(data, w, WeightMeaning.Precision, WeightedPriorGram.Unweighted)));
        double[] nudged = [.. w];
        nudged[7] = Math.BitIncrement(nudged[7]);
        Assert.NotEqual(weighted, Bass.Fingerprint(config, Weighted(data, nudged, WeightMeaning.Precision)));
        Assert.NotEqual(Bass.Fingerprint(config, Weighted(data, [.. Enumerable.Repeat(1.0, n)], WeightMeaning.Precision)), plain);

        // A resume continues on its own weights exactly as the uninterrupted run, and refuses others.
        BassFit first = Bass.Fit(config, Weighted(data, w, WeightMeaning.Precision), new BassRunOptions { Checkpoint = true });
        Assert.Throws<ArgumentException>(() => Bass.Resume(first.Checkpoint!, config, Weighted(data, nudged, WeightMeaning.Precision)));
        Assert.Throws<ArgumentException>(() => Bass.Resume(first.Checkpoint!, config, data));
        BassFit resumed = Bass.Resume(first.Checkpoint!, config, Weighted(data, w, WeightMeaning.Precision));
        BassFit whole = Bass.Fit(config with { Sampling = config.Sampling with { Samples = 20 } }, Weighted(data, w, WeightMeaning.Precision));
        for (int c = 0; c < whole.Draws.Count; c++)
        {
            Assert.Equal(whole.Draws[c].Count, resumed.Draws[c].Count);
            for (int d = 0; d < whole.Draws[c].Count; d++)
            {
                Assert.Equal(whole.Draws[c][d].Sigma2, resumed.Draws[c][d].Sigma2);
                Assert.Equal(whole.Draws[c][d].Coefficients, resumed.Draws[c][d].Coefficients);
            }
        }
    }

    // ---- Refusals ----

    [Fact]
    public void Refuses_WeightsThatAreNotOnePositiveFiniteValuePerCase()
    {
        BassData data = BassWeightedCoreTests.Scattered(20, 3);
        Assert.Throws<ArgumentException>(() => new CaseWeights([], WeightMeaning.Precision));
        Assert.Throws<ArgumentNullException>(() => new CaseWeights(null!, WeightMeaning.Precision));
        foreach (double bad in new[] { 0.0, -1.0, double.NaN, double.PositiveInfinity })
        {
            double[] w = [.. Enumerable.Repeat(1.5, 20)];
            w[4] = bad;
            ArgumentException e = Assert.Throws<ArgumentException>(() => new CaseWeights(w, WeightMeaning.Replication));
            Assert.Contains("Case 4", e.Message);
        }
        Assert.Throws<ArgumentException>(() => new CaseWeights([1.0, 2.0], WeightMeaning.Precision, WeightedPriorGram.Frozen));
        Assert.Throws<ArgumentOutOfRangeException>(() => new CaseWeights([1.0, 2.0], (WeightMeaning)7));
        Assert.Throws<ArgumentException>(() => Weighted(data, [.. Enumerable.Repeat(1.0, 19)], WeightMeaning.Precision));
        Assert.Throws<ArgumentException>(() => Weighted(data, [.. Enumerable.Repeat(1.0, 21)], WeightMeaning.Precision));
        // The weights are copied.
        double[] values = [.. Enumerable.Repeat(2.0, 20)];
        var copied = new CaseWeights(values, WeightMeaning.Precision);
        values[0] = 9.0;
        Assert.Equal(2.0, copied.Values.Span[0]);
    }

    [Fact]
    public void Refuses_WhatCaseWeightsAreNotBuiltBeside()
    {
        BassData data = BassStudentTTests.Data(n: 40);
        double[] w = BassWeightedCoreTests.Weights(40);
        BassData weighted = Weighted(data, w, WeightMeaning.Precision);
        // A noise, whose latents hand the core weights of their own.
        Assert.Throws<NotSupportedException>(() => BassModel.Create(BassStudentTTests.Config(BassNoise.StudentT(3, 0.7)), weighted));
        // Another observation model, which hands its own working response.
        Assert.Throws<NotSupportedException>(() => BassModel.Create(BassStudentTTests.Config(null), weighted,
            new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r))));
        // Random effects.
        var grouping = new CategoricalPredictor("g", [.. Enumerable.Range(0, 40).Select(i => i % 4)], ["a", "b", "c", "d"]);
        var grouped = new BassData(data.Predictors, data.Response, [grouping], new CaseWeights(w, WeightMeaning.Precision));
        Assert.Throws<NotSupportedException>(() => BassModel.Create(BassStudentTTests.Config(null) with { RandomIntercepts = [new BassRandomIntercept("g", 1.3)] }, grouped));
        // A curve.
        var t = new double[] { 0.0, 0.5, 1.0 };
        var x = new double[12];
        var curves = new double[12, 3];
        for (int i = 0; i < 12; i++) { x[i] = i / 11.0; for (int r = 0; r < 3; r++) curves[i, r] = x[i] + t[r] + 0.1 * Math.Sin(7 * i + r); }
        var functional = new BassData([new ContinuousPredictor("x", x)], new FunctionalResponse(curves, [new ContinuousPredictor("t", t)]), null,
                                      new CaseWeights([.. Enumerable.Repeat(2.0, 12)], WeightMeaning.Replication));
        Assert.Throws<NotSupportedException>(() => BassModel.Create(BassSweepTests.Config(), functional));
        // Gaussian BASS on the same data takes them.
        Assert.False(BassModel.Create(BassStudentTTests.Config(null), weighted).Working.IsUnweighted);
    }
}
