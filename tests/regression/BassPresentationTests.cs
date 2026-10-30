using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// A direct fit's presentation data (stage 9, O9; P09, P10, P12, P14) against direct computations: the traces against the
/// draw stores and their thinning; the posterior summaries' means, sds, extremes and quantiles from the draws under equal and
/// unequal weights, the quantiles against independent type-7 and weighted type-7 (Akinshin's) computations; fitted and residual data
/// against <see cref="BassFit.Predict"/> at the training inputs, the intervals the predictive mixture's (O-O1's construction,
/// where R's per-draw endpoints differ); the metadata round-tripping the fit's config, resolution, identity, account and
/// evidence; and the refusals.
/// </summary>
public sealed class BassPresentationTests
{
    private readonly ITestOutputHelper _out;
    public BassPresentationTests(ITestOutputHelper output) => _out = output;

    // ── Summaries ─────────────────────────────────────────────────────────────────────────────────────────────────

    // The type-7 quantile of equal-weight draws, by hand.
    internal static double Type7(IEnumerable<double> values, double p)
    {
        double[] sorted = [.. values.Order()];
        double position = p * (sorted.Length - 1);
        int lo = (int)Math.Floor(position);
        return sorted[lo] + (position - lo) * (sorted[Math.Min(lo + 1, sorted.Length - 1)] - sorted[lo]);
    }

    // Akinshin's weighted type 7 by hand, written apart from DrawSample's: the draws unpooled in value order, each holding
    // its slice [c₋, c] of the normalized cumulative weight, and taking its overlap with the window [(h − 1)/n*, h/n*],
    // h = (n* − 1)p + 1, over the window's width 1/n*.
    internal static double WeightedType7(IReadOnlyList<double> values, IReadOnlyList<double> weights, double p)
    {
        var kept = values.Zip(weights).Where(x => x.Second > 0).OrderBy(x => x.First).ToArray();
        double total = kept.Sum(x => x.Second), effective = total * total / kept.Sum(x => x.Second * x.Second);
        double h = (effective - 1) * p + 1, lo = (h - 1) / effective, hi = h / effective, below = 0.0, quantile = 0.0;
        foreach ((double value, double weight) in kept)
        {
            double from = below / total, to = (below + weight) / total;
            below += weight;
            quantile += Math.Max(0.0, Math.Min(to, hi) - Math.Max(from, lo)) * effective * value;
        }
        return quantile;
    }

    // A summary against its draws: the mean, sd and extremes under the weights, and each quantile by `quantile`.
    internal static void AssertSummary(PosteriorSummary s, IReadOnlyList<double> values, IReadOnlyList<double>? weights,
                                       Func<IReadOnlyList<double>, IReadOnlyList<double>?, double, double> quantile, int digits = 12)
    {
        IReadOnlyList<double> w = weights ?? Enumerable.Repeat(1.0, values.Count).ToArray();
        double mass = w.Sum(), mean = values.Zip(w).Sum(x => x.First * x.Second) / mass;
        Assert.Equal(mean, s.Mean, digits);
        Assert.Equal(Math.Sqrt(values.Zip(w).Sum(x => x.Second * (x.First - mean) * (x.First - mean)) / mass), s.Sd, digits);
        Assert.Equal(values.Where((_, d) => w[d] > 0).Min(), s.Min);
        Assert.Equal(values.Where((_, d) => w[d] > 0).Max(), s.Max);
        Assert.Equal(w.Count(x => x > 0), s.Draws);
        for (int q = 0; q < s.Probabilities.Count; q++) Assert.Equal(quantile(values, weights, s.Probabilities[q]), s.Quantiles[q], digits);
    }

    internal static double Reference(IReadOnlyList<double> values, IReadOnlyList<double>? weights, double p)
        => weights is null || weights.Where(w => w > 0).Distinct().Count() == 1
            ? Type7(values.Where((_, d) => weights is null || weights[d] > 0), p)
            : WeightedType7(values, weights, p);

    [Fact]
    public void PosteriorSummary_MatchesItsDraws_UnderEqualAndUnequalWeights()
    {
        double[] values = [3.0, -1.0, 2.5, 7.0, 0.5, 2.5, 4.0];
        AssertSummary(PosteriorSummary.Of(values), values, null, Reference);
        // Normalized equal weights are equal weights: the quantiles stay type 7.
        double[] normalized = [.. Enumerable.Repeat(1.0 / values.Length, values.Length)];
        PosteriorSummary equal = PosteriorSummary.Of(values, normalized);
        Assert.Equal(PosteriorSummary.Of(values).Quantiles, equal.Quantiles);
        double[] unequal = [0.5, 2.0, 0.0, 1.0, 3.0, 0.25, 1.0];
        AssertSummary(PosteriorSummary.Of(values, unequal, [0.0, 0.1, 0.5, 0.9, 1.0]), values, unequal, Reference);
        Assert.Equal(PosteriorSummary.DefaultProbabilities, PosteriorSummary.Of(values).Probabilities);

        // No mass: NaN throughout. Refusals: a mismatched or bad weight, a non-finite value, a probability outside [0, 1].
        PosteriorSummary none = PosteriorSummary.Of([1.0, 2.0], [0.0, 0.0]);
        Assert.True(double.IsNaN(none.Mean) && double.IsNaN(none.Max) && none.Quantiles.All(double.IsNaN) && none.Draws == 0);
        Assert.Throws<ArgumentException>(() => PosteriorSummary.Of(values, [1.0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => PosteriorSummary.Of([1.0, 2.0], [1.0, -1.0]));
        Assert.Throws<ArgumentException>(() => PosteriorSummary.Of([1.0, double.NaN]));
        Assert.Throws<ArgumentOutOfRangeException>(() => PosteriorSummary.Of(values, null, [1.5]));
    }

    [Fact]
    public void Summary_IsTheDrawsPooledOverTheChains_UnderTheFitsWeights()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        foreach (double[]? weights in new[] { null, Weights(draws.Length) })
        {
            BassFit weighted = fit with { DrawWeights = weights };
            BassSummary summary = weighted.Summary();
            AssertSummary(summary.Posterior.Count, [.. draws.Select(d => (double)d.Structure.TermIds.Length)], weights, Reference);
            AssertSummary(summary.Posterior.Sigma2, [.. draws.Select(d => d.Sigma2)], weights, Reference);
            AssertSummary(summary.Posterior.Tau!, [.. draws.Select(d => d.Tau)], weights, Reference);
            AssertSummary(summary.Posterior.Lambda, [.. draws.Select(d => d.Lambda)], weights, Reference);
            // R's summary: the range of the basis count and the posterior mean error variance.
            Assert.Equal(draws.Min(d => d.Structure.TermIds.Length), summary.Posterior.Count.Min);
            Assert.Equal(draws.Max(d => d.Structure.TermIds.Length), summary.Posterior.Count.Max);
        }
        Assert.Null((fit with { Config = fit.Config with { Coefficients = new FixedG() } }).Summary().Posterior.Tau);
        _out.WriteLine($"M in [{fit.Summary().Posterior.Count.Min}, {fit.Summary().Posterior.Count.Max}], E[σ²] = {fit.Summary().Posterior.Sigma2.Mean:G4}");
    }

    // Unequal weights that sum to anything: the summaries normalize them.
    internal static double[] Weights(int count) => [.. Enumerable.Range(0, count).Select(d => 0.2 + (d * 37 % 11) / 4.0)];

    // ── Traces ────────────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Traces_AreTheDrawStores_WithTheirSweeps()
    {
        BassFit fit = SensitivityFixture.FitScalar();   // thin 2
        BassTraces traces = fit.Traces();
        Assert.Equal(fit.Identity.Chains, traces.Chains.Count);
        Assert.Equal(2, traces.Thin);
        Assert.True(traces.TauSampled);
        for (int c = 0; c < fit.Draws.Count; c++)
        {
            BassChainTrace t = traces.Chains[c];
            DrawStore<BassDraw> store = fit.Draws[c];
            Assert.Equal(c, t.Chain);
            Assert.Equal(store.Count, t.Count.Length);
            Assert.Equal(store.Select(d => d.Structure.TermIds.Length).ToArray(), t.Count);
            Assert.Equal(store.Select(d => d.Sigma2).ToArray(), t.Sigma2);
            Assert.Equal(store.Select(d => d.Tau).ToArray(), t.Tau);
            Assert.Equal(store.Select(d => d.Lambda).ToArray(), t.Lambda);
            Assert.Equal(Enumerable.Range(1, store.Count).Select(i => 2L * i).ToArray(), t.Iterations);
        }
        // A store keeping one in three: draw i is the (3i + 1)-th retained step, each two sweeps.
        var sparse = new DrawStore<BassDraw>(every: 3);
        foreach (BassDraw d in fit.Draws[0]) sparse.Add(d);
        BassChainTrace thinned = (fit with { Draws = [sparse] }).Traces().Chains[0];
        Assert.Equal(Enumerable.Range(0, sparse.Count).Select(i => (3L * i + 1) * 2).ToArray(), thinned.Iterations);
        Assert.False((fit with { Config = fit.Config with { Coefficients = new FixedG() } }).Traces().TauSampled);
    }

    // ── Fitted against observed ───────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(BassPointwiseTests.Kind.Scalar)]
    [InlineData(BassPointwiseTests.Kind.Categorical)]
    [InlineData(BassPointwiseTests.Kind.Functional)]
    public void TrainingFit_IsThePredictionAtTheTrainingInputs(BassPointwiseTests.Kind kind)
    {
        (BassData data, BassFit fit) = BassPredictTests.FitOf(kind);
        foreach (double[]? weights in new[] { null, Weights(fit.Draws.Sum(c => c.Count)) })
        {
            BassFit weighted = fit with { DrawWeights = weights };
            BassTrainingFit t = weighted.TrainingFit(data);
            BassPrediction p = weighted.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Predictive = false });
            int n = data.Response.Cases, m = data.Response.Locations;
            Assert.Equal((n, m), (t.Cases, t.Locations));
            Assert.Equal(p.Mean(), t.Mean);
            Assert.Equal(IntervalSource.PredictiveMixture, t.Intervals);
            Assert.Equal(PosteriorSummary.DefaultProbabilities, t.Probabilities);
            for (int q = 0; q < t.Probabilities.Count; q++)
            {
                double[,] at = p.Quantile(t.Probabilities[q]);
                for (int i = 0; i < n; i++)
                    for (int l = 0; l < m; l++) Assert.Equal(at[i, l], t.Quantiles[q, i, l]);
            }
            for (int i = 0; i < n; i++)
                for (int l = 0; l < m; l++)
                {
                    double y = data.Response is FunctionalResponse f ? f[i, l] : ((ScalarResponse)data.Response).Y.Span[i];
                    Assert.Equal(y, t.Observed[i, l]);
                    Assert.Equal(y - t.Mean[i, l], t.Residuals[i, l]);
                }
            double sd = Enumerable.Range(0, p.DrawCount).Sum(d => p.Weights[d] * Math.Sqrt(p.NoiseVariances[d]));
            Assert.Equal(sd, t.NoiseSd!.Value, 13);
            if (data.Response is FunctionalResponse functional)
            {
                ContinuousColumn coordinates = Assert.Single(t.LocationCoordinates!);
                Assert.Equal("t", coordinates.Name);
                Assert.Equal(functional.OutputVariables[0].Values.ToArray(), coordinates.Values.ToArray());
            }
            else Assert.Null(t.LocationCoordinates);
        }
    }

    [Fact]
    public void TrainingFit_Intervals_AreThePredictiveMixtures_NotRsPerDrawEndpoints()
    {
        // O-O1's construction: 1000 draws whose means spread as Φ⁻¹((k − ½)/1000) at sd 1, or whose sds spread as e^{u/2}
        // at mean 0. Intercept-only draws make every case's latent the draw's mean. R's plot takes each quantile of the
        // per-draw endpoints μ + Φ⁻¹(p)·σ: ±3.91 and ±5.20 in the fixture, against the mixture's ±2.77 and ±2.64.
        (BassData data, BassFit fit) = BassPredictTests.FitOf(BassPointwiseTests.Kind.Scalar);
        const int draws = 1000;
        double[] u = [.. Enumerable.Range(1, draws).Select(k => Normal.StandardQuantile((k - 0.5) / draws))];
        (double[] Mu, double[] Sd, double Mixture, double R)[] spreads =
            [(u, [.. Enumerable.Repeat(1.0, draws)], 2.77, 3.91), (new double[draws], [.. u.Select(v => Math.Exp(v / 2))], 2.64, 5.20)];
        double[] probabilities = [0.025, 0.975];
        double[] equal = [.. Enumerable.Repeat(1.0, draws)];
        foreach ((double[] mu, double[] sd, double mixture, double r) in spreads)
        {
            var store = new DrawStore<BassDraw>();
            for (int k = 0; k < draws; k++) store.Add(new BassDraw(new BassStructure([]), 1.0, sd[k] * sd[k], [mu[k]], 1.0));
            BassTrainingFit t = (fit with { Draws = [store], DrawWeights = null }).TrainingFit(data, probabilities);
            for (int q = 0; q < 2; q++)
            {
                double exact = PredictiveMixture.Quantile(mu, sd, equal, probabilities[q]);
                double endpoints = Type7(mu.Zip(sd, (m, s) => m + Normal.StandardQuantile(probabilities[q]) * s), probabilities[q]);
                for (int i = 0; i < t.Cases; i++) Assert.Equal(exact, t.Quantiles[q, i, 0], 12);
                Assert.Equal((q == 0 ? -1 : 1) * mixture, exact, 0.006);
                Assert.Equal((q == 0 ? -1 : 1) * r, endpoints, 0.006);
            }
            _out.WriteLine($"mixture ({t.Quantiles[0, 0, 0]:F3}, {t.Quantiles[1, 0, 0]:F3}); R's endpoints ±{r}");
        }
    }

    [Fact]
    public void TrainingFit_RefusesOtherData_AndProbabilitiesOutsideTheOpenInterval()
    {
        (BassData data, BassFit fit) = BassPredictTests.FitOf(BassPointwiseTests.Kind.Scalar);
        Assert.Throws<ArgumentException>(() => fit.TrainingFit(BassPointwiseTests.Data(BassPointwiseTests.Kind.Scalar, seed: 6)));
        Assert.Throws<ArgumentOutOfRangeException>(() => fit.TrainingFit(data, [0.0, 0.5]));
        Assert.Throws<ArgumentOutOfRangeException>(() => fit.TrainingFit(data, [0.5, 1.0]));
        Assert.Equal(new[] { 0.1, 0.9 }, fit.TrainingFit(data, [0.1, 0.9]).Probabilities);
    }

    // ── Metadata ──────────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Metadata_RoundTripsWhatTheFitHolds()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassMetadata meta = fit.Metadata();
        Assert.Same(fit.Config, meta.Config);
        Assert.Same(fit.Resolved, meta.Resolved);
        Assert.Same(fit.Identity, meta.Identity);
        Assert.Same(fit.Diagnostics, meta.Diagnostics);
        Assert.Null(meta.Evidence);
        Assert.Null(meta.Diagnostics.Tempering);

        // The data's shape, R's print: three variables, one categorical, 60 cases, a scalar response.
        Assert.Equal((60, 1), (meta.Data.Cases, meta.Data.Locations));
        Assert.Equal(new[] { "x0", "x1", "z" }, meta.Data.Predictors.Select(p => p.Name));
        Assert.Equal(1, meta.Data.CategoricalPredictors);
        Assert.Equal(new[] { "a", "b", "c" }, ((TrainingLevels)meta.Data.Predictors[2]).Levels);
        Assert.Equal((0.0, 1.0), (((RangeScale)meta.Data.Predictors[0]).Min, ((RangeScale)meta.Data.Predictors[0]).Max));
        Assert.Empty(meta.Data.OutputVariables);

        // The run's lengths.
        BassSampling s = fit.Config.Sampling;
        Assert.Equal(new int?[] { s.Chains, s.BurnIn, s.Samples, s.Thin }, new[] { meta.Lengths.Chains, meta.Lengths.BurnIn, meta.Lengths.Samples, meta.Lengths.Thin });
        Assert.Null(meta.Lengths.Smc);
        Assert.Equal(fit.Draws.Select(c => c.Count).ToArray(), meta.Lengths.DrawsPerChain);
        Assert.Equal(fit.Draws.Sum(c => c.Count), meta.Lengths.Draws);
        Assert.Same(meta.Config, fit.Summary().Metadata.Config);

        // A direct functional fit's output variables and grid size.
        BassMetadata functional = SensitivityFixture.FitFunctional().Metadata();
        Assert.Equal((40, 7), (functional.Data.Cases, functional.Data.Locations));
        RangeScale t = Assert.Single(functional.Data.OutputVariables);
        Assert.Equal(("t", 0.0, 1.0), (t.Name, t.Min, t.Max));
        Assert.Equal(0, functional.Data.CategoricalPredictors);
    }

    [Fact]
    public void Metadata_CarriesTheEvidenceAndTheTemperingReadout()
    {
        (BassConfig config, BassData data) = BassTemperingTests.RunCase("reference", samples: 30);
        config = config with { Sampling = config.Sampling with { BurnIn = 60, Chains = 2 } };
        BassFit fit = Bass.Fit(config, data);
        BassMetadata meta = fit.Metadata();
        Assert.NotNull(fit.Evidence);
        Assert.Same(fit.Evidence, meta.Evidence);
        Assert.Same(fit.Diagnostics.Tempering, meta.Diagnostics.Tempering);
        Assert.NotNull(meta.Diagnostics.Tempering);
        Assert.Same(fit.Resolved.AutoPicks, meta.Resolved.AutoPicks);
        _out.WriteLine($"ln p(y) = {meta.Evidence!.LogEvidence:F3} ± {meta.Evidence.MonteCarloError:F3}, levels [{string.Join(", ", meta.Diagnostics.Tempering!.Levels)}]");
    }
}
