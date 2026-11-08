using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The conformal predictive distribution on BASS fits (E16, slice 4). BASS's log CDF seam under every noise law: the units' at the
/// training cases is the candidates' there, each draw's CDF is its density's integral, and what does not order one value is refused.
/// With the qualification switch (<c>BASS_QUALIFICATION</c> = 1), the reading the slice turns on: whether the in-sample predictive CDF
/// stays monotone in the candidate's y under the add-one-in weights, on data from outside every model (<see cref="BassCandidateTests.Draw"/>),
/// for a healthy Gaussian fit, a Student t fit, a chain stopped early and a fit that never sees x₀, one fact per fit.
/// </summary>
public sealed class BassConformalDistributionTests(ITestOutputHelper output)
{
    private const int Train = 99, Tests = 10, Replications = 5;

    // The laws whose CDF is a table of tails built per draw (NormalVarianceMeanMixture), its log tails held to 10⁻⁵ relative, so two
    // reads of F to 2·10⁻⁵.
    private static bool Tabled(string law) => law is "normal-wald" or "skewed t" or "generalized hyperbolic";

    [Theory]
    [MemberData(nameof(BassCandidateTests.Laws), MemberType = typeof(BassCandidateTests))]
    public void LogCdf_AtTrainingCases_TheUnitsAreTheCandidates(string law)
    {
        // CB-4 for the CDF: the units' log CDF (the training surface) and the candidates' at the training cases and responses (the
        // new-case surface) read one distribution, to the surfaces' rounding.
        BassData data = BassStudentTTests.Data();
        BassFit fit = Bass.Fit(BassStudentTTests.Config(BassCandidateTests.Noise(law)), data);
        PointwiseLogDensity training = fit.PointwiseLogDensity(data, cdf: true);
        double[] y = data.Response.Values.ToArray();
        CandidateLogDensity candidates = fit.CandidateLogDensity(data, PredictionCases.Of(data), [.. y.Select(v => new[] { v })], cdf: true);
        Assert.NotNull(training.LogCdf);
        Assert.Null(fit.PointwiseLogDensity(data).LogCdf);
        double worst = 0.0;
        int upper = 0;
        for (int s = 0; s < training.Draws; s++)
            for (int i = 0; i < y.Length; i++)
            {
                double a = training.LogCdf![s, i], b = candidates.LogCdf![s, i];
                // Above the median ln F holds 1 − F's digits, so the gap is read relative to it there.
                bool high = a > -Math.Log(2.0);
                if (high) upper++;
                worst = Math.Max(worst, Math.Abs(a - b) / (high ? Math.Max(-a, 1e-300) : Math.Max(1.0, -a)));
            }
        output.WriteLine($"{law}: {training.Draws} draws × {y.Length} units, {upper} entries above the median; largest relative gap {worst:E2}");
        Assert.True(worst <= 1e-9, $"{law}: {worst:E2}");
    }

    [Theory]
    [MemberData(nameof(BassCandidateTests.Laws), MemberType = typeof(BassCandidateTests))]
    public void LogCdf_IsEachDrawsDensityIntegrated(string law)
    {
        // At a new case (the first of two; the second, on one value, keeps the predictors from being constant), F(y | θ_s) − F(y₀ | θ_s)
        // against Simpson's rule on the same draw's density over a fine grid, for each of one chain's 20 draws: the CDF and the density the
        // seam carries are one law, the latent treated alike.
        BassData data = BassStudentTTests.Data();
        BassFit fit = Bass.Fit(BassStudentTTests.Config(BassCandidateTests.Noise(law), samples: 20, chains: 1), data);
        const double h = 0.0005;
        double[] grid = [.. Enumerable.Range(0, 12001).Select(g => -1.5 + h * g)];
        PredictionCases two = PredictionCases.Of([new ContinuousPredictor("x0", [0.3, 0.7]), new ContinuousPredictor("x1", [0.6, 0.2])]);
        CandidateLogDensity c = fit.CandidateLogDensity(data, two, [grid, [0.0]], cdf: true);
        double worst = 0.0;
        for (int s = 0; s < c.Draws; s++)
        {
            double integral = 0.0, f0 = Math.Exp(c.LogCdf![s, 0]);
            for (int g = 2; g < grid.Length; g += 2)
            {
                integral += h / 3 * (Math.Exp(c.Values[s, g - 2]) + 4 * Math.Exp(c.Values[s, g - 1]) + Math.Exp(c.Values[s, g]));
                worst = Math.Max(worst, Math.Abs(Math.Exp(c.LogCdf[s, g]) - f0 - integral));
            }
        }
        output.WriteLine($"{law}: largest gap between the CDF's rise and the density's integral {worst:E2}");
        Assert.True(worst <= (Tabled(law) ? 2e-5 : 1e-6), $"{law}: {worst:E2}");
    }

    [Theory]
    [InlineData("t")]
    [InlineData("quantile")]
    [InlineData("normal-wald")]
    [InlineData("skewed t")]
    [InlineData("generalized hyperbolic")]
    public void LogCdf_HoldsTheUpperTailsDigits(string law)
    {
        // Above the median the seam reads ln F from the survival, so it stays at or below 0 and −expm1(ln F) gives back 1 − F to its
        // relative precision, past the mixtures' tables too (they end at the 10⁻¹² tails); below it, the law's own ln F.
        IContinuousDistribution noise = new VarianceMeanMixture(BassCandidateTests.Noise(law)!).At(1.0, null);
        foreach (double e in new[] { -30.0, -5.0, 0.5, 5.0, 18.0, 32.0 })
        {
            double l = BassCandidates.LogCdf(noise, e), survival = Math.Exp(noise.LogSurvival(e));
            Assert.True(l <= 0.0, $"{law} at {e}: ln F {l:E3}");
            if (e < 0.0) Assert.Equal(noise.LogCdf(e), l);
            else Assert.True(Math.Abs(-SpecialFunctions.Expm1(l) - survival) <= 1e-12 * survival, $"{law} at {e}: 1 − F {-SpecialFunctions.Expm1(l):E6} against S {survival:E6}");
        }
    }

    [Fact]
    public void LogCdf_Refuses_WhatDoesNotOrderOneValue()
    {
        BassData data = BassStudentTTests.Data();
        BassFit fit = Bass.Fit(BassStudentTTests.Config(null), data);
        // A unit of two observations has no CDF; nor does the collapsed content, whose terms leave each unit out.
        int[] pairs = [.. Enumerable.Range(0, data.Response.Observations).Select(i => i / 2)];
        Assert.Throws<ArgumentException>(() => fit.PointwiseLogDensity(data, new UnitPartition(pairs), cdf: true));
        Assert.Throws<ArgumentException>(() => fit.PointwiseLogDensity(data, content: PointwiseContent.CollapsedPredictive, cdf: true));
        // The matrices: the log terms' shape, values in [−∞, 0].
        PointwiseLogDensity training = fit.PointwiseLogDensity(data);
        int s = training.Draws, n = data.Response.Observations;
        Assert.Throws<ArgumentException>(() => training with { LogCdf = new double[s, n + 1] });
        var positive = new double[s, n];
        positive[0, 0] = 1e-12;
        Assert.Throws<ArgumentException>(() => training with { LogCdf = positive });
        var nan = new double[s, n];
        nan[1, 1] = double.NaN;
        Assert.Throws<ArgumentException>(() => training with { LogCdf = nan });
        var lowest = new double[s, n];
        lowest[0, 0] = double.NegativeInfinity;
        Assert.NotNull((training with { LogCdf = lowest }).LogCdf);
        Assert.Throws<ArgumentException>(() => new CandidateLogDensity(new double[2, 1], [[0.0]], []) { LogCdf = new double[2, 1] { { 0.5 }, { -1.0 } } });
        // A weighted fit's units read their weights, which a candidate does not carry.
        BassConfig config = BassStudentTTests.Config(null);
        double[] y = data.Response.Values.ToArray();
        BassModel weighted = BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, weights: Enumerable.Repeat(2.0, y.Length).ToArray())));
        Assert.Throws<NotSupportedException>(() => BassCandidates.UnitLogCdf(weighted, null!, null!, UnitPartition.EachObservation(y.Length)));
    }

    [Fact]
    public void Monotonicity_AHealthyFit() => Monotonicity("a healthy fit");

    [Fact]
    public void Monotonicity_AStudentTFit() => Monotonicity("a Student t fit");

    [Fact]
    public void Monotonicity_AChainStoppedEarly() => Monotonicity("a chain stopped early");

    [Fact]
    public void Monotonicity_AWrongFit() => Monotonicity("a wrong fit");

    internal static BassConfig ConfigFor(string kind) => kind switch
    {
        "a chain stopped early" => BassStudentTTests.Config(null, burnIn: 0, samples: 20, chains: 1),
        "a Student t fit" => BassStudentTTests.Config(BassNoise.StudentT(3, 1.0 / 3), burnIn: 1000, samples: 500),
        _ => BassStudentTTests.Config(null, burnIn: 1000, samples: 500),
    };

    private void Monotonicity(string kind)
    {
        if (!BassQualification.Enabled) return;
        double[] grid = [.. Enumerable.Range(0, 901).Select(g => -4.0 + 0.01 * g)];
        int[] seeds = SeedTree.Derive(2026_1109, Replications);
        var readings = new List<ConformalDistributionTests.Monotonicity>();
        for (int r = 0; r < Replications; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            (double[] x0, double[] x1, double[] y) = BassCandidateTests.Draw(Train, rng);
            (double[] t0, double[] t1, _) = BassCandidateTests.Draw(Tests, rng);
            bool blind = kind == "a wrong fit";
            BassPredictor[] Predictors(double[] a, double[] b)
                => blind ? [new ContinuousPredictor("x1", b)] : [new ContinuousPredictor("x0", a), new ContinuousPredictor("x1", b)];
            var data = new BassData(Predictors(x0, x1), new ScalarResponse(y));
            BassFit fit = Bass.Fit(ConfigFor(kind), data);
            PointwiseLogDensity training = fit.PointwiseLogDensity(data, cdf: true);
            CandidateLogDensity candidates = fit.CandidateLogDensity(data, PredictionCases.Of(Predictors(t0, t1)), [.. Enumerable.Repeat(grid, Tests)], cdf: true);
            int[] members = [.. Enumerable.Range(0, Train)];
            for (int t = 0; t < Tests; t++)
            {
                (ConformalDistributionTests.TailScore[][] units, ConformalDistributionTests.TailScore[] own) = ConformalDistributionTests.Reference(training.LogCdf!, members, candidates.Values, candidates.LogCdf!,
                                                                                       candidates.Offset(t), grid.Length);
                readings.Add(ConformalDistributionTests.Monotonicity.Read(units, own));
            }
        }
        ConformalDistributionTests.Report(output, kind, readings, Train + 1);
    }
}
