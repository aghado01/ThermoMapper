using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// CB-4(c), the latent-posterior reading (E16, the conformal layer's slice 12b): every training unit reads its retained observation latents
/// and the candidate reads its own drawn from their conditional given its value, under the add-one-in weights of its integrated density.
/// On a finite model with exact weighted draws the layer computes the latent-posterior score of every unit (<c>latentScore_some_draws</c>,
/// <c>latentScore_none_draws</c>) and CB-1 holds by enumeration. BARS's draws of a candidate's latent average to the conditional mean,
/// against quadrature (<c>drawnLatent_mean_var</c>), and its seams pair under the option and are refused where they would mix. With the
/// qualification switch (<c>BASS_QUALIFICATION</c> = 1): coverage and width of the two readings on BARS t-weight fits under a right and a
/// wrong noise, and the cost row.
/// </summary>
public sealed class LatentPosteriorReadingTests(ITestOutputHelper output)
{
    // ---- The finite model: θ ∈ {0, 1}, a latent bit per unit, y ∈ {0, 1, 2} ----

    private static readonly double[] PriorTheta = [0.5, 0.5];
    private static readonly double[,] Q = { { 0.7, 0.3 }, { 0.4, 0.6 } };   // q(v | θ)
    private static readonly double[,,] F =                                     // f(y | θ, v)
    {
        { { 0.7, 0.2, 0.1 }, { 0.1, 0.3, 0.6 } },
        { { 0.5, 0.4, 0.1 }, { 0.2, 0.2, 0.6 } },
    };

    private static double G(int theta, int y) => Q[theta, 0] * F[theta, 0, y] + Q[theta, 1] * F[theta, 1, y];

    private static double ConditionalMean(int theta, int y)
        => (Q[theta, 0] * F[theta, 0, y] * F[theta, 0, y] + Q[theta, 1] * F[theta, 1, y] * F[theta, 1, y]) / G(theta, y);

    private static readonly LatentReport[] Conditioned = [new("θ", LatentScope.Global, LatentTreatment.Conditioned),
                                                          new("v", LatentScope.Observation, LatentTreatment.Conditioned)];

    private static readonly LatentReport[] Integrated = [new("θ", LatentScope.Global, LatentTreatment.Conditioned),
                                                         new("v", LatentScope.Observation, LatentTreatment.Integrated)];

    // The exact joint posterior of (θ, v₁…vₙ) given the units, as weighted draws: every combination once, weighted by p(θ) Πᵢ q(vᵢ | θ) f(yᵢ | θ, vᵢ).
    // The training density reads each unit's retained latent, or each unit integrated; the candidates read g for the weights and, conditioned,
    // the conditional mean E[f | θ, y].
    private static (PointwiseLogDensity Conditioned, PointwiseLogDensity Integrated, CandidateLogDensity Candidates, CandidateLogDensity Plain) Exact(int[] y)
    {
        int n = y.Length, draws = 2 << n;
        var cond = new double[draws, n];
        var integ = new double[draws, n];
        var weights = new double[draws];
        var values = new double[draws, 3];
        var means = new double[draws, 3];
        for (int s = 0; s < draws; s++)
        {
            int theta = s >> n;
            double w = PriorTheta[theta];
            for (int i = 0; i < n; i++)
            {
                int v = (s >> i) & 1;
                w *= Q[theta, v] * F[theta, v, y[i]];
                cond[s, i] = Math.Log(F[theta, v, y[i]]);
                integ[s, i] = Math.Log(G(theta, y[i]));
            }
            weights[s] = w;
            for (int g = 0; g < 3; g++)
            {
                values[s, g] = Math.Log(G(theta, g));
                means[s, g] = Math.Log(ConditionalMean(theta, g));
            }
        }
        double total = weights.Sum();
        double[] normalized = [.. weights.Select(w => w / total)];
        var units = UnitPartition.EachObservation(n);
        return (new PointwiseLogDensity(cond, units, Conditioned) { Weights = normalized },
                new PointwiseLogDensity(integ, units, Integrated) { Weights = normalized },
                new CandidateLogDensity(values, [[0.0, 1.0, 2.0]], Conditioned) { Weights = normalized, ConditionedValues = means },
                new CandidateLogDensity(values, [[0.0, 1.0, 2.0]], Integrated) { Weights = normalized });
    }

    // CB-4(c) of unit z among all units, in closed form: Σ_θ p(θ | all) Σ_v q f²/g (latentScore_eq).
    private static double LatentScore(int[] all, int z)
    {
        var post = new double[2];
        for (int theta = 0; theta < 2; theta++)
        {
            post[theta] = PriorTheta[theta];
            foreach (int a in all) post[theta] *= G(theta, a);
        }
        double total = post.Sum();
        return (post[0] * ConditionalMean(0, z) + post[1] * ConditionalMean(1, z)) / total;
    }

    [Fact]
    public void OnExactDraws_TheLayerReadsTheLatentPosteriorScore_OfEveryUnit()
    {
        foreach (int[] y in new[] { new[] { 0, 2, 1 }, new[] { 2, 2, 0 }, new[] { 1, 1, 1 } })
        {
            (PointwiseLogDensity training, _, CandidateLogDensity candidates, _) = Exact(y);
            for (int g = 0; g < 3; g++)
            {
                ConformalScores scores = ConformalBayes.Scores(training, candidates, g, latents: LatentReading.Conditioned);
                int[] all = [.. y, g];
                for (int i = 0; i < y.Length; i++) Assert.Equal(Math.Log(LatentScore(all, y[i])), scores.Scores[i], 12);
                Assert.Equal(Math.Log(LatentScore(all, g)), scores.Candidate, 12);
            }
        }
    }

    [Fact]
    public void Cb1_ByEnumeration_WithTheLatentPosteriorScore()
    {
        // Every data set of three units and a candidate, at k 1 of m 4: coverage at least 3/4 under an iid law outside the model and
        // under the model's own law, as the integrated score's is.
        double[] outside = [0.5, 0.3, 0.2];
        double Own(int[] all)
        {
            double sum = 0.0;
            for (int theta = 0; theta < 2; theta++)
            {
                double p = PriorTheta[theta];
                foreach (int a in all) p *= G(theta, a);
                sum += p;
            }
            return sum;
        }
        double[] coveredOutside = new double[2], coveredOwn = new double[2];
        var options = new ConformalOptions { Alpha = 0.3, ObservationLatents = LatentReading.Conditioned };
        for (int code = 0; code < 27; code++)
        {
            int[] y = [code % 3, code / 3 % 3, code / 9];
            (PointwiseLogDensity conditioned, PointwiseLogDensity integrated, CandidateLogDensity candidates, CandidateLogDensity plain) = Exact(y);
            ConformalSet c = ConformalBayes.FullStep(conditioned, candidates, options).Sets[0];
            ConformalSet b = ConformalBayes.FullStep(integrated, plain, options with { ObservationLatents = LatentReading.Integrated }).Sets[0];
            Assert.Equal(1, c.Level);
            for (int g = 0; g < 3; g++)
            {
                int[] all = [.. y, g];
                double pOutside = all.Aggregate(1.0, (p, a) => p * outside[a]), pOwn = Own(all);
                if (c.AcceptsPlain(g)) { coveredOutside[1] += pOutside; coveredOwn[1] += pOwn; }
                if (b.AcceptsPlain(g)) { coveredOutside[0] += pOutside; coveredOwn[0] += pOwn; }
            }
        }
        output.WriteLine($"coverage, bound 0.75: integrated {coveredOutside[0]:F4} outside the model and {coveredOwn[0]:F4} under its own law; "
                         + $"latent-posterior {coveredOutside[1]:F4} and {coveredOwn[1]:F4}");
        Assert.All(coveredOutside.Concat(coveredOwn), v => Assert.True(v >= 0.75 - 1e-12, $"coverage {v:F4}"));
    }

    [Fact]
    public void WhatTheReadingDoesNotServe_IsRefused()
    {
        (PointwiseLogDensity conditioned, PointwiseLogDensity integrated, CandidateLogDensity candidates, CandidateLogDensity plain) = Exact([0, 2, 1]);
        var on = new ConformalOptions { Alpha = 0.3, ObservationLatents = LatentReading.Conditioned };
        // The default refuses conditioned training latents, naming the option.
        Assert.Contains("ObservationLatents", Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(conditioned, candidates)).Message);
        // Under the option: integrated training latents; a seam with no conditioned reading; a mixed report; the margin; Hoff's score.
        Assert.Contains("integrated", Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(integrated, candidates, on)).Message);
        Assert.Contains("ConditionedValues", Assert.Throws<ArgumentException>(
            () => ConformalBayes.FullStep(conditioned, candidates with { ConditionedValues = null }, on)).Message);
        var mixed = new CandidateLogDensity(candidates.Values, candidates.Grids, Integrated) { Weights = candidates.Weights, ConditionedValues = candidates.ConditionedValues };
        Assert.Contains("not treated as the training units'", Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(conditioned, mixed, on)).Message);
        Assert.Contains("ConditionedValues", Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(conditioned, plain, on)).Message);
        Assert.Contains("carried bounds", Assert.Throws<ArgumentException>(
            () => ConformalBayes.FullStep(conditioned, candidates, on with { MarginBounds = MarginBounds.Carried })).Message);
        Assert.Contains("Hoff", Assert.Throws<ArgumentException>(
            () => ConformalBayes.FullStep(conditioned, candidates, on with { Score = FullStepScore.LeaveOneOut })).Message);
        // The predictive distribution keeps refusing conditioned latents, and its options overload refuses the reading.
        Assert.Throws<ArgumentException>(() => ConformalBayes.Distribution(conditioned, candidates));
        Assert.Contains("CDF score", Assert.Throws<ArgumentException>(() => ConformalBayes.Distribution(integrated, plain, on)).Message);
        // A conditioned matrix of another shape.
        Assert.Throws<ArgumentException>(() => candidates with { ConditionedValues = new double[3, 3] });
    }

    // ---- BARS ----

    private static BarsEnsemble Bars(IWeightResampler? resampler = null, BarsNoise? noise = null) =>
        new(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0), resampler: resampler, noise: noise);

    private static BarsDraw Draw(double sigma2, double skew) => new(new KnotConfig([0.5]), [0.0], sigma2, null, new PeakDraw(double.NaN, double.NaN, []), skew);

    // ∫ h(u) du over ln v by Simpson's rule on [−30, 15].
    private static double Simpson(Func<double, double> h)
    {
        const int steps = 6000;
        double a = -30.0, b = 15.0, step = (b - a) / steps, sum = h(a) + h(b);
        for (int k = 1; k < steps; k++) sum += (k % 2 == 1 ? 4.0 : 2.0) * h(a + k * step);
        return sum * step / 3.0;
    }

    [Theory]
    [InlineData("t weights")]
    [InlineData("t noise")]
    [InlineData("skewed gig")]
    public void TheDrawnLatent_AveragesToTheConditionalMean(string law)
    {
        // At one draw and one residual, f at 40,000 latents drawn from their conditional averages to E[f | θ, y] = ∫ q f²/∫ q f, read by
        // quadrature over ln v from the latent's prior, which the draw does not use.
        (BarsEnsemble ensemble, Func<double, double> logPrior, double c) = law switch
        {
            "t weights" => (Bars(new StudentTWeights(4.0)), (Func<double, double>)(v => new StudentTMixing(4.0).LogPrior(v)), 1.0),
            "t noise" => (Bars(noise: BarsNoise.StudentT(5.0)), v => BarsNoise.StudentT(5.0).Mixing.LogPrior(v), BarsNoise.StudentT(5.0).VarianceScale),
            _ => (Bars(noise: BarsNoise.Skewed(new GigMixing(new GigParameters(0.7, 1.2, 0.8)), Normal.MeanSd(0.3, 0.5))),
                  v => new GigMixing(new GigParameters(0.7, 1.2, 0.8)).LogPrior(v), 1.0),
        };
        const double sigma2 = 0.5, residual = 1.3;
        double skew = law == "skewed gig" ? 0.4 : 0.0;
        BarsDraw d = Draw(sigma2, skew);
        // The latent v, with weight w = 1/(c·v): N(e; m·v, σ²·c·v), m = √σ²·β.
        double m = Math.Sqrt(sigma2) * skew;
        double F(double v) => Math.Exp(-0.5 * Math.Log(2.0 * Math.PI * sigma2 * c * v) - 0.5 * (residual - m * v) * (residual - m * v) / (sigma2 * c * v));
        double first = Simpson(u => Math.Exp(logPrior(Math.Exp(u)) + u) * F(Math.Exp(u)));
        double second = Simpson(u => Math.Exp(logPrior(Math.Exp(u)) + u) * F(Math.Exp(u)) * F(Math.Exp(u)));
        double reference = second / first;

        const int draws = 40000;
        double[] weights = ensemble.LatentDraws()(d, [.. Enumerable.Repeat(residual, draws)], new Xoshiro256PlusPlus(7));
        double[] f = [.. weights.Select(w => Math.Exp(ensemble.ConditionedTerm(d, residual, w)))];
        double mean = f.Average(), sd = Math.Sqrt(f.Sum(v => (v - mean) * (v - mean)) / (draws - 1));
        output.WriteLine($"{law}: the drawn reading's mean {mean:F6} (se {sd / Math.Sqrt(draws):F6}) against the conditional mean {reference:F6}");
        Assert.True(Math.Abs(mean - reference) <= 4.5 * sd / Math.Sqrt(draws), $"{law}: {mean} against {reference}");
    }

    private static (double[] X, double[] Y) Data(int n, int seed, bool contaminated)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] x = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble()).Order()];
        double[] y = [.. x.Select(v => contaminated ? Truth(v, rng) : Math.Sin(2.0 * Math.PI * v) + 0.5 * v + 0.2 * StudentT(4.0, rng))];
        return (x, y);
    }

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double StudentT(double nu, Xoshiro256PlusPlus rng)
    {
        double chi = 0.0;
        for (int k = 0; k < (int)nu; k++) { double z = Gaussian(rng); chi += z * z; }
        return Gaussian(rng) / Math.Sqrt(chi / nu);
    }

    // Slice 5's truth outside every model: sin(2πx) + x/2 with noise 0.2ζ plus, one time in seven, 0.8·Exp(1).
    private static double Truth(double x, Xoshiro256PlusPlus rng)
        => Math.Sin(2.0 * Math.PI * x) + 0.5 * x + 0.2 * Gaussian(rng) + (rng.NextDouble() < 1.0 / 7.0 ? -0.8 * Math.Log(1.0 - rng.NextDouble()) : 0.0);

    [Theory]
    [InlineData("t weights")]
    [InlineData("t noise")]
    public void BarsSeams_PairUnderTheOption_AndAreRefusedWhereTheyWouldMix(string law)
    {
        BarsEnsemble ensemble = law == "t weights" ? Bars(new StudentTWeights(4.0)) : Bars(noise: BarsNoise.StudentT(5.0));
        (double[] x, double[] y) = Data(40, seed: 3, contaminated: false);
        BarsResult r = ensemble.Run(x, y, x, chains: 2, masterSeed: 5, burn: 200, samples: 100);
        PointwiseLogDensity conditioned = ensemble.PointwiseLogDensity(r, x, y, observationLatents: LatentTreatment.Conditioned);
        PointwiseLogDensity integrated = ensemble.PointwiseLogDensity(r, x, y);
        double[] grid = [.. Enumerable.Range(0, 61).Select(g => -2.0 + 0.07 * g)];
        CandidateLogDensity drawn = ensemble.CandidateLogDensity(r, [0.31, 0.77], [grid, grid], observationLatents: LatentTreatment.Conditioned, seed: 9);
        CandidateLogDensity plain = ensemble.CandidateLogDensity(r, [0.31, 0.77], [grid, grid]);

        // Each training term is its draw's Normal at its retained latent, written out; the candidates' values stay the integrated terms.
        BarsDraw d0 = r.Draws[0][17];
        double[,] z = new SplineBasis(3).Design(d0.Knots, x);
        double eta = Enumerable.Range(0, d0.Coefficients.Length).Sum(j => z[4, j] * d0.Coefficients[j]);
        double variance = d0.Sigma2 / d0.Weights![4];
        Assert.Equal(-0.5 * Math.Log(2.0 * Math.PI * variance) - 0.5 * (y[4] - eta) * (y[4] - eta) / variance, conditioned.Values[17, 4], 10);
        for (int s = 0; s < drawn.Draws; s += 13)
            for (int k = 0; k < grid.Length; k++) Assert.Equal(plain.Values[s, k], drawn.Values[s, k]);
        Assert.NotNull(drawn.ConditionedValues);
        Assert.Null(plain.ConditionedValues);

        // The full step takes the conditioned pair under the option, and the integrated pair by default.
        var on = new ConformalOptions { ObservationLatents = LatentReading.Conditioned };
        ConformalPrediction latentSets = ConformalBayes.FullStep(conditioned, drawn, on), integratedSets = ConformalBayes.FullStep(integrated, plain);
        for (int t = 0; t < 2; t++)
        {
            Assert.True(latentSets.Sets[t].Width > 0.0 && latentSets.Sets[t].Width < 4.2);
            output.WriteLine($"{law}, case {t}: widths {integratedSets.Sets[t].Width:F3} integrated, {latentSets.Sets[t].Width:F3} latent-posterior");
        }
        // Mixes are refused: conditioned training by default; conditioned training beside integrated candidates; the CDF under the reading.
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(conditioned, plain));
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(conditioned, plain, on));
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(integrated, drawn, on));
        Assert.Throws<ArgumentException>(() => ensemble.PointwiseLogDensity(r, x, y, cdf: true, observationLatents: LatentTreatment.Conditioned));
        Assert.Throws<ArgumentException>(() => ensemble.CandidateLogDensity(r, [0.5], [grid], cdf: true, observationLatents: LatentTreatment.Conditioned));
        Assert.Throws<ArgumentOutOfRangeException>(() => ensemble.PointwiseLogDensity(r, x, y, observationLatents: LatentTreatment.Mixed));
        // One seed repeats the candidates' draws.
        CandidateLogDensity again = ensemble.CandidateLogDensity(r, [0.31, 0.77], [grid, grid], observationLatents: LatentTreatment.Conditioned, seed: 9);
        Assert.Equal(drawn.ConditionedValues![11, 20], again.ConditionedValues![11, 20]);
    }

    [Fact]
    public void ANoiseWithoutAnExactConditionalDraw_IsRefused_AndAModelWithoutLatents_ReadsAlike()
    {
        BarsEnsemble gbp = Bars(noise: BarsNoise.Skewed(new GbpMixing(new GbpParameters(1.0, 0.5, 1.5)), Normal.MeanSd(0.0, 0.5)));
        Assert.Throws<NotSupportedException>(() => gbp.LatentDraws());

        BarsEnsemble normal = Bars();
        (double[] x, double[] y) = Data(30, seed: 4, contaminated: false);
        BarsResult r = normal.Run(x, y, x, chains: 1, masterSeed: 2, burn: 100, samples: 60);
        double[] grid = [.. Enumerable.Range(0, 41).Select(g => -2.0 + 0.1 * g)];
        PointwiseLogDensity training = normal.PointwiseLogDensity(r, x, y, observationLatents: LatentTreatment.Conditioned);
        CandidateLogDensity candidates = normal.CandidateLogDensity(r, [0.4], [grid], observationLatents: LatentTreatment.Conditioned);
        Assert.Null(candidates.ConditionedValues);
        ConformalSet a = ConformalBayes.FullStep(training, candidates).Sets[0];
        ConformalSet b = ConformalBayes.FullStep(training, candidates, new ConformalOptions { ObservationLatents = LatentReading.Conditioned }).Sets[0];
        Assert.Equal(a.PValues, b.PValues);
    }

    // ---- BASS ----

    private static BassPredictor[] Predictors(double[] a, double[] b) => [new ContinuousPredictor("x0", a), new ContinuousPredictor("x1", b)];

    [Fact]
    public void BassSeams_UnderAStudentTNoise_ReadEachUnitAtItsLatent()
    {
        var rng = new Xoshiro256PlusPlus(11);
        (double[] x0, double[] x1, double[] y) = BassCandidateTests.Draw(60, rng);
        var data = new BassData(Predictors(x0, x1), new ScalarResponse(y));
        BassConfig config = BassStudentTTests.Config(BassNoise.StudentT(3, 1.0 / 3), burnIn: 200, samples: 60);
        BassFit fit = Bass.Fit(config, data);
        PointwiseLogDensity conditioned = fit.PointwiseLogDensity(data, observationLatents: LatentTreatment.Conditioned);
        PointwiseLogDensity integrated = fit.PointwiseLogDensity(data);
        Assert.False(conditioned.LeavesUnitsOut);

        // The cases: the first three training inputs, where each draw's surface is its training fit's; grids of seven values each.
        double[][] grids = [.. Enumerable.Range(0, 3).Select(i => Enumerable.Range(0, 7).Select(g => y[i] - 0.6 + 0.2 * g).ToArray())];
        var cases = PredictionCases.Of(Predictors(x0[..3], x1[..3]));
        CandidateLogDensity drawn = fit.CandidateLogDensity(data, cases, grids, observationLatents: LatentTreatment.Conditioned, seed: 21);
        CandidateLogDensity plain = fit.CandidateLogDensity(data, cases, grids);

        // Each training term at its draw's retained latent, written out: N(e | σβv, σ²c·v), e net of the draw's fit. Each candidate's
        // conditioned term, its latent drawn from the noise's conditional on the seam's own stream, in draw, case and grid order.
        BassModel model = BassModel.Create(config, data);
        var at = new BassDrawFit(model, fit.Terms);
        var mixture = (VarianceMeanMixture)model.Observation;
        BassDraw[] pooled = [.. fit.Draws.SelectMany(chain => chain)];
        var f = new double[y.Length];
        var stream = new Xoshiro256PlusPlus(21);
        double Term(double e, double v, double unit, double drift)
        {
            double variance = unit * v, d = e - drift * v;
            return -0.5 * Math.Log(2.0 * Math.PI * variance) - 0.5 * d * d / variance;
        }
        for (int s = 0; s < pooled.Length; s++)
        {
            at.Evaluate(pooled[s], f);
            ObservationLatents latents = pooled[s].LatentsOn(model)!;
            double unit = pooled[s].Sigma2 * mixture.VarianceScale, drift = Math.Sqrt(pooled[s].Sigma2) * latents.Skew;
            if (s % 7 == 0)
                for (int i = 0; i < y.Length; i++)
                    Assert.Equal(Term(y[i] - f[i], latents.Values[i], unit, drift), conditioned.Values[s, i], 10);
            var mixing = (StudentTMixing)mixture.MixingAt(latents);
            for (int t = 0; t < 3; t++)
                for (int g = 0; g < 7; g++)
                {
                    double e = grids[t][g] - f[t];
                    double v = mixing.SampleConditional(e, unit, drift, 1.0, stream);
                    Assert.Equal(Term(e, v, unit, drift), drawn.ConditionedValues![s, drawn.Offset(t) + g], 10);
                    Assert.Equal(plain.Values[s, plain.Offset(t) + g], drawn.Values[s, drawn.Offset(t) + g]);
                }
        }

        // The full step takes the conditioned pair under the option; mixes and the CDF are refused.
        var on = new ConformalOptions { ObservationLatents = LatentReading.Conditioned };
        Assert.Equal(3, ConformalBayes.FullStep(conditioned, drawn, on).Sets.Count);
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(conditioned, plain, on));
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(integrated, drawn, on));
        Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(conditioned, drawn));
        Assert.Throws<ArgumentException>(() => fit.PointwiseLogDensity(data, cdf: true, observationLatents: LatentTreatment.Conditioned));
        Assert.Throws<ArgumentException>(() => fit.CandidateLogDensity(data, cases, grids, cdf: true, observationLatents: LatentTreatment.Conditioned));
    }

    // ---- Behind the qualification switch ----

    [Theory]
    [InlineData("a Student t fit")]
    public void Measurement_OnBassFits_TheLatentPosteriorReading_AgainstTheIntegrated(string kind)
    {
        // Slice 4's Student t fits from 12a's measurement: n 99 (m 100, k 10), 20 fits, 50 held-out cases each from a truth outside every
        // model and 20 more on a grid of step 0.04 for the width, both readings on the same fit and cases.
        if (!BassQualification.Enabled) return;
        const int fits = 20, tests = 50, widths = 20;
        double[] grid = [.. Enumerable.Range(0, 201).Select(g => -4.0 + 0.04 * g)];
        int[] seeds = SeedTree.Derive(2026_1112, fits);
        var covB = new List<double>();
        var covC = new List<double>();
        var diff = new List<double>();
        double widthB = 0.0, widthC = 0.0, least = double.PositiveInfinity;
        int fit = 0;
        foreach (int seed in seeds)
        {
            var rng = new Xoshiro256PlusPlus(seed);
            (double[] x0, double[] x1, double[] y) = BassCandidateTests.Draw(99, rng);
            (double[] t0, double[] t1, double[] yt) = BassCandidateTests.Draw(tests + widths, rng);
            var data = new BassData(Predictors(x0, x1), new ScalarResponse(y));
            BassFit bass = Bass.Fit(BassConformalDistributionTests.ConfigFor(kind), data);
            double[][] grids = [.. yt[..tests].Select(v => new[] { v }), .. Enumerable.Repeat(grid, widths)];
            var cases = PredictionCases.Of(Predictors(t0, t1));
            ConformalPrediction b = ConformalBayes.FullStep(bass.PointwiseLogDensity(data), bass.CandidateLogDensity(data, cases, grids));
            ConformalPrediction c = ConformalBayes.FullStep(bass.PointwiseLogDensity(data, observationLatents: LatentTreatment.Conditioned),
                                                            bass.CandidateLogDensity(data, cases, grids, observationLatents: LatentTreatment.Conditioned, seed: fit++),
                                                            new ConformalOptions { ObservationLatents = LatentReading.Conditioned });
            covB.Add(Enumerable.Range(0, tests).Count(t => b.Sets[t].Accepts(0)) / (double)tests);
            covC.Add(Enumerable.Range(0, tests).Count(t => c.Sets[t].Accepts(0)) / (double)tests);
            for (int t = 0; t < tests; t++) least = Math.Min(least, b.Sets[t].EffectiveSampleSizes![0]);
            for (int t = tests; t < tests + widths; t++)
            {
                widthB += b.Sets[t].Width;
                widthC += c.Sets[t].Width;
                diff.Add(c.Sets[t].Width - b.Sets[t].Width);
            }
        }
        (double Mean, double Se) Of(List<double> v)
        {
            double mean = v.Average();
            return (mean, Math.Sqrt(v.Sum(a => (a - mean) * (a - mean)) / (v.Count - 1) / v.Count));
        }
        (double mb, double sb) = Of(covB);
        (double mc, double sc) = Of(covC);
        (double wd, double wse) = Of(diff);
        double wb = widthB / (fits * widths), wc = widthC / (fits * widths);
        output.WriteLine($"BASS, {kind}: coverage integrated {mb:F4} (se {sb:F4}), latent-posterior {mc:F4} (se {sc:F4}), nominal 0.9000; width {wb:F3} and "
                         + $"{wc:F3}, latent-posterior less integrated {wd:+0.000;-0.000} (se {wse:F3}), {100.0 * wd / wb:+0.0;-0.0}%; least ESS at the truth {least:F1}");
        int count = fits * tests;
        Assert.True(Math.Abs(mb - 0.9) <= 3.5 * Math.Max(sb, Math.Sqrt(0.09 / count)), $"integrated {mb:F4}");
        Assert.True(Math.Abs(mc - 0.9) <= 3.5 * Math.Max(sc, Math.Sqrt(0.09 / count)), $"latent-posterior {mc:F4}");
    }


    [Theory]
    [InlineData("t noise, the model's own")]
    [InlineData("contaminated, outside the model")]
    public void Measurement_OnBarsTWeightFits_TheLatentPosteriorReading_AgainstTheIntegrated(string truth)
    {
        // 20 fits of 99 points under Student t₄ weights, m 100 and k 10; 50 held-out cases each at their true value, 20 on a grid for the
        // width, the same fits and cases for both readings; the se is the fits' spread, the widths' paired.
        if (!BassQualification.Enabled) return;
        bool contaminated = truth.StartsWith("contaminated");
        const int fits = 20, n = 99, tests = 50, widths = 20;
        double[] grid = [.. Enumerable.Range(0, 101).Select(g => -3.0 + 0.08 * g)];
        var covB = new double[fits];
        var covC = new double[fits];
        var widthB = new List<double>();
        var widthC = new List<double>();
        double leastB = double.PositiveInfinity, leastC = double.PositiveInfinity;
        int untrustedB = 0, untrustedC = 0;
        double seconds = 0.0;
        for (int fit = 0; fit < fits; fit++)
        {
            long t0 = System.Diagnostics.Stopwatch.GetTimestamp();
            var rng = new Xoshiro256PlusPlus(20261112 + 97 * fit);
            (double[] x, double[] y) = Data(n, 1000 + fit, contaminated);
            double[] xt = [.. Enumerable.Range(0, tests + widths).Select(_ => rng.NextDouble())];
            double[] yt = [.. xt[..tests].Select(v => contaminated ? Truth(v, rng) : Math.Sin(2.0 * Math.PI * v) + 0.5 * v + 0.2 * StudentT(4.0, rng))];
            double[][] grids = [.. yt.Select(v => new[] { v }), .. Enumerable.Repeat(grid, widths)];
            BarsEnsemble ensemble = Bars(new StudentTWeights(4.0));
            BarsResult r = ensemble.Run(x, y, x, chains: 2, masterSeed: fit + 1, burn: 500, samples: 500);
            ConformalPrediction b = ConformalBayes.FullStep(ensemble.PointwiseLogDensity(r, x, y), ensemble.CandidateLogDensity(r, xt, grids));
            ConformalPrediction c = ConformalBayes.FullStep(ensemble.PointwiseLogDensity(r, x, y, observationLatents: LatentTreatment.Conditioned),
                                                            ensemble.CandidateLogDensity(r, xt, grids, observationLatents: LatentTreatment.Conditioned, seed: fit),
                                                            new ConformalOptions { ObservationLatents = LatentReading.Conditioned });
            int cb = 0, cc = 0;
            for (int t = 0; t < tests; t++)
            {
                if (b.Sets[t].Accepts(0)) cb++;
                if (c.Sets[t].Accepts(0)) cc++;
                leastB = Math.Min(leastB, b.Sets[t].EffectiveSampleSizes![0]);
                leastC = Math.Min(leastC, c.Sets[t].EffectiveSampleSizes![0]);
                untrustedB += b.Sets[t].Margin!.Untrusted(b.Sets[t].Level);
                untrustedC += c.Sets[t].Margin!.Untrusted(c.Sets[t].Level);
            }
            for (int t = tests; t < tests + widths; t++)
            {
                widthB.Add(b.Sets[t].Width);
                widthC.Add(c.Sets[t].Width);
            }
            covB[fit] = cb / (double)tests;
            covC[fit] = cc / (double)tests;
            seconds += (System.Diagnostics.Stopwatch.GetTimestamp() - t0) / (double)System.Diagnostics.Stopwatch.Frequency;
        }
        (double Mean, double Se) Of(IReadOnlyList<double> v)
        {
            double mean = v.Average();
            return (mean, Math.Sqrt(v.Sum(a => (a - mean) * (a - mean)) / (v.Count - 1) / v.Count));
        }
        (double mb, double sb) = Of(covB);
        (double mc, double sc) = Of(covC);
        (double wd, double wse) = Of([.. widthC.Zip(widthB, (cw, bw) => cw - bw)]);
        double wb = widthB.Average(), wc = widthC.Average();
        output.WriteLine($"{truth}: coverage integrated {mb:F4} (se {sb:F4}), latent-posterior {mc:F4} (se {sc:F4}), nominal 0.9000; "
                         + $"width {wb:F3} and {wc:F3}, latent-posterior less integrated {wd:+0.000;-0.000} (se {wse:F3}), {100.0 * wd / wb:+0.0;-0.0}%; "
                         + $"least ESS at the truth {leastB:F1} and {leastC:F1}; values with k̂ above its threshold where the inflated rule accepts, "
                         + $"at the truth, {untrustedB} and {untrustedC}; {seconds:F0} s");
        int count = fits * tests;
        Assert.True(Math.Abs(mb - 0.9) <= 3.5 * Math.Max(sb, Math.Sqrt(0.09 / count)), $"integrated {mb:F4}");
        Assert.True(Math.Abs(mc - 0.9) <= 3.5 * Math.Max(sc, Math.Sqrt(0.09 / count)), $"latent-posterior {mc:F4}");
    }

    [Fact]
    public void CostRow()
    {
        // n 400, S 1,000, one grid of 200 at 20 cases, under t₄ weights: the candidate seam and the full step, each reading.
        if (!BassQualification.Enabled) return;
        (double[] x, double[] y) = Data(400, 77, contaminated: false);
        BarsEnsemble ensemble = Bars(new StudentTWeights(4.0));
        BarsResult r = ensemble.Run(x, y, x, chains: 2, masterSeed: 3, burn: 500, samples: 500);
        double[] xt = [.. Enumerable.Range(0, 20).Select(k => (k + 0.5) / 20.0)];
        double[] grid = [.. Enumerable.Range(0, 200).Select(g => -3.0 + 0.03 * g)];
        double[][] grids = [.. Enumerable.Repeat(grid, 20)];
        double Ms(long from) => (BassCost.Tick() - from) / BassCost.TicksPerSecond * 1e3;
        long t0 = BassCost.Tick();
        PointwiseLogDensity tb = ensemble.PointwiseLogDensity(r, x, y);
        double trainB = Ms(t0);
        t0 = BassCost.Tick();
        PointwiseLogDensity tc = ensemble.PointwiseLogDensity(r, x, y, observationLatents: LatentTreatment.Conditioned);
        double trainC = Ms(t0);
        t0 = BassCost.Tick();
        CandidateLogDensity cb = ensemble.CandidateLogDensity(r, xt, grids);
        double seamB = Ms(t0);
        t0 = BassCost.Tick();
        CandidateLogDensity cc = ensemble.CandidateLogDensity(r, xt, grids, observationLatents: LatentTreatment.Conditioned, seed: 1);
        double seamC = Ms(t0);
        t0 = BassCost.Tick();
        ConformalBayes.FullStep(tb, cb);
        double stepB = Ms(t0);
        t0 = BassCost.Tick();
        ConformalBayes.FullStep(tc, cc, new ConformalOptions { ObservationLatents = LatentReading.Conditioned });
        double stepC = Ms(t0);
        output.WriteLine($"n 400, S {r.Draws.Sum(ch => ch.Count)}, 20 cases of 200 ({(BassCost.Cycles ? "cycles" : "wall")} clock, indicative): "
                         + $"the training density {trainB:F0} ms integrated, {trainC:F0} ms conditioned; the candidates {seamB:F0} ms and {seamC:F0} ms; "
                         + $"the full step {stepB:F0} ms and {stepC:F0} ms");
    }
}
