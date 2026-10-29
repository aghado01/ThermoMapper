using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Basis prediction (stage 5, O3; D07; <see cref="BassBasisFit.Predict"/>), against each component's own direct
/// prediction and dense computations: the latent reconstruction y_m + y_s ∘ K w(x) over the joint draws; the component
/// noise ε drawn for each draw and case and the truncation residual drawn whole and uniformly over the training cases,
/// recovered from the predictive draws through the in-sample residuals' orthogonality to K; draw order and singleton
/// axes; the seed; the covariance in physical units against the draws' own covariance and a dense assembly; the
/// quantiles; the mean surface; and new output locations refused.
/// </summary>
public sealed class BassPredictBasisTests
{
    private readonly ITestOutputHelper _out;
    public BassPredictBasisTests(ITestOutputHelper output) => _out = output;

    private static PredictionCases Cases(double[] z0, double[] z1) => new([new ContinuousColumn("z0", z0), new ContinuousColumn("z1", z1)]);

    private static readonly PredictionCases Four = Cases([0.1, 0.45, 0.9, 1.3], [0.2, 0.7, 0.5, 0.4]);

    [Fact]
    public void Latent_ReconstructsEachComponentsOwnPrediction()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        BassDecomposition d = fit.Decomposition;
        BassBasisPrediction p = fit.Predict(Four, new BassPredictOptions { Predictive = false });
        int s = p.DrawCount, r = d.Rank, m = d.Locations;
        Assert.Equal((fit.JointDraws.Count, 4, m), (s, p.Cases, p.Locations));
        Assert.Equal([false, false, false, true], p.CaseOutsideRange);
        for (int k = 0; k < r; k++)
        {
            int[] draws = [.. Enumerable.Range(0, s).Select(j => fit.JointDraws.Draw(k, j))];
            BassPrediction own = fit.Components[k].Predict(fit.ComponentData(k), Four, new BassPredictOptions { Draws = draws, Predictive = false });
            for (int j = 0; j < s; j++)
                for (int i = 0; i < 4; i++)
                {
                    Assert.Equal(own.Latent[j, i, 0], p.Coefficients[j, i, k]);
                    Assert.Equal(own.NoiseVariances[j], p.NoiseVariances[j, k]);
                }
        }
        for (int j = 0; j < s; j++)
            for (int i = 0; i < 4; i++)
                for (int l = 0; l < m; l++)
                {
                    double kw = Enumerable.Range(0, r).Sum(k => d.Map.Basis(l, k) * p.Coefficients[j, i, k]);
                    Assert.Equal(d.Map.Offsets[l] + d.Map.Scales[l] * kw, p.Latent[j, i, l], 12);
                }
    }

    [Fact]
    public void Predictive_DrawsComponentNoisePerCase_AndWholeResidualsOverEveryCase()
    {
        // Under in-sample residuals Kᵀe = 0 (each is off the span of K, whose columns are orthogonal with norms d_k), so
        // the predictive draw's standardized deviation (y − latent)/y_s = Kε + e splits into ε = D⁻²Kᵀ(·) and e.
        BassBasisFit fit = BassOutputComponentsTests.FitOf(TruncationResidual.InSample);
        BassDecomposition d = fit.Decomposition;
        int r = d.Rank, m = d.Locations, n = d.Cases;
        var cases = Cases([.. Enumerable.Repeat(0.3, 6)], [.. Enumerable.Repeat(0.6, 6)]);   // six identical cases
        BassBasisPrediction p = fit.Predict(cases, new BassPredictOptions { Seed = 17 });
        int s = p.DrawCount;
        var epsilon = new double[s, 6, r];
        var drawn = new int[s, 6];
        double worstResidual = 0.0;
        for (int j = 0; j < s; j++)
            for (int i = 0; i < 6; i++)
            {
                double[] u = [.. Enumerable.Range(0, m).Select(l => (p.Predictive![j, i, l] - p.Latent[j, i, l]) / d.Map.Scales[l])];
                for (int k = 0; k < r; k++)
                {
                    double kk = Enumerable.Range(0, m).Sum(l => d.Map.Basis(l, k) * d.Map.Basis(l, k));
                    epsilon[j, i, k] = Enumerable.Range(0, m).Sum(l => d.Map.Basis(l, k) * u[l]) / kk / Math.Sqrt(p.NoiseVariances[j, k]);
                }
                double[] e = [.. Enumerable.Range(0, m).Select(l => u[l] - Enumerable.Range(0, r).Sum(k => d.Map.Basis(l, k) * epsilon[j, i, k] * Math.Sqrt(p.NoiseVariances[j, k])))];
                // e is one training case's whole residual vector.
                (int match, double gap) = Enumerable.Range(0, n)
                    .Select(c => (c, Enumerable.Range(0, m).Max(l => Math.Abs(e[l] - d.Truncation.Values[c, l]))))
                    .MinBy(t => t.Item2);
                drawn[j, i] = match;
                worstResidual = Math.Max(worstResidual, gap);
            }
        _out.WriteLine($"largest gap of a residual to its training case: {worstResidual:E1}");
        Assert.InRange(worstResidual, 0, 1e-10);

        // ε: standard normal, independent across the cases within a draw (R recycles one per draw).
        double[] all = epsilon.Cast<double>().ToArray();
        Assert.InRange(all.Average(), -4 / Math.Sqrt(all.Length), 4 / Math.Sqrt(all.Length));
        Assert.InRange(all.Average(v => v * v), 1 - 6 * Math.Sqrt(2.0 / all.Length), 1 + 6 * Math.Sqrt(2.0 / all.Length));
        double correlation = Enumerable.Range(0, s).Sum(j => epsilon[j, 0, 0] * epsilon[j, 1, 0])
                           / Math.Sqrt(Enumerable.Range(0, s).Sum(j => epsilon[j, 0, 0] * epsilon[j, 0, 0]) * Enumerable.Range(0, s).Sum(j => epsilon[j, 1, 0] * epsilon[j, 1, 0]));
        Assert.InRange(Math.Abs(correlation), 0, 4 / Math.Sqrt(s));

        // e: every training case is drawn, about uniformly, and independently across the cases of a draw.
        var counts = new int[n];
        foreach (int c in drawn) counts[c]++;
        double expected = s * 6.0 / n;
        _out.WriteLine($"residual cases drawn: {string.Join(" ", counts)} (expected {expected:F0} each)");
        Assert.All(counts, c => Assert.InRange(c, expected - 5 * Math.Sqrt(expected), expected + 5 * Math.Sqrt(expected)));
        int same = Enumerable.Range(0, s).Count(j => drawn[j, 0] == drawn[j, 1]);
        Assert.InRange(same, 0, s / (double)n + 5 * Math.Sqrt(s / (double)n));
    }

    [Fact]
    public void DrawOrder_SingletonAxes_AndTheSeed()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        int m = fit.Decomposition.Locations;
        BassBasisPrediction all = fit.Predict(Four, new BassPredictOptions { Predictive = false });
        int[] order = [30, 2, 17];
        BassBasisPrediction p = fit.Predict(Four, new BassPredictOptions { Draws = order, Seed = 3 });
        Assert.Equal(order, p.Draws);
        Assert.Equal((3, 4, m), (p.Latent.GetLength(0), p.Latent.GetLength(1), p.Latent.GetLength(2)));
        for (int j = 0; j < 3; j++)
            for (int i = 0; i < 4; i++)
                for (int l = 0; l < m; l++) Assert.Equal(all.Latent[order[j], i, l], p.Latent[j, i, l]);
        BassBasisPrediction single = fit.Predict(Cases([0.5], [0.5]), new BassPredictOptions { Draws = [4] });
        Assert.Equal((1, 1, m), (single.Predictive!.GetLength(0), single.Predictive.GetLength(1), single.Predictive.GetLength(2)));
        Assert.Equal((1, 1, fit.Decomposition.Rank), (single.Coefficients.GetLength(0), single.Coefficients.GetLength(1), single.Coefficients.GetLength(2)));

        BassBasisPrediction again = fit.Predict(Four, new BassPredictOptions { Draws = order, Seed = 3 });
        BassBasisPrediction other = fit.Predict(Four, new BassPredictOptions { Draws = order, Seed = 4 });
        Assert.Equal(p.Predictive!.Cast<double>(), again.Predictive!.Cast<double>());
        Assert.NotEqual(p.Predictive.Cast<double>(), other.Predictive!.Cast<double>());
    }

    [Fact]
    public void OneLocation_KeepsItsAxis()
    {
        (double[,] y, BassPredictor[] predictors) = BassOutputComponentsTests.Problem(m: 1);
        BassDecomposition d = BassDecomposition.FromPca(y);
        BassBasisFit fit = BassComponents.Fit(d, predictors, BassOutputComponentsTests.Config());
        BassBasisPrediction p = fit.Predict(Cases([0.5], [0.5]));
        Assert.Equal((fit.JointDraws.Count, 1, 1), (p.Predictive!.GetLength(0), p.Predictive.GetLength(1), p.Predictive.GetLength(2)));
    }

    [Fact]
    public void Covariance_InPhysicalUnits_MatchesDenseComputations()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        BassDecomposition d = fit.Decomposition;
        int r = d.Rank, m = d.Locations;
        BassBasisPrediction p = fit.Predict(Four, new BassPredictOptions { Seed = 2 });
        int s = p.DrawCount;
        for (int i = 0; i < 4; i++)
        {
            // Latent: exactly the weighted covariance of the latent draws, which are linear in w.
            double[,] latent = p.Covariance(i, predictive: false);
            double[] mean = [.. Enumerable.Range(0, m).Select(l => Enumerable.Range(0, s).Sum(j => p.Weights[j] * p.Latent[j, i, l]))];
            for (int l = 0; l < m; l++)
                for (int q = 0; q < m; q++)
                {
                    double dense = Enumerable.Range(0, s).Sum(j => p.Weights[j] * (p.Latent[j, i, l] - mean[l]) * (p.Latent[j, i, q] - mean[q]));
                    Assert.Equal(dense, latent[l, q], 1e-11 * (1 + Math.Abs(dense)));
                }
            // Predictive: the latent part plus y_s,l·y_s,q·(Σ_k K_lk K_qk E[σ_k²] + Cov[e]_lq), assembled by hand.
            double[,] predictive = p.Covariance(i);
            double[,] truncation = d.Truncation.Covariance();
            double[] noise = [.. Enumerable.Range(0, r).Select(k => Enumerable.Range(0, s).Sum(j => p.Weights[j] * p.NoiseVariances[j, k]))];
            for (int l = 0; l < m; l++)
                for (int q = 0; q < m; q++)
                {
                    double extra = d.Map.Scales[l] * d.Map.Scales[q]
                                 * (Enumerable.Range(0, r).Sum(k => d.Map.Basis(l, k) * d.Map.Basis(q, k) * noise[k]) + truncation[l, q]);
                    Assert.Equal(latent[l, q] + extra, predictive[l, q], 1e-11 * (1 + Math.Abs(predictive[l, q])));
                }
        }
        // And the predictive draws' own covariance agrees within Monte Carlo error.
        BassBasisPrediction many = fit.Predict(Cases([.. Enumerable.Repeat(0.45, 50)], [.. Enumerable.Repeat(0.7, 50)]), new BassPredictOptions { Seed = 5 });
        double[,] analytic = many.Covariance(0);
        int count = many.DrawCount * 50;
        double worst = 0.0;
        foreach ((int l, int q) in new[] { (0, 0), (3, 3), (7, 7), (2, 9) })
        {
            double ml = 0, mq = 0, lq = 0;
            for (int j = 0; j < many.DrawCount; j++)
                for (int i = 0; i < 50; i++)
                {
                    ml += many.Predictive![j, i, l];
                    mq += many.Predictive[j, i, q];
                    lq += many.Predictive[j, i, l] * many.Predictive[j, i, q];
                }
            double empirical = lq / count - ml / count * (mq / count);
            double scale = Math.Sqrt(analytic[l, l] * analytic[q, q]);
            worst = Math.Max(worst, Math.Abs(empirical - analytic[l, q]) / scale);
        }
        _out.WriteLine($"largest Monte Carlo gap in the covariance, relative: {worst:F3}");
        Assert.InRange(worst, 0, 0.1);
    }

    [Fact]
    public void Quantile_IsTheWeightedDrawsQuantile()
    {
        // The joint draws weigh equally, so the shared weighted quantile is type 7 over the predictive draws.
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        BassBasisPrediction p = fit.Predict(Four, new BassPredictOptions { Seed = 8 });
        double[,] q = p.Quantile(0.9);
        for (int i = 0; i < 4; i++)
            for (int l = 0; l < p.Locations; l += 4)
            {
                double[] sorted = [.. Enumerable.Range(0, p.DrawCount).Select(j => p.Predictive![j, i, l]).Order()];
                double position = 0.9 * (p.DrawCount - 1);
                int lo = (int)Math.Floor(position);
                double type7 = sorted[lo] + (position - lo) * (sorted[Math.Min(lo + 1, sorted.Length - 1)] - sorted[lo]);
                Assert.Equal(type7, q[i, l], 12);
            }
        Assert.Throws<InvalidOperationException>(() => fit.Predict(Four, new BassPredictOptions { Predictive = false }).Quantile(0.9));
    }

    [Fact]
    public void MeanSurface_ReconstructsEachComponentsMeanSurface()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        BassDecomposition d = fit.Decomposition;
        BassBasisPrediction p = fit.Predict(Four, new BassPredictOptions { Predictive = false, MeanSurface = true });
        int s = p.DrawCount, r = d.Rank;
        var w = new double[4, r];
        for (int k = 0; k < r; k++)
        {
            int[] draws = [.. Enumerable.Range(0, s).Select(j => fit.JointDraws.Draw(k, j))];
            BassPrediction own = fit.Components[k].Predict(fit.ComponentData(k), Four, new BassPredictOptions { Draws = draws, Predictive = false, MeanSurface = true });
            for (int i = 0; i < 4; i++) w[i, k] = own.MeanSurface![i, 0];
        }
        for (int i = 0; i < 4; i++)
            for (int l = 0; l < d.Locations; l++)
                Assert.Equal(d.Map.Offsets[l] + d.Map.Scales[l] * Enumerable.Range(0, r).Sum(k => d.Map.Basis(l, k) * w[i, k]), p.MeanSurface![i, l], 12);
    }

    [Fact]
    public void Refusals_SayWhy()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        Assert.Contains("training locations", Assert.Throws<ArgumentException>(() =>
            fit.Predict(Four, new BassPredictOptions { OutputLocations = [new ContinuousColumn("t", [0.5])] })).Message);
        Assert.Contains("not among", Assert.Throws<ArgumentException>(() => fit.Predict(Four, new BassPredictOptions { Draws = [100000] })).Message);
        Assert.Contains("'z1'", Assert.Throws<ArgumentException>(() => fit.Predict(new PredictionCases([new ContinuousColumn("z0", [0.5])]))).Message);
        Assert.Throws<ArgumentOutOfRangeException>(() => fit.Predict(Four).Covariance(4));
    }
}
