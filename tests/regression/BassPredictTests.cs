using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Direct prediction (stage 5, O3; D07; <see cref="BassFit.Predict"/>), against independent computations: the stage-4
/// fitted values at the training inputs bit for bit (scalar, categorical, direct functional); the hinges evaluated by hand
/// at new cases and new output locations, extrapolation flagged per case and per location; the requested draw order and
/// singleton axes; noise independent across draws, cases and locations at each draw's σ²; the call's seed; unseen levels;
/// the predictive mixture's quantile; the Rao–Blackwellized mean surface against the g-prior's closed form; the draw
/// weights; and the refusals.
/// </summary>
public sealed class BassPredictTests
{
    private readonly ITestOutputHelper _out;
    public BassPredictTests(ITestOutputHelper output) => _out = output;

    private static readonly Dictionary<BassPointwiseTests.Kind, (BassData Data, BassFit Fit)> Fits = [];

    internal static (BassData Data, BassFit Fit) FitOf(BassPointwiseTests.Kind kind)
    {
        lock (Fits)
        {
            if (Fits.TryGetValue(kind, out var cached)) return cached;
            BassData data = BassPointwiseTests.Data(kind);
            BassFit fit = Bass.Fit(BassRunTests.Config(chains: 3, burn: 300, samples: 150, maxBasis: 6), data);
            Fits.Add(kind, (data, fit));
            return (data, fit);
        }
    }

    private static List<BassDraw> Pooled(BassFit fit) => [.. fit.Draws.SelectMany(c => c)];

    // A term's case side by hand: Πⱼ (sⱼ(x̃ⱼ − tⱼ))₊^q / c · Πₗ 1(zₗ ∈ Dₗ), c = Πⱼ((sⱼ + 1)/2 − sⱼtⱼ)^q or 1 where that is 0.
    private static double CaseSide(BassTerm term, BassInputs inputs, int degree, double[] scaled, int[] levels)
    {
        double value = 1.0, peak = 1.0;
        foreach (HingeFactor h in term.Continuous)
        {
            double t = inputs.ContinuousKnots(h.Variable)[h.Knot];
            value *= Math.Pow(Math.Max(h.Sign * (scaled[h.Variable] - t), 0.0), degree);
            peak *= (h.Sign + 1) / 2.0 - h.Sign * t;
        }
        double c = Math.Pow(peak, degree);
        if (c != 0.0) value /= c;
        foreach (CategoricalFactor f in term.Categorical)
            if (!f.Levels.Contains(levels[f.Variable])) value = 0.0;
        return value;
    }

    private static double OutputSide(BassTerm term, BassInputs inputs, int degree, double scaledT)
    {
        double value = 1.0, peak = 1.0;
        foreach (HingeFactor h in term.Functional)
        {
            double t = inputs.FunctionalKnots(h.Variable)[h.Knot];
            value *= Math.Pow(Math.Max(h.Sign * (scaledT - t), 0.0), degree);
            peak *= (h.Sign + 1) / 2.0 - h.Sign * t;
        }
        double c = Math.Pow(peak, degree);
        return c != 0.0 ? value / c : value;
    }

    [Theory]
    [InlineData(BassPointwiseTests.Kind.Scalar)]
    [InlineData(BassPointwiseTests.Kind.Categorical)]
    [InlineData(BassPointwiseTests.Kind.Functional)]
    public void AtTheTrainingInputs_ReproducesTheStageFourFits(BassPointwiseTests.Kind kind)
    {
        (BassData data, BassFit fit) = FitOf(kind);
        BassPrediction p = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Predictive = false });
        var drawFit = new BassDrawFit(BassModel.Create(fit.Config, data), fit.Terms);
        List<BassDraw> pooled = Pooled(fit);
        int n = data.Response.Cases, m = data.Response.Locations;
        Assert.Equal((pooled.Count, n, m), (p.DrawCount, p.Cases, p.Locations));
        Assert.Null(p.Predictive);
        Assert.All(p.CaseOutsideRange, Assert.False);
        Assert.All(p.LocationOutsideRange, Assert.False);
        var f = new double[n * m];
        for (int d = 0; d < pooled.Count; d++)
        {
            drawFit.Evaluate(pooled[d], f);
            for (int i = 0; i < n; i++)
                for (int l = 0; l < m; l++) Assert.Equal(f[i * m + l], p.Latent[d, i, l]);
        }
    }

    [Fact]
    public void NewCases_EvaluateTheHinges_AndFlagExtrapolation()
    {
        (BassData data, BassFit fit) = FitOf(BassPointwiseTests.Kind.Scalar);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        double[] x0 = [-0.2, 0.3, 0.55, 1.4], x1 = [0.5, 0.1, 1.0, 0.9];
        var cases = new PredictionCases([new ContinuousColumn("x1", x1), new ContinuousColumn("x0", x0)]);   // by name, any order
        BassPrediction p = fit.Predict(data, cases, new BassPredictOptions { Predictive = false });
        Assert.Equal([true, false, false, true], p.CaseOutsideRange);
        var s0 = (RangeScale)fit.Transform.Predictors[0];
        var s1 = (RangeScale)fit.Transform.Predictors[1];
        List<BassDraw> pooled = Pooled(fit);
        double worst = 0.0;
        for (int d = 0; d < pooled.Count; d++)
            for (int i = 0; i < 4; i++)
            {
                double[] scaled = [s0.Scale(x0[i]), s1.Scale(x1[i])];
                double f = pooled[d].Coefficients[0];
                int[] ids = pooled[d].Structure.TermIds;
                for (int j = 0; j < ids.Length; j++) f += pooled[d].Coefficients[j + 1] * CaseSide(fit.Terms[ids[j]], inputs, 1, scaled, []);
                worst = Math.Max(worst, Math.Abs(f - p.Latent[d, i, 0]));
            }
        _out.WriteLine($"largest gap to the hand evaluation: {worst:E1}");
        Assert.InRange(worst, 0, 1e-12);
    }

    [Fact]
    public void NewOutputLocations_ForADirectFunctionalFit()
    {
        (BassData data, BassFit fit) = FitOf(BassPointwiseTests.Kind.Functional);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        double[] t = [-0.1, 0.25, 0.7, 1.2];
        var cases = new PredictionCases([new ContinuousColumn("x0", [0.4, 1.1]), new ContinuousColumn("x1", [0.5, 0.2])]);
        BassPrediction p = fit.Predict(data, cases, new BassPredictOptions { Predictive = false, OutputLocations = [new ContinuousColumn("t", t)] });
        Assert.Equal((2, 4), (p.Cases, p.Locations));
        Assert.Equal([true, false, false, true], p.LocationOutsideRange);
        Assert.Equal([false, true], p.CaseOutsideRange);
        var s0 = (RangeScale)fit.Transform.Predictors[0];
        var s1 = (RangeScale)fit.Transform.Predictors[1];
        RangeScale st = fit.Transform.OutputVariables[0];
        List<BassDraw> pooled = Pooled(fit);
        double worst = 0.0;
        for (int d = 0; d < pooled.Count; d++)
            for (int i = 0; i < 2; i++)
                for (int l = 0; l < 4; l++)
                {
                    double[] scaled = [s0.Scale(i == 0 ? 0.4 : 1.1), s1.Scale(i == 0 ? 0.5 : 0.2)];
                    double f = pooled[d].Coefficients[0];
                    int[] ids = pooled[d].Structure.TermIds;
                    for (int j = 0; j < ids.Length; j++)
                    {
                        BassTerm term = fit.Terms[ids[j]];
                        f += pooled[d].Coefficients[j + 1] * CaseSide(term, inputs, 1, scaled, []) * OutputSide(term, inputs, 1, st.Scale(t[l]));
                    }
                    worst = Math.Max(worst, Math.Abs(f - p.Latent[d, i, l]));
                }
        Assert.InRange(worst, 0, 1e-12);
    }

    [Fact]
    public void DrawOrder_AndSingletonAxes_AreKept()
    {
        (BassData data, BassFit fit) = FitOf(BassPointwiseTests.Kind.Scalar);
        var three = new PredictionCases([new ContinuousColumn("x0", [0.2, 0.5, 0.7]), new ContinuousColumn("x1", [0.3, 0.3, 0.3])]);
        BassPrediction all = fit.Predict(data, three, new BassPredictOptions { Predictive = false });
        int[] order = [4, 1, 5, 2];
        BassPrediction p = fit.Predict(data, three, new BassPredictOptions { Draws = order, Predictive = false });
        Assert.Equal(order, p.Draws);
        Assert.Equal((4, 3, 1), (p.Latent.GetLength(0), p.Latent.GetLength(1), p.Latent.GetLength(2)));
        for (int d = 0; d < 4; d++)
            for (int i = 0; i < 3; i++) Assert.Equal(all.Latent[order[d], i, 0], p.Latent[d, i, 0]);

        var one = new PredictionCases([new ContinuousColumn("x0", [0.2]), new ContinuousColumn("x1", [0.3])]);
        BassPrediction single = fit.Predict(data, one, new BassPredictOptions { Draws = [4] });
        Assert.Equal((1, 1, 1), (single.Latent.GetLength(0), single.Latent.GetLength(1), single.Latent.GetLength(2)));
        Assert.Equal((1, 1, 1), (single.Predictive!.GetLength(0), single.Predictive.GetLength(1), single.Predictive.GetLength(2)));
        // Repeats are allowed and keep their places.
        BassPrediction repeated = fit.Predict(data, one, new BassPredictOptions { Draws = [2, 2, 0], Predictive = false });
        Assert.Equal(repeated.Latent[0, 0, 0], repeated.Latent[1, 0, 0]);

        (BassData fdata, BassFit ffit) = FitOf(BassPointwiseTests.Kind.Functional);
        var fone = new PredictionCases([new ContinuousColumn("x0", [0.4]), new ContinuousColumn("x1", [0.5])]);
        BassPrediction functional = ffit.Predict(fdata, fone, new BassPredictOptions { Draws = [3] });
        Assert.Equal((1, 1, 6), (functional.Latent.GetLength(0), functional.Latent.GetLength(1), functional.Latent.GetLength(2)));
        BassPrediction oneLocation = ffit.Predict(fdata, fone, new BassPredictOptions { Draws = [3], OutputLocations = [new ContinuousColumn("t", [0.5])] });
        Assert.Equal((1, 1, 1), (oneLocation.Latent.GetLength(0), oneLocation.Latent.GetLength(1), oneLocation.Latent.GetLength(2)));
    }

    [Fact]
    public void Noise_IsIndependentAcrossDrawsCasesAndLocations()
    {
        (BassData data, BassFit fit) = FitOf(BassPointwiseTests.Kind.Functional);
        BassPrediction p = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Seed = 21 });
        int s = p.DrawCount, n = p.Cases, m = p.Locations;
        var z = new double[s, n, m];
        var all = new List<double>();
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++)
                for (int l = 0; l < m; l++)
                {
                    z[d, i, l] = (p.Predictive![d, i, l] - p.Latent[d, i, l]) / Math.Sqrt(p.NoiseVariances[d]);
                    all.Add(z[d, i, l]);
                }
        double mean = all.Average(), variance = all.Average(v => (v - mean) * (v - mean));
        double Correlation(Func<int, double> a, Func<int, double> b)
        {
            double ma = Enumerable.Range(0, s).Average(a), mb = Enumerable.Range(0, s).Average(b);
            double ab = Enumerable.Range(0, s).Sum(d => (a(d) - ma) * (b(d) - mb));
            return ab / Math.Sqrt(Enumerable.Range(0, s).Sum(d => (a(d) - ma) * (a(d) - ma)) * Enumerable.Range(0, s).Sum(d => (b(d) - mb) * (b(d) - mb)));
        }
        double acrossCases = Correlation(d => z[d, 0, 2], d => z[d, 1, 2]);
        double acrossLocations = Correlation(d => z[d, 0, 2], d => z[d, 0, 3]);
        double bound = 4.0 / Math.Sqrt(s);
        _out.WriteLine($"standardized noise: mean {mean:F4}, variance {variance:F4} over {all.Count}; correlations {acrossCases:F3}, {acrossLocations:F3} (bound {bound:F3})");
        Assert.InRange(mean, -4.0 / Math.Sqrt(all.Count), 4.0 / Math.Sqrt(all.Count));
        Assert.InRange(variance, 1 - 6 * Math.Sqrt(2.0 / all.Count), 1 + 6 * Math.Sqrt(2.0 / all.Count));
        Assert.InRange(Math.Abs(acrossCases), 0, bound);
        Assert.InRange(Math.Abs(acrossLocations), 0, bound);
        // Within one draw the noise differs across cases (R recycles one value per draw).
        Assert.True(Enumerable.Range(0, n).Select(i => z[0, i, 0]).Distinct().Count() == n);
    }

    [Fact]
    public void Seed_ReproducesTheDraws()
    {
        (BassData data, BassFit fit) = FitOf(BassPointwiseTests.Kind.Scalar);
        PredictionCases cases = PredictionCases.Of(data);
        BassPrediction a = fit.Predict(data, cases, new BassPredictOptions { Seed = 5 });
        BassPrediction b = fit.Predict(data, cases, new BassPredictOptions { Seed = 5 });
        BassPrediction c = fit.Predict(data, cases, new BassPredictOptions { Seed = 6 });
        Assert.Equal(a.Predictive!.Cast<double>(), b.Predictive!.Cast<double>());
        Assert.NotEqual(a.Predictive.Cast<double>(), c.Predictive!.Cast<double>());
        Assert.Equal(a.Latent.Cast<double>(), c.Latent.Cast<double>());
    }

    [Fact]
    public void UnseenLevel_IsRefusedByDefault_OrPlacedOutsideEverySubset()
    {
        (BassData data, BassFit fit) = FitOf(BassPointwiseTests.Kind.Categorical);
        var cases = new PredictionCases([new ContinuousColumn("x0", [0.3, 0.6]), new CategoricalColumn("z", ["b", "d"])]);
        var e = Assert.Throws<ArgumentException>(() => fit.Predict(data, cases));
        Assert.Contains("'d'", e.Message);
        Assert.Contains("OutsideEverySubset", e.Message);

        BassPrediction p = fit.Predict(data, cases, new BassPredictOptions { Predictive = false, UnseenLevels = UnseenLevel.OutsideEverySubset });
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        var s0 = (RangeScale)fit.Transform.Predictors[0];
        List<BassDraw> pooled = Pooled(fit);
        bool sawCategorical = false;
        for (int d = 0; d < pooled.Count; d++)
            for (int i = 0; i < 2; i++)
            {
                int level = i == 0 ? 1 : -1;
                double f = pooled[d].Coefficients[0];
                int[] ids = pooled[d].Structure.TermIds;
                for (int j = 0; j < ids.Length; j++)
                {
                    sawCategorical |= fit.Terms[ids[j]].Categorical.Count > 0;
                    f += pooled[d].Coefficients[j + 1] * CaseSide(fit.Terms[ids[j]], inputs, 1, [s0.Scale(i == 0 ? 0.3 : 0.6)], [level]);
                }
                Assert.Equal(f, p.Latent[d, i, 0], 12);
            }
        Assert.True(sawCategorical);
    }

    [Fact]
    public void Quantile_IsThePredictiveMixtures()
    {
        (BassData data, BassFit fit) = FitOf(BassPointwiseTests.Kind.Scalar);
        var cases = new PredictionCases([new ContinuousColumn("x0", [0.1, 0.9]), new ContinuousColumn("x1", [0.2, 0.8])]);
        BassPrediction p = fit.Predict(data, cases, new BassPredictOptions { Predictive = false });
        foreach (double probability in new[] { 0.025, 0.5, 0.975 })
        {
            double[,] q = p.Quantile(probability);
            for (int i = 0; i < 2; i++)
            {
                double mass = Enumerable.Range(0, p.DrawCount).Sum(d =>
                    p.Weights[d] * Normal.MeanSd(p.Latent[d, i, 0], Math.Sqrt(p.NoiseVariances[d])).Cdf(q[i, 0]));
                Assert.Equal(probability, mass, 12);
            }
        }
        // One component: its own quantile. Two symmetric ones: the median between them.
        Assert.Equal(2 + 3 * Normal.StandardQuantile(0.9), PredictiveMixture.Quantile([2.0], [3.0], [1.0], 0.9), 13);
        Assert.Equal(0.0, PredictiveMixture.Quantile([-1.0, 1.0], [1.0, 1.0], [1.0, 1.0], 0.5), 15);
        Assert.Equal(-1.0 + Normal.StandardQuantile(0.3), PredictiveMixture.Quantile([-1.0, 5.0], [1.0, 1.0], [1.0, 0.0], 0.3), 13);
        // The empirical form: the least value whose cumulative weight reaches p.
        Assert.Equal(2.0, PredictiveMixture.Quantile([3.0, 1.0, 2.0], [1.0, 1.0, 1.0], 0.5));
        Assert.Equal(1.0, PredictiveMixture.Quantile([3.0, 1.0, 2.0], [0.1, 0.1, 0.8], 0.05));
        Assert.Equal(2.0, PredictiveMixture.Quantile([3.0, 1.0, 2.0], [0.1, 0.1, 0.8], 0.25));
        Assert.Equal(3.0, PredictiveMixture.Quantile([3.0, 1.0, 2.0], [0.1, 0.1, 0.8], 0.95));
    }

    [Fact]
    public void MeanSurface_IsTheGPriorsConditionalMean()
    {
        // Under the mean centre the conditional mean at the training inputs is ȳ + (ŷ_OLS − ȳ)/(1 + τ) (the block's prior
        // β_c | σ² ~ N(0, σ²(τG)⁻¹) on the centred columns, the intercept at ȳ).
        (BassData data, BassFit fit) = FitOf(BassPointwiseTests.Kind.Scalar);
        BassModel model = BassModel.Create(fit.Config, data);
        var drawFit = new BassDrawFit(model, fit.Terms);
        double[] y = ((ScalarResponse)data.Response).Y.ToArray();
        double yBar = y.Average();
        List<BassDraw> pooled = Pooled(fit);
        foreach (int d in new[] { 0, 7, pooled.Count - 1 })
        {
            BassDraw draw = pooled[d];
            double[,] x = drawFit.CentredDesign(draw.Structure.TermIds);
            double[] fitted = LeastSquares(x, y);
            BassPrediction p = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Draws = [d], Predictive = false, MeanSurface = true });
            for (int i = 0; i < y.Length; i++) Assert.Equal(yBar + (fitted[i] - yBar) / (1 + draw.Tau), p.MeanSurface![i, 0], 10);
        }
        // Over every draw it is a lower-variance estimate of the posterior mean the draws estimate.
        BassPrediction whole = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Predictive = false, MeanSurface = true });
        double[,] mean = whole.Mean();
        double gap = Enumerable.Range(0, y.Length).Max(i => Math.Abs(mean[i, 0] - whole.MeanSurface![i, 0]));
        _out.WriteLine($"largest gap between the draws' mean and the mean surface: {gap:F4}");
        Assert.InRange(gap, 0, 0.1);
    }

    // The least-squares fit X(XᵀX)⁻¹Xᵀy.
    private static double[] LeastSquares(double[,] x, double[] y)
    {
        int n = x.GetLength(0), p = x.GetLength(1);
        var gram = new double[p, p];
        var rhs = new double[p];
        for (int a = 0; a < p; a++)
        {
            for (int b = 0; b < p; b++)
                for (int i = 0; i < n; i++) gram[a, b] += x[i, a] * x[i, b];
            for (int i = 0; i < n; i++) rhs[a] += x[i, a] * y[i];
        }
        var factor = new DenseCholesky(p);
        factor.Decompose(gram);
        double[] beta = factor.Solve(rhs);
        var fitted = new double[n];
        for (int i = 0; i < n; i++)
            for (int a = 0; a < p; a++) fitted[i] += x[i, a] * beta[a];
        return fitted;
    }

    [Fact]
    public void Weights_AreTaken()
    {
        (BassData data, BassFit fit) = FitOf(BassPointwiseTests.Kind.Scalar);
        int s = fit.Draws.Sum(c => c.Count);
        var w = new double[s];
        w[3] = 3.0;
        w[10] = 1.0;
        BassFit weighted = fit with { DrawWeights = w };
        var one = new PredictionCases([new ContinuousColumn("x0", [0.5]), new ContinuousColumn("x1", [0.5])]);
        BassPrediction p = weighted.Predict(data, one, new BassPredictOptions { Predictive = false });
        Assert.Equal(0.75, p.Weights[3], 15);
        Assert.Equal(0.25, p.Weights[10], 15);
        Assert.Equal(0.75 * p.Latent[3, 0, 0] + 0.25 * p.Latent[10, 0, 0], p.Mean()[0, 0], 13);
        // All the weight on one draw: its own Gaussian's quantile.
        w = new double[s];
        w[3] = 1.0;
        BassPrediction single = (fit with { DrawWeights = w }).Predict(data, one, new BassPredictOptions { Predictive = false });
        Assert.Equal(single.Latent[3, 0, 0] + Math.Sqrt(single.NoiseVariances[3]) * Normal.StandardQuantile(0.9), single.Quantile(0.9)[0, 0], 12);
        Assert.Contains("no weight", Assert.Throws<ArgumentException>(() =>
            (fit with { DrawWeights = w }).Predict(data, one, new BassPredictOptions { Draws = [0, 1] })).Message);
    }

    [Fact]
    public void Refusals_SayWhy()
    {
        (BassData data, BassFit fit) = FitOf(BassPointwiseTests.Kind.Scalar);
        var good = new PredictionCases([new ContinuousColumn("x0", [0.5]), new ContinuousColumn("x1", [0.5])]);
        Assert.Contains("other data", Assert.Throws<ArgumentException>(() => fit.Predict(BassPointwiseTests.Data(BassPointwiseTests.Kind.Scalar, seed: 6), good)).Message);
        Assert.Contains("'x1'", Assert.Throws<ArgumentException>(() =>
            fit.Predict(data, new PredictionCases([new ContinuousColumn("x0", [0.5])]))).Message);
        Assert.Contains("names no predictor", Assert.Throws<ArgumentException>(() =>
            fit.Predict(data, new PredictionCases([.. good.Columns, new ContinuousColumn("x9", [0.5])]))).Message);
        Assert.Contains("kind", Assert.Throws<ArgumentException>(() =>
            fit.Predict(data, new PredictionCases([new ContinuousColumn("x0", [0.5]), new CategoricalColumn("x1", ["a"])]))).Message);
        Assert.Contains("no output locations", Assert.Throws<ArgumentException>(() =>
            fit.Predict(data, good, new BassPredictOptions { OutputLocations = [new ContinuousColumn("t", [0.5])] })).Message);
        Assert.Contains("not among", Assert.Throws<ArgumentException>(() => fit.Predict(data, good, new BassPredictOptions { Draws = [100000] })).Message);
        Assert.Contains("at least one draw", Assert.Throws<ArgumentException>(() => fit.Predict(data, good, new BassPredictOptions { Draws = [] })).Message);
        Assert.Contains("given twice", Assert.Throws<ArgumentException>(() =>
            new PredictionCases([new ContinuousColumn("x0", [0.5]), new ContinuousColumn("x0", [0.5])])).Message);
        Assert.Contains("non-finite", Assert.Throws<ArgumentException>(() => new ContinuousColumn("x0", [double.NaN])).Message);

        (BassData fdata, BassFit ffit) = FitOf(BassPointwiseTests.Kind.Functional);
        var fcases = new PredictionCases([new ContinuousColumn("x0", [0.5]), new ContinuousColumn("x1", [0.5])]);
        Assert.Contains("'t'", Assert.Throws<ArgumentException>(() =>
            ffit.Predict(fdata, fcases, new BassPredictOptions { OutputLocations = [new ContinuousColumn("u", [0.5])] })).Message);
    }
}
