using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// BASS's output path (stage 5, O1 and O3; D07) against the O-O1 defect fixtures (<c>bass_output_defects_oracle.R</c>),
/// which run the pinned R code as written beside D07's corrected expectations, on 8 curves at 40 locations with spikes
/// that skew the residuals. Each section asserts R's defect as the fixture records it and the corrected value on this
/// side: the truncation sampler over every training case and whole residual vectors; component noise per case; the
/// residual's orientation and units; 100% retention; constant columns; draw order and axes, and each row's noise from
/// its own draw; interval quantiles of the predictive mixture. Opt-in: it runs where the R toolchain and the BASS
/// checkout are present.
/// </summary>
public sealed class BassOutputDefectsParityTests
{
    private readonly ITestOutputHelper _out;
    public BassOutputDefectsParityTests(ITestOutputHelper output) => _out = output;

    private static bool Available => BassTargetParityTests.Available;

    // The stage-1 checks' recipe: sin(2π(t + 0.1i)) with a spike on every third curve, noise of sd 0.02 (this side's).
    private static double[,] Curves()
    {
        var rng = new Xoshiro256PlusPlus(27);
        var y = new double[8, 40];
        for (int i = 1; i <= 8; i++)
            for (int l = 0; l < 40; l++)
            {
                double t = l / 39.0;
                y[i - 1, l] = Math.Sin(2 * Math.PI * (t + 0.1 * i)) + (i % 3 == 0 ? 2 * Math.Exp(-Math.Pow((t - 0.2 - 0.07 * i) / 0.04, 2)) : 0.0)
                            + 0.02 * Normal.SampleStandard(rng);
            }
        return y;
    }

    private static readonly Lazy<JsonElement> Oracle = new(() => BassOutputTruncationParityTests.Run("oracles/bass_output_defects_oracle.R", Curves()));

    private static readonly Lazy<BassBasisFit> BasisFit = new(() =>
    {
        double[,] y = Curves();
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 0.9, Residual = TruncationResidual.InSample });
        var config = new BassConfig { Sampling = new BassSampling { Chains = 2, MasterSeed = 3, BurnIn = 200, Samples = 100 } };
        return BassComponents.Fit(d, [new ContinuousPredictor("x", [.. Enumerable.Range(1, 8).Select(i => i / 8.0)])], config);
    });

    private static int[] Ints(JsonElement e) => e.ValueKind == JsonValueKind.Array ? [.. e.EnumerateArray().Select(v => v.GetInt32())] : [e.GetInt32()];

    // Each predictive draw's component noise (standardized by its σ_k) and the training case its residual came from, at
    // `cases` identical prediction cases: under in-sample residuals Kᵀe = 0, so (y − latent)/y_s = Kε + e splits.
    private static (double[,,] Epsilon, int[,] Drawn) Recover(BassBasisFit fit, int cases, int seed)
    {
        BassDecomposition d = fit.Decomposition;
        int r = d.Rank, m = d.Locations, n = d.Cases;
        BassBasisPrediction p = fit.Predict(new PredictionCases([new ContinuousColumn("x", [.. Enumerable.Repeat(0.5, cases)])]),
                                            new BassPredictOptions { Seed = seed });
        var epsilon = new double[p.DrawCount, cases, r];
        var drawn = new int[p.DrawCount, cases];
        for (int j = 0; j < p.DrawCount; j++)
            for (int i = 0; i < cases; i++)
            {
                double[] u = [.. Enumerable.Range(0, m).Select(l => (p.Predictive![j, i, l] - p.Latent[j, i, l]) / d.Map.Scales[l])];
                var kEpsilon = new double[m];
                for (int k = 0; k < r; k++)
                {
                    double kk = Enumerable.Range(0, m).Sum(l => d.Map.Basis(l, k) * d.Map.Basis(l, k));
                    double raw = Enumerable.Range(0, m).Sum(l => d.Map.Basis(l, k) * u[l]) / kk;
                    epsilon[j, i, k] = raw / Math.Sqrt(p.NoiseVariances[j, k]);
                    for (int l = 0; l < m; l++) kEpsilon[l] += d.Map.Basis(l, k) * raw;
                }
                (int c, double gap) = Enumerable.Range(0, n).Select(c => (c, Enumerable.Range(0, m).Max(l => Math.Abs(u[l] - kEpsilon[l] - d.Truncation.Values[c, l]))))
                                                .MinBy(t => t.Item2);
                Assert.InRange(gap, 0, 1e-10);   // a whole training residual vector
                drawn[j, i] = c;
            }
        return (epsilon, drawn);
    }

    [Fact]
    public void TruncationSampler_DrawsWholeResidualsOverEveryCase()
    {
        if (!Available) return;
        JsonElement ts = Oracle.Value.GetProperty("truncation_sampler");
        Assert.Equal(JsonValueKind.String, ts.GetProperty("wide").GetProperty("r_error").ValueKind);   // R stops at m > n
        Assert.Equal([1, 2, 3], Ints(ts.GetProperty("narrow").GetProperty("r_cases_drawn")));         // R draws 3 of 8
        Assert.Equal(8, Ints(ts.GetProperty("narrow").GetProperty("corrected_cases_drawn")).Length);
        (_, int[,] drawn) = Recover(BasisFit.Value, cases: 10, seed: 4);
        int[] seen = [.. drawn.Cast<int>().Distinct().Order()];
        _out.WriteLine($"R stops: \"{ts.GetProperty("wide").GetProperty("r_error").GetString()}\"; this side drew cases {string.Join(",", seen)}");
        Assert.Equal(Enumerable.Range(0, 8), seen);
    }

    [Fact]
    public void ComponentNoise_IsDrawnPerCase()
    {
        if (!Available) return;
        JsonElement nr = Oracle.Value.GetProperty("noise_recycling");
        Assert.Equal(0.0, nr.GetProperty("r_max_spread_across_cases").GetDouble());
        Assert.Equal(1.0, nr.GetProperty("r_correlation_cases_1_2").GetDouble(), 12);
        (double[,,] epsilon, _) = Recover(BasisFit.Value, cases: 5, seed: 6);
        int s = epsilon.GetLength(0);
        double spread = Enumerable.Range(0, s).Max(j => Enumerable.Range(0, 5).Max(i => epsilon[j, i, 0]) - Enumerable.Range(0, 5).Min(i => epsilon[j, i, 0]));
        double correlation = Enumerable.Range(0, s).Sum(j => epsilon[j, 0, 0] * epsilon[j, 1, 0])
                           / Math.Sqrt(Enumerable.Range(0, s).Sum(j => epsilon[j, 0, 0] * epsilon[j, 0, 0]) * Enumerable.Range(0, s).Sum(j => epsilon[j, 1, 0] * epsilon[j, 1, 0]));
        _out.WriteLine($"R: spread 0, correlation 1; this side: spread {spread:F2}, correlation {correlation:F3} over {s} draws");
        Assert.True(spread > 1.0);
        Assert.InRange(Math.Abs(correlation), 0, 4 / Math.Sqrt(s));
    }

    [Fact]
    public void Residual_IsObservedMinusReconstructed_AddedInStandardizedUnits()
    {
        if (!Available) return;
        JsonElement rs = Oracle.Value.GetProperty("residual");
        double[,] y = Curves();
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 0.9, ScaleByColumnSd = true, Residual = TruncationResidual.InSample });
        Assert.Equal(rs.GetProperty("r").GetInt32(), d.Rank);
        double[] range = [.. rs.GetProperty("y_s_range").EnumerateArray().Select(v => v.GetDouble())];
        Assert.Equal(range[0], Enumerable.Range(0, 40).Min(l => d.Map.Scales[l]), 12);
        Assert.Equal(range[1], Enumerable.Range(0, 40).Max(l => d.Map.Scales[l]), 12);
        // R stores −e and adds it after restoring physical units: its reconstruction misses Y.
        Assert.Equal(0.0, rs.GetProperty("r_stored_is_negative_of_corrected").GetDouble());
        Assert.True(rs.GetProperty("r_reconstruction_max_error").GetDouble() > 0.1);
        double worst = 0.0;
        var w = new double[d.Rank];
        var e = new double[40];
        var back = new double[40];
        for (int i = 0; i < 8; i++)
        {
            for (int k = 0; k < d.Rank; k++) w[k] = d.Coefficient(k, i);
            for (int l = 0; l < 40; l++) e[l] = d.Truncation.Values[i, l];
            d.Map.Reconstruct(w, e, back);
            for (int l = 0; l < 40; l++) worst = Math.Max(worst, Math.Abs(back[l] - y[i, l]));
        }
        double[] all = d.Truncation.Values.Cast<double>().ToArray();
        double mean = all.Average();
        double skew = all.Average(v => Math.Pow(v - mean, 3)) / Math.Pow(all.Average(v => (v - mean) * (v - mean)), 1.5);
        _out.WriteLine($"reconstruction: R {rs.GetProperty("r_reconstruction_max_error").GetDouble():F2}, this side {worst:E1}; " +
                       $"skewness: R {rs.GetProperty("skewness_r").GetDouble():F3}, corrected {rs.GetProperty("skewness_corrected").GetDouble():F3}, this side {skew:F3}");
        Assert.InRange(worst, 0, 1e-12 * rs.GetProperty("scale").GetDouble());
        Assert.Equal(rs.GetProperty("skewness_corrected").GetDouble(), skew, 9);
        Assert.Equal(-rs.GetProperty("skewness_r").GetDouble(), skew, 9);
        Assert.True(skew > 0);
    }

    [Fact]
    public void FullRetention_IsTheNumericalRank()
    {
        if (!Available) return;
        JsonElement fr = Oracle.Value.GetProperty("full_retention");
        Assert.Equal(JsonValueKind.String, fr.GetProperty("r_error").ValueKind);
        Assert.Equal(JsonValueKind.Null, fr.GetProperty("r_n_pc").ValueKind);
        BassDecomposition d = BassDecomposition.FromPca(Curves(), new BassPcaSpec { VarianceShare = 1.0 });
        Assert.Equal(fr.GetProperty("corrected_rank").GetInt32(), d.Rank);
        Assert.Equal(7, d.Rank);
        Assert.InRange(d.Truncation.InSample!.Cast<double>().Max(Math.Abs), 0, d.Pca!.RankTolerance);
    }

    [Fact]
    public void ConstantColumns_AreRefusedAsInputs_AndScaledByOneAsOutputs()
    {
        if (!Available) return;
        JsonElement cc = Oracle.Value.GetProperty("constant_columns");
        JsonElement input = cc.GetProperty("input");
        Assert.Equal(0.0, input.GetProperty("r_training").GetDouble());
        Assert.Equal(["NaN", "Inf"], input.GetProperty("r_prediction").EnumerateArray().Select(v => v.GetString()));
        // The fixture's corrected mapping to 0 predates D07's amendment (2026-10-29): a constant training column is refused.
        Assert.Contains("constant", Assert.Throws<ArgumentException>(() => new ContinuousPredictor("x", [.. Enumerable.Repeat(0.3, 8)])).Message);

        JsonElement output = cc.GetProperty("output");
        Assert.Equal(JsonValueKind.String, output.GetProperty("r_error").ValueKind);
        double[,] y = Curves();
        for (int i = 0; i < 8; i++) y[i, 1] = 7.0;
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 0.9, ScaleByColumnSd = true });
        Assert.Equal(output.GetProperty("corrected_y_s").GetDouble(), d.Map.Scales[1]);
        Assert.Equal([1], d.Pca!.ConstantColumns);
    }

    [Fact]
    public void DrawOrder_AxesAndNoise_FollowTheRequest()
    {
        if (!Available) return;
        JsonElement order = Oracle.Value.GetProperty("draw_order");
        int[] use = Ints(order.GetProperty("requested"));
        Assert.Equal([4, 1, 2, 5], Ints(order.GetProperty("r_rows")));   // R reorders by model group
        Assert.NotEqual(use, Ints(order.GetProperty("r_rows")));

        // A mock fit whose draw d (1-based) is the constant d, with σ² = (0, 1, 0, 0, 0) as R's mock (0 as 10⁻³⁰⁰, below
        // what the draw's level resolves).
        var x = new ContinuousPredictor("x", [0.0, 0.25, 0.5, 0.75, 1.0]);
        var data = new BassData([x], new ScalarResponse([0.1, 0.5, 0.2, 0.9, 0.4]));
        BassFit fit = Bass.Fit(new BassConfig { Sampling = new BassSampling { Chains = 1, BurnIn = 0, Samples = 1 } }, data);
        double[] s2 = [1e-300, 1.0, 1e-300, 1e-300, 1e-300];
        var store = new DrawStore<BassDraw>();
        for (int d = 1; d <= 5; d++) store.Add(new BassDraw(new BassStructure([]), 1.0, s2[d - 1], [d], 1.0));
        BassFit mock = fit with { Draws = [store] };
        int[] request = [.. use.Select(u => u - 1)];
        PredictionCases Cases(params double[] values) => new([new ContinuousColumn("x", values)]);

        BassPrediction p43 = mock.Predict(data, Cases(0.2, 0.5, 0.7), new BassPredictOptions { Draws = request, Seed = 13 });
        Assert.Equal(Ints(order.GetProperty("corrected_rows")), Enumerable.Range(0, 4).Select(d => (int)p43.Latent[d, 0, 0]));
        int[] noisy = [.. Enumerable.Range(0, 4).Where(d => Enumerable.Range(0, 3).Any(i => p43.Predictive![d, i, 0] != p43.Latent[d, i, 0])).Select(d => use[d])];
        Assert.Equal(Ints(order.GetProperty("corrected_rows_with_noise")), noisy);

        JsonElement[] dims = [.. order.GetProperty("corrected_dims").EnumerateArray()];
        BassPrediction p41 = mock.Predict(data, Cases(0.2), new BassPredictOptions { Draws = request });
        BassPrediction p11 = mock.Predict(data, Cases(0.2), new BassPredictOptions { Draws = [3] });
        foreach ((BassPrediction p, JsonElement dim) in new[] { (p43, dims[0]), (p41, dims[1]), (p11, dims[2]) })
            Assert.Equal(Ints(dim), new[] { p.Predictive!.GetLength(0), p.Predictive.GetLength(1), p.Predictive.GetLength(2) });
        _out.WriteLine($"R: rows {string.Join(",", Ints(order.GetProperty("r_rows")))}, dims 4x1 {string.Join("x", Ints(order.GetProperty("r_dim_4x1")))}, " +
                       $"1x1 {string.Join("x", Ints(order.GetProperty("r_dim_1x1")))}, noise on {string.Join(",", Ints(order.GetProperty("r_rows_with_noise")))}");
    }

    [Fact]
    public void Intervals_AreThePredictiveMixturesQuantiles()
    {
        if (!Available) return;
        JsonElement[] cases = [.. Oracle.Value.GetProperty("interval_quantiles").EnumerateArray()];
        const int draws = 1000;
        double[] u = [.. Enumerable.Range(1, draws).Select(k => Normal.StandardQuantile((k - 0.5) / draws))];
        double[] equal = [.. Enumerable.Repeat(1.0, draws)];
        (double[] Mu, double[] Sd)[] spreads = [(u, [.. Enumerable.Repeat(1.0, draws)]), (new double[draws], [.. u.Select(v => Math.Exp(v / 2))])];
        double[] probabilities = [0.025, 0.975];
        for (int c = 0; c < 2; c++)
        {
            double[] mixture = [.. cases[c].GetProperty("mixture").EnumerateArray().Select(v => v.GetDouble())];
            double[] r = [.. cases[c].GetProperty("r").EnumerateArray().Select(v => v.GetDouble())];
            double[] ours = [.. probabilities.Select(p => PredictiveMixture.Quantile(spreads[c].Mu, spreads[c].Sd, equal, p))];
            _out.WriteLine($"case {c + 1}: R's interval ({r[0]:F3}, {r[1]:F3}), the mixture's ({mixture[0]:F3}, {mixture[1]:F3}), this side ({ours[0]:F3}, {ours[1]:F3})");
            for (int q = 0; q < 2; q++)
            {
                Assert.Equal(mixture[q], ours[q], 1e-10);
                Assert.True(Math.Abs(r[q] - mixture[q]) > 0.5);
            }
        }
    }
}
