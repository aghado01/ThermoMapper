using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 8b gates on the comparison (<see cref="BassProjectionComparison"/>):
/// <list type="bullet">
/// <item>C, the held-out log predictive density: against a direct computation from each draw's prediction and σ² under Gaussian noise,
/// on both sides; under a Student t noise against the seam's own values; weighted draws against the draws replicated by their
/// weights; and the cross-fitted folds' densities summing to the report's;</item>
/// <item>D, cross-fitting with groupings and case weights: integer replication weights scoring as duplicated rows in the RMSE, the sd
/// and the density; precision weights reading σ²/wᵢ in the density and each case once in the RMSE; group folds disjoint by level and
/// covering every case, each fold's new groups read through the integrated density under a random intercept, and case folds spreading
/// a grouping's levels; each fold's frame read from the other folds alone under group folds; and the refusals.</item>
/// </list>
/// </summary>
public sealed class ProjectionComparisonHeldOutTests(ITestOutputHelper output)
{
    private static BassConfig Light(int seed) => BassProjectedFitTests.Config() with
    {
        Sampling = new BassSampling { Chains = 1, MasterSeed = seed, BurnIn = 600, Samples = 30, Thin = 5 },
    };

    private sealed record Pipeline(BassData Data, BassFit Native, ActiveSubspaceAnalysis Analysis, BassProjectedFit Refit);

    private static Pipeline Build(BassData data, BassConfig config)
    {
        BassFit native = Bass.Fit(config, data);
        ActiveSubspaceAnalysis analysis = native.ActiveSubspace(data);
        BassProjection projection = BassProjection.Of(native.Transform, analysis, 1);
        return new Pipeline(data, native, analysis, BassProjectedFit.Of(projection, data, config));
    }

    // The ridge's held-out cases and responses, n of them from the seed.
    private static (PredictionCases Cases, double[,] Response, double[][] X) HeldOut(int n, int seed)
    {
        (_, double[][] x) = BassProjectedFitTests.RidgeData(n, seed);
        double[] w = BassProjectedFitTests.W;
        var response = new double[n, 1];
        for (int r = 0; r < n; r++) response[r, 0] = BassProjectedFitTests.Ridge(Enumerable.Range(0, w.Length).Sum(i => w[i] * x[i][r]));
        return (BassProjectedFitTests.Cases(x), response, x);
    }

    // Σ_r ln Σ_d w_d N(y_r; f_d(x_r), σ_d²/c_r) from the prediction's draws, c_r a precision or 1.
    private static double Direct(BassPrediction prediction, double[,] response, double[]? precision = null)
    {
        double total = 0.0;
        for (int r = 0; r < response.GetLength(0); r++)
        {
            double top = double.NegativeInfinity;
            var terms = new double[prediction.Draws.Count];
            for (int d = 0; d < terms.Length; d++)
            {
                double variance = prediction.NoiseVariances[d] / (precision?[r] ?? 1.0), e = response[r, 0] - prediction.Latent[d, r, 0];
                terms[d] = Math.Log(prediction.Weights[d]) - 0.5 * Math.Log(2.0 * Math.PI * variance) - 0.5 * e * e / variance;
                top = Math.Max(top, terms[d]);
            }
            total += top + Math.Log(terms.Sum(t => Math.Exp(t - top)));
        }
        return total;
    }

    /// <summary>
    /// On slice 7b's ridge with 300 held-out cases, Gaussian noise: each side's held-out density against Σ_r ln Σ_d w_d N(y_r; f_d, σ_d²)
    /// from its prediction's draws, to 10⁻¹² relative. Under weights 1, 2 and 3 on the native draws in turn, the density equals the
    /// one from the draws replicated by those weights.
    /// </summary>
    [Fact]
    public void HeldOutDensity_AgainstPerDrawPredictions_AndWeightedDraws()
    {
        BassProjectedFitTests.RidgeFits f = BassProjectedFitTests.Fixture.Value;
        (PredictionCases cases, double[,] response, _) = HeldOut(300, 41);
        var spec = new ProjectionComparisonSpec { TestCases = cases, TestResponse = response };
        BassProjectionComparison c = f.Refit.Compare(f.Native, f.Data, f.Analysis, spec);
        double native = Direct(f.Native.Predict(f.Data, cases), response), refit = Direct(f.Refit.Predict(cases), response);
        double gapNative = Math.Abs(c.Native.HeldOutLogDensity!.Value - native) / Math.Abs(native);
        double gapRefit = Math.Abs(c.Refit.HeldOutLogDensity!.Value - refit) / Math.Abs(refit);

        BassDraw[] pooled = [.. f.Native.Draws.SelectMany(d => d)];
        double[] weights = [.. Enumerable.Range(0, pooled.Length).Select(s => 1.0 + s % 3)];
        var replicated = new DrawStore<BassDraw>();
        for (int s = 0; s < pooled.Length; s++)
            for (int copy = 0; copy < (int)weights[s]; copy++) replicated.Add(pooled[s]);
        double weighted = f.Refit.Compare(f.Native with { DrawWeights = weights }, f.Data, f.Analysis, spec).Native.HeldOutLogDensity!.Value;
        double copies = f.Refit.Compare(f.Native with { Draws = [replicated], DrawWeights = null }, f.Data, f.Analysis, spec).Native.HeldOutLogDensity!.Value;
        double gapWeights = Math.Abs(weighted - copies) / Math.Abs(copies);
        output.WriteLine($"held-out density: native {c.Native.HeldOutLogDensity:F3} (direct {native:F3}, {gapNative:E2}), refit {c.Refit.HeldOutLogDensity:F3} " +
                         $"(direct {refit:F3}, {gapRefit:E2}); weighted draws against replicated {gapWeights:E2}");
        Assert.True(gapNative <= 1e-12, $"{gapNative:E2}");
        Assert.True(gapRefit <= 1e-12, $"{gapRefit:E2}");
        Assert.True(gapWeights <= 1e-12, $"{gapWeights:E2}");
    }

    /// <summary>Under a Student t noise (ν = 5) the native side's density is the seam's own values, each case's one candidate, averaged
    /// over the draws in logs, to rounding.</summary>
    [Fact]
    public void HeldOutDensity_UnderANoiseLaw_IsTheSeams()
    {
        BassProjectedFitTests.RidgeFits f = BassProjectedFitTests.Fixture.Value;
        BassFit noisy = Bass.Fit(Light(43) with { Noise = BassNoise.StudentT(5.0) }, f.Data);
        (PredictionCases cases, double[,] response, _) = HeldOut(120, 47);
        double got = f.Refit.Compare(noisy, f.Data, f.Analysis, new ProjectionComparisonSpec { TestCases = cases, TestResponse = response }).Native.HeldOutLogDensity!.Value;
        CandidateLogDensity seam = noisy.CandidateLogDensity(f.Data, cases, [.. Enumerable.Range(0, 120).Select(r => new[] { response[r, 0] })], quadrature: true);
        double expected = 0.0;
        for (int t = 0; t < 120; t++)
        {
            double[] v = [.. Enumerable.Range(0, seam.Draws).Select(s => seam.Values[s, seam.Offset(t)])];
            double top = v.Max();
            expected += top + Math.Log(v.Sum(x => Math.Exp(x - top)) / seam.Draws);
        }
        double gap = Math.Abs(got - expected) / Math.Abs(expected);
        output.WriteLine($"Student t (ν = 5): {got:F4} against the seam's {expected:F4}, {gap:E2}");
        Assert.True(gap <= 1e-12, $"{gap:E2}");
    }

    /// <summary>
    /// Replication weights 1 to 3 on the ridge's 120 training cases: held out, 60 cases with weights 1 to 3 score as the same cases
    /// duplicated by their weights at weight 1, the RMSE, the percentage and the density alike, on both sides. Precision weights 0.5 to
    /// 2: the density is Σ_r ln Σ_d w_d N(y_r; f_d, σ_d²/w_r), and the RMSE counts each case once.
    /// </summary>
    [Fact]
    public void CaseWeights_ScoreAtTheirMeaning()
    {
        (BassData plain, _) = BassProjectedFitTests.RidgeData(120, 53);
        double[] replication = [.. Enumerable.Range(0, 120).Select(r => 1.0 + r % 3)];
        Pipeline p = Build(new BassData(plain.Predictors, plain.Response, null, new CaseWeights(replication, WeightMeaning.Replication)), Light(59));
        (PredictionCases cases, double[,] response, double[][] x) = HeldOut(60, 61);
        double[] w = [.. Enumerable.Range(0, 60).Select(r => 1.0 + r % 3)];
        int[] rows = [.. Enumerable.Range(0, 60).SelectMany(r => Enumerable.Repeat(r, (int)w[r]))];
        var duplicatedResponse = new double[rows.Length, 1];
        for (int i = 0; i < rows.Length; i++) duplicatedResponse[i, 0] = response[rows[i], 0];
        PredictionCases duplicated = BassProjectedFitTests.Cases([.. x.Select(column => rows.Select(r => column[r]).ToArray())]);
        BassProjectionComparison once = p.Refit.Compare(p.Native, p.Data, p.Analysis, new ProjectionComparisonSpec { TestCases = cases, TestResponse = response, TestWeights = w });
        BassProjectionComparison copies = p.Refit.Compare(p.Native, p.Data, p.Analysis, new ProjectionComparisonSpec
        {
            TestCases = duplicated, TestResponse = duplicatedResponse, TestWeights = [.. rows.Select(_ => 1.0)],
        });
        double worst = 0.0;
        foreach ((ProjectionFitReport a, ProjectionFitReport b) in new[] { (once.Native, copies.Native), (once.Refit, copies.Refit) })
        {
            worst = Math.Max(worst, Math.Abs(a.HeldOutRmse!.Value - b.HeldOutRmse!.Value) / b.HeldOutRmse.Value);
            worst = Math.Max(worst, Math.Abs(a.HeldOutExplained!.Value - b.HeldOutExplained!.Value) / Math.Abs(b.HeldOutExplained.Value));
            worst = Math.Max(worst, Math.Abs(a.HeldOutLogDensity!.Value - b.HeldOutLogDensity!.Value) / Math.Abs(b.HeldOutLogDensity.Value));
        }

        double[] precision = [.. Enumerable.Range(0, 120).Select(r => 0.5 + 1.5 * (r % 4) / 3.0)];
        Pipeline q = Build(new BassData(plain.Predictors, plain.Response, null, new CaseWeights(precision, WeightMeaning.Precision)), Light(67));
        double[] wp = [.. Enumerable.Range(0, 60).Select(r => 0.5 + 1.5 * (r % 4) / 3.0)];
        BassProjectionComparison prec = q.Refit.Compare(q.Native, q.Data, q.Analysis, new ProjectionComparisonSpec { TestCases = cases, TestResponse = response, TestWeights = wp });
        double direct = Direct(q.Native.Predict(q.Data, cases), response, wp);
        double gapPrecision = Math.Abs(prec.Native.HeldOutLogDensity!.Value - direct) / Math.Abs(direct);
        double[,] mean = q.Native.Predict(q.Data, cases).Mean();
        double rmse = Math.Sqrt(Enumerable.Range(0, 60).Average(r => Math.Pow(mean[r, 0] - response[r, 0], 2)));
        double gapRmse = Math.Abs(prec.Native.HeldOutRmse!.Value - rmse) / rmse;
        output.WriteLine($"replication: weighted against duplicated, worst {worst:E2} (native density {once.Native.HeldOutLogDensity:F3}); " +
                         $"precision: density against σ²/w {gapPrecision:E2}, RMSE once per case {gapRmse:E2}");
        Assert.True(worst <= 1e-12, $"{worst:E2}");
        Assert.True(gapPrecision <= 1e-12, $"{gapPrecision:E2}");
        Assert.True(gapRmse <= 1e-12, $"{gapRmse:E2}");
    }

    // The ridge's n cases with a grouping g of `levels` levels, case r in level r mod levels; under `effects` y gains a level's offset.
    private static BassData Grouped(int n, int levels, int seed, bool effects)
    {
        (BassData plain, _) = BassProjectedFitTests.RidgeData(n, seed);
        int[] codes = [.. Enumerable.Range(0, n).Select(r => r % levels)];
        double[] y = ((ScalarResponse)plain.Response).Y.ToArray();
        if (effects)
            for (int r = 0; r < n; r++) y[r] += 0.2 * Math.Sin(1.7 * codes[r]);
        var g = new CategoricalPredictor("g", codes, [.. Enumerable.Range(0, levels).Select(l => $"L{l}")]);
        return new BassData(plain.Predictors, new ScalarResponse(y), [g]);
    }

    /// <summary>
    /// A random intercept on a grouping of 8 levels over 120 cases. Folds by the grouping, 4 of them: disjoint, covering every case,
    /// each level's cases in one fold, two levels a fold; every held-out case is in a new group, read through the integrated density,
    /// finite; the folds' densities sum to the report's totals. Folds by cases on the same data spread the levels across folds.
    /// </summary>
    [Fact]
    public void GroupFolds_ByLevel_NewGroupsThroughTheIntegratedDensity()
    {
        BassData data = Grouped(120, 8, 71, effects: true);
        BassConfig config = Light(73) with { RandomIntercepts = [new BassRandomIntercept("g", 0.8, ScalePrior: new BassGroupScalePrior())] };
        Pipeline p = Build(data, config);
        CrossFitReport byLevel = p.Refit.Compare(p.Native, p.Data, p.Analysis, new ProjectionComparisonSpec { Folds = 4, FoldSeed = 3, FoldGrouping = "g" }).CrossFitted!;
        int[] codes = ((CategoricalPredictor)data.Groupings[0]).Codes.ToArray();
        Assert.Equal(Enumerable.Range(0, 120), byLevel.Folds.SelectMany(fold => fold).Order());
        int[][] levelsOf = [.. byLevel.Folds.Select(fold => fold.Select(r => codes[r]).Distinct().Order().ToArray())];
        Assert.Equal(8, levelsOf.SelectMany(l => l).Distinct().Count());
        Assert.Equal(8, levelsOf.Sum(l => l.Length));
        Assert.All(levelsOf, l => Assert.Equal(2, l.Length));
        Assert.Equal("g", byLevel.FoldGrouping);
        Assert.All(byLevel.NativeLogDensity!.Concat(byLevel.RefitLogDensity!), d => Assert.True(double.IsFinite(d)));
        Assert.Equal(byLevel.NativeLogDensity!.Sum(), byLevel.NativeLogDensityTotal);
        Assert.Equal(byLevel.RefitLogDensity!.Sum(), byLevel.RefitLogDensityTotal);

        CrossFitReport byCase = p.Refit.Compare(p.Native, p.Data, p.Analysis, new ProjectionComparisonSpec { Folds = 4, FoldSeed = 3 }).CrossFitted!;
        int spread = Enumerable.Range(0, 8).Count(l => byCase.Folds.Count(fold => fold.Any(r => codes[r] == l)) > 1);
        Assert.Null(byCase.FoldGrouping);
        output.WriteLine($"by level: folds of {string.Join(", ", byLevel.Folds.Select(fold => fold.Length))}, densities native {byLevel.NativeLogDensityTotal:F2}, refit " +
                         $"{byLevel.RefitLogDensityTotal:F2}, RMSE native {byLevel.NativePooled:F4}, refit {byLevel.RefitPooled:F4}; by case: {spread} of 8 levels across folds, " +
                         $"densities native {byCase.NativeLogDensityTotal:F2}, refit {byCase.RefitLogDensityTotal:F2}");
        Assert.True(spread >= 7, $"{spread}");
    }

    /// <summary>Folds by a grouping's levels, under a random intercept (a grouping is read only by random effects): with fold 0's
    /// responses raised by 50, fold 0's frame is the same bits, read from the other folds alone, and the others' move, as slice 7c's case
    /// folds.</summary>
    [Fact]
    public void GroupFolds_EachFoldsFrameFromTheOthers()
    {
        BassData data = Grouped(90, 6, 79, effects: true);
        BassConfig config = Light(83) with { RandomIntercepts = [new BassRandomIntercept("g", 0.8, ScalePrior: new BassGroupScalePrior())] };
        Pipeline p = Build(data, config);
        var spec = new ProjectionComparisonSpec { Folds = 3, FoldSeed = 5, FoldGrouping = "g" };
        CrossFitReport cross = p.Refit.Compare(p.Native, p.Data, p.Analysis, spec).CrossFitted!;
        double[] y = ((ScalarResponse)data.Response).Y.ToArray();
        foreach (int r in cross.Folds[0]) y[r] += 50.0;
        Pipeline raised = Build(new BassData(data.Predictors, new ScalarResponse(y), data.Groupings), config);
        CrossFitReport moved = raised.Refit.Compare(raised.Native, raised.Data, raised.Analysis, spec).CrossFitted!;
        Assert.Equal(cross.Folds.Select(fold => fold.ToArray()), moved.Folds.Select(fold => fold.ToArray()));
        Assert.Equal(Flatten(cross.Frames[0].Frame), Flatten(moved.Frames[0].Frame));
        Assert.NotEqual(Flatten(cross.Frames[1].Frame), Flatten(moved.Frames[1].Frame));
        Assert.NotEqual(Flatten(cross.Frames[2].Frame), Flatten(moved.Frames[2].Frame));
    }

    /// <summary>Refused: held-out cases without weights on weighted data, weights on unweighted data or without cases, a weight count
    /// off, a fold grouping the data do not hold, more folds than its observed levels, and a fold grouping without folds.</summary>
    [Fact]
    public void Refusals()
    {
        BassProjectedFitTests.RidgeFits f = BassProjectedFitTests.Fixture.Value;
        (PredictionCases cases, double[,] response, _) = HeldOut(10, 89);
        Assert.Throws<ArgumentException>(() => f.Refit.Compare(f.Native, f.Data, f.Analysis, new ProjectionComparisonSpec { TestCases = cases, TestResponse = response, TestWeights = new double[10] }));
        Assert.Throws<ArgumentException>(() => f.Refit.Compare(f.Native, f.Data, f.Analysis, new ProjectionComparisonSpec { TestWeights = [1.0] }));
        Assert.Throws<ArgumentException>(() => f.Refit.Compare(f.Native, f.Data, f.Analysis, new ProjectionComparisonSpec { Folds = 2, FoldGrouping = "g" }));
        Assert.Throws<ArgumentException>(() => f.Refit.Compare(f.Native, f.Data, f.Analysis, new ProjectionComparisonSpec { FoldGrouping = "g" }));

        (BassData plain, _) = BassProjectedFitTests.RidgeData(60, 97);
        Pipeline p = Build(new BassData(plain.Predictors, plain.Response, null,
                                        new CaseWeights([.. Enumerable.Range(0, 60).Select(r => 1.0 + r % 2)], WeightMeaning.Replication)), Light(101));
        Assert.Throws<ArgumentException>(() => p.Refit.Compare(p.Native, p.Data, p.Analysis, new ProjectionComparisonSpec { TestCases = cases, TestResponse = response }));
        Assert.Throws<ArgumentException>(() => p.Refit.Compare(p.Native, p.Data, p.Analysis, new ProjectionComparisonSpec { TestCases = cases, TestResponse = response, TestWeights = [1.0, 2.0] }));
        Pipeline g = Build(Grouped(60, 3, 103, effects: true), Light(107) with { RandomIntercepts = [new BassRandomIntercept("g", 0.8, ScalePrior: new BassGroupScalePrior())] });
        Assert.Throws<ArgumentOutOfRangeException>(() => g.Refit.Compare(g.Native, g.Data, g.Analysis, new ProjectionComparisonSpec { Folds = 4, FoldGrouping = "g" }));
    }

    private static double[] Flatten(double[,] m) => [.. Enumerable.Range(0, m.GetLength(0)).SelectMany(i => Enumerable.Range(0, m.GetLength(1)).Select(j => m[i, j]))];
}
