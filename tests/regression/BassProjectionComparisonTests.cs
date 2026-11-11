using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 7c gates (the comparison, <see cref="BassProjectionComparison"/>), on slice 7b's six-input ridge and on an additive
/// function with no low-dimensional structure:
/// <list type="bullet">
/// <item>each criterion against its direct computation, the caveat present, and the refusals;</item>
/// <item>cross-fitting's folds disjoint and covering every case, and each fold's frame read from the other folds alone: a fold's
/// responses perturbed leave its own frame the same bits and move the others';</item>
/// <item>out of sample the refit wins on the ridge and does not win where there is no low-dimensional structure.</item>
/// </list>
/// </summary>
public sealed class BassProjectionComparisonTests(ITestOutputHelper output)
{
    /// <summary>
    /// On the ridge with 500 held-out cases: the basis functions per draw, PSIS-LOO, the held-out RMSE and percentage against their
    /// direct computation; no evidence from untempered runs; the refit's mapped subspace on the frame at k = d; each native draw's
    /// spread about the frame against its principal angles; the caveat. Refused: an analysis the frame was not read from, held-out
    /// cases without responses, and one fold.
    /// </summary>
    [Fact]
    public void Criteria_AgainstTheirDirectComputation()
    {
        BassProjectedFitTests.RidgeFits f = BassProjectedFitTests.Fixture.Value;
        (_, double[][] test) = BassProjectedFitTests.RidgeData(500, 31);
        double[] w = BassProjectedFitTests.W;
        var response = new double[500, 1];
        for (int r = 0; r < 500; r++) response[r, 0] = BassProjectedFitTests.Ridge(Enumerable.Range(0, w.Length).Sum(i => w[i] * test[i][r]));
        PredictionCases cases = BassProjectedFitTests.Cases(test);
        BassProjectionComparison c = f.Refit.Compare(f.Native, f.Data, f.Analysis, new ProjectionComparisonSpec { TestCases = cases, TestResponse = response });

        foreach ((ProjectionFitReport side, BassFit fit, BassData data, double[,] mean) in new[]
        {
            (c.Native, f.Native, f.Data, f.Native.Predict(f.Data, cases).Mean()),
            (c.Refit, f.Refit.Fit, f.Refit.Data, f.Refit.Predict(cases).Mean()),
        })
        {
            Assert.Equal(fit.Draws.SelectMany(d => d).Average(d => d.Structure.TermIds.Length), side.BasisFunctions, 12);
            Assert.Equal(PredictiveChecks.PsisLoo(fit.PointwiseLogDensity(data)).Elpd, side.Loo.Elpd);
            Assert.Null(side.Evidence);
            double rmse = Math.Sqrt(Enumerable.Range(0, 500).Average(r => Math.Pow(mean[r, 0] - response[r, 0], 2)));
            double mu = Enumerable.Range(0, 500).Average(r => response[r, 0]);
            double sd = Math.Sqrt(Enumerable.Range(0, 500).Sum(r => Math.Pow(response[r, 0] - mu, 2)) / 499);
            Assert.True(Math.Abs(side.HeldOutRmse!.Value - rmse) <= 1e-12 * rmse, $"{side.HeldOutRmse} against {rmse}");
            Assert.True(Math.Abs(side.HeldOutExplained!.Value - 100.0 * (1.0 - rmse / sd)) <= 1e-9);
        }
        Assert.True(c.RefitToFrame[^1] <= 1e-12, $"{c.RefitToFrame[^1]:E2}");
        for (int s = 0; s < c.FrameSpread.Length; s++)
        {
            double[,] v = f.Analysis.Posterior.Draws[s].Vectors;
            var lead = new double[v.GetLength(0), 1];
            for (int i = 0; i < lead.GetLength(0); i++) lead[i, 0] = v[i, 0];
            Assert.Equal(SubspacePosterior.PrincipalAngles(lead, f.Projection.Frame), c.FrameSpread[s]);
        }
        Assert.Contains("favour the refit", c.Caveat);
        double[] spread = [.. c.FrameSpread.Select(a => a[0] * 180 / Math.PI).Order()];
        output.WriteLine($"native: {c.Native.BasisFunctions:F1} basis functions, LOO elpd {c.Native.Loo.Elpd:F2}, held-out RMSE {c.Native.HeldOutRmse:F4} ({c.Native.HeldOutExplained:F2}%); " +
                         $"refit: {c.Refit.BasisFunctions:F1}, {c.Refit.Loo.Elpd:F2}, {c.Refit.HeldOutRmse:F4} ({c.Refit.HeldOutExplained:F2}%); " +
                         $"native draws about the frame: median {spread[spread.Length / 2]:F2}°, largest {spread[^1]:F2}°");

        var otherLaw = new ActiveSubspaceSpec { Measures = new Dictionary<string, InputMeasure> { ["x0"] = new MixtureMeasure([new NormalComponent(1.0, 0.2, 0.1)]) } };
        Assert.Throws<ArgumentException>(() => f.Refit.Compare(f.Native, f.Data, f.Native.ActiveSubspace(f.Data, otherLaw)));
        Assert.Throws<ArgumentException>(() => f.Refit.Compare(f.Native, f.Data, f.Analysis, new ProjectionComparisonSpec { TestCases = cases }));
        Assert.Throws<ArgumentOutOfRangeException>(() => f.Refit.Compare(f.Native, f.Data, f.Analysis, new ProjectionComparisonSpec { Folds = 1 }));
    }

    /// <summary>
    /// Three folds on the ridge's 200 cases: disjoint, covering every case, of sizes 66 to 67. With fold 0's responses raised by 50,
    /// fold 0's frame is the same bits, read from the other folds alone, folds 1 and 2's frames move, and fold 0's scores move. Out of
    /// sample the refit's pooled RMSE is below the native fit's.
    /// </summary>
    [Fact]
    public void CrossFitting_FoldsAndFrames()
    {
        BassProjectedFitTests.RidgeFits f = BassProjectedFitTests.Fixture.Value;
        var spec = new ProjectionComparisonSpec { Folds = 3, FoldSeed = 7 };
        CrossFitReport cross = f.Refit.Compare(f.Native, f.Data, f.Analysis, spec).CrossFitted!;
        int n = f.Data.Response.Cases;
        int[] all = [.. cross.Folds.SelectMany(fold => fold).Order()];
        Assert.Equal(Enumerable.Range(0, n), all);
        Assert.InRange(cross.Folds.Max(fold => fold.Length) - cross.Folds.Min(fold => fold.Length), 0, 1);

        // Fold 0's responses raised by 50: the data the other folds' fits see change, fold 0's own fits do not.
        double[] y = ((ScalarResponse)f.Data.Response).Y.ToArray();
        foreach (int r in cross.Folds[0]) y[r] += 50.0;
        var raised = new BassData(f.Data.Predictors, new ScalarResponse(y));
        BassFit nativeRaised = Bass.Fit(f.Native.Config, raised);
        ActiveSubspaceAnalysis analysisRaised = nativeRaised.ActiveSubspace(raised);
        BassProjection projectionRaised = BassProjection.Of(nativeRaised.Transform, analysisRaised, 1);
        BassProjectedFit refitRaised = BassProjectedFit.Of(projectionRaised, raised, f.Refit.Fit.Config);
        CrossFitReport moved = refitRaised.Compare(nativeRaised, raised, analysisRaised, spec).CrossFitted!;
        Assert.Equal(cross.Folds.Select(fold => fold.ToArray()), moved.Folds.Select(fold => fold.ToArray()));
        Assert.Equal(Flatten(cross.Frames[0].Frame), Flatten(moved.Frames[0].Frame));
        Assert.NotEqual(Flatten(cross.Frames[1].Frame), Flatten(moved.Frames[1].Frame));
        Assert.NotEqual(Flatten(cross.Frames[2].Frame), Flatten(moved.Frames[2].Frame));
        Assert.NotEqual(cross.NativeRmse[0], moved.NativeRmse[0]);
        output.WriteLine($"folds of {string.Join(", ", cross.Folds.Select(fold => fold.Length))}; fold 0's frame the same bits with its responses raised, " +
                         $"folds 1 and 2's moved; out of sample: native RMSE {cross.NativePooled:F4} ({cross.NativeExplained:F2}%), refit {cross.RefitPooled:F4} ({cross.RefitExplained:F2}%)");
        Assert.True(cross.RefitPooled < cross.NativePooled, $"{cross.RefitPooled:F4} against {cross.NativePooled:F4}");
    }

    /// <summary>
    /// An additive function of four inputs with no low-dimensional structure, Σᵢ sin(2πxᵢ), 200 cases, noise sd 0.05: E[C] is near a
    /// multiple of the identity, and at d = 1 the cross-fitted refit does not beat the native fit.
    /// </summary>
    [Fact]
    public void NoLowDimensionalStructure_RefitDoesNotWin()
    {
        var rng = new Xoshiro256PlusPlus(20261112);
        const int n = 200, p = 4;
        double[][] x = [.. Enumerable.Range(0, p).Select(_ => Enumerable.Range(0, n).Select(_ => rng.NextDouble()).ToArray())];
        double[] y = [.. Enumerable.Range(0, n).Select(r => Enumerable.Range(0, p).Sum(i => Math.Sin(2 * Math.PI * x[i][r])) + 0.05 * Normal.SampleStandard(rng))];
        var data = new BassData([.. Enumerable.Range(0, p).Select(i => new ContinuousPredictor($"x{i}", x[i], new InputDomain(0.0, 1.0)))], new ScalarResponse(y));
        BassConfig config = BassProjectedFitTests.Config();
        BassFit native = Bass.Fit(config, data);
        ActiveSubspaceAnalysis analysis = native.ActiveSubspace(data);
        BassProjectedFit refit = BassProjectedFit.Of(BassProjection.Of(native.Transform, analysis, 1), data, config);
        CrossFitReport cross = refit.Compare(native, data, analysis, new ProjectionComparisonSpec { Folds = 3, FoldSeed = 7 }).CrossFitted!;
        double[] lambda = analysis.Posterior.MeanCEigenvalues;
        output.WriteLine($"E[C]'s eigenvalues over the first: {string.Join(", ", lambda.Select(l => (l / lambda[0]).ToString("F3")))}; out of sample: native RMSE " +
                         $"{cross.NativePooled:F4} ({cross.NativeExplained:F2}%), refit {cross.RefitPooled:F4} ({cross.RefitExplained:F2}%)");
        Assert.True(cross.RefitPooled > cross.NativePooled, $"{cross.RefitPooled:F4} against {cross.NativePooled:F4}");
    }

    private static double[] Flatten(double[,] m) => [.. m.Cast<double>()];
}
