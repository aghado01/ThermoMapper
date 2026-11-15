using System;
using System.Collections.Generic;
using System.Diagnostics;
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
/// E8's slice 8c gates (<see cref="BassProjectedBasisFit"/>), on a functional ridge g(wᵀx)·φ(t) over five inputs and twenty output
/// locations, from outside the model:
/// <list type="bullet">
/// <item>the decomposition reused, the same object; the component mode carried over, apart and joint; the training-range rule;</item>
/// <item>the reconstruction's gradient in the original inputs against central differences of each joint draw's prediction;</item>
/// <item>the refit's mapped C in span(W₁);</item>
/// <item>recovery at d = 1: the frame on w, and the refit with fewer basis functions and a lower held-out error than the native fit;</item>
/// <item>the comparison's criteria against their direct computation, and the refusals.</item>
/// </list>
/// </summary>
public sealed class BassProjectedBasisFitTests(ITestOutputHelper output)
{
    private static readonly double[] W = Normalize([0.7, 0.5, -0.3, 0.4, 0.1]);

    private static double[] Normalize(double[] v)
    {
        double n = Math.Sqrt(v.Sum(x => x * x));
        return [.. v.Select(x => x / n)];
    }

    private static double G(double u) => Math.Sin(3.0 * u) + u * u;

    private static double Phi(double t) => 1.0 + Math.Sin(Math.PI * t) + 0.5 * t;

    private static readonly double[] T = [.. Enumerable.Range(0, 20).Select(l => l / 19.0)];

    // n cases uniform in [0, 1]^5, domains declared [0, 1], and y(x, t) = g(wᵀx)·φ(t) plus noise of sd `noise`.
    internal static (IReadOnlyList<BassPredictor> Predictors, double[,] Y, double[][] X) Data(int n, int seed, double noise)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[][] x = [.. Enumerable.Range(0, W.Length).Select(_ => Enumerable.Range(0, n).Select(_ => rng.NextDouble()).ToArray())];
        var y = new double[n, T.Length];
        for (int r = 0; r < n; r++)
        {
            double g = G(Enumerable.Range(0, W.Length).Sum(i => W[i] * x[i][r]));
            for (int l = 0; l < T.Length; l++) y[r, l] = g * Phi(T[l]) + noise * Normal.SampleStandard(rng);
        }
        BassPredictor[] predictors = [.. Enumerable.Range(0, W.Length).Select(i => new ContinuousPredictor($"x{i}", x[i], new InputDomain(0.0, 1.0)))];
        return (predictors, y, x);
    }

    private static BassConfig Config() => SensitivityFixture.Config(maxBasis: 30) with
    {
        Sampling = new BassSampling { Chains = 2, MasterSeed = 37, BurnIn = 3_000, Samples = 30, Thin = 10 },
    };

    private sealed record Fits(BassBasisFit Native, ActiveSubspaceAnalysis Analysis, BassProjection Projection, BassProjectedBasisFit Refit, double NativeSeconds,
                               double RefitSeconds);

    private static readonly Lazy<Fits> Fixture = new(() =>
    {
        (IReadOnlyList<BassPredictor> predictors, double[,] y, _) = Data(300, 20261115, 0.02);
        BassDecomposition decomposition = BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 2 });
        long start = Stopwatch.GetTimestamp();
        BassBasisFit native = BassComponents.Fit(decomposition, predictors, Config());
        double nativeSeconds = Stopwatch.GetElapsedTime(start).TotalSeconds;
        ActiveSubspaceAnalysis analysis = native.ActiveSubspace();
        BassProjection projection = BassProjection.Of(native.Components[0].Transform, analysis, 1);
        start = Stopwatch.GetTimestamp();
        BassProjectedBasisFit refit = BassProjectedBasisFit.Of(projection, native);
        return new Fits(native, analysis, projection, refit, nativeSeconds, Stopwatch.GetElapsedTime(start).TotalSeconds);
    });

    /// <summary>The refit's decomposition is the native fit's, the same object; components fitted apart are refitted apart under the
    /// training-range rule; a joint native fit is refitted jointly, its decomposition reused too.</summary>
    [Fact]
    public void Decomposition_Reused_AndTheModeCarriedOver()
    {
        Fits f = Fixture.Value;
        Assert.Same(f.Native.Decomposition, f.Refit.Fit.Decomposition);
        Assert.Null(f.Refit.Fit.Joint);
        Assert.All(f.Refit.Fit.Components, c => Assert.Equal(HingeNormalization.TrainingRange, c.Config.Structure.HingeNormalization));
        Assert.Equal(f.Projection.Names, f.Refit.Fit.ComponentData(0).Predictors.Select(p => p.Name));

        (IReadOnlyList<BassPredictor> predictors, double[,] y, _) = Data(80, 7, 0.02);
        BassDecomposition decomposition = BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 2 });
        BassBasisFit joint = BassComponents.FitJoint(decomposition, predictors, Config() with { Sampling = Config().Sampling with { BurnIn = 300, Samples = 10 } });
        ActiveSubspaceAnalysis analysis = joint.ActiveSubspace();
        BassProjectedBasisFit refit = BassProjectedBasisFit.Of(BassProjection.Of(joint.Components[0].Transform, analysis, 1), joint);
        Assert.NotNull(refit.Fit.Joint);
        Assert.Same(joint.Decomposition, refit.Fit.Decomposition);
        Assert.Equal(joint.MasterSeed, refit.Fit.MasterSeed);
    }

    /// <summary>At 15 interior points and 4 joint draws, each partial of the reconstruction in the original inputs against the central
    /// difference of that draw's prediction at h = 10⁻⁷, within 10⁻⁶ of the gradient's size (degree-one hinges, piecewise linear).</summary>
    [Fact]
    public void Gradient_IsTheChainRuleThroughTheReconstruction()
    {
        Fits f = Fixture.Value;
        var rng = new Xoshiro256PlusPlus(11);
        const int n = 15, p = 5;
        double[][] x = [.. Enumerable.Range(0, p).Select(_ => Enumerable.Range(0, n).Select(_ => 0.1 + 0.8 * rng.NextDouble()).ToArray())];
        int[] draws = [0, 17, 33, 59];
        BassGradient g = f.Refit.Gradient(BassProjectedFitTests.Cases(x), draws);
        const double h = 1e-7;
        double worst = 0.0, size = 0.0;
        for (int b = 0; b < draws.Length; b++)
            for (int i = 0; i < p; i++)
            {
                double[][] up = [.. x.Select(c => c.ToArray())], down = [.. x.Select(c => c.ToArray())];
                for (int r = 0; r < n; r++)
                {
                    up[i][r] += h;
                    down[i][r] -= h;
                }
                var options = new BassPredictOptions { Draws = [draws[b]] };
                double[,,] yu = f.Refit.Predict(BassProjectedFitTests.Cases(up), options).Latent, yd = f.Refit.Predict(BassProjectedFitTests.Cases(down), options).Latent;
                for (int r = 0; r < n; r++)
                    for (int l = 0; l < T.Length; l++)
                    {
                        double difference = (yu[0, r, l] - yd[0, r, l]) / (2.0 * h), partial = g.Partials[i][b, r, l];
                        worst = Math.Max(worst, Math.Abs(partial - difference));
                        size = Math.Max(size, Math.Abs(partial));
                    }
            }
        output.WriteLine($"worst |∂y/∂x − central difference| {worst:E2} against the gradient's size {size:E2}");
        Assert.True(worst <= 1e-6 * size, $"{worst:E2} against {size:E2}");
    }

    /// <summary>The refit's mapped C, every draw, lies in span(W₁): (I − W₁W₁ᵀ)C(I − W₁W₁ᵀ) within 10⁻¹² of ‖C‖.</summary>
    [Fact]
    public void MappedC_InTheFramesSpan()
    {
        Fits f = Fixture.Value;
        ProjectedActiveSubspace read = f.Refit.ActiveSubspace(new ActiveSubspaceSpec { Points = 2048, Seed = 3 });
        double[,] w = f.Projection.Frame;
        int p = w.GetLength(0);
        var complement = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++) complement[i, j] = (i == j ? 1.0 : 0.0) - w[i, 0] * w[j, 0];
        double worst = 0.0;
        foreach (DrawSubspace draw in read.InputSpace.Draws)
        {
            double[,] c = draw.C, off = Multiply(Multiply(complement, c), complement);
            worst = Math.Max(worst, Norm(off) / Norm(c));
        }
        output.WriteLine($"worst ‖P⊥CP⊥‖/‖C‖ {worst:E2} over {read.InputSpace.Draws.Count} draws");
        Assert.True(worst <= 1e-12, $"{worst:E2}");
    }

    /// <summary>
    /// Recovery at d = 1 and the comparison: the native frame within 5° of w; with 300 held-out cases, the refit uses fewer basis functions
    /// per joint draw and has a lower held-out RMSE than the native fit; each criterion against its direct computation (the basis
    /// functions, PSIS-LOO on the coefficients, the RMSE), no evidence and no density, the refit's mapped subspace on the frame.
    /// </summary>
    [Fact]
    public void Recovery_AndTheComparisonsCriteria()
    {
        Fits f = Fixture.Value;
        double cosine = Math.Abs(Enumerable.Range(0, W.Length).Sum(i => f.Projection.Frame[i, 0] * W[i]));
        double angle = Math.Acos(Math.Min(1.0, cosine)) * 180.0 / Math.PI;
        (_, double[,] response, double[][] x) = Data(300, 101, 0.0);
        PredictionCases cases = BassProjectedFitTests.Cases(x);
        BassProjectionComparison c = f.Refit.Compare(f.Native, f.Analysis, new ProjectionComparisonSpec { TestCases = cases, TestResponse = response });
        foreach ((ProjectionFitReport side, BassBasisFit fit, double[,] mean) in new[]
        {
            (c.Native, f.Native, f.Native.Predict(cases).Mean()),
            (c.Refit, f.Refit.Fit, f.Refit.Predict(cases).Mean()),
        })
        {
            BassDraw[][] pooled = [.. fit.Components.Select(k => k.Draws.SelectMany(chain => chain).ToArray())];
            double terms = Enumerable.Range(0, fit.JointDraws.Count).Average(b => Enumerable.Range(0, pooled.Length).Sum(k => pooled[k][fit.JointDraws.Draw(k, b)].Structure.TermIds.Length));
            Assert.Equal(terms, side.BasisFunctions, 12);
            Assert.Equal(PredictiveChecks.PsisLoo(fit.PointwiseLogDensity()).Elpd, side.Loo.Elpd);
            double rmse = Math.Sqrt(Enumerable.Range(0, 300).Sum(r => Enumerable.Range(0, T.Length).Sum(l => Math.Pow(mean[r, l] - response[r, l], 2))) / (300.0 * T.Length));
            Assert.True(Math.Abs(side.HeldOutRmse!.Value - rmse) <= 1e-12 * rmse, $"{side.HeldOutRmse} against {rmse}");
            Assert.Null(side.Evidence);
            Assert.Null(side.HeldOutLogDensity);
        }
        Assert.True(c.RefitToFrame[^1] <= 1e-12, $"{c.RefitToFrame[^1]:E2}");
        Assert.Null(c.CrossFitted);
        output.WriteLine($"frame against w {angle:F2}°; native {c.Native.BasisFunctions:F1} basis functions per joint draw, held-out RMSE {c.Native.HeldOutRmse:F4} " +
                         $"({c.Native.HeldOutExplained:F2}%), LOO {c.Native.Loo.Elpd:F1}; refit {c.Refit.BasisFunctions:F1}, {c.Refit.HeldOutRmse:F4} ({c.Refit.HeldOutExplained:F2}%), LOO {c.Refit.Loo.Elpd:F1}; " +
                         $"indicative fit times: native {f.NativeSeconds:F2} s, refit {f.RefitSeconds:F2} s");
        Assert.True(angle <= 5.0, $"{angle:F2}°");
        Assert.True(c.Refit.BasisFunctions < c.Native.BasisFunctions, $"{c.Refit.BasisFunctions:F1} against {c.Native.BasisFunctions:F1}");
        Assert.True(c.Refit.HeldOutRmse < c.Native.HeldOutRmse, $"{c.Refit.HeldOutRmse:F4} against {c.Native.HeldOutRmse:F4}");
    }

    /// <summary>Refused: cross-fitting a basis fit, held-out weights, an analysis the frame was not read from, a projection on other
    /// inputs, and a declared product law for the refit's readout.</summary>
    [Fact]
    public void Refusals()
    {
        Fits f = Fixture.Value;
        Assert.Throws<ArgumentException>(() => f.Refit.Compare(f.Native, f.Analysis, new ProjectionComparisonSpec { Folds = 3 }));
        (_, double[,] response, double[][] x) = Data(5, 3, 0.0);
        Assert.Throws<ArgumentException>(() => f.Refit.Compare(f.Native, f.Analysis, new ProjectionComparisonSpec
        {
            TestCases = BassProjectedFitTests.Cases(x), TestResponse = response, TestWeights = [1.0, 1.0, 1.0, 1.0, 1.0],
        }));
        ActiveSubspaceAnalysis other = f.Native.ActiveSubspace(new ActiveSubspaceSpec { Measures = new Dictionary<string, InputMeasure> { ["x0"] = new MixtureMeasure([new NormalComponent(1.0, 0.2, 0.1)]) } });
        Assert.Throws<ArgumentException>(() => f.Refit.Compare(f.Native, other));
        BassProjectedFitTests.RidgeFits ridge = BassProjectedFitTests.Fixture.Value;
        Assert.Throws<ArgumentException>(() => BassProjectedBasisFit.Of(ridge.Projection, f.Native));
        Assert.Throws<ArgumentException>(() => f.Refit.ActiveSubspace(new ActiveSubspaceSpec { Measures = new Dictionary<string, InputMeasure> { ["u0"] = new UniformMeasure() } }));
    }

    private static double[,] Multiply(double[,] a, double[,] b)
    {
        int n = a.GetLength(0), k = a.GetLength(1), m = b.GetLength(1);
        var c = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
            {
                double s = 0.0;
                for (int q = 0; q < k; q++) s += a[i, q] * b[q, j];
                c[i, j] = s;
            }
        return c;
    }

    private static double Norm(double[,] a)
    {
        double s = 0.0;
        foreach (double x in a) s += x * x;
        return Math.Sqrt(s);
    }
}
