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
/// E8's slice 7b gates (the projected refit, <see cref="BassProjectedFit"/>), on a ridge function of six inputs from outside the
/// model, g(wᵀx) with g(u) = sin(4u) + u²:
/// <list type="bullet">
/// <item>the gradient in the original inputs against central differences of the prediction, per draw;</item>
/// <item>the mapped C lies in span(W₁), and equals C read directly from the original-input gradient at the native points the
/// pushforward came from;</item>
/// <item>the ridge recovered at d = 1 with fewer basis functions and a lower held-out error than the native fit;</item>
/// <item>the readout's refusals;</item>
/// <item>under the qualification switch, 2307.11241 §4.1.3's function refitted at d = 1 to 6 against its table (exploration).</item>
/// </list>
/// </summary>
public sealed class BassProjectedFitTests(ITestOutputHelper output)
{
    private static readonly double[] W = Normalize([3.0, -2.0, 1.0, 0.5, 0.0, 0.0]);

    private static double Ridge(double u) => Math.Sin(4.0 * u) + u * u;

    private static double[] Normalize(double[] v)
    {
        double n = Math.Sqrt(v.Sum(x => x * x));
        return [.. v.Select(x => x / n)];
    }

    private static BassConfig Config() => SensitivityFixture.Config(maxBasis: 30) with
    {
        Sampling = new BassSampling { Chains = 2, MasterSeed = 29, BurnIn = 2_000, Samples = 40, Thin = 10 },
    };

    // n cases uniform in [0, 1]^p, each input's domain declared [0, 1], and y = g(wᵀx) + noise of sd 0.05.
    private static (BassData Data, double[][] X) RidgeData(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[][] x = [.. Enumerable.Range(0, W.Length).Select(_ => Enumerable.Range(0, n).Select(_ => rng.NextDouble()).ToArray())];
        double[] y = [.. Enumerable.Range(0, n).Select(r => Ridge(Enumerable.Range(0, W.Length).Sum(i => W[i] * x[i][r])) + 0.05 * Normal.SampleStandard(rng))];
        BassPredictor[] predictors = [.. Enumerable.Range(0, W.Length).Select(i => new ContinuousPredictor($"x{i}", x[i], new InputDomain(0.0, 1.0)))];
        return (new BassData(predictors, new ScalarResponse(y)), x);
    }

    private sealed record RidgeFits(BassData Data, BassFit Native, ActiveSubspaceAnalysis Analysis, BassProjection Projection, BassProjectedFit Refit);

    private static readonly Lazy<RidgeFits> Fixture = new(() =>
    {
        (BassData data, _) = RidgeData(200, 20261111);
        BassFit native = Bass.Fit(Config(), data);
        ActiveSubspaceAnalysis analysis = native.ActiveSubspace(data);
        BassProjection projection = BassProjection.Of(native.Transform, analysis, 1);
        return new RidgeFits(data, native, analysis, projection, BassProjectedFit.Of(projection, data, Config()));
    });

    private static PredictionCases Cases(double[][] x) => new([.. x.Select((column, i) => new ContinuousColumn($"x{i}", column))]);

    // ── The gradient ─────────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>At 20 random interior points and 5 draws, each partial in the original inputs against the central difference of the
    /// draw's prediction at h = 10⁻⁷, within 10⁻⁶ of the gradient's size (degree-one hinges, piecewise linear between knots).</summary>
    [Fact]
    public void Gradient_IsTheChainRule()
    {
        RidgeFits f = Fixture.Value;
        var rng = new Xoshiro256PlusPlus(5);
        const int n = 20, p = 6;
        double[][] x = [.. Enumerable.Range(0, p).Select(_ => Enumerable.Range(0, n).Select(_ => 0.05 + 0.9 * rng.NextDouble()).ToArray())];
        int total = f.Refit.Fit.Draws.Sum(c => c.Count);
        int[] draws = [0, total / 4, total / 2, 3 * total / 4, total - 1];
        BassGradient g = f.Refit.Gradient(Cases(x), draws);
        const double h = 1e-7;
        double worst = 0.0;
        for (int i = 0; i < p; i++)
        {
            double[][] up = [.. x.Select(c => (double[])c.Clone())], down = [.. x.Select(c => (double[])c.Clone())];
            for (int r = 0; r < n; r++)
            {
                up[i][r] += h;
                down[i][r] -= h;
            }
            var options = new BassPredictOptions { Draws = draws };
            double[,,] lu = f.Refit.Predict(Cases(up), options).Latent, ld = f.Refit.Predict(Cases(down), options).Latent;
            for (int s = 0; s < draws.Length; s++)
                for (int r = 0; r < n; r++)
                {
                    double difference = (lu[s, r, 0] - ld[s, r, 0]) / (2.0 * h), size = 0.0;
                    for (int k = 0; k < p; k++) size = Math.Max(size, Math.Abs(g.Partials[k][s, r, 0]));
                    worst = Math.Max(worst, Math.Abs(g.Partials[i][s, r, 0] - difference) / Math.Max(size, 1e-12));
                }
        }
        output.WriteLine($"{n} points × {draws.Length} draws × {p} inputs: worst |∂f/∂x − central difference| {worst:E2} of the gradient's size");
        Assert.True(worst <= 1e-6, $"{worst:E2}");
        Assert.Equal(f.Projection.Inputs, g.Inputs);
    }

    // ── The readout ──────────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// The readout at 2,048 pushforward points: every mapped C lies in span(W₁) to 10⁻¹² of its size, and equals C read directly
    /// from the original-input gradient at the native points the pushforward came from, in the box metric, to 10⁻¹⁰.
    /// </summary>
    [Fact]
    public void Readout_MapsIntoTheFrame()
    {
        RidgeFits f = Fixture.Value;
        const int points = 2048, seed = 3;
        ProjectedActiveSubspace read = f.Refit.ActiveSubspace(new ActiveSubspaceSpec { Points = points, Seed = seed });
        double[,] w1 = f.Projection.Frame;
        int p = f.Projection.Inputs.Count;
        PredictionCases native = BassActiveSubspaces.Points(f.Native.Transform, f.Analysis.Measure!, points, seed);
        int[] draws = [.. read.Active.DrawIndex];
        BassGradient g = f.Refit.Gradient(native, draws);
        IReadOnlyList<double> widths = f.Projection.Widths;
        double outside = 0.0, direct = 0.0;
        for (int s = 0; s < draws.Length; s++)
        {
            double[,] c = read.InputSpace.Draws[s].C;
            double size = 0.0;
            foreach (double v in c) size = Math.Max(size, Math.Abs(v));
            // (I − W₁W₁ᵀ)C.
            for (int i = 0; i < p; i++)
                for (int k = 0; k < p; k++)
                {
                    double projected = 0.0;
                    for (int l = 0; l < p; l++) projected += w1[i, 0] * w1[l, 0] * c[l, k];
                    outside = Math.Max(outside, Math.Abs(c[i, k] - projected) / size);
                }
            var reference = new double[p, p];
            for (int r = 0; r < points; r++)
                for (int i = 0; i < p; i++)
                    for (int k = 0; k < p; k++)
                        reference[i, k] += g.Partials[i][s, r, 0] * widths[i] * g.Partials[k][s, r, 0] * widths[k] / points;
            direct = Math.Max(direct, FunctionalActiveSubspaceTests.RelativeGap(c, reference));
        }
        output.WriteLine($"{draws.Length} draws: C off span(W₁) {outside:E2}; against the direct reading {direct:E2}; " +
                         $"E[C]'s leading eigenvalue {read.InputSpace.MeanCEigenvalues[0]:F4}, the native analysis's {f.Analysis.Posterior.MeanCEigenvalues[0]:F4}");
        Assert.True(outside <= 1e-12, $"{outside:E2}");
        Assert.True(direct <= 1e-10, $"{direct:E2}");
    }

    /// <summary>Refused: a declared product law on the active variables, another metric, the exact route, and k above d.</summary>
    [Fact]
    public void Readout_RefusesWhatItDoesNotRead()
    {
        RidgeFits f = Fixture.Value;
        Assert.Throws<ArgumentException>(() => f.Refit.ActiveSubspace(new ActiveSubspaceSpec { Measures = new Dictionary<string, InputMeasure> { ["u1"] = new UniformMeasure() } }));
        Assert.Throws<ArgumentException>(() => f.Refit.ActiveSubspace(new ActiveSubspaceSpec { Metric = InputMetric.Raw }));
        Assert.Throws<ArgumentException>(() => f.Refit.ActiveSubspace(new ActiveSubspaceSpec { Route = ActiveSubspaceRoute.Exact }));
        Assert.Throws<ArgumentOutOfRangeException>(() => f.Refit.ActiveSubspace(new ActiveSubspaceSpec { Dimension = 2 }));
    }

    // ── Recovery ─────────────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// The ridge recovered at d = 1: the frame's angle to w; the refit's posterior mean against g at 1,000 held-out points, with a lower
    /// RMSE than the native fit's; and fewer basis functions per draw on average.
    /// </summary>
    [Fact]
    public void Ridge_RecoveredWithFewerTermsAndLowerError()
    {
        RidgeFits f = Fixture.Value;
        (_, double[][] test) = RidgeData(1000, 77);
        double[] truth = [.. Enumerable.Range(0, 1000).Select(r => Ridge(Enumerable.Range(0, W.Length).Sum(i => W[i] * test[i][r])))];
        double[,] native = f.Native.Predict(f.Data, Cases(test)).Mean(), refit = f.Refit.Predict(Cases(test)).Mean();
        double rmseNative = Math.Sqrt(Enumerable.Range(0, 1000).Average(r => Math.Pow(native[r, 0] - truth[r], 2)));
        double rmseRefit = Math.Sqrt(Enumerable.Range(0, 1000).Average(r => Math.Pow(refit[r, 0] - truth[r], 2)));
        double termsNative = f.Native.Draws.SelectMany(c => c).Average(d => d.Structure.TermIds.Length);
        double termsRefit = f.Refit.Fit.Draws.SelectMany(c => c).Average(d => d.Structure.TermIds.Length);
        var w = new double[W.Length, 1];
        for (int i = 0; i < W.Length; i++) w[i, 0] = W[i];
        double angle = SubspacePosterior.PrincipalAngles(f.Projection.Frame, w)[0];
        output.WriteLine($"frame to w {angle * 180 / Math.PI:F3}°; held-out RMSE against g: native {rmseNative:F4}, refit {rmseRefit:F4}; " +
                         $"basis functions per draw: native {termsNative:F1}, refit {termsRefit:F1}; eigengap at d {f.Projection.Gaps[0]:F4}");
        Assert.True(rmseRefit < rmseNative, $"{rmseRefit:F4} against {rmseNative:F4}");
        Assert.True(termsRefit < termsNative, $"{termsRefit:F1} against {termsNative:F1}");
    }

    // ── The paper's function (exploration) ───────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// 2307.11241 §4.1.3: g(x) = f(x·t₁, x·t₂) with f = x₁² + x₁x₂ + x₂³/9, t₁ = (12, 9, 6, 1, 1, 1)/3, t₂ = (1, 1, 1, 15, 6, 6)/3, a
    /// Latin hypercube of n = 500 on [0, 1]⁶. A native fit, its E[C], and refits at d = 1 to 6, each read as the percentage of g's sd it
    /// accounts for at 2,000 held-out points, 100·(1 − RMSE/sd), against the paper's 91.06, 98.37, 99.03, 99.39, 99.45 and 99.45.
    /// Exploration, not validation; under the qualification switch.
    /// </summary>
    [Fact]
    public void PapersFunction_RefitsAtEachDimension()
    {
        if (!BassQualification.Enabled) return;
        double[] t1 = [4.0, 3.0, 2.0, 1.0 / 3, 1.0 / 3, 1.0 / 3], t2 = [1.0 / 3, 1.0 / 3, 1.0 / 3, 5.0, 2.0, 2.0];
        double G(double[] x)
        {
            double a = x.Zip(t1, (u, v) => u * v).Sum(), b = x.Zip(t2, (u, v) => u * v).Sum();
            return a * a + a * b + b * b * b / 9.0;
        }
        var rng = new Xoshiro256PlusPlus(20260714);
        const int n = 500, p = 6;
        double[][] x = [.. Enumerable.Range(0, p).Select(_ => Latin(rng, n))];
        double[] y = [.. Enumerable.Range(0, n).Select(r => G([.. x.Select(c => c[r])]))];
        var data = new BassData([.. Enumerable.Range(0, p).Select(i => new ContinuousPredictor($"x{i}", x[i], new InputDomain(0.0, 1.0)))], new ScalarResponse(y));
        BassConfig config = Config();
        BassFit native = Bass.Fit(config, data);
        ActiveSubspaceAnalysis analysis = native.ActiveSubspace(data);
        double[][] test = [.. Enumerable.Range(0, p).Select(_ => Enumerable.Range(0, 2000).Select(_ => rng.NextDouble()).ToArray())];
        double[] truth = [.. Enumerable.Range(0, 2000).Select(r => G([.. test.Select(c => c[r])]))];
        double sd = Math.Sqrt(truth.Select(v => v - truth.Average()).Sum(v => v * v) / (truth.Length - 1));
        double[] paper = [91.057, 98.368, 99.025, 99.393, 99.448, 99.452];
        var lines = new List<string>();
        double[,] mean = native.Predict(data, Cases(test)).Mean();
        lines.Add($"native: {100.0 * (1.0 - Rmse(mean, truth) / sd):F3}");
        for (int d = 1; d <= p; d++)
        {
            BassProjectedFit refit = BassProjectedFit.Of(BassProjection.Of(native.Transform, analysis, d), data, config);
            double explained = 100.0 * (1.0 - Rmse(refit.Predict(Cases(test)).Mean(), truth) / sd);
            lines.Add($"d = {d}: {explained:F3} (paper {paper[d - 1]:F3})");
        }
        output.WriteLine("percent of g's sd accounted for: " + string.Join("; ", lines) +
                         $"; √eigenvalues over the first {string.Join(", ", analysis.Posterior.MeanCEigenvalues.Select(v => Math.Sqrt(v / analysis.Posterior.MeanCEigenvalues[0]).ToString("F3")))}");
    }

    /// <summary>
    /// The cost row (indicative; the qualification switch): on the ridge's data (n = 200, p = 6, the same configuration), the native
    /// fit's and the refit's run times and basis functions per draw, and the projected readout per draw at 4,096 pushforward points,
    /// the best of three.
    /// </summary>
    [Fact]
    public void CostRow()
    {
        if (!BassQualification.Enabled) return;
        RidgeFits f = Fixture.Value;
        double Seconds(Action run)
        {
            double best = double.PositiveInfinity;
            for (int k = 0; k < 3; k++)
            {
                long t0 = BassCost.Tick();
                run();
                best = Math.Min(best, (BassCost.Tick() - t0) / BassCost.TicksPerSecond);
            }
            return best;
        }
        double native = Seconds(() => Bass.Fit(Config(), f.Data)), refit = Seconds(() => BassProjectedFit.Of(f.Projection, f.Data, Config()));
        int draws = f.Refit.Fit.Draws.Sum(c => c.Count);
        double readout = Seconds(() => f.Refit.ActiveSubspace(new ActiveSubspaceSpec { Points = 4096 })) * 1e6 / draws;
        output.WriteLine($"clock: {(BassCost.Cycles ? "thread cycles" : "wall")}; native fit {native:F2} s, {f.Native.Draws.SelectMany(c => c).Average(d => d.Structure.TermIds.Length):F1} " +
                         $"basis functions per draw; refit {refit:F2} s, {f.Refit.Fit.Draws.SelectMany(c => c).Average(d => d.Structure.TermIds.Length):F1}; " +
                         $"the projected readout at 4,096 points {readout:F0} µs per draw");
    }

    private static double Rmse(double[,] mean, double[] truth) => Math.Sqrt(Enumerable.Range(0, truth.Length).Average(r => Math.Pow(mean[r, 0] - truth[r], 2)));

    private static double[] Latin(Xoshiro256PlusPlus rng, int n)
    {
        int[] order = [.. Enumerable.Range(0, n)];
        for (int i = n - 1; i > 0; i--)
        {
            int j = rng.NextInt(i + 1);
            (order[i], order[j]) = (order[j], order[i]);
        }
        return [.. order.Select(k => (k + rng.NextDouble()) / n)];
    }
}
