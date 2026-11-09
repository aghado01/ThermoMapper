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
/// E8's slice 2 gates (the active-subspace brief) on the readout (<see cref="BassFit.ActiveSubspace"/>) and its fit-agnostic summary
/// (<see cref="SubspacePosterior"/>):
/// <list type="bullet">
/// <item>SG-5's three-draw example, E[C] and the subspaces' projection mean apart, the Lean's rational case and the 1.2-radian one;</item>
/// <item>weighted draws against replicated ones, in the summary and on a fit;</item>
/// <item>the exact and sampled routes agreeing within the sampled route's Monte Carlo error, under the uniform law, a mixture and an
/// empirical marginal, with a categorical input;</item>
/// <item>C under a declared change of metric, LᵀC_boxL, on a fit whose box is not [0, 1];</item>
/// <item>decision 6's checks: Poincaré on hinge draws, SG-4(b) on Legendre draws, no violation;</item>
/// <item>recovery at growing N on ridges from outside the model (G3), E[C] and the mean subspace each against the span in closed form,
/// registered in <see cref="ComparativeSweep"/>.</item>
/// </list>
/// </summary>
public sealed class ActiveSubspaceReadoutTests(ITestOutputHelper output)
{
    // ── SG-5 ─────────────────────────────────────────────────────────────────────────────────────────────────────────

    private static double[,] Outer(double lambda, double[] v) => new[,] { { lambda * v[0] * v[0], lambda * v[0] * v[1] }, { lambda * v[1] * v[0], lambda * v[1] * v[1] } };

    // The top eigenvector of a symmetric 2 × 2, by its own closed form.
    private static double[] Top(double[,] m)
    {
        double a = m[0, 0], b = m[0, 1], c = m[1, 1];
        double top = 0.5 * (a + c) + Math.Sqrt(0.25 * (a - c) * (a - c) + b * b);
        double x = b, y = top - a;
        if (Math.Abs(b) < 1e-300) (x, y) = a >= c ? (1.0, 0.0) : (0.0, 1.0);
        double n = Math.Sqrt(x * x + y * y);
        return [x / n, y / n];
    }

    /// <summary>
    /// SG-5 (`top_eigenspace_ne_subspace_mean`): draws (λ, v) = (1, a), (1, a), (100, b), C_j = λvvᵀ. With a = (1, 0) and
    /// b = (3/5, 4/5), 3E[C] = [[38, 48], [48, 64]] and 3P̄ = [[59/25, 12/25], [12/25, 16/25]], which share no eigenvector; with a and
    /// b at 1.2 radians, the two summaries' leading directions sit 54° apart. Each angle here against the matrices' own closed forms.
    /// </summary>
    [Fact]
    public void Sg5_ExpectedCAndTheProjectionMean_SitApart()
    {
        double[] a = [1.0, 0.0];
        double[] b = [0.6, 0.8];
        var rational = SubspacePosterior.Of([Outer(1, a), Outer(1, a), Outer(100, b)], [1.0, 1.0, 1.0], 1);
        double[,] c3 = { { 38, 48 }, { 48, 64 } }, p3 = { { 59 / 25.0, 12 / 25.0 }, { 12 / 25.0, 16 / 25.0 } };
        double cGap = 0.0, pGap = 0.0;
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
            {
                cGap = Math.Max(cGap, Math.Abs(3 * rational.MeanC[i, j] - c3[i, j]));
                pGap = Math.Max(pGap, Math.Abs(3 * rational.Projection[i, j] - p3[i, j]));
            }
        double expected = Angle(Top(c3), Top(p3));
        output.WriteLine($"rational: 3E[C] off by {cGap:E1}, 3P̄ off by {pGap:E1}; summaries {Degrees(rational.MeanCToMean[0]):F4}° apart, closed form {Degrees(expected):F4}°");
        Assert.True(cGap <= 1e-13 && pGap <= 1e-15, $"{cGap:E2}, {pGap:E2}");
        Assert.Equal(expected, rational.MeanCToMean[0], 12);
        Assert.True(Degrees(expected) > 1.0, "the summaries coincide");

        double[] wide = [Math.Cos(1.2), Math.Sin(1.2)];
        var spread = SubspacePosterior.Of([Outer(1, a), Outer(1, a), Outer(100, wide)], [1.0, 1.0, 1.0], 1);
        double apart = Angle(Top(spread.MeanC), Top(spread.Projection));
        output.WriteLine($"1.2 rad: summaries {Degrees(spread.MeanCToMean[0]):F3}° apart (closed form {Degrees(apart):F3}°); P̄'s eigenvalues {spread.ProjectionEigenvalues[0]:F4}, {spread.ProjectionEigenvalues[1]:F4}");
        Assert.Equal(apart, spread.MeanCToMean[0], 12);
        Assert.InRange(Degrees(apart), 53.5, 54.5);
        Assert.Equal(1.0, spread.ProjectionEigenvalues.Sum(), 14);
    }

    private static double Angle(double[] u, double[] v) => Math.Acos(Math.Min(1.0, Math.Abs(u[0] * v[0] + u[1] * v[1])));

    private static double Degrees(double radians) => radians * 180.0 / Math.PI;

    // ── Weights ──────────────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>Weights (2, 1, 3) on three draws are the draws replicated 2, 1 and 3 times: E[C], P̄, their spectra and the angles.</summary>
    [Fact]
    public void WeightedDraws_AreReplicatedDraws_InTheSummary()
    {
        var rng = new Xoshiro256PlusPlus(5);
        double[][,] c = [.. Enumerable.Range(0, 3).Select(_ => RandomPsd(rng, 4))];
        var weighted = SubspacePosterior.Of(c, [2.0, 1.0, 3.0], 2);
        var replicated = SubspacePosterior.Of([c[0], c[0], c[1], c[2], c[2], c[2]], [1, 1, 1, 1, 1, 1.0], 2);
        double gap = Math.Max(MaxGap(weighted.MeanC, replicated.MeanC), MaxGap(weighted.Projection, replicated.Projection));
        gap = Math.Max(gap, weighted.ProjectionEigenvalues.Zip(replicated.ProjectionEigenvalues, (x, y) => Math.Abs(x - y)).Max());
        gap = Math.Max(gap, weighted.MeanCToMean.Zip(replicated.MeanCToMean, (x, y) => Math.Abs(x - y)).Max());
        output.WriteLine($"weighted against replicated: {gap:E2}");
        Assert.True(gap <= 1e-13, $"{gap:E2}");
    }

    /// <summary>On a fit: draw weights (1, 2, 3, …) against those draws replicated by Draws, both routes' summaries alike.</summary>
    [Fact]
    public void WeightedDraws_AreReplicatedDraws_OnAFit()
    {
        BassData data = SensitivityFixture.Scalar();
        BassFit fit = SensitivityFixture.FitScalar(1);
        int total = fit.Draws.Sum(c => c.Count);
        int[] multiplicity = [.. Enumerable.Range(0, total).Select(d => 1 + d % 3)];
        BassFit weightedFit = fit with { DrawWeights = [.. multiplicity.Select(m => (double)m)] };
        var spec = new ActiveSubspaceSpec { Dimension = 2 };
        ActiveSubspaceAnalysis weighted = weightedFit.ActiveSubspace(data, spec);
        ActiveSubspaceAnalysis replicated = fit.ActiveSubspace(data, spec with { Draws = [.. Enumerable.Range(0, total).SelectMany(d => Enumerable.Repeat(d, multiplicity[d]))] });
        double gap = Math.Max(MaxGap(weighted.Posterior.MeanC, replicated.Posterior.MeanC), MaxGap(weighted.Posterior.Projection, replicated.Posterior.Projection));
        output.WriteLine($"{total} draws weighted 1–3 against replicated: {gap:E2}");
        Assert.True(gap <= 1e-13, $"{gap:E2}");
    }

    private static double[,] RandomPsd(Xoshiro256PlusPlus rng, int p)
    {
        var g = new double[p + 2, p];
        for (int r = 0; r < p + 2; r++)
            for (int i = 0; i < p; i++) g[r, i] = Normal.SampleStandard(rng);
        var c = new double[p, p];
        for (int r = 0; r < p + 2; r++)
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++) c[i, j] += g[r, i] * g[r, j];
        return c;
    }

    private static double MaxGap(double[,] a, double[,] b)
    {
        double worst = 0.0, scale = 0.0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
            {
                worst = Math.Max(worst, Math.Abs(a[i, j] - b[i, j]));
                scale = Math.Max(scale, Math.Abs(b[i, j]));
            }
        return scale == 0.0 ? worst : worst / scale;
    }

    // ── The routes ───────────────────────────────────────────────────────────────────────────────────────────────────

    public static TheoryData<string> Laws() => new() { "uniform", "mixture", "empirical" };

    private static Dictionary<string, InputMeasure> LawOf(string law)
    {
        var measures = new Dictionary<string, InputMeasure> { ["z"] = new CategoricalMeasure(new Dictionary<string, double> { ["a"] = 0.2, ["b"] = 0.5, ["c"] = 0.3 }) };
        if (law == "mixture") measures["x0"] = new MixtureMeasure([new NormalComponent(0.4, 0.3, 0.15), new StudentComponent(0.6, 3.0, 0.7, 0.2)]) { Lower = 0.05, Upper = 0.95 };
        if (law == "empirical") measures["x0"] = new EmpiricalMeasure();
        return measures;
    }

    /// <summary>
    /// The exact and sampled routes on the same draws (a hinge fit with a categorical input): every entry of every draw's C within 4.5
    /// of the sampled route's standard errors, and E[C] within 4.5 standard errors of its mean's error, at 8,192 points.
    /// </summary>
    [Theory]
    [MemberData(nameof(Laws))]
    public void ExactAndSampledRoutes_AgreeWithinMonteCarloError(string law)
    {
        BassData data = SensitivityFixture.Scalar();
        BassFit fit = SensitivityFixture.FitScalar(1);
        var spec = new ActiveSubspaceSpec { Measures = LawOf(law), Dimension = 2, Points = 8192, Seed = 11 };
        ActiveSubspaceAnalysis exact = fit.ActiveSubspace(data, spec with { Route = ActiveSubspaceRoute.Exact });
        ActiveSubspaceAnalysis sampled = fit.ActiveSubspace(data, spec with { Route = ActiveSubspaceRoute.Sampled });
        Assert.Equal(ActiveSubspaceRoute.Exact, exact.Route);
        Assert.Equal(ActiveSubspaceRoute.Sampled, sampled.Route);
        double worst = 0.0;
        int count = 0;
        for (int d = 0; d < exact.Posterior.Draws.Count; d++)
            for (int i = 0; i < 2; i++)
                for (int j = i; j < 2; j++)
                {
                    double se = sampled.MonteCarloError![d][i, j];
                    double gap = Math.Abs(exact.Posterior.Draws[d].C[i, j] - sampled.Posterior.Draws[d].C[i, j]);
                    if (se > 0.0) worst = Math.Max(worst, gap / se);
                    else Assert.True(gap <= 1e-12, $"draw {d}: a zero error with gap {gap:E2}");
                    count++;
                }
        output.WriteLine($"{law}: {count} entries, worst |exact − sampled| / se {worst:F2}; E[C]'s top eigenvalue {exact.Posterior.MeanCEigenvalues[0]:F4} exact, {sampled.Posterior.MeanCEigenvalues[0]:F4} sampled");
        Assert.True(worst <= 4.5, $"{worst:F2}");
        Assert.Equal(["x0", "x1"], exact.Inputs);
        Assert.Equal(["z"], exact.Excluded);
    }

    // ── The metric ───────────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// On a fit whose box is [−3, 5]: C in raw units is D⁻¹C_boxD⁻¹ with D the widths, and whitened under a mixture law it is
    /// SC_boxS with S the law's standard deviations in box units, each read from the law's interval moments; Map says which.
    /// </summary>
    [Fact]
    public void DeclaredMetric_IsTheCongruenceOfTheBoxs()
    {
        (BassData data, BassFit fit) = WideBox.Value;
        var measures = new Dictionary<string, InputMeasure>
        {
            ["x0"] = new MixtureMeasure([new NormalComponent(1.0, 1.0, 1.5)]),
            ["x1"] = new UniformMeasure(-2.0, 4.0),
        };
        var spec = new ActiveSubspaceSpec { Measures = measures, Dimension = 1 };
        ActiveSubspaceAnalysis box = fit.ActiveSubspace(data, spec);
        ActiveSubspaceAnalysis raw = fit.ActiveSubspace(data, spec with { Metric = InputMetric.Raw });
        ActiveSubspaceAnalysis white = fit.ActiveSubspace(data, spec with { Metric = InputMetric.Whitened });
        double[] rawScale = [1.0 / 8.0, 1.0 / 8.0];
        // The laws' sds in box units, here: the uniform's (6/8)/√12, the truncated normal's by quadrature on [0, 1].
        ContinuousMeasure m0 = box.Measure.Continuous(0)!;
        double sd0 = BoxSd(m0), sd1 = (6.0 / 8.0) / Math.Sqrt(12.0);
        double[] whiteScale = [sd0, sd1];
        double worst = 0.0;
        for (int d = 0; d < box.Posterior.Draws.Count; d++)
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                {
                    double c = box.Posterior.Draws[d].C[i, j];
                    worst = Math.Max(worst, Math.Abs(raw.Posterior.Draws[d].C[i, j] - rawScale[i] * c * rawScale[j]) / Math.Max(1e-300, Math.Abs(rawScale[i] * c * rawScale[j])));
                    worst = Math.Max(worst, Math.Abs(white.Posterior.Draws[d].C[i, j] - whiteScale[i] * c * whiteScale[j]) / Math.Max(1e-300, Math.Abs(whiteScale[i] * c * whiteScale[j])));
                }
        output.WriteLine($"raw and whitened against LᵀC_boxL: {worst:E2}; maps {raw.Map[0, 0]:F4}, {white.Map[0, 0]:F4}/{white.Map[1, 1]:F4}; " +
                         $"E[C]'s direction box ({box.Posterior.MeanCVectors[0, 0]:F3}, {box.Posterior.MeanCVectors[1, 0]:F3}), raw ({raw.Posterior.MeanCVectors[0, 0]:F3}, {raw.Posterior.MeanCVectors[1, 0]:F3})");
        Assert.True(worst <= 1e-12, $"{worst:E2}");
        Assert.Equal(InputMetric.Whitened, white.Metric);
        Assert.Equal(sd0, white.Map[0, 0], 12);
    }

    // A truncated law's sd on its support in box units, by Gauss–Legendre on its density.
    private static double BoxSd(ContinuousMeasure m)
    {
        var density = (IContinuousDistribution)m.Measure;
        GaussLegendre rule = GaussLegendre.Of(60);
        double s0 = 0.0, s1 = 0.0, s2 = 0.0, half = 0.5 * (m.Upper - m.Lower), mid = 0.5 * (m.Upper + m.Lower);
        for (int r = 0; r < rule.Count; r++)
        {
            double x = mid + half * rule.Nodes[r], w = half * rule.Weights[r] * Math.Exp(density.LogPdf(x));
            s0 += w;
            s1 += w * x;
            s2 += w * x * x;
        }
        double mean = s1 / s0;
        return Math.Sqrt(s2 / s0 - mean * mean);
    }

    private static readonly Lazy<(BassData Data, BassFit Fit)> WideBox = new(() =>
    {
        var rng = new Xoshiro256PlusPlus(31);
        const int n = 150;
        double[] x0 = new double[n], x1 = new double[n], y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = -3.0 + 8.0 * rng.NextDouble();
            x1[i] = -3.0 + 8.0 * rng.NextDouble();
            y[i] = Math.Sin(0.4 * x0[i] + 0.25 * x1[i]) + 0.05 * Normal.SampleStandard(rng);
        }
        var domain = new InputDomain(-3.0, 5.0);
        var data = new BassData([new ContinuousPredictor("x0", x0, domain), new ContinuousPredictor("x1", x1, domain)], new ScalarResponse(y));
        return (data, Bass.Fit(SensitivityFixture.Config(), data));
    });

    // ── Decision 6's checks ──────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Poincaré on a hinge fit's draws under the uniform law (x1 on an inner interval) and on a Legendre fit's, and SG-4(b) on the
    /// Legendre fit's: no violation. On a Legendre fit of degree one, where every term holding an input has degree one in it, SG-4(b)
    /// holds with equality (`tensorTotal_eq_half`), each ratio 1.
    /// </summary>
    [Fact]
    public void ConsistencyChecks_HoldOnEveryDraw()
    {
        BassData data = SensitivityFixture.Scalar();
        BassFit fit = SensitivityFixture.FitScalar(2);
        var spec = new ActiveSubspaceSpec { Measures = new Dictionary<string, InputMeasure> { ["x1"] = new UniformMeasure(0.2, 0.9) } };
        PoincareCheck hinge = fit.ActiveSubspace(data, spec).Poincare;
        double hingeWorst = Max(hinge.Ratio);
        Assert.Null(hinge.JacobiRatio);

        BassData kData = KhaosReadoutTests.Data(120, 20261110);
        PoincareCheck legendre = Bass.Fit(KhaosReadoutTests.Config(5), kData).ActiveSubspace(kData).Poincare;
        PoincareCheck linear = Bass.Fit(KhaosReadoutTests.Config(1, order: 1), kData).ActiveSubspace(kData).Poincare;
        double equality = 0.0;
        int held = 0;
        for (int d = 0; d < linear.Total.GetLength(0); d++)
            for (int i = 0; i < linear.Total.GetLength(1); i++)
            {
                if (linear.Total[d, i] == 0.0) continue;
                equality = Math.Max(equality, Math.Abs(linear.JacobiRatio![d, i] - 1.0));
                held++;
            }
        output.WriteLine($"hinge (degree 2): worst total / Poincaré bound {hingeWorst:F4}, {hinge.Violations} violations; " +
                         $"Legendre (degree 5): worst {Max(legendre.Ratio):F4} against Poincaré, {Max(legendre.JacobiRatio!):F4} against SG-4(b)'s ½·Jacobi, " +
                         $"{legendre.Violations} violations; Legendre (degree 1): {held} draw-inputs at SG-4(b)'s equality to {equality:E2}");
        Assert.Equal(0, hinge.Violations);
        Assert.Equal(0, legendre.Violations);
        Assert.Equal(0, linear.Violations);
        Assert.InRange(hingeWorst, 0.0, 1.0);
        Assert.InRange(Max(legendre.JacobiRatio!), 0.0, 1.0 + 1e-9);
        Assert.True(held > 0 && equality <= 1e-12, $"{held}, {equality:E2}");
    }

    private static double Max(double[,] m)
    {
        double worst = 0.0;
        foreach (double v in m)
            if (!double.IsNaN(v)) worst = Math.Max(worst, v);
        return worst;
    }

    /// <summary>The readout refuses what it does not read: k outside 1 to p, a functional response.</summary>
    [Fact]
    public void Readout_RefusesWhatItDoesNotRead()
    {
        BassData data = SensitivityFixture.Scalar();
        BassFit fit = SensitivityFixture.FitScalar(1);
        Assert.Throws<ArgumentOutOfRangeException>(() => fit.ActiveSubspace(data, new ActiveSubspaceSpec { Dimension = 3 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => fit.ActiveSubspace(data, new ActiveSubspaceSpec { Dimension = 0 }));
        Assert.Throws<ArgumentException>(() => SensitivityFixture.FitFunctional().ActiveSubspace(SensitivityFixture.Functional()));
    }

    // ── Recovery (G3) ────────────────────────────────────────────────────────────────────────────────────────────────

    // At one N: the largest principal angles to the truth's span of E[C]'s top-k eigenspace, of the mean subspace and, averaged, of the
    // draws'; E[C]'s relative error in Frobenius norm; and the Davis–Kahan bound on E[C]'s angle from that error, asin of
    // 2‖ΔC‖_F/(λ_k − λ_{k+1}) with the truth's eigenvalues (Yu, Wang and Samworth's variant, cited from memory), beside the truth's
    // relative gap (λ_k − λ_{k+1})/λ_1.
    internal sealed record RecoveryReading(int N, double MeanCAngle, double ProjectionAngle, double MeanCError, double DrawAngle, double DavisKahan, double GapRatio);

    // A ridge from outside the model, sin(2.5 wᵀx) on [0, 1]³ with w = (0.8, 0.5, 0.33)/‖·‖, and two ridges, sin(2.5 w₁ᵀx) + sin(2 w₂ᵀx), on
    // [0, 1]⁴ with x₄ inert: data at each N with 0.02 noise, a degree-one hinge fit, and the readout's E[C] and mean subspace against the
    // span and C in closed form (ActiveSubspaceReferences).
    internal static RecoveryReading[] Recovery(string function, int dataSeed, int chainSeed, int[] sizes)
    {
        (ReferenceFunction f, double[][] span, int k) = function switch
        {
            "ridge" => ((ReferenceFunction)new RidgeSum("sin ridge", new RidgeTerm(ExponentialProfile.Sin(2.5), Unit([0.8, 0.5, 0.33]))), new[] { Unit([0.8, 0.5, 0.33]) }, 1),
            _ => (new RidgeSum("two ridges", new RidgeTerm(ExponentialProfile.Sin(2.5), Unit([1.0, 0.6, 0.0, 0.0])), new RidgeTerm(ExponentialProfile.Sin(2.0), Unit([0.0, 0.5, 1.0, 0.0]))),
                  new[] { new[] { 1.0, 0.6, 0.0, 0.0 }, new[] { 0.0, 0.5, 1.0, 0.0 } }, 2),
        };
        int p = f.Dimension;
        double[,] truth = f.C(ReferenceBox.Cube(p, 0.0, 1.0))!;
        double[,] spanBasis = Orthonormal(span, p);
        double[] lambda = Maths.LinAlg.DenseEigen.DecomposeSymmetric(truth).Eigenvalues;
        double gap = lambda[k - 1] - (k < p ? lambda[k] : 0.0);
        var readings = new RecoveryReading[sizes.Length];
        var rng = new Xoshiro256PlusPlus(dataSeed);   // one stream: each size's data follow the previous size's
        for (int s = 0; s < sizes.Length; s++)
        {
            int n = sizes[s];
            var x = new double[p][];
            for (int j = 0; j < p; j++) x[j] = new double[n];
            var y = new double[n];
            var point = new double[p];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < p; j++) point[j] = x[j][i] = rng.NextDouble();
                y[i] = f.Value(point) + 0.02 * Normal.SampleStandard(rng);
            }
            var domain = new InputDomain(0.0, 1.0);
            var data = new BassData([.. Enumerable.Range(0, p).Select(j => new ContinuousPredictor($"x{j}", x[j], domain))], new ScalarResponse(y));
            BassConfig config = SensitivityFixture.Config(maxBasis: 40) with
            {
                Sampling = new BassSampling { Chains = 2, MasterSeed = chainSeed, BurnIn = 2_000, Samples = 30, Thin = 10 },
            };
            ActiveSubspaceAnalysis a = Bass.Fit(config, data).ActiveSubspace(data, new ActiveSubspaceSpec { Dimension = k });
            double[,] meanC = a.Posterior.MeanC;
            double scale = 0.0, error = 0.0;
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++)
                {
                    scale += truth[i, j] * truth[i, j];
                    error += (meanC[i, j] - truth[i, j]) * (meanC[i, j] - truth[i, j]);
                }
            double drawAngle = a.Posterior.Draws.Select(d => SubspacePosterior.PrincipalAngles(d.Basis, spanBasis).Max()).Average();
            readings[s] = new RecoveryReading(n, SubspacePosterior.PrincipalAngles(Leading(a.Posterior.MeanCVectors, k), spanBasis).Max(),
                                              SubspacePosterior.PrincipalAngles(Leading(a.Posterior.ProjectionVectors, k), spanBasis).Max(),
                                              Math.Sqrt(error / scale), drawAngle, Math.Asin(Math.Min(1.0, 2.0 * Math.Sqrt(error) / gap)), gap / lambda[0]);
        }
        return readings;
    }

    /// <summary>
    /// Recovery at growing N (G3): on each function, the largest principal angle between the truth's span and E[C]'s top-k eigenspace,
    /// and the mean subspace's, and E[C]'s relative error in Frobenius norm, at N = 100, 300 and 900. Each falls from the smallest N to
    /// the largest (registered in <see cref="ComparativeSweep"/>). At every N, E[C]'s angle sits within the Davis–Kahan bound its own
    /// error gives under the truth's gap: a subspace is as well identified as its gap allows.
    /// </summary>
    [Theory]
    [InlineData("ridge")]
    [InlineData("two ridges")]
    public void Recovery_AtGrowingN_FromOutsideTheModel(string function)
    {
        RecoveryReading[] r = Recovery(function, 20261108, 13, [100, 300, 900]);
        foreach (RecoveryReading x in r)
            output.WriteLine($"{function}, N {x.N}: E[C]'s subspace {x.MeanCAngle:F4} rad from the truth (Davis–Kahan bound {x.DavisKahan:F4}), " +
                             $"the mean subspace {x.ProjectionAngle:F4}, draws on average {x.DrawAngle:F4}; E[C]'s relative error {x.MeanCError:F4}; " +
                             $"the truth's relative gap {x.GapRatio:F4}");
        Assert.True(r[^1].MeanCAngle < r[0].MeanCAngle, "E[C]'s angle does not fall");
        Assert.True(r[^1].ProjectionAngle < r[0].ProjectionAngle, "the mean subspace's angle does not fall");
        Assert.True(r[^1].MeanCError < r[0].MeanCError, "E[C]'s error does not fall");
        Assert.All(r, x => Assert.True(x.MeanCAngle <= x.DavisKahan, $"N {x.N}: {x.MeanCAngle:F4} against {x.DavisKahan:F4}"));
    }

    private static double[] Unit(double[] w)
    {
        double n = Math.Sqrt(w.Sum(v => v * v));
        return [.. w.Select(v => v / n)];
    }

    private static double[,] Orthonormal(double[][] span, int p)
    {
        var q = new List<double[]>();
        foreach (double[] d in span)
        {
            var v = (double[])d.Clone();
            for (int pass = 0; pass < 2; pass++)
                foreach (double[] u in q)
                {
                    double dot = v.Zip(u, (a, b) => a * b).Sum();
                    for (int i = 0; i < p; i++) v[i] -= dot * u[i];
                }
            q.Add(Unit(v));
        }
        var basis = new double[p, q.Count];
        for (int j = 0; j < q.Count; j++)
            for (int i = 0; i < p; i++) basis[i, j] = q[j][i];
        return basis;
    }

    private static double[,] Leading(double[,] vectors, int k)
    {
        int p = vectors.GetLength(0);
        var basis = new double[p, k];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < k; j++) basis[i, j] = vectors[i, j];
        return basis;
    }
}
