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
/// E8's slice 3 gates (the active-subspace brief) on co-active pairs, C_fg = E[∇f∇gᵀ] (Rumsey, Hardy, Ahrens and Vander Wiel, arXiv
/// 2311.18146), and their posterior (<see cref="CoActivePosterior"/>):
/// <list type="bullet">
/// <item>2311.18146's simple polynomial example in closed form, f₁ = x₁² + x₁x₂ and f₂ = f₁ + βx₂³: C₁, C₂, C₁₂ against the paper's
/// matrices and against quadrature, the concordance, V's directions ordered by |λ|, the contributions and the co-activity scores; and
/// the paper's quoted readings, three of which its own definitions do not give (each pinned to what it does read);</item>
/// <item>the same example from BASS fits at the paper's n = 200 and 1,000, E[C₁₂] against the closed form;</item>
/// <item>two fits on the same inputs: swapping them transposes C_fg, a fit paired with itself gives its own C and concordance 1, C_fg
/// against quadrature of the two gradient readouts, the exact and sampled routes within Monte Carlo error;</item>
/// <item>two outputs of one fit, direct and basis, and a functional fit paired with a scalar one: the same reductions, against
/// quadrature of the output functionals' gradients;</item>
/// <item>the refusals.</item>
/// </list>
/// </summary>
public sealed class CoActiveSubspaceTests(ITestOutputHelper output)
{
    // ── 2311.18146's simple polynomial example (§3.1.1) ──────────────────────────────────────────────────────────────

    private static readonly ReferenceBox UnitSquare = ReferenceBox.Cube(2, 0.0, 1.0);

    private static PolynomialReference F1() => new("f1", [new Monomial(1.0, [2, 0]), new Monomial(1.0, [1, 1])]);

    private static PolynomialReference F2(double beta) => new("f2", [new Monomial(1.0, [2, 0]), new Monomial(1.0, [1, 1]), new Monomial(beta, [0, 3])]);

    // The paper's eq. (simple_results1): a 2 × 2 matrix over 180.
    private static double[,] Over180(double a, double b, double c, double d) => new[,] { { a / 180.0, b / 180.0 }, { c / 180.0, d / 180.0 } };

    private static CoActivePosterior SimpleExample(double beta)
    {
        PolynomialReference f1 = F1(), f2 = F2(beta);
        return CoActivePosterior.Of([new CoActivePair(f1.C(UnitSquare)!, f2.C(UnitSquare)!, f1.Cross(f2, UnitSquare)!)], [1.0]);
    }

    // conc(f₁, f₂) = (3 + β/2)/√(9 + 3β(9β/5 + 1)), the paper's closed form.
    private static double Concordance(double beta) => (3.0 + beta / 2.0) / Math.Sqrt(9.0 + 3.0 * beta * (1.8 * beta + 1.0));

    /// <summary>
    /// C₁, C₂ and C₁₂ by the polynomial closed form against the paper's matrices and against tensor quadrature; C₁₂ at β = 0 is C₁; the
    /// traces t₁ = 3, t₂ = 3 + β(9β/5 + 1), t₁₂ = 3 + β/2; the concordance, the contributions summing to it, V symmetric with its
    /// eigenvalues ordered by magnitude, and the co-activity scores: at q = 1 signed and unsigned agree up to sign, and at q = p the signed
    /// scores sum to t₁₂.
    /// </summary>
    [Theory]
    [InlineData(0.5)]
    [InlineData(3.0)]
    [InlineData(-12.0)]
    [InlineData(0.0)]
    public void SimpleExample_ClosedForms(double beta)
    {
        PolynomialReference f1 = F1(), f2 = F2(beta);
        double[,] c1 = f1.C(UnitSquare)!, c2 = f2.C(UnitSquare)!, c12 = f1.Cross(f2, UnitSquare)!;
        AssertClose(Over180(480, 165, 165, 60), c1, 1e-14);
        AssertClose(Over180(480, 165 + 315 * beta, 165 + 315 * beta, 60 + beta * (324 * beta + 180)), c2, 1e-14);
        AssertClose(Over180(480, 165 + 315 * beta, 165, 60 + 90 * beta), c12, 1e-14);
        QuadratureRule[] rules = [QuadratureRule.Uniform(0.0, 1.0, 4), QuadratureRule.Uniform(0.0, 1.0, 4)];
        AssertClose(c12, GradientQuadrature.Cross(f1.Gradient, f2.Gradient, rules), 1e-14);
        AssertClose(c1, f1.Cross(f1, UnitSquare)!, 1e-15);

        CoActivePosterior post = SimpleExample(beta);
        (double t1, double t2, double t12) = post.MeanTraces;
        Assert.Equal(3.0, t1, 13);
        Assert.Equal(3.0 + beta * (1.8 * beta + 1.0), t2, 12);
        Assert.Equal(3.0 + beta / 2.0, t12, 13);
        Assert.Equal(Concordance(beta), post.Concordance, 13);
        Assert.Equal(post.Concordance, post.Contributions.Sum(), 13);
        Assert.Equal(post.V[0, 1], post.V[1, 0]);
        Assert.True(Math.Abs(post.Values[0]) >= Math.Abs(post.Values[1]));
        double[] one = post.CoActivityScores(1), oneUnsigned = post.CoActivityScores(1, signed: false);
        for (int i = 0; i < 2; i++) Assert.Equal(Math.Abs(one[i]), oneUnsigned[i], 14);
        Assert.Equal(t12, post.CoActivityScores(2).Sum(), 12);
        if (beta == 0.0) Assert.Equal(1.0, post.Concordance, 14);
        output.WriteLine($"β = {beta}: conc {post.Concordance:F5}, λ ({post.Values[0]:F4}, {post.Values[1]:F4}), " +
                         $"w₁ ({post.Directions[0, 0]:F4}, {post.Directions[1, 0]:F4}), π ({post.Contributions[0]:F4}, {post.Contributions[1]:F4})");
    }

    /// <summary>
    /// The paper's quoted readings (2311.18146 v2, §3.1.1), against its own definitions. Its concordances 0.944 at β = ½ and 0.551 at
    /// β = 3 hold. Three do not, and each is pinned to what it reads: −0.131 "for β = −12" is the concordance at β = −15 (at −12 it is
    /// −0.1095); w₁⁽¹²⁾ ≈ (0.907, 0.422) at β = ½ is V's leading eigenvector at β = 0.4 (at ½ it is (0.8971, 0.4418)); and
    /// π⁽¹²⁾ = (0.9518, −0.0077) at β = ½ are the unsymmetrized C₁₂'s eigenvalues over √(t₁t₂), where V's, which its eq. (contrib)
    /// defines, give (0.9684, −0.0242). Both pairs sum to the concordance, since tr V = tr C₁₂. Its single-function directions
    /// w₁⁽¹⁾ ≈ (0.945, 0.327) and w₁⁽²⁾ ≈ (0.825, 0.566) hold. At β = −12 the leading co-active value is negative: ordering by
    /// magnitude matters.
    /// </summary>
    [Fact]
    public void SimpleExample_ThePapersReadings()
    {
        Assert.Equal(0.944, SimpleExample(0.5).Concordance, 3);
        Assert.Equal(0.551, SimpleExample(3.0).Concordance, 3);
        Assert.Equal(-0.1095, SimpleExample(-12.0).Concordance, 4);
        Assert.Equal(-0.131, SimpleExample(-15.0).Concordance, 3);

        CoActivePosterior half = SimpleExample(0.5);
        double sign = Math.Sign(half.Directions[0, 0]);
        Assert.Equal(0.8971, sign * half.Directions[0, 0], 4);
        Assert.Equal(0.4418, sign * half.Directions[1, 0], 4);
        Assert.Equal(0.9684, half.Contributions[0], 4);
        Assert.Equal(-0.0242, half.Contributions[1], 4);

        CoActivePosterior near = SimpleExample(0.4);
        double s4 = Math.Sign(near.Directions[0, 0]);
        Assert.Equal(0.907, s4 * near.Directions[0, 0], 3);
        Assert.Equal(0.422, s4 * near.Directions[1, 0], 3);

        // The unsymmetrized C₁₂'s eigenvalues, by the 2 × 2 closed form, over √(t₁t₂).
        double[,] c12 = half.MeanCfg;
        double tr = c12[0, 0] + c12[1, 1], det = c12[0, 0] * c12[1, 1] - c12[0, 1] * c12[1, 0], root = Math.Sqrt(tr * tr - 4.0 * det);
        (double t1, double t2, _) = half.MeanTraces;
        Assert.Equal(0.9518, (tr + root) / 2.0 / Math.Sqrt(t1 * t2), 4);
        Assert.Equal(-0.0077, (tr - root) / 2.0 / Math.Sqrt(t1 * t2), 4);

        Assert.Equal(0.945, Leading(half.MeanCf)[0], 3);
        Assert.Equal(0.327, Leading(half.MeanCf)[1], 3);
        Assert.Equal(0.825, Leading(half.MeanCg)[0], 3);
        Assert.Equal(0.566, Leading(half.MeanCg)[1], 3);

        CoActivePosterior negative = SimpleExample(-12.0);
        Assert.True(negative.Values[0] < 0.0 && Math.Abs(negative.Values[0]) > Math.Abs(negative.Values[1]), $"{negative.Values[0]}, {negative.Values[1]}");
    }

    /// <summary>
    /// 2311.18146's estimate (§3.2): BASS fits of f₁ and f₂ at β = 3 from n = 200 and 1,000 Latin hypercube points, E[C₁₂] in Frobenius
    /// distance from the closed form. The paper reports 0.050 and 0.0086 from one fit each through <c>concordance</c>'s <c>Cfg_bass</c>;
    /// a fit's error depends on its draws, so the gate is ours: the distance falls with n (registered in <see cref="ComparativeSweep"/>),
    /// and at n = 1,000 it is under 2% of ‖C₁₂‖_F.
    /// </summary>
    [Fact]
    public void SimpleExample_FromBassFits()
    {
        (double[] distances, double norm) = SimpleExampleFromFits(20231130, 23, output);
        Assert.True(distances[1] < distances[0], $"{distances[0]:F4} → {distances[1]:F4}");
        Assert.True(distances[1] < 0.02 * norm, $"{distances[1]:F4} against ‖C₁₂‖_F {norm:F4}");
    }

    // At n = 200 and then 1,000 Latin hypercube points from one stream: E[C₁₂]'s Frobenius distance from the closed form, and ‖C₁₂‖_F.
    internal static (double[] Distances, double Norm) SimpleExampleFromFits(int dataSeed, int chainSeed, ITestOutputHelper? output = null)
    {
        PolynomialReference f1 = F1(), f2 = F2(3.0);
        double[,] truth = f1.Cross(f2, UnitSquare)!;
        double norm = Frobenius(truth, new double[2, 2]);
        var rng = new Xoshiro256PlusPlus(dataSeed);
        var distances = new List<double>();
        foreach (int n in new[] { 200, 1000 })
        {
            double[][] x = [LatinColumn(rng, n), LatinColumn(rng, n)];
            var point = new double[2];
            var y1 = new double[n];
            var y2 = new double[n];
            for (int i = 0; i < n; i++)
            {
                point[0] = x[0][i];
                point[1] = x[1][i];
                y1[i] = f1.Value(point);
                y2[i] = f2.Value(point);
            }
            var domain = new InputDomain(0.0, 1.0);
            BassPredictor[] predictors = [new ContinuousPredictor("x1", x[0], domain), new ContinuousPredictor("x2", x[1], domain)];
            BassData d1 = new(predictors, new ScalarResponse(y1)), d2 = new(predictors, new ScalarResponse(y2));
            BassConfig config = SensitivityFixture.Config(maxBasis: 40) with
            {
                Sampling = new BassSampling { Chains = 2, MasterSeed = chainSeed, BurnIn = 2_000, Samples = 30, Thin = 10 },
            };
            CoActiveAnalysis a = Bass.Fit(config, d1).CoActiveSubspace(d1, Bass.Fit(config, d2), d2);
            double distance = Frobenius(a.Posterior.MeanCfg, truth);
            distances.Add(distance);
            output?.WriteLine($"n {n}: E[C₁₂] = [[{a.Posterior.MeanCfg[0, 0]:F3}, {a.Posterior.MeanCfg[0, 1]:F3}], [{a.Posterior.MeanCfg[1, 0]:F3}, {a.Posterior.MeanCfg[1, 1]:F3}]], " +
                             $"Frobenius distance {distance:F4} ({distance / norm:P2} of ‖C₁₂‖_F; the paper's single fits: {(n == 200 ? 0.050 : 0.0086)}); " +
                             $"concordance of the means {a.Posterior.Concordance:F4} (truth {Concordance(3.0):F4}), draws' concordance " +
                             $"{a.Posterior.Concordances.Min():F4} to {a.Posterior.Concordances.Max():F4}");
        }
        return ([.. distances], norm);
    }

    // One Latin hypercube column on [0, 1]: a point uniform in each of n strata, the strata in random order.
    private static double[] LatinColumn(Xoshiro256PlusPlus rng, int n)
    {
        int[] order = [.. Enumerable.Range(0, n)];
        for (int i = n - 1; i > 0; i--)
        {
            int j = rng.NextInt(i + 1);
            (order[i], order[j]) = (order[j], order[i]);
        }
        return [.. order.Select(k => (k + rng.NextDouble()) / n)];
    }

    // ── Two fits on the same inputs ──────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Swapping the fits transposes each pair's C_fg and swaps C_f and C_g; a fit paired with itself, draw for draw, gives C_fg = C_f,
    /// its own C from the active-subspace readout, and concordance 1; the posterior's weights are the normalized products.
    /// </summary>
    [Fact]
    public void TwoFits_SwapTransposes_AndASelfPairIsC()
    {
        (BassData df, BassData dg) = CoActiveFixture.Data();
        (BassFit f, BassFit g) = CoActiveFixture.Fits();
        int[] sel = Spread(f.Draws.Sum(c => c.Count), 8), other = [.. sel.Reverse()];
        CoActiveAnalysis fg = f.CoActiveSubspace(df, g, dg, new CoActiveSpec { Draws = sel, OtherDraws = other });
        CoActiveAnalysis gf = g.CoActiveSubspace(dg, f, df, new CoActiveSpec { Draws = other, OtherDraws = sel });
        double worst = 0.0;
        for (int s = 0; s < sel.Length; s++)
        {
            CoActivePair a = fg.Posterior.Pairs[s], b = gf.Posterior.Pairs[s];
            Assert.Equal(a.Cf, b.Cg);
            Assert.Equal(a.Cg, b.Cf);
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(b.Cfg, Transpose(a.Cfg)));
        }
        Assert.True(worst <= 1e-13, $"swap: {worst:E2}");
        Assert.Equal(fg.Posterior.Concordance, gf.Posterior.Concordance, 13);
        Assert.Equal(sel, fg.DrawIndex);
        Assert.Equal(other, fg.OtherDrawIndex!);

        CoActiveAnalysis ff = f.CoActiveSubspace(df, f, df, new CoActiveSpec { Draws = sel, OtherDraws = sel });
        ActiveSubspaceAnalysis af = f.ActiveSubspace(df, new ActiveSubspaceSpec { Draws = sel });
        double self = 0.0;
        for (int s = 0; s < sel.Length; s++)
        {
            CoActivePair p = ff.Posterior.Pairs[s];
            Assert.Equal(af.Posterior.Draws[s].C, p.Cf);
            self = Math.Max(self, FunctionalActiveSubspaceTests.RelativeGap(p.Cfg, p.Cf));
            Assert.Equal(1.0, p.Concordance, 13);
        }
        output.WriteLine($"swap {worst:E2}; self-pair C_fg against C_f {self:E2}; concordance of the means f, g {fg.Posterior.Concordance:F4}");
        Assert.True(self <= 1e-13, $"self: {self:E2}");
        Assert.All(ff.Posterior.Weights, w => Assert.Equal(1.0 / sel.Length, w, 15));
    }

    /// <summary>Each pair's C_f, C_g and C_fg against tensor quadrature of the two fits' gradient readouts, split at the knots.</summary>
    [Fact]
    public void TwoFits_C_IsTheQuadratureOfTheirGradients()
    {
        (BassData df, BassData dg) = CoActiveFixture.Data();
        (BassFit f, BassFit g) = CoActiveFixture.Fits();
        int[] sel = Spread(f.Draws.Sum(c => c.Count), 6), other = [.. sel.Select(d => (d + 7) % f.Draws.Sum(c => c.Count))];
        CoActiveAnalysis a = f.CoActiveSubspace(df, g, dg, new CoActiveSpec { Draws = sel, OtherDraws = other });
        (PredictionCases cases, double[] w) = FunctionalActiveSubspaceTests.Tensor(BassInputs.Prepare(df, f.Transform), a.Measure);
        BassGradient gf = f.Gradient(df, cases, sel), gg = g.Gradient(dg, cases, other);
        double worst = 0.0;
        for (int s = 0; s < sel.Length; s++)
        {
            int ss = s;
            CoActivePair p = a.Posterior.Pairs[s];
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(p.Cfg, Cross(w, (r, i) => gf.Partials[i][ss, r, 0], (r, i) => gg.Partials[i][ss, r, 0])));
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(p.Cf, Cross(w, (r, i) => gf.Partials[i][ss, r, 0], (r, i) => gf.Partials[i][ss, r, 0])));
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(p.Cg, Cross(w, (r, i) => gg.Partials[i][ss, r, 0], (r, i) => gg.Partials[i][ss, r, 0])));
        }
        output.WriteLine($"worst gap to quadrature {worst:E2}");
        Assert.True(worst <= 1e-10, $"{worst:E2}");
    }

    /// <summary>The exact and sampled routes for two fits, each entry of C_f, C_g and C_fg within 4.5 of the sampled route's standard errors.</summary>
    [Fact]
    public void TwoFits_ExactAndSampledRoutes_AgreeWithinMonteCarloError()
    {
        (BassData df, BassData dg) = CoActiveFixture.Data();
        (BassFit f, BassFit g) = CoActiveFixture.Fits();
        int[] sel = Spread(f.Draws.Sum(c => c.Count), 10);
        var spec = new CoActiveSpec { Draws = sel, OtherDraws = sel, Points = 8192, Seed = 3 };
        CoActiveAnalysis exact = f.CoActiveSubspace(df, g, dg, spec with { Route = ActiveSubspaceRoute.Exact });
        CoActiveAnalysis sampled = f.CoActiveSubspace(df, g, dg, spec with { Route = ActiveSubspaceRoute.Sampled });
        Assert.Null(exact.MonteCarloError);
        double worst = WorstInErrors(exact, sampled);
        output.WriteLine($"worst |exact − sampled| / se {worst:F2} over {sel.Length} pairs, three matrices each");
        Assert.True(worst <= 4.5, $"{worst:F2}");
        Assert.Equal(8192, sampled.Points);
    }

    // ── Output functionals ───────────────────────────────────────────────────────────────────────────────────────────

    private static double[] OutputA(int m) => [.. Enumerable.Range(0, m).Select(l => 0.2 + 0.3 * Math.Sin(1.3 * l))];

    private static double[] OutputB(int m) => [.. Enumerable.Range(0, m).Select(l => l % 3 == 0 ? 1.0 : -0.4 + 0.05 * l)];

    /// <summary>
    /// Two outputs aᵀy and bᵀy of a direct functional fit: swapping them transposes C_fg; a = b gives C_fg = C_f, the readout's C under
    /// W = aaᵀ, and concordance 1; C_f, C_g and C_fg against quadrature of aᵀJ and bᵀJ by the gradient readout; the two-fit route with
    /// the fit paired with itself draw for draw gives the same pairs; and the exact and sampled routes within Monte Carlo error.
    /// </summary>
    [Fact]
    public void TwoOutputs_OfADirectFit()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        int m = data.Response.Locations;
        double[] a = OutputA(m), b = OutputB(m);
        int[] sel = Spread(fit.Draws.Sum(c => c.Count), 6);
        var spec = new CoActiveSpec { Draws = sel, Output = a, OtherOutput = b };
        CoActiveAnalysis ab = fit.CoActiveSubspace(data, spec);
        CoActiveAnalysis ba = fit.CoActiveSubspace(data, spec with { Output = b, OtherOutput = a });
        CoActiveAnalysis aa = fit.CoActiveSubspace(data, spec with { OtherOutput = a });
        ActiveSubspaceAnalysis waa = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Draws = sel, OutputMetric = new OutputMetric(Outer(a, a)) });
        CoActiveAnalysis twoFit = fit.CoActiveSubspace(data, fit, data, spec with { OtherDraws = sel });
        (PredictionCases cases, double[] w) = FunctionalActiveSubspaceTests.Tensor(BassInputs.Prepare(data, fit.Transform), ab.Measure);
        BassGradient g = fit.Gradient(data, cases, sel);
        double swap = 0.0, reduce = 0.0, quadrature = 0.0;
        for (int s = 0; s < sel.Length; s++)
        {
            int ss = s;
            swap = Math.Max(swap, FunctionalActiveSubspaceTests.RelativeGap(ba.Posterior.Pairs[s].Cfg, Transpose(ab.Posterior.Pairs[s].Cfg)));
            reduce = Math.Max(reduce, FunctionalActiveSubspaceTests.RelativeGap(aa.Posterior.Pairs[s].Cfg, aa.Posterior.Pairs[s].Cf));
            reduce = Math.Max(reduce, FunctionalActiveSubspaceTests.RelativeGap(aa.Posterior.Pairs[s].Cf, waa.Posterior.Draws[s].C));
            Assert.Equal(1.0, aa.Posterior.Pairs[s].Concordance, 12);
            Func<int, int, double> fa = (r, i) => Along(a, g, ss, r, i), fb = (r, i) => Along(b, g, ss, r, i);
            quadrature = Math.Max(quadrature, FunctionalActiveSubspaceTests.RelativeGap(ab.Posterior.Pairs[s].Cfg, Cross(w, fa, fb)));
            quadrature = Math.Max(quadrature, FunctionalActiveSubspaceTests.RelativeGap(ab.Posterior.Pairs[s].Cf, Cross(w, fa, fa)));
            quadrature = Math.Max(quadrature, FunctionalActiveSubspaceTests.RelativeGap(ab.Posterior.Pairs[s].Cg, Cross(w, fb, fb)));
            Assert.Equal(ab.Posterior.Pairs[s].Cfg, twoFit.Posterior.Pairs[s].Cfg);
        }
        output.WriteLine($"swap {swap:E2}, a = b against C under aaᵀ {reduce:E2}, quadrature {quadrature:E2}; concordance of the means {ab.Posterior.Concordance:F4}");
        Assert.True(swap <= 1e-13, $"swap {swap:E2}");
        Assert.True(reduce <= 1e-12, $"reduce {reduce:E2}");
        Assert.True(quadrature <= 1e-10, $"quadrature {quadrature:E2}");
        Assert.Null(ab.OtherDrawIndex);

        var routes = spec with { Draws = Spread(fit.Draws.Sum(c => c.Count), 10), Points = 8192, Seed = 9 };
        double worst = WorstInErrors(fit.CoActiveSubspace(data, routes with { Route = ActiveSubspaceRoute.Exact }),
                                     fit.CoActiveSubspace(data, routes with { Route = ActiveSubspaceRoute.Sampled }));
        output.WriteLine($"routes: worst |exact − sampled| / se {worst:F2}");
        Assert.True(worst <= 4.5, $"{worst:F2}");
    }

    /// <summary>A functional fit read through an output functional, paired with a scalar fit on the same inputs, against quadrature.</summary>
    [Fact]
    public void FunctionalAndScalarFits_Pair()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        (BassData scalarData, BassFit scalar) = CoActiveFixture.OnFunctionalInputs();
        int m = data.Response.Locations;
        double[] a = OutputA(m);
        int[] sel = Spread(fit.Draws.Sum(c => c.Count), 6);
        CoActiveAnalysis pair = fit.CoActiveSubspace(data, scalar, scalarData, new CoActiveSpec { Draws = sel, OtherDraws = sel, Output = a });
        (PredictionCases cases, double[] w) = FunctionalActiveSubspaceTests.Tensor(BassInputs.Prepare(data, fit.Transform), pair.Measure);
        BassGradient gf = fit.Gradient(data, cases, sel), gs = scalar.Gradient(scalarData, cases, sel);
        double worst = 0.0;
        for (int s = 0; s < sel.Length; s++)
        {
            int ss = s;
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(pair.Posterior.Pairs[s].Cfg, Cross(w, (r, i) => Along(a, gf, ss, r, i), (r, i) => gs.Partials[i][ss, r, 0])));
        }
        output.WriteLine($"functional against scalar: worst gap to quadrature {worst:E2}; concordance of the means {pair.Posterior.Concordance:F4}");
        Assert.True(worst <= 1e-10, $"{worst:E2}");
    }

    /// <summary>
    /// Two outputs of a basis fit over its joint draws: a = b gives C_fg = C_f, the basis readout's C under W = aaᵀ; C_fg against quadrature
    /// of aᵀJ and bᵀJ with J the reconstructed Jacobian.
    /// </summary>
    [Fact]
    public void TwoOutputs_OfABasisFit()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        int m = fit.Decomposition.Locations;
        double[] a = OutputA(m), b = OutputB(m);
        int[] sel = Spread(fit.JointDraws.Count, 5);
        var spec = new CoActiveSpec { Draws = sel, Output = a, OtherOutput = b };
        CoActiveAnalysis ab = fit.CoActiveSubspace(spec);
        CoActiveAnalysis aa = fit.CoActiveSubspace(spec with { OtherOutput = a });
        ActiveSubspaceAnalysis waa = fit.ActiveSubspace(new ActiveSubspaceSpec { Draws = sel, OutputMetric = new OutputMetric(Outer(a, a)) });
        (PredictionCases cases, double[] w) = FunctionalActiveSubspaceTests.Tensor(BassInputs.Prepare(fit.ComponentData(0), fit.Components[0].Transform), ab.Measure);
        double reduce = 0.0, quadrature = 0.0;
        for (int s = 0; s < sel.Length; s++)
        {
            reduce = Math.Max(reduce, FunctionalActiveSubspaceTests.RelativeGap(aa.Posterior.Pairs[s].Cfg, aa.Posterior.Pairs[s].Cf));
            reduce = Math.Max(reduce, FunctionalActiveSubspaceTests.RelativeGap(aa.Posterior.Pairs[s].Cf, waa.Posterior.Draws[s].C));
            double[][,] j = FunctionalActiveSubspaceTests.BasisJacobians(fit, cases, sel[s]);
            quadrature = Math.Max(quadrature, FunctionalActiveSubspaceTests.RelativeGap(ab.Posterior.Pairs[s].Cfg,
                Cross(w, (r, i) => Contract(a, j[r], i), (r, i) => Contract(b, j[r], i))));
        }
        output.WriteLine($"basis: a = b against C under aaᵀ {reduce:E2}, quadrature {quadrature:E2}");
        Assert.True(reduce <= 1e-12, $"reduce {reduce:E2}");
        Assert.True(quadrature <= 1e-10, $"quadrature {quadrature:E2}");
        Assert.Null(ab.MonteCarloError);
    }

    // ── Refusals ─────────────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// The readout refuses what it does not read: fits on other inputs, or of other hinge degrees; selections of different lengths; an
    /// output functional on a scalar fit or of the wrong length; a functional fit without one; OtherDraws for two outputs of one fit;
    /// two outputs of a scalar fit; a basis fit's sampled route or a missing output.
    /// </summary>
    [Fact]
    public void CoActive_RefusesWhatItDoesNotRead()
    {
        (BassData df, BassData dg) = CoActiveFixture.Data();
        (BassFit f, BassFit g) = CoActiveFixture.Fits();
        BassData scalar = SensitivityFixture.Scalar();
        Assert.Throws<ArgumentException>(() => f.CoActiveSubspace(df, SensitivityFixture.FitScalar(1), scalar));
        Assert.Throws<ArgumentException>(() => SensitivityFixture.FitScalar(1).CoActiveSubspace(scalar, SensitivityFixture.FitScalar(2), scalar));
        Assert.Throws<ArgumentException>(() => f.CoActiveSubspace(df, g, dg, new CoActiveSpec { Draws = [0, 1], OtherDraws = [0] }));
        Assert.Throws<ArgumentException>(() => f.CoActiveSubspace(df, g, dg, new CoActiveSpec { Draws = [0], OtherDraws = [0], Output = [1.0] }));
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        int m = data.Response.Locations;
        Assert.Throws<ArgumentException>(() => fit.CoActiveSubspace(data, fit, data, new CoActiveSpec { Draws = [0], OtherDraws = [0], Output = OutputA(m) }));
        Assert.Throws<ArgumentException>(() => fit.CoActiveSubspace(data, new CoActiveSpec { Draws = [0], Output = OutputA(m) }));
        Assert.Throws<ArgumentException>(() => fit.CoActiveSubspace(data, new CoActiveSpec { Draws = [0], Output = OutputA(m + 1), OtherOutput = OutputB(m + 1) }));
        Assert.Throws<ArgumentException>(() => fit.CoActiveSubspace(data, new CoActiveSpec { Draws = [0], OtherDraws = [0], Output = OutputA(m), OtherOutput = OutputB(m) }));
        Assert.Throws<ArgumentException>(() => f.CoActiveSubspace(df, new CoActiveSpec { Draws = [0], Output = [1.0], OtherOutput = [1.0] }));
        BassBasisFit basis = BassOutputComponentsTests.FitOf();
        int mb = basis.Decomposition.Locations;
        Assert.Throws<ArgumentException>(() => basis.CoActiveSubspace(new CoActiveSpec { Draws = [0], Output = OutputA(mb) }));
        Assert.Throws<ArgumentException>(() => basis.CoActiveSubspace(new CoActiveSpec { Draws = [0], Output = OutputA(mb), OtherOutput = OutputB(mb), Route = ActiveSubspaceRoute.Sampled }));
    }

    /// <summary>The posterior's means under weights: E[C_fg] is the weighted mean, the weights normalized; V its symmetric part.</summary>
    [Fact]
    public void CoActivePosterior_WeightsThePairs()
    {
        double[,] i2 = { { 1.0, 0.0 }, { 0.0, 1.0 } };
        var first = new CoActivePair(i2, i2, new[,] { { 1.0, 2.0 }, { 0.0, 1.0 } });
        var second = new CoActivePair(i2, i2, new[,] { { -1.0, 0.0 }, { 4.0, 1.0 } });
        CoActivePosterior post = CoActivePosterior.Of([first, second], [1.0, 3.0]);
        Assert.Equal(new[] { 0.25, 0.75 }, post.Weights);
        AssertClose(new[,] { { -0.5, 0.5 }, { 3.0, 1.0 } }, post.MeanCfg, 1e-15);
        AssertClose(new[,] { { -0.5, 1.75 }, { 1.75, 1.0 } }, post.V, 1e-15);
        Assert.Equal(0.25, post.Concordance, 15);
        Assert.Equal(new[] { 1.0, 0.0 }, post.Concordances);
        Assert.Throws<ArgumentException>(() => CoActivePosterior.Of([first], [-1.0]));
        Assert.Throws<ArgumentException>(() => CoActivePosterior.Of([first, second], [1.0]));
    }

    // ── Helpers ──────────────────────────────────────────────────────────────────────────────────────────────────────

    private static int[] Spread(int count, int k) => [.. Enumerable.Range(0, k).Select(i => (int)((long)i * (count - 1) / (k - 1)))];

    // Σ_r w_r ∂_i f(x_r) ∂_j g(x_r) over the rule's points, two inputs.
    private static double[,] Cross(double[] w, Func<int, int, double> gradF, Func<int, int, double> gradG)
    {
        var c = new double[2, 2];
        for (int r = 0; r < w.Length; r++)
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++) c[i, j] += w[r] * gradF(r, i) * gradG(r, j);
        return c;
    }

    // aᵀ∂y/∂x_i at point r of draw s.
    private static double Along(double[] a, BassGradient g, int s, int r, int i)
    {
        double v = 0.0;
        for (int l = 0; l < a.Length; l++) v += a[l] * g.Partials[i][s, r, l];
        return v;
    }

    private static double Contract(double[] a, double[,] j, int i)
    {
        double v = 0.0;
        for (int l = 0; l < a.Length; l++) v += a[l] * j[l, i];
        return v;
    }

    private static double[,] Outer(double[] a, double[] b)
    {
        var o = new double[a.Length, b.Length];
        for (int i = 0; i < a.Length; i++)
            for (int j = 0; j < b.Length; j++) o[i, j] = a[i] * b[j];
        return o;
    }

    private static double[,] Transpose(double[,] m)
    {
        var t = new double[m.GetLength(1), m.GetLength(0)];
        for (int i = 0; i < m.GetLength(0); i++)
            for (int j = 0; j < m.GetLength(1); j++) t[j, i] = m[i, j];
        return t;
    }

    private static double Frobenius(double[,] a, double[,] b)
    {
        double s = 0.0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) s += (a[i, j] - b[i, j]) * (a[i, j] - b[i, j]);
        return Math.Sqrt(s);
    }

    // The leading eigenvector of a symmetric 2 × 2 matrix, first entry nonnegative.
    private static double[] Leading(double[,] c)
    {
        double tr = c[0, 0] + c[1, 1], det = c[0, 0] * c[1, 1] - c[0, 1] * c[1, 0];
        double lambda = 0.5 * (tr + Math.Sqrt(tr * tr - 4.0 * det));
        double x = c[0, 1], y = lambda - c[0, 0], n = Math.Sqrt(x * x + y * y);
        return x < 0.0 ? [-x / n, -y / n] : [x / n, y / n];
    }

    // The worst |exact − sampled|/se over each pair's three matrices, where the error is positive.
    private static double WorstInErrors(CoActiveAnalysis exact, CoActiveAnalysis sampled)
    {
        double worst = 0.0;
        for (int s = 0; s < exact.Posterior.Pairs.Count; s++)
        {
            CoActivePair e = exact.Posterior.Pairs[s], q = sampled.Posterior.Pairs[s];
            CoActiveError se = sampled.MonteCarloError![s];
            foreach ((double[,] x, double[,] y, double[,] z) in new[] { (e.Cf, q.Cf, se.Cf), (e.Cg, q.Cg, se.Cg), (e.Cfg, q.Cfg, se.Cfg) })
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++)
                    {
                        double gap = Math.Abs(x[i, j] - y[i, j]);
                        if (z[i, j] > 0.0) worst = Math.Max(worst, gap / z[i, j]);
                        else Assert.True(gap <= 1e-10 * Math.Max(1.0, Math.Abs(x[i, j])), $"pair {s}: a zero error with gap {gap:E2}");
                    }
        }
        return worst;
    }

    private static void AssertClose(double[,] expected, double[,] actual, double tolerance)
    {
        double gap = FunctionalActiveSubspaceTests.RelativeGap(actual, expected);
        Assert.True(gap <= tolerance, $"relative gap {gap:E2}");
    }
}

/// <summary>Two scalar responses on one set of inputs, and a scalar fit on the functional fixture's inputs (E8's slice 3).</summary>
internal static class CoActiveFixture
{
    // x0 and x1 on the scalar fixture's grids, ends 0 and 1; f a hinge plus an interaction, g a curve, an opposed hinge and an interaction.
    public static (BassData F, BassData G) Data()
    {
        var rng = new Xoshiro256PlusPlus(17);
        const int n = 60;
        var x0 = new double[n];
        var x1 = new double[n];
        var yf = new double[n];
        var yg = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = (i % 10) / 9.0;
            x1[i] = (i * 7 % 11) / 10.0;
            yf[i] = 1 + 2 * Math.Max(x0[i] - 0.4, 0) + 3 * x0[i] * x1[i] + 0.1 * (rng.NextDouble() - 0.5);
            yg[i] = 0.5 + 2 * x1[i] * x1[i] - 1.5 * Math.Max(0.6 - x0[i], 0) + 2 * x0[i] * x1[i] + 0.1 * (rng.NextDouble() - 0.5);
        }
        BassPredictor[] predictors = [new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)];
        return (new BassData(predictors, new ScalarResponse(yf)), new BassData(predictors, new ScalarResponse(yg)));
    }

    private static readonly Lazy<(BassFit, BassFit)> TwoFits = new(() =>
    {
        (BassData f, BassData g) = Data();
        return (Bass.Fit(SensitivityFixture.Config(), f), Bass.Fit(SensitivityFixture.Config(), g));
    });

    public static (BassFit F, BassFit G) Fits() => TwoFits.Value;

    // The functional fixture's inputs with its response's mean over the output locations as a scalar response.
    private static readonly Lazy<(BassData, BassFit)> Scalar = new(() =>
    {
        BassData functional = SensitivityFixture.Functional();
        var response = (FunctionalResponse)functional.Response;
        int n = response.Cases, m = response.Locations;
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            for (int l = 0; l < m; l++) y[i] += response[i, l];
            y[i] /= m;
        }
        var data = new BassData(functional.Predictors, new ScalarResponse(y));
        return (data, Bass.Fit(SensitivityFixture.Config(), data));
    });

    public static (BassData Data, BassFit Fit) OnFunctionalInputs() => Scalar.Value;
}
