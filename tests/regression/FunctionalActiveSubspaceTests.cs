using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 3 gates (the active-subspace brief) on functional outputs, C = E[JᵀWJ] under a declared output metric W:
/// <list type="bullet">
/// <item>a basis fit's C against tensor quadrature of its reconstructed Jacobian, built from each component's gradient readout and the
/// reconstruction's map G = diag(y_s)K, under the identity and a full W;</item>
/// <item>a direct functional fit's C against tensor quadrature of the gradient readout at the training output locations;</item>
/// <item>the direct fit's exact and sampled routes within the sampled route's Monte Carlo error under a full W;</item>
/// <item>W's readings: a location's own C under e_ℓe_ℓᵀ, linearity in W, the default identity, the checks without a violation, and the
/// refusals.</item>
/// </list>
/// The quadrature splits each input at its knots, so it integrates the piecewise polynomial gradients to rounding; it shares no code
/// with the moment route.
/// </summary>
public sealed class FunctionalActiveSubspaceTests(ITestOutputHelper output)
{
    // ── Quadrature references ────────────────────────────────────────────────────────────────────────────────────────

    // The tensor rule over two continuous inputs under the law, each input split at its knots (BassGradientMomentTests.Rule), as cases
    // in the inputs' own units (the fixtures' boxes are [0, 1], so box and own units agree).
    internal static (PredictionCases Cases, double[] W) Tensor(BassInputs inputs, ProductMeasure measure, int nodes = 4)
    {
        var rules = new (double[] X, double[] W)[2];
        var names = new string[2];
        for (int k = 0; k < 2; k++)
        {
            int v = measure.OfContinuous(k);
            rules[k] = BassGradientMomentTests.Rule(measure.Continuous(v)!, inputs.ContinuousKnots(k).ToArray(), nodes, 1.0);
            names[k] = measure.Names[v];
        }
        int n = rules[0].X.Length * rules[1].X.Length, at = 0;
        var x0 = new double[n];
        var x1 = new double[n];
        var w = new double[n];
        for (int i = 0; i < rules[0].X.Length; i++)
            for (int j = 0; j < rules[1].X.Length; j++, at++)
            {
                x0[at] = rules[0].X[i];
                x1[at] = rules[1].X[j];
                w[at] = rules[0].W[i] * rules[1].W[j];
            }
        return (new PredictionCases([new ContinuousColumn(names[0], x0), new ContinuousColumn(names[1], x1)]), w);
    }

    // Each point's Jacobian J[ℓ, i] of a direct fit's draw b, from the gradient readout.
    internal static double[][,] DirectJacobians(BassGradient g, int b)
    {
        int n = g.Partials[0].GetLength(1), m = g.Partials[0].GetLength(2), p = g.Partials.Count;
        var j = new double[n][,];
        for (int r = 0; r < n; r++)
        {
            j[r] = new double[m, p];
            for (int l = 0; l < m; l++)
                for (int i = 0; i < p; i++) j[r][l, i] = g.Partials[i][b, r, l];
        }
        return j;
    }

    // Each point's Jacobian of a basis fit's joint draw: Σ_k G_ℓk ∂w_k/∂x_i, each component's partials from its own gradient readout.
    internal static double[][,] BasisJacobians(BassBasisFit fit, PredictionCases cases, int joint)
    {
        int r = fit.Components.Count;
        var map = fit.Decomposition.Map;
        int m = map.Locations;
        var component = new BassGradient[r];
        for (int k = 0; k < r; k++) component[k] = fit.Components[k].Gradient(fit.ComponentData(k), cases, [fit.JointDraws.Draw(k, joint)]);
        int n = component[0].Partials[0].GetLength(1), p = component[0].Partials.Count;
        var j = new double[n][,];
        for (int q = 0; q < n; q++)
        {
            j[q] = new double[m, p];
            for (int l = 0; l < m; l++)
                for (int k = 0; k < r; k++)
                {
                    double g = map.Scales[l] * map.Basis(l, k);
                    for (int i = 0; i < p; i++) j[q][l, i] += g * component[k].Partials[i][0, q, 0];
                }
        }
        return j;
    }

    // Σ_r w_r J_rᵀ W J_r.
    internal static double[,] Quadrature(double[][,] jacobians, double[] w, double[,] metric)
    {
        int m = jacobians[0].GetLength(0), p = jacobians[0].GetLength(1);
        var c = new double[p, p];
        for (int r = 0; r < jacobians.Length; r++)
        {
            double[,] j = jacobians[r];
            for (int a = 0; a < p; a++)
                for (int b = 0; b < p; b++)
                {
                    double s = 0.0;
                    for (int l = 0; l < m; l++)
                    {
                        double jl = j[l, a];
                        if (jl == 0.0) continue;
                        for (int q = 0; q < m; q++) s += jl * metric[l, q] * j[q, b];
                    }
                    c[a, b] += w[r] * s;
                }
        }
        return c;
    }

    internal static double[,] Identity(int m)
    {
        var w = new double[m, m];
        for (int l = 0; l < m; l++) w[l, l] = 1.0;
        return w;
    }

    // AᵀA/m for a standard normal A: full, positive definite almost surely.
    internal static double[,] RandomMetric(int m, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var a = new double[m, m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++) a[i, j] = Normal.SampleStandard(rng);
        var w = new double[m, m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
            {
                double s = 0.0;
                for (int k = 0; k < m; k++) s += a[k, i] * a[k, j];
                w[i, j] = s / m;
            }
        return w;
    }

    internal static double RelativeGap(double[,] actual, double[,] expected)
    {
        double worst = 0.0, scale = 0.0;
        foreach (double x in expected) scale = Math.Max(scale, Math.Abs(x));
        for (int i = 0; i < expected.GetLength(0); i++)
            for (int j = 0; j < expected.GetLength(1); j++) worst = Math.Max(worst, Math.Abs(actual[i, j] - expected[i, j]));
        return scale > 0.0 ? worst / scale : worst;
    }

    private static int[] Spread(int count, int k) => [.. Enumerable.Range(0, k).Select(i => (int)((long)i * (count - 1) / (k - 1)))];

    // ── Basis fits ───────────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// The brief's gate: a basis fit's C through the components' cross moments, Σ_kk′ (GᵀWG)_kk′ E[∇w_k∇w_k′ᵀ], against the direct
    /// computation on its reconstructed outputs, J = GJ_w by each component's gradient readout and Σ w JᵀWJ by quadrature, per joint draw,
    /// under the identity and a full W; the checks find no violation.
    /// </summary>
    [Theory]
    [InlineData("identity")]
    [InlineData("full")]
    public void BasisFit_C_IsTheQuadratureOfItsReconstructedJacobian(string metric)
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        int m = fit.Decomposition.Locations;
        double[,] w = metric == "identity" ? Identity(m) : RandomMetric(m, 41);
        int[] draws = Spread(fit.JointDraws.Count, 5);
        ActiveSubspaceAnalysis a = fit.ActiveSubspace(new ActiveSubspaceSpec { Draws = draws, OutputMetric = metric == "identity" ? null : new OutputMetric(w) });
        Assert.True(fit.Components.Count >= 2, $"the fixture has {fit.Components.Count} components");
        BassInputs inputs = BassInputs.Prepare(fit.ComponentData(0), fit.Components[0].Transform);
        (PredictionCases cases, double[] weights) = Tensor(inputs, a.Measure);
        double worst = 0.0;
        for (int d = 0; d < draws.Length; d++)
        {
            double[,] reference = Quadrature(BasisJacobians(fit, cases, draws[d]), weights, w);
            worst = Math.Max(worst, RelativeGap(a.Posterior.Draws[d].C, reference));
        }
        output.WriteLine($"{metric}: {fit.Components.Count} components, {m} locations, {draws.Length} joint draws; worst gap to quadrature {worst:E2}; " +
                         $"Poincaré ratio up to {MaxFinite(a.Poincare.Ratio):F3}");
        Assert.True(worst <= 1e-10, $"{worst:E2}");
        Assert.Equal(ActiveSubspaceRoute.Exact, a.Route);
        Assert.Equal(0, a.Poincare.Violations);
        Assert.Equal(draws, a.DrawIndex);
    }

    // ── Direct functional fits ───────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// A direct functional fit's C, Σ_ab (HᵀWH)_ab E[∇B_a∇B_bᵀ] with H's columns the terms' output-side factors at the training locations
    /// times their coefficients, against quadrature of the gradient readout's Jacobian at those locations, under the identity, a
    /// diagonal of trapezoid weights and a full W.
    /// </summary>
    [Theory]
    [InlineData("identity")]
    [InlineData("trapezoid")]
    [InlineData("full")]
    public void DirectFunctional_C_IsTheQuadratureOfItsJacobian(string metric)
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        int m = data.Response.Locations;
        double[,] w = metric switch
        {
            "identity" => Identity(m),
            "trapezoid" => OutputMetric.Diagonal([.. Enumerable.Range(0, m).Select(l => (l == 0 || l == m - 1 ? 0.5 : 1.0) / (m - 1))]).Matrix,
            _ => RandomMetric(m, 43),
        };
        int[] draws = Spread(fit.Draws.Sum(c => c.Count), 6);
        var spec = new ActiveSubspaceSpec { Draws = draws, Dimension = 2, OutputMetric = metric == "identity" ? null : new OutputMetric(w) };
        ActiveSubspaceAnalysis a = fit.ActiveSubspace(data, spec);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        (PredictionCases cases, double[] weights) = Tensor(inputs, a.Measure);
        BassGradient g = fit.Gradient(data, cases, draws);
        double worst = 0.0;
        for (int d = 0; d < draws.Length; d++) worst = Math.Max(worst, RelativeGap(a.Posterior.Draws[d].C, Quadrature(DirectJacobians(g, d), weights, w)));
        output.WriteLine($"{metric}: {m} locations, {draws.Length} draws; worst gap to quadrature {worst:E2}; E[C]'s eigenvalues " +
                         $"{a.Posterior.MeanCEigenvalues[0]:F4}, {a.Posterior.MeanCEigenvalues[1]:F4}; Poincaré ratio up to {MaxFinite(a.Poincare.Ratio):F3}");
        Assert.True(worst <= 1e-10, $"{worst:E2}");
        Assert.Equal(0, a.Poincare.Violations);
        Assert.Null(a.Poincare.JacobiRatio);
    }

    /// <summary>The direct fit's exact and sampled routes under a full W, each entry within 4.5 of the sampled route's standard errors.</summary>
    [Fact]
    public void DirectFunctional_ExactAndSampledRoutes_AgreeWithinMonteCarloError()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        int m = data.Response.Locations;
        var spec = new ActiveSubspaceSpec { Draws = Spread(fit.Draws.Sum(c => c.Count), 12), OutputMetric = new OutputMetric(RandomMetric(m, 47)), Points = 8192, Seed = 5 };
        ActiveSubspaceAnalysis exact = fit.ActiveSubspace(data, spec with { Route = ActiveSubspaceRoute.Exact });
        ActiveSubspaceAnalysis sampled = fit.ActiveSubspace(data, spec with { Route = ActiveSubspaceRoute.Sampled });
        double worst = 0.0;
        for (int d = 0; d < exact.Posterior.Draws.Count; d++)
            for (int i = 0; i < 2; i++)
                for (int j = i; j < 2; j++)
                {
                    double se = sampled.MonteCarloError![d][i, j];
                    double gap = Math.Abs(exact.Posterior.Draws[d].C[i, j] - sampled.Posterior.Draws[d].C[i, j]);
                    if (se > 0.0) worst = Math.Max(worst, gap / se);
                    else Assert.True(gap <= 1e-10 * Math.Max(1.0, Math.Abs(exact.Posterior.Draws[d].C[i, j])), $"draw {d}: a zero error with gap {gap:E2}");
                }
        output.WriteLine($"worst |exact − sampled| / se {worst:F2} over {exact.Posterior.Draws.Count} draws; sampled violations {sampled.Poincare.Violations}");
        Assert.True(worst <= 4.5, $"{worst:F2}");
        Assert.Equal(0, sampled.Poincare.Violations);
    }

    /// <summary>
    /// W's readings: e_ℓe_ℓᵀ reads location ℓ's own C, the scalar C of that location's surface; C is linear in W, so a diagonal W gives
    /// Σ_ℓ w_ℓ C_ℓ; and no metric is the identity.
    /// </summary>
    [Fact]
    public void OutputMetric_ReadsLocationsLinearly()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        int m = data.Response.Locations;
        int[] draws = Spread(fit.Draws.Sum(c => c.Count), 5);
        var spec = new ActiveSubspaceSpec { Draws = draws };
        var perLocation = new ActiveSubspaceAnalysis[m];
        for (int l = 0; l < m; l++)
        {
            var e = new double[m];
            e[l] = 1.0;
            perLocation[l] = fit.ActiveSubspace(data, spec with { OutputMetric = OutputMetric.Diagonal(e) });
        }
        double[] weights = [.. Enumerable.Range(0, m).Select(l => 0.3 + 0.1 * l)];
        ActiveSubspaceAnalysis diagonal = fit.ActiveSubspace(data, spec with { OutputMetric = OutputMetric.Diagonal(weights) });
        ActiveSubspaceAnalysis identity = fit.ActiveSubspace(data, spec);
        ActiveSubspaceAnalysis declared = fit.ActiveSubspace(data, spec with { OutputMetric = new OutputMetric(Identity(m)) });
        double worst = 0.0;
        for (int d = 0; d < draws.Length; d++)
        {
            var sum = new double[2, 2];
            var plain = new double[2, 2];
            for (int l = 0; l < m; l++)
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++)
                    {
                        sum[i, j] += weights[l] * perLocation[l].Posterior.Draws[d].C[i, j];
                        plain[i, j] += perLocation[l].Posterior.Draws[d].C[i, j];
                    }
            worst = Math.Max(worst, RelativeGap(diagonal.Posterior.Draws[d].C, sum));
            worst = Math.Max(worst, RelativeGap(identity.Posterior.Draws[d].C, plain));
            Assert.Equal(identity.Posterior.Draws[d].C, declared.Posterior.Draws[d].C);
        }
        output.WriteLine($"linearity in W: worst gap {worst:E2}");
        Assert.True(worst <= 1e-12, $"{worst:E2}");

        // Location ℓ's C is its surface's: quadrature of that location's partials alone.
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        (PredictionCases cases, double[] w) = Tensor(inputs, identity.Measure);
        BassGradient g = fit.Gradient(data, cases, draws);
        var single = new double[m, m];
        single[3, 3] = 1.0;
        for (int d = 0; d < draws.Length; d++)
            Assert.True(RelativeGap(perLocation[3].Posterior.Draws[d].C, Quadrature(DirectJacobians(g, d), w, single)) <= 1e-10);
    }

    /// <summary>The readout refuses a metric it cannot read: on a scalar fit, of the wrong size, not symmetric, not finite, with a
    /// negative eigenvalue on either route; a basis fit's sampled route waits for slice 4.</summary>
    [Fact]
    public void OutputMetric_RefusesWhatItDoesNotRead()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        int m = data.Response.Locations;
        var spec = new ActiveSubspaceSpec { Draws = [0, 1] };
        Assert.Throws<ArgumentException>(() => SensitivityFixture.FitScalar(1).ActiveSubspace(SensitivityFixture.Scalar(), spec with { OutputMetric = new OutputMetric(Identity(1)) }));
        Assert.Throws<ArgumentException>(() => fit.ActiveSubspace(data, spec with { OutputMetric = new OutputMetric(Identity(m + 1)) }));
        double[,] skew = Identity(m);
        skew[0, 1] = 0.5;
        Assert.Throws<ArgumentException>(() => fit.ActiveSubspace(data, spec with { OutputMetric = new OutputMetric(skew) }));
        double[,] infinite = Identity(m);
        infinite[2, 2] = double.PositiveInfinity;
        Assert.Throws<ArgumentException>(() => fit.ActiveSubspace(data, spec with { OutputMetric = new OutputMetric(infinite) }));
        double[,] indefinite = Identity(m);
        indefinite[0, 0] = -1.0;
        Assert.Throws<ArgumentException>(() => fit.ActiveSubspace(data, spec with { OutputMetric = new OutputMetric(indefinite) }));
        Assert.Throws<ArgumentException>(() => fit.ActiveSubspace(data, spec with { OutputMetric = new OutputMetric(indefinite), Route = ActiveSubspaceRoute.Sampled, Points = 64 }));
        BassBasisFit basis = BassOutputComponentsTests.FitOf();
        Assert.Throws<ArgumentException>(() => basis.ActiveSubspace(spec with { Route = ActiveSubspaceRoute.Sampled }));
        Assert.Throws<ArgumentException>(() => basis.ActiveSubspace(spec with { OutputMetric = new OutputMetric(Identity(basis.Decomposition.Locations + 1)) }));
    }

    private static double MaxFinite(double[,] m)
    {
        double worst = 0.0;
        foreach (double v in m)
            if (double.IsFinite(v)) worst = Math.Max(worst, v);
        return worst;
    }
}
