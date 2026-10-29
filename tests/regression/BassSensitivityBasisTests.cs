using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// Basis sensitivity (O7, <see cref="BassBasisSensitivity.Of"/>) on basis fits:
/// <list type="bullet">
/// <item>a one-component fit is the direct decomposition of its component, scaled by y_s(ℓ)²·K_ℓ1² at each location;</item>
/// <item>each location's and the integrated decomposition against a tensor quadrature of the reconstructed output under a
/// non-uniform product measure, on a fit of several components whose cross terms matter;</item>
/// <item>the Sobol identities in both forms, and the law of total variance between them;</item>
/// <item>the pairing of draws across components (thinned and resampled), singleton locations, zero variance (undefined,
/// never 0 or 1), the order, the subset budget and the refusals.</item>
/// </list>
/// </summary>
public sealed class BassSensitivityBasisTests
{
    private static readonly Dictionary<string, InputMeasure> Measures = new()
    {
        ["z0"] = new MixtureMeasure([new NormalComponent(0.5, 0.35, 0.15), new NormalComponent(0.5, 0.75, 0.25)]),
        ["z1"] = new UniformMeasure(0.1, 0.9),
    };

    [Fact]
    public void OneComponent_IsTheDirectDecompositionScaled()
    {
        BassBasisFit full = BassOutputComponentsTests.FitOf();
        (double[,] y, _) = BassOutputComponentsTests.Problem();
        // Rank 1 keeps the first right singular vector, component 0's coefficient response, so its fit carries over.
        BassDecomposition one = BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 1, ScaleByColumnSd = true });
        BassBasisFit fit = BassComponents.Assemble(one, full.Predictors, [full.Components[0]], 3);
        var spec = new SensitivitySpec { Measures = Measures };
        BasisSobolAnalysis basis = fit.Sensitivity(spec);
        SobolDecomposition direct = BassSensitivity.Of(fit.Components[0], fit.ComponentData(0), spec).Decompositions[0];
        Assert.Equal(direct.Count, basis.Count);
        Assert.Equal(direct.Subsets, basis.Subsets);
        for (int d = 0; d < direct.Count; d++)
        {
            // The ℓ-free state is the direct decomposition itself, to the bit, before the direct one clamps.
            Assert.Equal(direct.Variance[d], basis.ComponentCovariance(d)[0, 0]);
            for (int u = 0; u < direct.Subsets.Count; u++) AssertSameOrClamped(direct.MainVariance[d, u], basis.ComponentMainCovariance(d, u)[0, 0], direct.Variance[d]);
            for (int v = 0; v < direct.Inputs.Count; v++) AssertSameOrClamped(direct.TotalVariance[d, v], basis.ComponentTotalCovariance(d, v)[0, 0], direct.Variance[d]);
        }
        LinearReconstruction map = one.Map;
        foreach (int l in new[] { 0, 3, 7, map.Locations - 1 })
        {
            SobolDecomposition at = basis.AtLocation(l);
            double scale = map.Scales[l] * map.Scales[l] * map.Basis(l, 0) * map.Basis(l, 0);
            Assert.True(scale > 0.0);
            for (int d = 0; d < direct.Count; d++)
            {
                Assert.True(Close(at.Variance[d], scale * direct.Variance[d]), $"location {l}, draw {d}");
                for (int u = 0; u < direct.Subsets.Count; u++) Assert.True(Close(at.MainVariance[d, u], scale * direct.MainVariance[d, u]));
                for (int v = 0; v < direct.Inputs.Count; v++) Assert.True(Close(at.TotalVariance[d, v], scale * direct.TotalVariance[d, v]));
                Assert.Equal(direct.Status[d], at.Status[d]);
                if (!direct.IsDefined(d)) continue;
                for (int u = 0; u < direct.Subsets.Count; u++) Assert.Equal(direct.MainIndex(d, u), at.MainIndex(d, u), 13);
            }
        }
    }

    [Fact]
    public void CrossComponentTerms_MatchTensorQuadratureOfTheReconstruction()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        BasisSobolAnalysis analysis = fit.Sensitivity(new SensitivitySpec { Measures = Measures });
        int r = analysis.Components, m = analysis.Locations;
        Assert.True(r >= 2, $"the fixture has {r} components");
        BassInputs inputs = BassInputs.Prepare(fit.ComponentData(0), fit.Components[0].Transform);
        var reference = new SensitivityReference(inputs, analysis.Measure, 1);
        (double[] x0, double[] w0) = reference.Nodes(0);
        (double[] x1, double[] w1) = reference.Nodes(1);
        int n0 = x0.Length, n1 = x1.Length, cells = n0 * n1;
        var cases = new CaseInputs([[.. x0.SelectMany(u => x1.Select(_ => u))], [.. x0.SelectMany(_ => x1)]], []);
        double mass0 = w0.Sum(), mass1 = w1.Sum();
        var basis = new BassBasis(inputs, 1);
        LinearReconstruction map = fit.Decomposition.Map;
        var pooled = fit.Components.Select(c => c.Draws.SelectMany(chain => chain).ToList()).ToArray();

        // The draw holding the most terms, an interaction among them, and two others.
        int richest = Enumerable.Range(0, analysis.Count).OrderByDescending(d => Enumerable.Range(0, r).Sum(k => pooled[k][fit.JointDraws.Draw(k, d)].Structure.TermIds.Length)).First();
        double crossShare = 0.0;
        foreach (int d in new[] { 0, analysis.Count / 2, richest })
        {
            // Each component's surface on the grid, then the output at every location.
            var w = new double[r][];
            var column = new double[cells];
            for (int k = 0; k < r; k++)
            {
                BassDraw draw = pooled[k][fit.JointDraws.Draw(k, d)];
                w[k] = Enumerable.Repeat(draw.Coefficients[0], cells).ToArray();
                for (int j = 0; j < draw.Structure.TermIds.Length; j++)
                {
                    basis.CaseColumn(fit.Components[k].Terms[draw.Structure.TermIds[j]], cases, column);
                    for (int i = 0; i < cells; i++) w[k][i] += draw.Coefficients[1 + j] * column[i];
                }
            }
            var f = new double[m][];
            var mean = new double[m];
            var cond0 = new double[m, n0];   // E[f_ℓ | z0]
            var cond1 = new double[m, n1];   // E[f_ℓ | z1]
            for (int l = 0; l < m; l++)
            {
                f[l] = new double[cells];
                for (int i = 0; i < cells; i++)
                {
                    double s = 0.0;
                    for (int k = 0; k < r; k++) s += map.Basis(l, k) * w[k][i];
                    f[l][i] = map.Offsets[l] + map.Scales[l] * s;
                    double weight = w0[i / n1] * w1[i % n1];
                    mean[l] += weight * f[l][i];
                    cond0[l, i / n1] += w1[i % n1] * f[l][i] / mass1;
                    cond1[l, i % n1] += w0[i / n1] * f[l][i] / mass0;
                }
                mean[l] /= mass0 * mass1;
            }

            // Per location: the variance, the two mains, and the interaction as the rest.
            for (int l = 0; l < m; l++)
            {
                SobolDecomposition at = analysis.AtLocation(l);
                double v = 0.0, m0 = 0.0, m1 = 0.0;
                for (int i = 0; i < cells; i++) v += w0[i / n1] * w1[i % n1] * Sq(f[l][i] - mean[l]) / (mass0 * mass1);
                for (int a = 0; a < n0; a++) m0 += w0[a] * Sq(cond0[l, a] - mean[l]) / mass0;
                for (int b = 0; b < n1; b++) m1 += w1[b] * Sq(cond1[l, b] - mean[l]) / mass1;
                Assert.True(Near(at.Variance[d], v, v), $"draw {d}, location {l}: variance {at.Variance[d]} against {v}");
                Assert.True(Near(MainOf(at, d, [0]), m0, v), $"draw {d}, location {l}: main z0 {MainOf(at, d, [0])} against {m0}");
                Assert.True(Near(MainOf(at, d, [1]), m1, v), $"draw {d}, location {l}: main z1 {MainOf(at, d, [1])} against {m1}");
                Assert.True(Near(MainOf(at, d, [0, 1]), v - m0 - m1, v));
                Assert.True(Near(at.TotalVariance[d, 0], v - m1, v) && Near(at.TotalVariance[d, 1], v - m0, v));

                // The cross terms' share: the variance less each component's own.
                double own = 0.0;
                double[,] c = analysis.ComponentCovariance(d);
                for (int k = 0; k < r; k++) own += Sq(map.Scales[l] * map.Basis(l, k)) * c[k, k];
                if (v > 0.0) crossShare = Math.Max(crossShare, Math.Abs(v - own) / v);
            }

            // Integrated: the location uniform over the m locations, one more input.
            SobolDecomposition whole = analysis.Integrated;
            double grand = mean.Average(), total = 0.0, location = 0.0, i0 = 0.0, i1 = 0.0, i0l = 0.0, i1l = 0.0;
            for (int l = 0; l < m; l++)
            {
                for (int i = 0; i < cells; i++) total += w0[i / n1] * w1[i % n1] * Sq(f[l][i] - grand) / (mass0 * mass1 * m);
                location += Sq(mean[l] - grand) / m;
            }
            for (int a = 0; a < n0; a++)
            {
                double avg = 0.0;
                for (int l = 0; l < m; l++) avg += cond0[l, a] / m;
                i0 += w0[a] * Sq(avg - grand) / mass0;
                for (int l = 0; l < m; l++) i0l += w0[a] * Sq(cond0[l, a] - grand) / (mass0 * m);
            }
            for (int b = 0; b < n1; b++)
            {
                double avg = 0.0;
                for (int l = 0; l < m; l++) avg += cond1[l, b] / m;
                i1 += w1[b] * Sq(avg - grand) / mass1;
                for (int l = 0; l < m; l++) i1l += w1[b] * Sq(cond1[l, b] - grand) / (mass1 * m);
            }
            Assert.True(Near(whole.Variance[d], total, total), $"draw {d}: integrated {whole.Variance[d]} against {total}");
            Assert.True(Near(MainOf(whole, d, [2]), location, total));
            Assert.True(Near(MainOf(whole, d, [0]), i0, total) && Near(MainOf(whole, d, [1]), i1, total));
            // Var(E[F | z, ℓ]) = V_z + V_ℓ + V_zℓ.
            Assert.True(Near(MainOf(whole, d, [0, 2]), i0l - i0 - location, total), $"draw {d}: z0 with the location");
            Assert.True(Near(MainOf(whole, d, [1, 2]), i1l - i1 - location, total), $"draw {d}: z1 with the location");
            Assert.True(Near(whole.TotalVariance[d, 2], total - Closed(whole, d, [0, 1]), total));
            // An input's total, V − Var(E[F | the other input, ℓ]).
            Assert.True(Near(whole.TotalVariance[d, 0], total - i1l, total) && Near(whole.TotalVariance[d, 1], total - i0l, total), $"draw {d}: totals");
        }
        Assert.True(crossShare > 1e-3, $"the cross terms' largest share is {crossShare}");
    }

    [Fact]
    public void Forms_ObeyTheSobolIdentities_AndTheLawOfTotalVariance()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        BasisSobolAnalysis analysis = fit.Sensitivity(new SensitivitySpec { Measures = Measures });
        IReadOnlyList<SobolDecomposition> locations = analysis.PerLocation();
        SobolDecomposition whole = analysis.Integrated;
        int m = analysis.Locations, p = analysis.Measure.Count;
        Assert.Equal(new[] { "z0", "z1", BasisSobolAnalysis.LocationInput }, whole.Inputs);
        Assert.Equal(analysis.MaxOrder + 1, analysis.IntegratedMaxOrder);
        Assert.Equal(analysis.IntegratedMaxOrder, whole.Subsets.Max(u => u.Length));
        LinearReconstruction map = fit.Decomposition.Map;
        for (int d = 0; d < analysis.Count; d++)
        {
            foreach (SobolDecomposition part in locations.Append(whole))
                Identities(part, d);
            // Within the locations plus the variance of their means; the inputs' totals are their locations' means, and each
            // subset's main splits into the location-averaged surface's and its interaction with the location.
            double[] mu = analysis.ComponentMeans(d);
            double[] means = [.. Enumerable.Range(0, m).Select(l => map.Offsets[l] + map.Scales[l] * Enumerable.Range(0, analysis.Components).Sum(k => map.Basis(l, k) * mu[k]))];
            double grand = means.Average(), between = means.Average(v => Sq(v - grand));
            double within = locations.Average(s => s.Variance[d]);
            Assert.True(Near(whole.Variance[d], within + between, whole.Variance[d]), $"draw {d}");
            Assert.True(Near(MainOf(whole, d, [p]), between, whole.Variance[d]));
            for (int v = 0; v < p; v++) Assert.True(Near(whole.TotalVariance[d, v], locations.Average(s => s.TotalVariance[d, v]), whole.Variance[d]));
            for (int u = 0; u < analysis.Subsets.Count; u++)
            {
                int[] subset = analysis.Subsets[u];
                double split = MainOf(whole, d, subset) + MainOf(whole, d, [.. subset, p]);
                Assert.True(Near(split, locations.Average(s => s.MainVariance[d, u]), whole.Variance[d]), $"draw {d}, subset {string.Join(",", subset)}");
            }
        }
        Assert.Equal(1.0, whole.MeanMainIndices().Sum(), 10);
        Assert.All(locations, s => Assert.Equal(1.0, s.MeanMainIndices().Sum(), 10));
    }

    [Fact]
    public void JointDraws_PairTheComponentsDraws_ThinnedOrResampled()
    {
        BassBasisFit full = BassOutputComponentsTests.FitOf();
        // Component 1 keeps half of each chain: the others are thinned evenly to it.
        BassFit[] components = [.. full.Components];
        components[1] = components[1] with { Draws = [.. components[1].Draws.Select(chain => Store(chain.Take(chain.Count / 2)))] };
        BassBasisFit thinned = BassComponents.Assemble(full.Decomposition, full.Predictors, components, 5);
        // Every component weighted: resampled, then equal-weight joint draws with no chains.
        BassFit[] weighted = [.. full.Components.Select(c => c with { DrawWeights = [.. Enumerable.Range(0, c.Draws.Sum(ch => ch.Count)).Select(i => 1.0 + i % 4)] })];
        BassBasisFit resampled = BassComponents.Assemble(full.Decomposition, full.Predictors, weighted, 5);
        Assert.True(resampled.JointDraws.Resampled);

        foreach (BassBasisFit fit in new[] { thinned, resampled })
        {
            BasisSobolAnalysis analysis = fit.Sensitivity(new SensitivitySpec { Measures = Measures });
            BassJointDraws joint = fit.JointDraws;
            Assert.Equal(joint.Count, analysis.Count);
            Assert.Null(analysis.Integrated.Weights);
            int at = 0;
            if (joint.ChainLengths is { } lengths)
                for (int c = 0; c < lengths.Count; c++)
                    for (int i = 0; i < lengths[c]; i++, at++)
                        Assert.True(analysis.Integrated.Chains[at] == c && analysis.Integrated.Draws[at] == i);
            else
                for (; at < joint.Count; at++) Assert.True(analysis.Integrated.Chains[at] == 0 && analysis.Integrated.Draws[at] == at);
            Assert.Equal(joint.Count, at);
            for (int k = 0; k < analysis.Components; k++)
            {
                BassFit component = fit.Components[k] with { DrawWeights = null };
                var spec = new SensitivitySpec { Measures = Measures };
                SobolDecomposition direct = BassSensitivity.Of(component, fit.ComponentData(k), spec).Decompositions[0];
                double[] f0 = BassSensitivity.Effects(component, fit.ComponentData(k), spec, new EffectSpec { Subsets = [] }).At[0].Mean;
                for (int d = 0; d < joint.Count; d++)
                {
                    int own = joint.Draw(k, d);
                    Assert.True(Close(analysis.ComponentCovariance(d)[k, k], direct.Variance[own]), $"component {k}, joint draw {d}");
                    Assert.Equal(f0[own], analysis.ComponentMeans(d)[k], 12);
                }
            }
        }
    }

    [Fact]
    public void CrossComponentCovariances_AreKept_AndTheComponentOrderIsFree()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        var spec = new SensitivitySpec { Measures = Measures };
        BasisSobolAnalysis analysis = fit.Sensitivity(spec);
        int r = analysis.Components, m = analysis.Locations;
        bool crossed = false;
        for (int d = 0; d < analysis.Count; d++)
        {
            double[,] c = analysis.ComponentCovariance(d);
            for (int k = 0; k < r; k++)
            {
                Assert.True(c[k, k] >= 0.0);
                for (int k2 = k + 1; k2 < r; k2++)
                {
                    Assert.True(c[k, k2] * c[k, k2] <= c[k, k] * c[k2, k2] * (1 + 1e-12), $"draw {d}: Cauchy–Schwarz for ({k}, {k2})");
                    crossed |= Math.Abs(c[k, k2]) > 1e-6 * Math.Sqrt(c[k, k] * c[k2, k2]);
                }
            }
            // The subsets' matrices sum to the whole.
            for (int k = 0; k < r; k++)
                for (int k2 = 0; k2 < r; k2++)
                {
                    double sum = Enumerable.Range(0, analysis.Subsets.Count).Sum(u => analysis.ComponentMainCovariance(d, u)[k, k2]);
                    Assert.True(Math.Abs(sum - c[k, k2]) <= 1e-12 * (Math.Abs(c[k, k]) + Math.Abs(c[k2, k2]) + 1e-300));
                }
        }
        Assert.True(crossed, "independently fitted components still covary under the input measure");

        // The components and K's columns reversed: the same output, so the same decompositions.
        BassDecomposition d0 = fit.Decomposition;
        int n = d0.Cases;
        double[,] basis = d0.Map.BasisMatrix();
        var kr = new double[m, r];
        var wr = new double[r, n];
        for (int q = 0; q < r; q++)
        {
            for (int l = 0; l < m; l++) kr[l, q] = basis[l, r - 1 - q];
            for (int i = 0; i < n; i++) wr[q, i] = d0.Coefficient(r - 1 - q, i);
        }
        BassDecomposition reversed = BassDecomposition.FromBasis(kr, wr, d0.Map.Offsets.ToArray(), d0.Map.Scales.ToArray(), residuals: d0.Truncation.Values);
        BasisSobolAnalysis flipped = BassComponents.Assemble(reversed, fit.Predictors, [.. fit.Components.Reverse()], 5).Sensitivity(spec);
        foreach (SobolDecomposition[] pair in Enumerable.Range(0, m).Select(l => new[] { analysis.AtLocation(l), flipped.AtLocation(l) })
                                                        .Append([analysis.Integrated, flipped.Integrated]))
            for (int d = 0; d < analysis.Count; d++)
            {
                double v = pair[0].Variance[d];
                Assert.True(Near(pair[1].Variance[d], v, v), $"draw {d}: {pair[1].Variance[d]} against {v}");
                for (int u = 0; u < pair[0].Subsets.Count; u++) Assert.True(Near(pair[1].MainVariance[d, u], pair[0].MainVariance[d, u], v));
                for (int i = 0; i < pair[0].Inputs.Count; i++) Assert.True(Near(pair[1].TotalVariance[d, i], pair[0].TotalVariance[d, i], v));
            }
    }

    [Fact]
    public void SingletonLocation_MakesTheFormsAgree()
    {
        BassBasisFit fit = BasisSensitivityFixture.Mock(new double[,] { { 0.6, -1.1 } }, [3.0], [2.5]);
        BasisSobolAnalysis analysis = fit.Sensitivity();
        Assert.Equal(1, analysis.Locations);
        SobolDecomposition at = analysis.AtLocation(0), whole = analysis.Integrated;
        int p = analysis.Measure.Count;
        Assert.Equal(at.Variance[0], whole.Variance[0], 14);
        Assert.Equal(0.0, MainOf(whole, 0, [p]));
        Assert.Equal(0.0, whole.TotalIndex(0, p));
        for (int u = 0; u < at.Subsets.Count; u++) Assert.Equal(at.MainIndex(0, u), MainOf(whole, 0, at.Subsets[u]) / whole.Variance[0], 14);
        Assert.Throws<ArgumentOutOfRangeException>(() => analysis.AtLocation(1));
    }

    [Fact]
    public void ZeroVariance_LeavesTheIndicesUndefined()
    {
        // Location 1's basis row is zero, as a constant output column's is: its output does not vary.
        BassBasisFit dead = BasisSensitivityFixture.Mock(new double[,] { { 1.0, 0.5 }, { 0.0, 0.0 }, { -0.4, 1.0 } }, [1.0, 2.0, 3.0], [2.0, 3.0, 4.0]);
        BasisSobolAnalysis a = dead.Sensitivity();
        SobolDecomposition zero = a.AtLocation(1);
        Assert.Equal(DecompositionStatus.ZeroVariance, Assert.Single(zero.Status));
        Assert.True(double.IsNaN(zero.MainIndex(0, 0)) && double.IsNaN(zero.TotalIndex(0, 0)));
        Assert.All(zero.MeanMainIndices(), v => Assert.True(double.IsNaN(v)));   // undefined, never 0
        Assert.True(a.AtLocation(0).IsDefined(0) && a.Integrated.IsDefined(0));

        // Intercept-only draws: no location varies in x, but the location means do, so the integrated form is all location.
        var flat = new (BassTerm[], double[][])[] { ([], [[0.7]]), ([], [[-0.3]]) };
        double[,] k = { { 1.0, 0.5 }, { 0.2, 0.8 }, { -0.4, 1.0 } };
        BasisSobolAnalysis b = BasisSensitivityFixture.Mock(k, [1.0, 2.0, 3.0], [2.0, 3.0, 4.0], flat).Sensitivity();
        Assert.Empty(b.Subsets);
        Assert.All(b.PerLocation(), s => Assert.Equal(DecompositionStatus.ZeroVariance, Assert.Single(s.Status)));
        SobolDecomposition whole = b.Integrated;
        Assert.Equal(new[] { 2 }, Assert.Single(whole.Subsets));
        Assert.Equal(1.0, whole.MainIndex(0, 0), 14);
        Assert.Equal(1.0, whole.TotalIndex(0, 2), 14);
        Assert.Equal(0.0, whole.TotalIndex(0, 0));

        // And with constant offsets and zero intercepts, nothing varies at all.
        var none = new (BassTerm[], double[][])[] { ([], [[0.0]]), ([], [[0.0]]) };
        SobolDecomposition still = BasisSensitivityFixture.Mock(k, [5.0, 5.0, 5.0], [2.0, 3.0, 4.0], none).Sensitivity().Integrated;
        Assert.Equal(DecompositionStatus.ZeroVariance, Assert.Single(still.Status));
        Assert.True(double.IsNaN(still.TotalIndex(0, 2)));

        // A small variation on a large offset is real: the location means' part is zero only when its standard deviation is
        // within the tolerance of the means' magnitude (here 10⁻⁹ against a tolerance of 10⁻¹⁰ squared times 10¹²).
        SobolDecomposition small = BasisSensitivityFixture.Mock(k, [1e6, 1e6 + 1e-3, 1e6 - 1e-3], [2.0, 3.0, 4.0], none).Sensitivity().Integrated;
        Assert.True(small.IsDefined(0));
        Assert.Equal(2e-6 / 3, small.Variance[0], 1e-12);
        Assert.Equal(1.0, small.MainIndex(0, 0), 14);
    }

    [Fact]
    public void OrderAndBudget_BoundTheSubsets()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        BasisSobolAnalysis full = fit.Sensitivity();
        Assert.Contains(full.Subsets, u => u.Length == 2);
        BasisSobolAnalysis first = fit.Sensitivity(new SensitivitySpec { MaxOrder = 1 });
        Assert.Equal(1, first.MaxOrder);
        Assert.Equal(1, first.IntegratedMaxOrder);
        Assert.All(first.Subsets, u => Assert.Single(u));
        Assert.Equal(new[] { new[] { 0 }, new[] { 1 }, new[] { 2 } }, first.Integrated.Subsets);
        // The totals do not depend on the reported order.
        for (int d = 0; d < full.Count; d++)
            for (int v = 0; v < 2; v++) Assert.Equal(full.ComponentTotalCovariance(d, v), first.ComponentTotalCovariance(d, v));
        Assert.Equal(full.Integrated.TotalVariance, first.Integrated.TotalVariance);
        // Both forms are views of one analysis, so the functional form is not read.
        Assert.Equal(full.Integrated.Variance, fit.Sensitivity(new SensitivitySpec { Functional = FunctionalForm.Condition }).Integrated.Variance);
        // A declared order counts the location in the integrated form.
        BasisSobolAnalysis second = fit.Sensitivity(new SensitivitySpec { MaxOrder = 2 });
        Assert.DoesNotContain(second.Integrated.Subsets, u => u.Length == 3);
        Assert.Contains(second.Integrated.Subsets, u => u.SequenceEqual(new[] { 0, 2 }));
        var ex = Assert.Throws<ArgumentException>(() => fit.Sensitivity(new SensitivitySpec { SubsetBudget = 1 }));
        Assert.Contains("budget", ex.Message);
    }

    [Fact]
    public void Refusals_NameTheirReason()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        Assert.Contains("training locations", Assert.Throws<ArgumentException>(() => fit.Sensitivity(new SensitivitySpec { Coordinates = [[0.5]] })).Message);
        Assert.Contains("equally weighted", Assert.Throws<ArgumentException>(() => fit.Sensitivity(new SensitivitySpec
        {
            OutputMeasures = new Dictionary<string, InputMeasure> { ["t"] = new UniformMeasure() },
        })).Message);
        Assert.Throws<ArgumentException>(() => fit.Sensitivity(new SensitivitySpec
        {
            Measures = new Dictionary<string, InputMeasure> { ["w"] = new UniformMeasure() },
        }));
        Assert.Throws<ArgumentOutOfRangeException>(() => fit.Sensitivity(new SensitivitySpec { RoundingTolerance = double.NaN }));
        Assert.Throws<ArgumentOutOfRangeException>(() => fit.Sensitivity(new SensitivitySpec { MaxOrder = 0 }));
        BassBasisFit mixed = BasisSensitivityFixture.Mock(new double[,] { { 1.0, 0.5 }, { 0.2, 0.8 } }, [0.0, 0.0], [1.0, 1.0], degrees: [1, 2]);
        Assert.Contains("degree", Assert.Throws<ArgumentException>(() => mixed.Sensitivity()).Message);
    }

    internal static double MainOf(SobolDecomposition s, int d, int[] subset)
    {
        for (int u = 0; u < s.Subsets.Count; u++)
            if (s.Subsets[u].SequenceEqual(subset)) return s.MainVariance[d, u];
        return 0.0;
    }

    // Var(E[F | the inputs in u]): the mains of u's subsets, summed.
    private static double Closed(SobolDecomposition s, int d, int[] subset)
        => Enumerable.Range(0, s.Subsets.Count).Where(u => s.Subsets[u].All(subset.Contains)).Sum(u => s.MainVariance[d, u]);

    private static void Identities(SobolDecomposition s, int d)
    {
        if (!s.IsDefined(d)) return;
        double v = s.Variance[d], sum = 0.0;
        for (int u = 0; u < s.Subsets.Count; u++) sum += s.MainVariance[d, u];
        Assert.True(Math.Abs(sum - v) <= 1e-12 * v, $"draw {d}: mains sum to {sum}, total {v}");
        for (int i = 0; i < s.Inputs.Count; i++)
        {
            double holding = Enumerable.Range(0, s.Subsets.Count).Where(u => s.Subsets[u].Contains(i)).Sum(u => s.MainVariance[d, u]);
            Assert.True(Math.Abs(s.TotalVariance[d, i] - holding) <= 1e-12 * v, $"draw {d}, input {i}: total {s.TotalVariance[d, i]}, subsets {holding}");
            Assert.True(s.TotalVariance[d, i] >= MainOf(s, d, [i]) - 1e-12 * v);
        }
        Assert.Equal(DecompositionStatus.Defined, s.Status[d] & ~DecompositionStatus.Clamped);
    }

    // The direct decomposition clamps a negative within the tolerance to 0; the basis keeps the raw value and clamps per view.
    private static void AssertSameOrClamped(double direct, double raw, double variance)
    {
        if (direct == raw) return;
        Assert.True(direct == 0.0 && raw < 0.0 && raw >= -1e-10 * variance, $"direct {direct}, raw {raw}");
    }

    private static DrawStore<BassDraw> Store(IEnumerable<BassDraw> draws)
    {
        var store = new DrawStore<BassDraw>();
        foreach (BassDraw draw in draws) store.Add(draw);
        return store;
    }

    private static double Sq(double x) => x * x;

    private static bool Close(double actual, double expected)
        => Math.Abs(actual - expected) <= 1e-13 * Math.Max(Math.Abs(expected), 1e-300) || actual == expected;

    private static bool Near(double actual, double expected, double scale)
        => Math.Abs(actual - expected) <= 1e-10 * Math.Max(Math.Abs(scale), 1e-12);
}

/// <summary>
/// Basis fits with declared terms and draws: two predictors x0, x1 on a grid holding 0, 0.2, 0.3, 0.4, 0.6 and 1 (so the
/// transform is the identity and those are knots), each component run briefly on its coefficient response so the
/// fingerprints hold, then given the declared terms and one chain of draws.
/// </summary>
internal static class BasisSensitivityFixture
{
    public static readonly double[] X = [0.0, 0.2, 0.3, 0.4, 0.6, 1.0];

    private static readonly Dictionary<(int, int), BassFit> Fits = [];

    // Two components by default: w₀ = 0.5 + 1.2·h(x0; +, 0.3) − 0.8·h(x1; −, 0.6) and
    // w₁ = −0.2 + 0.9·h(x0; +, 0.3)·h(x1; +, 0.4) + 0.5·h(x1; +, 0.2), one draw each.
    public static (BassTerm[] Terms, double[][] Draws)[] Default() =>
    [
        ([Hinge(0, 1, 0.3), Hinge(1, -1, 0.6)], [[0.5, 1.2, -0.8]]),
        ([new BassTerm([new HingeFactor(0, 1, Knot(0.3)), new HingeFactor(1, 1, Knot(0.4))]), Hinge(1, 1, 0.2)], [[-0.2, 0.9, 0.5]]),
    ];

    public static BassPredictor[] Predictors() => [new ContinuousPredictor("x0", X), new ContinuousPredictor("x1", [.. X.Reverse()])];

    public static BassBasisFit Mock(double[,] basis, double[] offsets, double[] scales,
                                    (BassTerm[] Terms, double[][] Draws)[]? components = null, int[]? degrees = null)
    {
        components ??= Default();
        int r = components.Length, n = X.Length;
        var w = new double[r, n];
        for (int k = 0; k < r; k++)
            for (int i = 0; i < n; i++) w[k, i] = Math.Sin(1.0 + i + 2.0 * k);
        BassDecomposition decomposition = BassDecomposition.FromBasis(basis, w, offsets, scales, residuals: new double[n, offsets.Length]);
        BassPredictor[] predictors = Predictors();
        var fits = new BassFit[r];
        for (int k = 0; k < r; k++)
        {
            int degree = degrees?[k] ?? 1;
            BassFit run;
            lock (Fits)
            {
                if (!Fits.TryGetValue((k, degree), out run!))
                {
                    var config = new BassConfig
                    {
                        Structure = new BassStructureSpec { Degree = degree, MaxBasis = 2 },
                        Sampling = new BassSampling { Chains = 1, MasterSeed = 5, BurnIn = 20, Samples = 2 },
                    };
                    run = Bass.Fit(config, new BassData(predictors, new ScalarResponse(decomposition.CoefficientResponse(k))));
                    Fits[(k, degree)] = run;
                }
            }
            var table = new BassTermTable();
            int[] ids = [.. components[k].Terms.Select(table.Intern)];
            var store = new DrawStore<BassDraw>();
            foreach (double[] a in components[k].Draws) store.Add(new BassDraw(new BassStructure(ids), 1.0, 1.0, a, 1.0));
            fits[k] = run with { Terms = table, Draws = [store], DrawWeights = null };
        }
        return BassComponents.Assemble(decomposition, predictors, fits, 7);
    }

    public static BassTerm Hinge(int variable, int sign, double knot) => new([new HingeFactor(variable, sign, Knot(knot))]);

    public static int Knot(double value) => Array.IndexOf(X, value);
}
