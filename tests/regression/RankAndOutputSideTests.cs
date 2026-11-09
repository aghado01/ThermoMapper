using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 5 gates (the active-subspace brief; the Grassmann discussion's items, settled by the user on 2026-11-08):
/// <list type="bullet">
/// <item>without a rule, the fixed-k readout unchanged by asking for ranks or output patterns;</item>
/// <item>per-draw k by a declared rule: the rules on chosen spectra; sums of r ridges recover r; tr P̄ against E[k_s]; the rank summary
/// unbeaten by random projections of every rank, and its loss and distances against direct Frobenius norms, on a seeded sweep;</item>
/// <item>the output side G = E[J W_x Jᵀ]: tr G = tr C under the same metrics, on both routes and every metric; a basis fit's G and a
/// direct fit's against quadrature of their Jacobians;</item>
/// <item>the matrix Poincaré check: equality on cos(πx), a deliberately mis-scaled D failing, mass outside the bound infinite; no
/// violation on fitted draws; absent where its law is not a product of uniforms over continuous inputs;</item>
/// <item>basis adequacy: zero inside the span and the known share outside it, exactly; on fits, a truncated basis against a direct fit
/// read against the share in closed form;</item>
/// <item>the angles through LinAlg's routine: within arccos's error of slice 2's above 10⁻⁶ rad, and accurate below where arccos is
/// not.</item>
/// </list>
/// </summary>
public sealed class RankAndOutputSideTests(ITestOutputHelper output)
{
    private static int[] Spread(int count, int k) => [.. Enumerable.Range(0, k).Select(i => (int)((long)i * (count - 1) / (k - 1)))];

    // ── Without a rule ───────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Asking for per-draw ranks and output patterns leaves the fixed-k readout as it is, bit for bit: each draw's C, eigenvalues and
    /// basis, E[C], P̄, their eigensystems and the angles, on both routes. Without them, no ranks and no patterns are read.
    /// </summary>
    [Fact]
    public void NoRule_LeavesTheFixedKReadout()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        foreach (ActiveSubspaceRoute route in new[] { ActiveSubspaceRoute.Exact, ActiveSubspaceRoute.Sampled })
        {
            var spec = new ActiveSubspaceSpec { Draws = Spread(fit.Draws.Sum(c => c.Count), 8), Route = route, Points = 2048 };
            ActiveSubspaceAnalysis plain = fit.ActiveSubspace(data, spec);
            ActiveSubspaceAnalysis asked = fit.ActiveSubspace(data, spec with { Rank = RankRule.TraceShare(0.9), OutputPatterns = 2 });
            Assert.Null(plain.Ranks);
            Assert.Null(plain.Outputs);
            Assert.NotNull(asked.Ranks);
            Assert.NotNull(asked.Outputs);
            for (int d = 0; d < plain.Posterior.Draws.Count; d++)
            {
                Assert.Equal(plain.Posterior.Draws[d].C, asked.Posterior.Draws[d].C);
                Assert.Equal(plain.Posterior.Draws[d].Eigenvalues, asked.Posterior.Draws[d].Eigenvalues);
                Assert.Equal(plain.Posterior.Draws[d].Basis, asked.Posterior.Draws[d].Basis);
                Assert.Equal(plain.Posterior.AnglesToMean[d], asked.Posterior.AnglesToMean[d]);
            }
            Assert.Equal(plain.Posterior.MeanC, asked.Posterior.MeanC);
            Assert.Equal(plain.Posterior.Projection, asked.Posterior.Projection);
            Assert.Equal(plain.Posterior.MeanCEigenvalues, asked.Posterior.MeanCEigenvalues);
            Assert.Equal(plain.Posterior.ProjectionVectors, asked.Posterior.ProjectionVectors);
        }
    }

    // ── Per-draw k ───────────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>The rules on chosen spectra: the trace share's smallest k, the largest gap against λ₁ up to k_max (k = p included, a
    /// trailing zero not a gap of its own), C = 0 at rank 0, and the refusals.</summary>
    [Fact]
    public void RankRules_OnChosenSpectra()
    {
        Assert.Equal(1, RankRule.TraceShare(0.9).Apply([10.0, 0.5, 0.1]));
        Assert.Equal(2, RankRule.TraceShare(0.97).Apply([10.0, 0.5, 0.1]));
        Assert.Equal(3, RankRule.TraceShare(1.0).Apply([10.0, 0.5, 0.1]));
        Assert.Equal(0, RankRule.TraceShare(0.9).Apply([0.0, 0.0, 0.0]));
        Assert.Equal(2, RankRule.Eigengap(3).Apply([5.0, 4.0, 0.01, 0.0]));
        Assert.Equal(1, RankRule.Eigengap(3).Apply([5.0, 0.4, 0.3, 0.0]));
        Assert.Equal(2, RankRule.Eigengap(2).Apply([5.0, 4.0, 0.01]));
        Assert.Equal(3, RankRule.Eigengap(3).Apply([1.0, 1.0, 1.0]));            // k = p, its gap λ_p/λ₁
        Assert.Equal(2, RankRule.Eigengap(4).Apply([5.0, 4.0, 1e-3, 0.0]));      // the trailing zero is no gap of its own
        Assert.Equal(1, RankRule.Eigengap(1).Apply([5.0, 4.0, 0.01]));
        Assert.Equal(1, RankRule.Eigengap(4).Apply([2.0]));
        Assert.Equal(0, RankRule.Eigengap(2).Apply([0.0, 0.0]));
        Assert.Equal(2, RankRule.Eigengap(3).Apply([3.0, 3.0, 1e-20, -1e-18]));   // a negative rounding read as 0
        Assert.Throws<ArgumentOutOfRangeException>(() => RankRule.TraceShare(0.0).Apply([1.0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => RankRule.TraceShare(1.5).Apply([1.0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => RankRule.Eigengap(0).Apply([1.0]));
    }

    /// <summary>
    /// Sums of r ridges recover r: BASS fits of a ridge on [0, 1]³ and of two ridges on [0, 1]⁴ with an inert input (slice 2's recovery
    /// functions, from outside the model) at N = 900, read with the trace-share rule at 0.95: the rank summary is r and the draws' ranks
    /// sit at r. At 0.99 the fits read one more, since each puts a few percent of its gradient energy off the true span (slice 2's
    /// recovery: E[C] within 2.4%); the gap rule against λ₁ reads the two ridges, whose second eigenvalue is 0.37 of the first, as rank
    /// 1, its documented bias. Both are printed.
    /// </summary>
    [Theory]
    [InlineData("ridge", 1)]
    [InlineData("two ridges", 2)]
    public void SumsOfRidges_RecoverTheirRank(string function, int r)
    {
        (BassData data, BassFit fit, int p) = RidgeFit(function);
        ActiveSubspaceAnalysis a = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Rank = RankRule.TraceShare(0.95) });
        RankPosterior ranks = a.Ranks!;
        RankPosterior gap = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Rank = RankRule.Eigengap(p) }).Ranks!;
        RankPosterior strict = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Rank = RankRule.TraceShare(0.99) }).Ranks!;
        output.WriteLine($"{function}: mean rank {ranks.MeanRank:F3}, the summary's rank {ranks.SummaryRank}, P̄'s eigenvalues " +
                         $"[{string.Join(", ", ranks.ProjectionEigenvalues.Select(x => x.ToString("F3")))}], share at r {ranks.Distribution[r]:P1}; " +
                         $"at 0.99 the summary's rank {strict.SummaryRank}, mean {strict.MeanRank:F3}; the gap rule's summary rank {gap.SummaryRank}, mean {gap.MeanRank:F3}");
        Assert.Equal(r, ranks.SummaryRank);
        Assert.True(ranks.Distribution[r] >= 0.9, $"{ranks.Distribution[r]:P1}");
        Assert.Equal(ranks.MeanRank, Trace(ranks.Projection), 12);
    }

    private static readonly Dictionary<string, Lazy<(BassData, BassFit, int)>> RidgeFits = new()
    {
        ["ridge"] = new(() => FitRidges("ridge")),
        ["two ridges"] = new(() => FitRidges("two ridges")),
    };

    private static (BassData Data, BassFit Fit, int P) RidgeFit(string function) => RidgeFits[function].Value;

    private static (BassData, BassFit, int) FitRidges(string function)
    {
        ReferenceFunction f = function == "ridge"
            ? new RidgeSum("sin ridge", new RidgeTerm(ExponentialProfile.Sin(2.5), Unit([0.8, 0.5, 0.33])))
            : new RidgeSum("two ridges", new RidgeTerm(ExponentialProfile.Sin(2.5), Unit([1.0, 0.6, 0.0, 0.0])), new RidgeTerm(ExponentialProfile.Sin(2.0), Unit([0.0, 0.5, 1.0, 0.0])));
        int p = f.Dimension;
        var rng = new Xoshiro256PlusPlus(20261109);
        const int n = 900;
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
            Sampling = new BassSampling { Chains = 2, MasterSeed = 13, BurnIn = 2_000, Samples = 30, Thin = 10 },
        };
        return (data, Bass.Fit(config, data), p);
    }

    private static double[] Unit(double[] w)
    {
        double n = Math.Sqrt(w.Sum(v => v * v));
        return [.. w.Select(v => v / n)];
    }

    /// <summary>
    /// A seeded sweep of 150 mixed-rank posteriors (p from 2 to 6, 3 to 40 draws of random rank and spectrum, random weights, both
    /// rules): tr P̄ = E[k_s]; the summary's loss equals Σ w_s‖P* − P_s‖²_F computed directly and is unbeaten by random projections of
    /// every rank (20 each) and by P̄'s leading k for every k; each draw's rank term plus containment defect is ½‖P_s − P*‖²_F.
    /// </summary>
    [Fact]
    public void RankSummary_IsUnbeaten_SeededSweep()
    {
        var rng = new Xoshiro256PlusPlus(20261110);
        double traceGap = 0.0, lossGap = 0.0, distanceGap = 0.0, worstMargin = double.PositiveInfinity;
        for (int c = 0; c < 150; c++)
        {
            int p = 2 + rng.NextInt(5), draws = 3 + rng.NextInt(38);
            var subspaces = new DrawSubspace[draws];
            for (int s = 0; s < draws; s++)
            {
                int rank = rng.NextInt(p + 1);
                double[,] v = Orthogonal(rng, p);
                var cm = new double[p, p];
                for (int k = 0; k < rank; k++)
                {
                    double lambda = Math.Pow(10.0, -3.0 * rng.NextDouble());
                    for (int i = 0; i < p; i++)
                        for (int j = 0; j < p; j++) cm[i, j] += lambda * v[i, k] * v[j, k];
                }
                subspaces[s] = DrawSubspace.FromC(cm, 1);
            }
            double[] w = [.. Enumerable.Range(0, draws).Select(_ => 0.1 + rng.NextDouble())];
            RankRule rule = c % 2 == 0 ? RankRule.TraceShare(0.5 + 0.5 * rng.NextDouble()) : RankRule.Eigengap(1 + rng.NextInt(p));
            RankPosterior post = RankPosterior.Of(subspaces, w, rule);
            traceGap = Math.Max(traceGap, Math.Abs(Trace(post.Projection) - post.MeanRank));
            double[][,] bases = [.. Enumerable.Range(0, draws).Select(s => Leading(subspaces[s].Vectors, post.Ranks[s]))];
            double[,] star = post.Summary;
            double direct = 0.0;
            for (int s = 0; s < draws; s++)
            {
                double f2 = Frobenius2(Projector(star, p), Projector(bases[s], p));
                direct += post.Weights[s] * f2;
                distanceGap = Math.Max(distanceGap, Math.Abs(post.RankTerms[s] + post.Containment[s] - 0.5 * f2));
            }
            lossGap = Math.Max(lossGap, Math.Abs(direct - post.SummaryLoss));
            for (int k = 0; k <= p; k++)
            {
                worstMargin = Math.Min(worstMargin, post.Loss(Leading(post.ProjectionVectors, k)) - post.SummaryLoss);
                for (int t = 0; t < 20; t++) worstMargin = Math.Min(worstMargin, post.Loss(Leading(Orthogonal(rng, p), k)) - post.SummaryLoss);
            }
        }
        output.WriteLine($"150 posteriors: tr P̄ against E[k] {traceGap:E2}; the summary's loss against Σw‖P* − P_s‖² {lossGap:E2}; rank term plus " +
                         $"containment against ½‖P_s − P*‖² {distanceGap:E2}; the least margin over random and leading projections {worstMargin:E2}");
        Assert.True(traceGap <= 1e-12, $"{traceGap:E2}");
        Assert.True(lossGap <= 1e-11, $"{lossGap:E2}");
        Assert.True(distanceGap <= 1e-11, $"{distanceGap:E2}");
        Assert.True(worstMargin >= -1e-11, $"{worstMargin:E2}");
    }

    // ── The output side ──────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// tr G_s = tr C_s for every draw under the same metrics: a direct functional fit and a basis fit, the box, raw and whitened input
    /// metrics, W_y = I and a full W_y, on both routes.
    /// </summary>
    [Fact]
    public void OutputSide_TraceIsTheTraceOfC()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        BassBasisFit basis = BassOutputComponentsTests.FitOf();
        double worst = 0.0;
        int read = 0;
        foreach (InputMetric metric in new[] { InputMetric.Box, InputMetric.Raw, InputMetric.Whitened })
            foreach (bool full in new[] { false, true })
                foreach (ActiveSubspaceRoute route in new[] { ActiveSubspaceRoute.Exact, ActiveSubspaceRoute.Sampled })
                {
                    var direct = new ActiveSubspaceSpec
                    {
                        Draws = Spread(fit.Draws.Sum(c => c.Count), 4), Metric = metric, Route = route, Points = 1024, OutputPatterns = 2,
                        OutputMetric = full ? new OutputMetric(FunctionalActiveSubspaceTests.RandomMetric(data.Response.Locations, 7)) : null,
                    };
                    var basisSpec = direct with
                    {
                        Draws = Spread(basis.JointDraws.Count, 4),
                        OutputMetric = full ? new OutputMetric(FunctionalActiveSubspaceTests.RandomMetric(basis.Decomposition.Locations, 11)) : null,
                    };
                    foreach (ActiveSubspaceAnalysis a in new[] { fit.ActiveSubspace(data, direct), basis.ActiveSubspace(basisSpec) })
                        for (int d = 0; d < a.Posterior.Draws.Count; d++)
                        {
                            double tc = Trace(a.Posterior.Draws[d].C);
                            worst = Math.Max(worst, Math.Abs(a.Outputs!.Traces[d] - tc) / tc);
                            read++;
                        }
                }
        output.WriteLine($"{read} draws: tr G against tr C {worst:E2}");
        Assert.True(worst <= 1e-12, $"{worst:E2}");
    }

    /// <summary>
    /// G per draw against quadrature of the Jacobian: a basis fit's reconstructed J = G_map·J_w by each component's gradient readout, and
    /// a direct fit's by its own readout, Σ_r w_r J_rJ_rᵀ in the box metric with W_y = I, against FFᵀ from the readout's factor; and E[G]'s
    /// eigenvalues against the weighted mean's.
    /// </summary>
    [Fact]
    public void OutputSide_IsTheQuadratureOfTheJacobian()
    {
        BassBasisFit basis = BassOutputComponentsTests.FitOf();
        int[] joint = Spread(basis.JointDraws.Count, 4);
        ActiveSubspaceAnalysis a = basis.ActiveSubspace(new ActiveSubspaceSpec { Draws = joint, OutputPatterns = 3 });
        BassInputs inputs = BassInputs.Prepare(basis.ComponentData(0), basis.Components[0].Transform);
        (PredictionCases cases, double[] w) = FunctionalActiveSubspaceTests.Tensor(inputs, a.Measure!);
        double worst = 0.0;
        var mean = new double[basis.Decomposition.Locations, basis.Decomposition.Locations];
        for (int d = 0; d < joint.Length; d++)
        {
            double[,] g = OuterSum(FunctionalActiveSubspaceTests.BasisJacobians(basis, cases, joint[d]), w);
            worst = Math.Max(worst, Math.Abs(Trace(g) - a.Outputs!.Traces[d]) / Trace(g));
            double[] values = Eigenvalues(g);
            for (int k = 0; k < a.Outputs.DrawEigenvalues[d].Length; k++) worst = Math.Max(worst, Math.Abs(values[k] - a.Outputs.DrawEigenvalues[d][k]) / values[0]);
            for (int i = 0; i < g.GetLength(0); i++)
                for (int j = 0; j < g.GetLength(1); j++) mean[i, j] += g[i, j] / joint.Length;
        }
        double[] meanValues = Eigenvalues(mean);
        for (int k = 0; k < a.Outputs!.MeanEigenvalues.Length; k++) worst = Math.Max(worst, Math.Abs(meanValues[k] - a.Outputs.MeanEigenvalues[k]) / meanValues[0]);
        Assert.True(a.Outputs.MeanEigenvalues.Length <= basis.Components.Count, "a basis fit's G has rank at most K");

        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        int[] draws = Spread(fit.Draws.Sum(c => c.Count), 4);
        ActiveSubspaceAnalysis direct = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Draws = draws, OutputPatterns = 2 });
        (PredictionCases dc, double[] dw) = FunctionalActiveSubspaceTests.Tensor(BassInputs.Prepare(data, fit.Transform), direct.Measure!);
        BassGradient gradient = fit.Gradient(data, dc, draws);
        for (int d = 0; d < draws.Length; d++)
        {
            double[,] g = OuterSum(FunctionalActiveSubspaceTests.DirectJacobians(gradient, d), dw);
            double[] values = Eigenvalues(g);
            for (int k = 0; k < 2; k++) worst = Math.Max(worst, Math.Abs(values[k] - direct.Outputs!.DrawEigenvalues[d][k]) / values[0]);
        }
        output.WriteLine($"basis and direct G against quadrature: worst {worst:E2}; E[G]'s leading eigenvalues {string.Join(", ", a.Outputs!.MeanEigenvalues.Take(3).Select(x => x.ToString("G4")))}");
        Assert.True(worst <= 1e-10, $"{worst:E2}");
    }

    // Σ_r w_r J_rJ_rᵀ.
    private static double[,] OuterSum(double[][,] jacobians, double[] w)
    {
        int m = jacobians[0].GetLength(0), p = jacobians[0].GetLength(1);
        var g = new double[m, m];
        for (int r = 0; r < jacobians.Length; r++)
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                {
                    double s = 0.0;
                    for (int k = 0; k < p; k++) s += jacobians[r][i, k] * jacobians[r][j, k];
                    g[i, j] += w[r] * s;
                }
        return g;
    }

    // ── The matrix Poincaré check ────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// The matrix check on reference fields by quadrature: f = (cos πx₁, cos πx₂) on [0, 1]², the Poincaré inequality's extremal case,
    /// at ratio 1; f = (cos πx₁ + x₂², x₁x₂) below 1; the bound's D halved, which fails at 2 on the cosines; and covariance outside the
    /// bound's range, infinite.
    /// </summary>
    [Fact]
    public void MatrixPoincare_EqualityOnCosines_AndAMisScaledBoundFails()
    {
        (double[,] cov, double[,] bound) Moments(Func<double, double, double[]> f, Func<double, double, double[,]> jac, double scale)
        {
            QuadratureRule rule = QuadratureRule.Uniform(0.0, 1.0, 20, 4);
            int m = f(0.5, 0.5).Length;
            var mean = new double[m];
            var second = new double[m, m];
            var b = new double[m, m];
            for (int i = 0; i < rule.Nodes.Length; i++)
                for (int j = 0; j < rule.Nodes.Length; j++)
                {
                    double x1 = rule.Nodes[i], x2 = rule.Nodes[j], w = rule.Weights[i] * rule.Weights[j];
                    double[] y = f(x1, x2);
                    double[,] jj = jac(x1, x2);
                    for (int a = 0; a < m; a++)
                    {
                        mean[a] += w * y[a];
                        for (int c = 0; c < m; c++)
                        {
                            second[a, c] += w * y[a] * y[c];
                            b[a, c] += w * scale / (Math.PI * Math.PI) * (jj[a, 0] * jj[c, 0] + jj[a, 1] * jj[c, 1]);
                        }
                    }
                }
            var cv = new double[m, m];
            for (int a = 0; a < m; a++)
                for (int c = 0; c < m; c++) cv[a, c] = second[a, c] - mean[a] * mean[c];
            return (cv, b);
        }
        double pi = Math.PI;
        Func<double, double, double[]> cosines = (x1, x2) => [Math.Cos(pi * x1), Math.Cos(pi * x2)];
        Func<double, double, double[,]> cosineJ = (x1, x2) => new[,] { { -pi * Math.Sin(pi * x1), 0.0 }, { 0.0, -pi * Math.Sin(pi * x2) } };
        (double[,] c1, double[,] b1) = Moments(cosines, cosineJ, 1.0);
        (double[,] c2, double[,] b2) = Moments(cosines, cosineJ, 0.5);
        Func<double, double, double[]> mixed = (x1, x2) => [Math.Cos(pi * x1) + x2 * x2, x1 * x2];
        Func<double, double, double[,]> mixedJ = (x1, x2) => new[,] { { -pi * Math.Sin(pi * x1), 2 * x2 }, { x2, x1 } };
        (double[,] c3, double[,] b3) = Moments(mixed, mixedJ, 1.0);
        double equal = MatrixPoincare.Ratio(c1, b1), halved = MatrixPoincare.Ratio(c2, b2), below = MatrixPoincare.Ratio(c3, b3);
        double outside = MatrixPoincare.Ratio(new[,] { { 1.0, 0.0 }, { 0.0, 1.0 } }, new[,] { { 1.0, 0.0 }, { 0.0, 0.0 } });
        output.WriteLine($"cosines {equal:F15}; D halved {halved:F12}; a mixed field {below:F6}; outside the bound's range {outside}");
        Assert.Equal(1.0, equal, 12);
        Assert.Equal(2.0, halved, 11);
        Assert.True(below < 1.0, $"{below}");
        Assert.True(double.IsPositiveInfinity(outside));
    }

    /// <summary>
    /// On fitted draws under the uniform law: no violation for a scalar fit, a direct functional fit under a full W_y and a basis fit;
    /// absent for a fit with a categorical input, under a non-uniform law and under a joint sample.
    /// </summary>
    [Fact]
    public void MatrixPoincare_HoldsOnFittedDraws()
    {
        (BassData df, _) = CoActiveFixture.Data();
        (BassFit f, _) = CoActiveFixture.Fits();
        ActiveSubspaceAnalysis scalar = f.ActiveSubspace(df);
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        ActiveSubspaceAnalysis direct = fit.ActiveSubspace(data, new ActiveSubspaceSpec { OutputMetric = new OutputMetric(FunctionalActiveSubspaceTests.RandomMetric(data.Response.Locations, 13)) });
        ActiveSubspaceAnalysis basis = BassOutputComponentsTests.FitOf().ActiveSubspace();
        double worst = new[] { scalar, direct, basis }.SelectMany(a => a.MatrixPoincare!.Ratios).Max();
        output.WriteLine($"worst ratio {worst:F4} over {scalar.MatrixPoincare!.Ratios.Count + direct.MatrixPoincare!.Ratios.Count + basis.MatrixPoincare!.Ratios.Count} draws");
        Assert.Equal(0, scalar.MatrixPoincare.Violations + direct.MatrixPoincare.Violations + basis.MatrixPoincare.Violations);
        Assert.Null(SensitivityFixture.FitScalar(1).ActiveSubspace(SensitivityFixture.Scalar()).MatrixPoincare);
        Assert.Null(f.ActiveSubspace(df, new ActiveSubspaceSpec { Measures = new Dictionary<string, InputMeasure> { ["x0"] = new MixtureMeasure([new NormalComponent(1.0, 0.5, 0.2)]) } }).MatrixPoincare);
        var rng = new Xoshiro256PlusPlus(3);
        double[] u = [.. Enumerable.Range(0, 500).Select(_ => rng.NextDouble())], v = [.. Enumerable.Range(0, 500).Select(_ => rng.NextDouble())];
        var sample = new JointSample(new PredictionCases([new ContinuousColumn("x0", u), new ContinuousColumn("x1", v)]));
        Assert.Null(f.ActiveSubspace(df, new ActiveSubspaceSpec { Sample = sample }).MatrixPoincare);
    }

    // ── Basis adequacy ───────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Adequacy on chosen factors: structure inside the span reads 0; F = [a·e₁, b·e₄] against span(e₁, e₂) reads the share b²/(a² + b²),
    /// the leading pattern contained, the second not (defect 1); a span given by a non-orthonormal basis reads the same.
    /// </summary>
    [Fact]
    public void BasisAdequacy_ZeroInside_TheKnownShareOutside()
    {
        const int m = 6;
        double[,] span = new double[m, 2], skewed = new double[m, 2];
        span[0, 0] = span[1, 1] = 1.0;
        skewed[0, 0] = 2.0;
        skewed[0, 1] = 1.0;
        skewed[1, 1] = 3.0;
        var inside = new double[m, 3];
        inside[0, 0] = 1.0;
        inside[1, 1] = 0.5;
        inside[0, 2] = 0.3;
        inside[1, 2] = -0.2;
        const double a = 2.0, b = 0.7;
        var mixed = new double[m, 2];
        mixed[0, 0] = a;
        mixed[3, 1] = b;
        BasisAdequacy zero = BasisAdequacy.Of([inside], [1.0], span, 2);
        BasisAdequacy known = BasisAdequacy.Of([mixed], [1.0], span, 1);
        BasisAdequacy second = BasisAdequacy.Of([mixed], [1.0], skewed, 2);
        output.WriteLine($"inside: share {zero.Shares[0]:E2}, defect {zero.Containment[0]:E2}; outside: share {known.Shares[0]:F15} against {b * b / (a * a + b * b):F15}");
        Assert.True(zero.Shares[0] <= 1e-15 && zero.Containment[0] <= 1e-15);
        Assert.Equal(b * b / (a * a + b * b), known.Shares[0], 15);
        Assert.Equal(0.0, known.Containment[0], 15);
        Assert.Equal(1.0, second.Containment[0], 14);
        Assert.Equal(known.MeanShare, second.MeanShare, 14);
        Assert.Equal(2, second.SpanRank);
    }

    /// <summary>
    /// On fits: curves y(x, t) = a(x)φ₁(t) + b(x)φ₂(t) + 0.01 noise, φ₁ and φ₂ orthonormal over 12 locations, a = 2x₁ + x₂ and
    /// b = 0.7x₁ − x₂ linear. A direct-mode fit beside basis fits of rank 1 and 2. Under the box metric and the uniform law the truth is
    /// G = 5φ₁φ₁ᵀ + 1.49φ₂φ₂ᵀ + 0.4(φ₁φ₂ᵀ + φ₂φ₁ᵀ), so its share outside a span is 1 − tr(P G)/tr G, here read for each fitted basis's own
    /// span (1.49/6.49 = 0.230 for span(φ₁) exactly). The readout meets it within 0.03 at both ranks, the direct fit's hinges in t
    /// approximating φ₁ and φ₂.
    /// </summary>
    [Fact]
    public void BasisAdequacy_OnFits_ReadsTheDroppedShare()
    {
        (BassData data, BassFit direct, BassPredictor[] predictors, double[,] y, double[] phi1, double[] phi2) = AdequacyFixture.Value;
        int m = phi1.Length;
        var gTrue = new double[m, m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++) gTrue[i, j] = 5.0 * phi1[i] * phi1[j] + 1.49 * phi2[i] * phi2[j] + 0.4 * (phi1[i] * phi2[j] + phi2[i] * phi1[j]);
        double worst = 0.0;
        for (int rank = 1; rank <= 2; rank++)
        {
            BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { Rank = rank });
            BassBasisFit basis = BassComponents.Fit(d, predictors, BassOutputComponentsTests.Config());
            BasisAdequacy adequacy = basis.Adequacy(direct, data, new ActiveSubspaceSpec { Draws = Spread(direct.Draws.Sum(c => c.Count), 10), OutputPatterns = 1 });
            var span = new double[m, rank];
            for (int l = 0; l < m; l++)
                for (int k = 0; k < rank; k++) span[l, k] = d.Map.Scales[l] * d.Map.Basis(l, k);
            double truth = BasisAdequacy.Of([Root(gTrue)], [1.0], span, 1).MeanShare;
            worst = Math.Max(worst, Math.Abs(adequacy.MeanShare - truth));
            output.WriteLine($"rank {rank}: E[G]'s share outside the span {adequacy.MeanShare:F4}, the truth's for the same span {truth:F4} " +
                             $"(for span(φ₁) {1.49 / 6.49:F4}); draws {adequacy.Shares.Min():F4} to {adequacy.Shares.Max():F4}; the leading pattern's defect {adequacy.MeanContainment:E2}");
        }
        Assert.True(worst <= 0.03, $"{worst:F4}");
        Assert.Throws<ArgumentException>(() => BassOutputComponentsTests.FitOf().Adequacy(direct, data));
    }

    // A factor F with FFᵀ = G for a symmetric positive semidefinite G.
    private static double[,] Root(double[,] g)
    {
        EigenResult e = DenseEigen.DecomposeSymmetric(g);
        int m = g.GetLength(0);
        var f = new double[m, m];
        for (int q = 0; q < m; q++)
        {
            double s = Math.Sqrt(Math.Max(0.0, e.Eigenvalues[q]));
            for (int l = 0; l < m; l++) f[l, q] = s * e.Eigenvectors[q][l];
        }
        return f;
    }

    internal static readonly Lazy<(BassData, BassFit, BassPredictor[], double[,], double[], double[])> AdequacyFixture = new(() =>
    {
        var rng = new Xoshiro256PlusPlus(20261111);
        const int n = 60, m = 12;
        double[] t = [.. Enumerable.Range(0, m).Select(l => l / (m - 1.0))];
        double[] phi1 = Normalize([.. t.Select(v => Math.Cos(Math.PI * v))]), phi2Raw = [.. t.Select(v => Math.Sin(Math.PI * v))];
        double dot = phi1.Zip(phi2Raw, (u, v) => u * v).Sum();
        double[] phi2 = Normalize([.. phi2Raw.Select((v, l) => v - dot * phi1[l])]);
        double[] x1 = new double[n], x2 = new double[n];
        var y = new double[n, m];
        for (int i = 0; i < n; i++)
        {
            x1[i] = (i % 10) / 9.0;
            x2[i] = (i * 7 % 11) / 10.0;
            double a = 2 * x1[i] + x2[i], b = 0.7 * x1[i] - x2[i];
            for (int l = 0; l < m; l++) y[i, l] = a * phi1[l] + b * phi2[l] + 0.01 * Normal.SampleStandard(rng);
        }
        BassPredictor[] predictors = [new ContinuousPredictor("x1", x1), new ContinuousPredictor("x2", x2)];
        var data = new BassData(predictors, new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
        BassConfig config = SensitivityFixture.Config(maxBasis: 30) with
        {
            Sampling = new BassSampling { Chains = 2, MasterSeed = 17, BurnIn = 2_000, Samples = 30, Thin = 10 },
        };
        return (data, Bass.Fit(config, data), predictors, y, phi1, phi2);
    });

    private static double[] Normalize(double[] v)
    {
        double n = Math.Sqrt(v.Sum(x => x * x));
        return [.. v.Select(x => x / n)];
    }

    // ── The angles ───────────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// The readout's angles through LinAlg's routine against slice 2's arccos of the cosines: on 300 random pairs of subspaces (p to 8,
    /// k to p, mixed ranks), every angle above 10⁻⁶ rad within arccos's own error, 4√(pk)·ε/sin θ; and on rotations by a known angle
    /// from 10⁻¹² to 10⁻⁶ rad, the new angles within 10⁻¹⁴ relative where arccos returns its floor.
    /// </summary>
    [Fact]
    public void Angles_ThroughLinAlg_WithinArccosError()
    {
        var rng = new Xoshiro256PlusPlus(20261112);
        double worst = 0.0;
        for (int c = 0; c < 300; c++)
        {
            int p = 2 + rng.NextInt(7), ka = 1 + rng.NextInt(p), kb = 1 + rng.NextInt(p);
            double[,] a = Leading(Orthogonal(rng, p), ka), b = Leading(Orthogonal(rng, p), kb);
            double[] angles = SubspacePosterior.PrincipalAngles(a, b), old = ArccosAngles(a, b);
            for (int i = 0; i < angles.Length; i++)
            {
                if (old[i] <= 1e-6) continue;
                double tolerance = 4.0 * Math.Sqrt(p * Math.Max(ka, kb)) * 2.220446049250313e-16 / Math.Sin(old[i]);
                worst = Math.Max(worst, Math.Abs(angles[i] - old[i]) / tolerance);
            }
        }
        double small = 0.0, arccosFloor = double.PositiveInfinity;
        foreach (double theta in new[] { 1e-12, 1e-10, 1e-8, 1e-6 })
        {
            double[,] e = { { 1.0 }, { 0.0 }, { 0.0 } }, r = { { Math.Cos(theta) }, { Math.Sin(theta) }, { 0.0 } };
            small = Math.Max(small, Math.Abs(SubspacePosterior.PrincipalAngles(e, r)[0] - theta) / theta);
            if (theta <= 1e-8) arccosFloor = Math.Min(arccosFloor, Math.Abs(ArccosAngles(e, r)[0] - theta) / theta);
        }
        output.WriteLine($"300 pairs: worst |new − arccos| over arccos's error {worst:F3}; known small angles, relative error {small:E2} (arccos's at most " +
                         $"10⁻⁸: {arccosFloor:E2} at least)");
        Assert.True(worst <= 1.0, $"{worst:F3}");
        Assert.True(small <= 1e-14, $"{small:E2}");
    }

    // Slice 2's angles: the arccos of the singular values of aᵀb.
    private static double[] ArccosAngles(double[,] a, double[,] b)
    {
        int p = a.GetLength(0), ka = a.GetLength(1), kb = b.GetLength(1);
        var m = new double[ka, kb];
        for (int i = 0; i < ka; i++)
            for (int j = 0; j < kb; j++)
                for (int r = 0; r < p; r++) m[i, j] += a[r, i] * b[r, j];
        double[] cosines = ThinSvd.Decompose(m).SingularValues;
        return [.. Enumerable.Range(0, Math.Min(ka, kb)).Select(i => Math.Acos(Math.Clamp(i < cosines.Length ? cosines[i] : 0.0, 0.0, 1.0)))];
    }

    // ── Helpers ──────────────────────────────────────────────────────────────────────────────────────────────────────

    private static double Trace(double[,] m)
    {
        double t = 0.0;
        for (int i = 0; i < m.GetLength(0); i++) t += m[i, i];
        return t;
    }

    private static double[] Eigenvalues(double[,] m) => [.. DenseEigen.DecomposeSymmetric(m).Eigenvalues.OrderByDescending(x => x)];

    private static double[,] Leading(double[,] vectors, int k)
    {
        int p = vectors.GetLength(0);
        var basis = new double[p, k];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < k; j++) basis[i, j] = vectors[i, j];
        return basis;
    }

    private static double[,] Projector(double[,] basis, int p)
    {
        var proj = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
                for (int k = 0; k < basis.GetLength(1); k++) proj[i, j] += basis[i, k] * basis[j, k];
        return proj;
    }

    private static double Frobenius2(double[,] a, double[,] b)
    {
        double s = 0.0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) s += (a[i, j] - b[i, j]) * (a[i, j] - b[i, j]);
        return s;
    }

    private static double[,] Orthogonal(Xoshiro256PlusPlus rng, int p)
    {
        var q = new double[p, p];
        for (int k = 0; k < p; k++)
        {
            var v = new double[p];
            for (int i = 0; i < p; i++) v[i] = Normal.SampleStandard(rng);
            for (int pass = 0; pass < 2; pass++)
                for (int j = 0; j < k; j++)
                {
                    double d = 0.0;
                    for (int i = 0; i < p; i++) d += v[i] * q[i, j];
                    for (int i = 0; i < p; i++) v[i] -= d * q[i, j];
                }
            double n = Math.Sqrt(v.Sum(x => x * x));
            for (int i = 0; i < p; i++) q[i, k] = v[i] / n;
        }
        return q;
    }
}
