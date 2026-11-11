using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E34's slice 5, the PCA mode's joint option (decision 8): a decomposition's scores fitted together as several outputs on one
/// structure (<see cref="BassComponents.FitJoint"/>), each component a view of the joint fit.
/// <list type="bullet">
/// <item><b>Rank 1</b> fits as a component fitted apart, draw for draw: structures, τ, σ² = Σ₀₀ and the coefficients, then the
/// prediction (latent, predictive and mean surface), the pointwise density in both spaces and the emulator's noise.</item>
/// <item><b>Equivariance on the scores:</b> a basis K and scores W against KL⁻ᵀ and LᵀW, the same span mixed. Under decision 2's Σ
/// prior the structure and τ chains are the same draw for draw, and so are what reads them alone: the Rao–Blackwellized mean
/// surface, and each draw's matrix-t conditional, whose location goes to ·L and output scale to Lᵀ·L.</item>
/// <item><b>Independent scores:</b> on surfaces outside the model with independent noise, the joint fit meets components fitted
/// apart: the components' posterior mean surfaces, their noise variances, and Σ's correlations near 0.</item>
/// <item><b>The readers of the noise</b> under correlated scores read Σ whole: the prediction's covariances and moments, its sampled
/// noise, the pointwise density in coefficient and y space against dense Gaussians, and the emulator's noise; a unit splitting a
/// case's components is refused.</item>
/// <item><b>Refusals:</b> random intercepts, and a joint fit run on other data.</item>
/// <item><b>The cost row</b> (qualification, indicative): the joint fit against the components fitted apart at r = 2, 3 and 5.</item>
/// </list>
/// </summary>
public sealed class MultiOutputBasisTests(ITestOutputHelper output)
{
    // ── Rank 1 ────────────────────────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void RankOne_TheJointFitIsTheComponentFit_DrawForDraw()
    {
        (double[,] y, BassPredictor[] predictors) = BassOutputComponentsTests.Problem();
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 1, ScaleByColumnSd = true });
        BassConfig config = BassOutputComponentsTests.Config();
        BassBasisFit apart = BassComponents.Fit(d, predictors, config);
        BassBasisFit joint = BassComponents.FitJoint(d, predictors, config);
        Assert.Null(apart.Joint);
        Assert.NotNull(joint.Joint);

        BassFit a = apart.Components[0], j = joint.Components[0];
        Assert.Equal(a.Draws.Count, j.Draws.Count);
        double worst = 0.0;
        for (int c = 0; c < a.Draws.Count; c++)
        {
            Assert.Equal(a.Draws[c].Count, j.Draws[c].Count);
            for (int t = 0; t < a.Draws[c].Count; t++)
            {
                BassDraw p = a.Draws[c][t], q = j.Draws[c][t];
                Assert.Equal(p.Structure.TermIds.Select(id => a.Terms[id]), q.Structure.TermIds.Select(id => j.Terms[id]));
                worst = Math.Max(worst, Relative(p.Tau, q.Tau));
                worst = Math.Max(worst, Relative(p.Sigma2, q.Sigma2));
                Assert.Equal(p.Coefficients.Length, q.Coefficients.Length);
                for (int k = 0; k < p.Coefficients.Length; k++) worst = Math.Max(worst, Math.Abs(p.Coefficients[k] - q.Coefficients[k]) / (1.0 + Math.Abs(p.Coefficients[k])));
            }
        }
        output.WriteLine($"rank 1, draws: worst relative difference {worst:G3}");
        Assert.True(worst < 1e-12, $"draws differ by {worst:G3}");

        var cases = new PredictionCases([new ContinuousColumn("z0", [0.1, 0.45, 0.9]), new ContinuousColumn("z1", [0.2, 0.7, 0.5])]);
        var options = new BassPredictOptions { Predictive = true, MeanSurface = true, Seed = 5 };
        BassBasisPrediction pa = apart.Predict(cases, options), pj = joint.Predict(cases, options);
        Assert.Null(pa.NoiseCovariances);
        Assert.NotNull(pj.NoiseCovariances);
        double latent = MaxDifference(pa.Latent, pj.Latent), predictive = MaxDifference(pa.Predictive!, pj.Predictive!);
        double mean = MaxDifference(pa.MeanSurface!, pj.MeanSurface!), covariance = MaxDifference(pa.Covariance(1), pj.Covariance(1));
        output.WriteLine($"rank 1, prediction: latent {latent:G3}, predictive {predictive:G3}, mean surface {mean:G3}, covariance {covariance:G3}");
        Assert.True(Math.Max(Math.Max(latent, predictive), Math.Max(mean, covariance)) < 1e-10);

        foreach (BasisDensitySpace space in new[] { BasisDensitySpace.Coefficients, BasisDensitySpace.ResponseWithDiagonalTruncation })
        {
            double density = MaxDifference(apart.PointwiseLogDensity(space: space).Values, joint.PointwiseLogDensity(space: space).Values);
            output.WriteLine($"rank 1, density in {space}: {density:G3}");
            Assert.True(density < 1e-9, $"{space}: {density:G3}");
        }

        ContinuousParameter z0 = new("z0", 0.0, 1.0), z1 = new("z1", 0.0, 1.0);
        var ea = new BassBasisEmulator(apart, [z0, z1]);
        var ej = new BassBasisEmulator(joint, [z0, z1]);
        for (int s = 0; s < ea.DrawCount; s += 17)
            Assert.True(MaxDifference(Dense((ComponentNoise)ea.Noise(s)), Dense((ComponentNoise)ej.Noise(s))) < 1e-12);
    }

    // ── Equivariance on the scores ────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Equivariance_TheScoresSpanDecides_NotTheBasisThatSpansIt()
    {
        // Continuous inputs: on a grid, hinges at different knots can span one column space, an exact tie in the target that
        // rounding breaks differently in each basis (a uniform drawn or not at log ratio ±10⁻¹⁴), which a draw-for-draw gate
        // cannot tell from a defect. Continuous inputs leave no such tie, and the target's own equivariance is read at every state.
        (BassDecomposition source, BassPredictor[] predictors, _, _) = Scores(60, Diagonal([0.1, 0.15, 0.12]), 2026_1113);
        int r = source.Rank, n = source.Cases, m = source.Locations;
        double[,] basis = source.Map.BasisMatrix(), residuals = source.Truncation.Values;
        double[] offsets = source.Map.Offsets.ToArray(), scales = source.Map.Scales.ToArray();
        var w = new double[r, n];
        for (int k = 0; k < r; k++)
            for (int i = 0; i < n; i++) w[k, i] = source.Coefficient(k, i);
        double[,] l = { { 1.3, 0.5, -0.2 }, { -0.4, 0.8, 0.3 }, { 0.7, 0.2, 2.1 } };
        double lnDet = Math.Log(Math.Abs(l[0, 0] * (l[1, 1] * l[2, 2] - l[1, 2] * l[2, 1]) - l[0, 1] * (l[1, 0] * l[2, 2] - l[1, 2] * l[2, 0])
                                         + l[0, 2] * (l[1, 0] * l[2, 1] - l[1, 1] * l[2, 0])));
        double[,] lInverse = Inverse3(l);
        var mixedBasis = new double[m, r];   // K L⁻ᵀ
        for (int j = 0; j < m; j++)
            for (int k = 0; k < r; k++)
                for (int q = 0; q < r; q++) mixedBasis[j, k] += basis[j, q] * lInverse[k, q];
        var mixedScores = new double[r, n];   // Lᵀ W
        for (int k = 0; k < r; k++)
            for (int i = 0; i < n; i++)
                for (int q = 0; q < r; q++) mixedScores[k, i] += l[q, k] * w[q, i];
        BassDecomposition plain = BassDecomposition.FromBasis(basis, w, offsets, scales, residuals: residuals);
        BassDecomposition mixed = BassDecomposition.FromBasis(mixedBasis, mixedScores, offsets, scales, residuals: residuals);

        BassConfig config = BassRunTests.Config(chains: 2, burn: 200, samples: 120, seed: 53, maxBasis: 6);
        BassBasisFit a = BassComponents.FitJoint(plain, predictors, config), b = BassComponents.FitJoint(mixed, predictors, config);
        BassFit ja = a.Joint!, jb = b.Joint!;
        BassModel pm = BassModel.Create(ja.Config, a.JointData()), qm = BassModel.Create(jb.Config, b.JointData());
        double LogMarginal(BassModel model, IReadOnlyList<BassTerm> terms, double tauValue)
        {
            BassDesign design = model.NewDesign();
            Assert.Null(design.Load(terms));
            return model.Target.Evaluate(BassStructureFit.Of(design, 0.0), Math.Log(tauValue)).LogMarginal;
        }

        // The target at every state the first fit visits: ln p(YL) − ln p(Y) = −N ln|det L| (level centring counts N).
        int draws = 0;
        double tau = 0.0, shift = 0.0;
        for (int c = 0; c < ja.Draws.Count; c++)
            for (int t = 0; t < ja.Draws[c].Count; t++, draws++)
            {
                BassDraw p = ja.Draws[c][t], q = jb.Draws[c][t];
                List<BassTerm> terms = [.. p.Structure.TermIds.Select(id => ja.Terms[id])], other = [.. q.Structure.TermIds.Select(id => jb.Terms[id])];
                double vp = LogMarginal(pm, terms, p.Tau), vq = LogMarginal(qm, terms, p.Tau);
                shift = Math.Max(shift, Math.Abs(vq - vp + n * lnDet) / (1.0 + Math.Abs(vp)));
                if (!terms.SequenceEqual(other))
                {
                    double tie = Math.Abs(LogMarginal(pm, other, q.Tau) - vp);
                    Assert.Fail($"chain {c}, draw {t}: the structures part, the marginals {(tie < 1e-9 * (1.0 + Math.Abs(vp)) ? "tied (rounding broke an exact tie)" : $"apart by {tie:G3}")}");
                }
                tau = Math.Max(tau, Relative(p.Tau, q.Tau));
            }
        output.WriteLine($"equivariance: {draws} draws, structures equal, τ to {tau:G3}; the target's shift −N ln|det L| at every state to {shift:G3}");
        Assert.True(tau < 1e-9, $"τ differs by {tau:G3}");
        Assert.True(shift < 1e-12, $"the target's shift misses −N ln|det L| by {shift:G3}");

        (PredictionCases cases, _) = HeldOut(6, 19);
        var options = new BassPredictOptions { MeanSurface = true };
        double surface = MaxDifference(a.Predict(cases, options).MeanSurface!, b.Predict(cases, options).MeanSurface!);
        output.WriteLine($"equivariance: Rao–Blackwellized mean surface to {surface:G3}");
        Assert.True(surface < 1e-8, $"mean surface differs by {surface:G3}");

        BassMatrixTPrediction ta = BassMultiOutput.Conditionals(ja, a.JointData(), cases), tb = BassMultiOutput.Conditionals(jb, b.JointData(), cases);
        Assert.Equal(ta.Components.Count, tb.Components.Count);
        double location = 0.0, columns = 0.0, rows = 0.0;
        for (int s = 0; s < ta.Components.Count; s++)
        {
            double[,] la = ta.Components[s].Location, lb = tb.Components[s].Location;
            double[,] sa = ta.Components[s].ColumnScale, sb = tb.Components[s].ColumnScale;
            location = Math.Max(location, MaxDifference(Times(la, l), lb) / (1.0 + MaxAbs(lb)));
            columns = Math.Max(columns, MaxDifference(Congruent(sa, l), sb) / MaxAbs(sb));
            rows = Math.Max(rows, MaxDifference(ta.Components[s].RowScale, tb.Components[s].RowScale));
            Assert.Equal(ta.Components[s].Nu, tb.Components[s].Nu, 12);
        }
        output.WriteLine($"equivariance: conditionals' location {location:G3}, output scale {columns:G3}, row scale {rows:G3}");
        Assert.True(Math.Max(location, columns) < 1e-8 && rows < 1e-9);
    }

    // ── Independent scores ────────────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void IndependentScores_TheJointFitMeetsComponentsFittedApart()
    {
        double[] sd = [0.10, 0.15, 0.12];
        (BassDecomposition d, BassPredictor[] predictors, _, _) = Scores(200, Diagonal(sd), 2026_1110);
        BassConfig config = BassRunTests.Config(chains: 2, burn: 400, samples: 300, seed: 41, maxBasis: 12);
        BassBasisFit apart = BassComponents.Fit(d, predictors, config), joint = BassComponents.FitJoint(d, predictors, config);

        (PredictionCases cases, double[,] truth) = HeldOut(120, 77);
        double[,] ma = MeanCoefficients(apart.Predict(cases)), mj = MeanCoefficients(joint.Predict(cases));
        BassOutputCovarianceSummary sigma = BassMultiOutput.Covariance(joint.Joint!, joint.JointData());
        for (int k = 0; k < 3; k++)
        {
            double between = Rmse(ma, mj, k), toTruthApart = Rmse(ma, truth, k), toTruthJoint = Rmse(mj, truth, k);
            double sigma2Apart = apart.Components[k].Draws.SelectMany(c => c).Average(x => x.Sigma2);
            double ratio = sigma.Mean![k, k] / sigma2Apart;
            output.WriteLine($"component {k}: RMSE joint−apart {between:F4}, to truth apart {toTruthApart:F4} joint {toTruthJoint:F4} (noise sd {sd[k]}); " +
                             $"E[Σ_kk] {sigma.Mean[k, k]:F5} / E[σ²] apart {sigma2Apart:F5} = {ratio:F3}");
            Assert.True(between < 0.5 * sd[k], $"component {k}: the two fits' surfaces differ by {between:F4}");
            Assert.True(Math.Abs(toTruthJoint - toTruthApart) < 0.25 * sd[k], $"component {k}: errors {toTruthJoint:F4} and {toTruthApart:F4}");
            Assert.InRange(ratio, 0.85, 1.18);
        }
        for (int k = 0; k < 3; k++)
            for (int q = 0; q < k; q++)
            {
                PosteriorSummary rho = sigma.CorrelationEntries[k, q];
                output.WriteLine($"ρ[{k},{q}] {rho.Mean:F4} ± {rho.Sd:F4}");
                Assert.True(Math.Abs(rho.Mean) < 3.0 * rho.Sd + 0.02, $"ρ[{k},{q}] = {rho.Mean:F4} ± {rho.Sd:F4} on independent noise");
            }
    }

    // ── The readers of the noise ──────────────────────────────────────────────────────────────────────────────────────────────

    private static readonly double[,] CorrelatedNoise = { { 0.0100, 0.0072, 0.0 }, { 0.0072, 0.0144, -0.0050 }, { 0.0, -0.0050, 0.0225 } };

    private static readonly Lazy<(BassBasisFit Fit, BassDecomposition D, double[][] X)> Correlated = new(() =>
    {
        (BassDecomposition d, BassPredictor[] predictors, _, double[][] x) = Scores(120, CorrelatedNoise, 2026_1111);
        return (BassComponents.FitJoint(d, predictors, BassRunTests.Config(chains: 2, burn: 300, samples: 150, seed: 43, maxBasis: 10)), d, x);
    });

    [Fact]
    public void CorrelatedScores_TheJointFitRecoversSigma()
    {
        BassBasisFit fit = Correlated.Value.Fit;
        BassOutputCovarianceSummary sigma = BassMultiOutput.Covariance(fit.Joint!, fit.JointData());
        for (int k = 0; k < 3; k++)
            for (int q = 0; q <= k; q++)
            {
                PosteriorSummary entry = sigma.CovarianceEntries[k, q];
                output.WriteLine($"Σ[{k},{q}] {entry.Mean:F5} ± {entry.Sd:F5} against {CorrelatedNoise[k, q]:F5}");
                Assert.True(Math.Abs(entry.Mean - CorrelatedNoise[k, q]) < 4.0 * entry.Sd + 0.002, $"Σ[{k},{q}]");
            }
    }

    [Fact]
    public void CorrelatedScores_ThePredictionReadsSigmaWhole()
    {
        BassBasisFit fit = Correlated.Value.Fit;
        (PredictionCases cases, _) = HeldOut(40, 5);
        BassBasisPrediction p = fit.Predict(cases, new BassPredictOptions { Predictive = true, BasisParts = BasisPredictiveParts.Noise, Seed = 9 });
        double[,,] sigma = p.NoiseCovariances!;
        BassDraw[] pooled = [.. fit.Joint!.Draws.SelectMany(c => c)];
        int s = p.DrawCount, r = 3, m = fit.Decomposition.Locations;
        for (int d = 0; d < s; d++)
        {
            double[,] own = pooled[fit.JointDraws.Draw(0, p.Draws[d])].Sigma!;
            for (int k = 0; k < r; k++)
            {
                Assert.Equal(own[k, k], p.NoiseVariances[d, k]);
                for (int q = 0; q < r; q++) Assert.Equal(own[k, q], sigma[d, k, q]);
            }
        }

        // The moments: diag(y_s)·K(Cov[w] + E[Σ])Kᵀ·diag(y_s), built here from the draws.
        double[,] basis = fit.Decomposition.Map.BasisMatrix();
        double[] scales = fit.Decomposition.Map.Scales.ToArray();
        for (int i = 0; i < 3; i++)
        {
            var inner = new double[r, r];
            var mean = new double[r];
            for (int d = 0; d < s; d++)
                for (int k = 0; k < r; k++) mean[k] += p.Weights[d] * p.Coefficients[d, i, k];
            for (int d = 0; d < s; d++)
                for (int k = 0; k < r; k++)
                    for (int q = 0; q < r; q++)
                        inner[k, q] += p.Weights[d] * ((p.Coefficients[d, i, k] - mean[k]) * (p.Coefficients[d, i, q] - mean[q]) + sigma[d, k, q]);
            var expected = new double[m, m];
            for (int a = 0; a < m; a++)
                for (int b = 0; b < m; b++)
                {
                    double v = 0.0;
                    for (int k = 0; k < r; k++)
                        for (int q = 0; q < r; q++) v += basis[a, k] * inner[k, q] * basis[b, q];
                    expected[a, b] = scales[a] * v * scales[b];
                }
            Assert.True(MaxDifference(expected, p.Covariance(i)) < 1e-12 * (1.0 + MaxAbs(expected)));
        }

        // The sampled noise: ε = Kᵀ(predictive − latent)/y_s (K orthonormal here), its second moment E[Σ] over the draws.
        var second = new double[r, r];
        var expectedSecond = new double[r, r];
        int samples = 0;
        for (int d = 0; d < s; d++)
        {
            for (int k = 0; k < r; k++)
                for (int q = 0; q < r; q++) expectedSecond[k, q] += sigma[d, k, q] / s;
            for (int i = 0; i < cases.Cases; i++, samples++)
            {
                var e = new double[r];
                for (int k = 0; k < r; k++)
                    for (int a = 0; a < m; a++) e[k] += basis[a, k] * (p.Predictive![d, i, a] - p.Latent[d, i, a]) / scales[a];
                for (int k = 0; k < r; k++)
                    for (int q = 0; q < r; q++) second[k, q] += e[k] * e[q];
            }
        }
        for (int k = 0; k < r; k++)
            for (int q = 0; q < r; q++)
            {
                second[k, q] /= samples;
                double se = Math.Sqrt((expectedSecond[k, k] * expectedSecond[q, q] + expectedSecond[k, q] * expectedSecond[k, q]) / samples);
                output.WriteLine($"sampled E[εεᵀ][{k},{q}] {second[k, q]:F5} against E[Σ] {expectedSecond[k, q]:F5} (se {se:F5})");
                Assert.True(Math.Abs(second[k, q] - expectedSecond[k, q]) < 5.0 * se, $"[{k},{q}]");
            }
    }

    [Fact]
    public void CorrelatedScores_TheDensityAndTheEmulatorReadSigmaWhole()
    {
        (BassBasisFit fit, BassDecomposition d, double[][] xs) = Correlated.Value;
        int n = d.Cases, r = d.Rank, m = d.Locations;
        // The components' fits at the training cases, each joint draw's, from the prediction at the training inputs.
        var training = new PredictionCases([new ContinuousColumn("x0", xs[0]), new ContinuousColumn("x1", xs[1]), new ContinuousColumn("x2", xs[2])]);
        BassBasisPrediction at = fit.Predict(training, new BassPredictOptions { Draws = [.. Enumerable.Range(0, fit.JointDraws.Count)] });
        BassDraw[] pooled = [.. fit.Joint!.Draws.SelectMany(c => c)];

        PointwiseLogDensity coefficients = fit.PointwiseLogDensity(space: BasisDensitySpace.Coefficients);
        PointwiseLogDensity response = fit.PointwiseLogDensity(space: BasisDensitySpace.ResponseWithDiagonalTruncation);
        double[,] basis = d.Map.BasisMatrix();
        double[] scales = d.Map.Scales.ToArray(), truncationMean = d.Truncation.Mean();
        double[,] truncationCovariance = d.Truncation.Covariance();
        double logScales = scales.Sum(Math.Log);
        double worstCoefficients = 0.0, worstResponse = 0.0;
        for (int s = 0; s < fit.JointDraws.Count; s += 7)
        {
            double[,] sigma = pooled[fit.JointDraws.Draw(0, s)].Sigma!;
            var covariance = new double[m, m];
            for (int a = 0; a < m; a++)
            {
                covariance[a, a] = truncationCovariance[a, a];
                for (int b = 0; b < m; b++)
                    for (int k = 0; k < r; k++)
                        for (int q = 0; q < r; q++) covariance[a, b] += basis[a, k] * sigma[k, q] * basis[b, q];
            }
            for (int i = 0; i < n; i++)
            {
                var e = new double[r];
                for (int k = 0; k < r; k++) e[k] = d.Coefficient(k, i) - at.Coefficients[s, i, k];
                worstCoefficients = Math.Max(worstCoefficients, Math.Abs(DenseLogNormal(e, sigma) - coefficients.Values[s, i]));
                var x = new double[m];
                for (int a = 0; a < m; a++)
                {
                    double y = d.Truncation.InSample![i, a] - truncationMean[a];
                    for (int k = 0; k < r; k++) y += basis[a, k] * (d.Coefficient(k, i) - at.Coefficients[s, i, k]);
                    x[a] = y;
                }
                worstResponse = Math.Max(worstResponse, Math.Abs(DenseLogNormal(x, covariance) - logScales - response.Values[s, i]));
            }
        }
        output.WriteLine($"density against dense Gaussians: coefficient space {worstCoefficients:G3}, y space {worstResponse:G3}");
        Assert.True(worstCoefficients < 1e-9 && worstResponse < 1e-8);

        // A unit splitting a case's components: refused jointly, read apart.
        UnitPartition split = UnitPartition.EachObservation(n * r);
        Assert.Contains("whole cases", Assert.Throws<ArgumentException>(() => fit.PointwiseLogDensity(split)).Message);

        // The emulator's noise: (y_s ∘ K)Σ(y_s ∘ K)ᵀ.
        var emulator = new BassBasisEmulator(fit, [.. fit.Predictors.Select(pr => new ContinuousParameter(pr.Name, 0.0, 1.0))]);
        for (int s = 0; s < emulator.DrawCount; s += 11)
        {
            double[,] sigma = pooled[fit.JointDraws.Draw(0, s)].Sigma!;
            var expected = new double[m, m];
            for (int a = 0; a < m; a++)
                for (int b = 0; b < m; b++)
                    for (int k = 0; k < r; k++)
                        for (int q = 0; q < r; q++) expected[a, b] += scales[a] * basis[a, k] * sigma[k, q] * basis[b, q] * scales[b];
            Assert.True(MaxDifference(expected, Dense((ComponentNoise)emulator.Noise(s))) < 1e-14);
        }
    }

    [Fact]
    public void TheViews_ReadAsComponentsFittedApart()
    {
        (BassBasisFit fit, BassDecomposition d, _) = Correlated.Value;
        Assert.Equal(d.Rank, fit.Components.Count);
        for (int k = 0; k < d.Rank; k++)
        {
            BassFit view = fit.Components[k];
            Assert.Null(view.Config.OutputCovariance);
            Assert.Null(view.Evidence);
            Assert.Null(view.Checkpoint);
            Assert.Equal(Bass.Fingerprint(view.Config, fit.ComponentData(k)), view.Identity.Fingerprint);
            BassDraw first = view.Draws[0][0], own = fit.Joint!.Draws[0][0];
            Assert.Equal(own.Sigma![k, k], first.Sigma2);
            Assert.Null(first.Sigma);
            Assert.Null(first.OutputCoefficients);
        }
        for (int s = 0; s < fit.JointDraws.Count; s++)
            for (int k = 1; k < d.Rank; k++) Assert.Equal(fit.JointDraws.Draw(0, s), fit.JointDraws.Draw(k, s));
        // The surface readers run on the views as on components fitted apart.
        Assert.Equal(d.Rank, fit.Traces().Count);
        Assert.Equal(d.Rank, fit.Metadata().ComponentFits.Count);
        Assert.NotNull(fit.Sensitivity(new SensitivitySpec { Draws = [0, 1, 2] }));
    }

    // ── Weighted draws ────────────────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void WeightedDraws_AreResampledOnce_ForEveryComponent()
    {
        (double[,] y, BassPredictor[] predictors) = BassOutputComponentsTests.Problem();
        BassConfig config = BassOutputComponentsTests.Config();
        config = config with { Sampling = config.Sampling with { Smc = new SmcSpec { Particles = 200, Replicates = 2 } } };

        // One resample of the joint fit's population, the same draw for every component; the views carry its weights. (A rank-1
        // joint fit under SMC is not a component fitted apart draw for draw: several outputs' SMC takes its own path at K = 1, so
        // the pairing is gated against the joint fit itself.)
        BassDecomposition two = BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 2, ScaleByColumnSd = true });
        BassBasisFit pair = BassComponents.FitJoint(two, predictors, config);
        BassFit joint = pair.Joint!;
        Assert.True(pair.JointDraws.Resampled);
        Assert.Null(pair.JointDraws.ChainLengths);
        int pooled = joint.Draws.Sum(c => c.Count);
        var used = new HashSet<int>();
        for (int s = 0; s < pair.JointDraws.Count; s++)
        {
            int d = pair.JointDraws.Draw(0, s);
            Assert.Equal(d, pair.JointDraws.Draw(1, s));
            Assert.InRange(d, 0, pooled - 1);
            Assert.True(joint.DrawWeights![d] > 0.0);
            used.Add(d);
        }
        foreach (BassFit view in pair.Components) Assert.Equal(joint.DrawWeights, view.DrawWeights);
        output.WriteLine($"SMC: {pair.JointDraws.Count} joint draws from {used.Count} distinct particles of {pooled}");
        var at = new PredictionCases([new ContinuousColumn("z0", [0.2, 0.8]), new ContinuousColumn("z1", [0.5, 0.1])]);
        BassBasisPrediction p = pair.Predict(at);
        Assert.NotNull(p.NoiseCovariances);
        for (int s = 0; s < p.DrawCount; s++)
        {
            double[,] own = joint.Draws.SelectMany(c => c).ElementAt(pair.JointDraws.Draw(0, p.Draws[s])).Sigma!;
            Assert.Equal(own[1, 0], p.NoiseCovariances![s, 1, 0]);
        }
    }

    // ── Refusals ──────────────────────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Refuses_RandomIntercepts_AndAJointFitOfOtherData()
    {
        (double[,] y, BassPredictor[] predictors) = BassOutputComponentsTests.Problem();
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 2, ScaleByColumnSd = true });
        BassConfig grouped = BassOutputComponentsTests.Config() with { RandomIntercepts = [new BassRandomIntercept("g", 1.0)] };
        Assert.Contains("random intercepts", Assert.Throws<ArgumentException>(() => BassComponents.FitJoint(d, predictors, grouped)).Message);

        BassBasisFit fit = BassComponents.FitJoint(d, predictors, BassOutputComponentsTests.Config(samples: 20));
        BassDecomposition other = BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 2 });
        Assert.Throws<ArgumentException>(() => BassComponents.AssembleJoint(other, predictors, fit.Joint!, 31));
        BassBasisFit again = BassComponents.AssembleJoint(d, predictors, fit.Joint!, BassOutputComponentsTests.Config().Sampling.MasterSeed);
        var at = new PredictionCases([new ContinuousColumn("z0", [0.3]), new ContinuousColumn("z1", [0.6])]);
        Assert.Equal(0.0, MaxDifference(fit.Predict(at).Latent, again.Predict(at).Latent));
    }

    // ── The cost row (qualification) ──────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Cost_TheJointFitAgainstComponentsFittedApart()
    {
        if (!BassQualification.Enabled) return;   // indicative, wall clock
        var one = new BassRunOptions { MaxDegreeOfParallelism = 1 };
        BassConfig config = BassRunTests.Config(chains: 1, burn: 300, samples: 600, seed: 47, maxBasis: 15);
        {   // a pass to warm both paths, discarded
            (BassDecomposition d, BassPredictor[] predictors, _, _) = Scores(300, Diagonal([0.1, 0.1]), 2026_1112, 2);
            BassComponents.Fit(d, predictors, config, one);
            BassComponents.FitJoint(d, predictors, config, one);
        }
        foreach (int r in new[] { 2, 3, 5 })
        {
            var noise = new double[r, r];
            for (int k = 0; k < r; k++) noise[k, k] = 0.01 * (1 + k);
            (BassDecomposition d, BassPredictor[] predictors, _, _) = Scores(300, noise, 2026_1112 + r, r);
            var apart = new List<double>();
            var joint = new List<double>();
            for (int repeat = 0; repeat < 3; repeat++)   // interleaved
            {
                var watch = System.Diagnostics.Stopwatch.StartNew();
                BassComponents.Fit(d, predictors, config, one);
                apart.Add(watch.Elapsed.TotalMilliseconds);
                watch.Restart();
                BassComponents.FitJoint(d, predictors, config, one);
                joint.Add(watch.Elapsed.TotalMilliseconds);
            }
            double a = apart.Order().ElementAt(1), j = joint.Order().ElementAt(1);
            output.WriteLine($"r = {r}: apart {a:F0} ms, joint {j:F0} ms, joint/apart {j / a:F2} (medians of 3; one thread, 900 sweeps, n = 300, M_max = 15)");
        }
    }

    // ── Fixtures ──────────────────────────────────────────────────────────────────────────────────────────────────────────────

    // r scores at n cases from surfaces outside the model, plus noise of covariance `noise`, on an orthonormal cosine basis of m = 8
    // locations with offsets 0, scales 1 and small supplied residuals.
    internal static (BassDecomposition D, BassPredictor[] Predictors, double[,] Truth, double[][] X) Scores(int n, double[,] noise, int seed, int r = 3)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        const int m = 8;
        var x = new double[3][];
        for (int v = 0; v < 3; v++) x[v] = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        double[,] root = BassRootForTests(noise);
        var truth = new double[n, r];
        var w = new double[r, n];
        for (int i = 0; i < n; i++)
        {
            double[] f = Surface(x[0][i], x[1][i], x[2][i], r);
            var z = new double[r];
            for (int k = 0; k < r; k++) z[k] = Normal.SampleStandard(rng);
            for (int k = 0; k < r; k++)
            {
                double e = 0.0;
                for (int q = 0; q <= k; q++) e += root[k, q] * z[q];
                truth[i, k] = f[k];
                w[k, i] = f[k] + e;
            }
        }
        var basis = new double[m, r];
        for (int l = 0; l < m; l++)
            for (int k = 0; k < r; k++) basis[l, k] = Math.Sqrt((k == 0 ? 1.0 : 2.0) / m) * Math.Cos(Math.PI * k * (l + 0.5) / m);
        // The training responses Kwᵢ + eᵢ, so the decomposition holds its in-sample residuals e (the y-space density reads them).
        var responses = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++)
            {
                double v = 0.02 * Normal.SampleStandard(rng);
                for (int k = 0; k < r; k++) v += basis[l, k] * w[k, i];
                responses[i, l] = v;
            }
        BassDecomposition d = BassDecomposition.FromBasis(basis, w, new double[m], [.. Enumerable.Repeat(1.0, m)], responses: responses);
        return (d, [new ContinuousPredictor("x0", x[0]), new ContinuousPredictor("x1", x[1]), new ContinuousPredictor("x2", x[2])], truth, x);
    }

    private static double[] Surface(double x0, double x1, double x2, int r)
    {
        double[] all =
        [
            Math.Sin(3.0 * x0) + 0.5 * x1,
            2.0 * Math.Exp(-x1) * x2,
            Math.Cos(2.0 * x2) + x0 * x1,
            Math.Tanh(2.0 * (x0 - x2)),
            0.8 * x1 * x1 - 0.3 * x0,
        ];
        return all[..r];
    }

    private static (PredictionCases Cases, double[,] Truth) HeldOut(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[3][];
        for (int v = 0; v < 3; v++) x[v] = [.. Enumerable.Range(0, n).Select(_ => 0.05 + 0.9 * rng.NextDouble())];
        var truth = new double[n, 3];
        for (int i = 0; i < n; i++)
        {
            double[] f = Surface(x[0][i], x[1][i], x[2][i], 3);
            for (int k = 0; k < 3; k++) truth[i, k] = f[k];
        }
        return (new PredictionCases([new ContinuousColumn("x0", x[0]), new ContinuousColumn("x1", x[1]), new ContinuousColumn("x2", x[2])]), truth);
    }

    private static double[,] MeanCoefficients(BassBasisPrediction p)
    {
        int s = p.DrawCount, n = p.Cases, r = p.Coefficients.GetLength(2);
        var mean = new double[n, r];
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++)
                for (int k = 0; k < r; k++) mean[i, k] += p.Weights[d] * p.Coefficients[d, i, k];
        return mean;
    }

    private static double Rmse(double[,] a, double[,] b, int k)
    {
        double sum = 0.0;
        int n = a.GetLength(0);
        for (int i = 0; i < n; i++) sum += (a[i, k] - b[i, k]) * (a[i, k] - b[i, k]);
        return Math.Sqrt(sum / n);
    }

    private static double[,] Diagonal(double[] sd)
    {
        var c = new double[sd.Length, sd.Length];
        for (int k = 0; k < sd.Length; k++) c[k, k] = sd[k] * sd[k];
        return c;
    }

    // The emulator's noise as a dense m × m matrix, F·diag(v)·Fᵀ.
    private static double[,] Dense(ComponentNoise noise)
    {
        int m = noise.Factor.GetLength(0), r = noise.Factor.GetLength(1);
        var c = new double[m, m];
        for (int a = 0; a < m; a++)
            for (int b = 0; b < m; b++)
                for (int k = 0; k < r; k++) c[a, b] += noise.Factor[a, k] * noise.Variances[k] * noise.Factor[b, k];
        return c;
    }

    // log N(x | 0, c) by a dense Cholesky.
    private static double DenseLogNormal(double[] x, double[,] c)
    {
        double[,] l = BassRootForTests(c);
        int n = x.Length;
        var u = new double[n];
        double logDet = 0.0, q = 0.0;
        for (int i = 0; i < n; i++)
        {
            double v = x[i];
            for (int j = 0; j < i; j++) v -= l[i, j] * u[j];
            u[i] = v / l[i, i];
            q += u[i] * u[i];
            logDet += Math.Log(l[i, i]);
        }
        return -0.5 * n * Math.Log(2.0 * Math.PI) - logDet - 0.5 * q;
    }

    private static double[,] BassRootForTests(double[,] c)
    {
        int n = c.GetLength(0);
        var l = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = c[i, j];
                for (int q = 0; q < j; q++) sum -= l[i, q] * l[j, q];
                l[i, j] = i == j ? Math.Sqrt(sum) : sum / l[j, j];
            }
        return l;
    }

    private static double[,] Inverse3(double[,] a)
    {
        double det = a[0, 0] * (a[1, 1] * a[2, 2] - a[1, 2] * a[2, 1]) - a[0, 1] * (a[1, 0] * a[2, 2] - a[1, 2] * a[2, 0])
                     + a[0, 2] * (a[1, 0] * a[2, 1] - a[1, 1] * a[2, 0]);
        var inv = new double[3, 3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                int i1 = (j + 1) % 3, i2 = (j + 2) % 3, j1 = (i + 1) % 3, j2 = (i + 2) % 3;
                inv[i, j] = (a[i1, j1] * a[i2, j2] - a[i1, j2] * a[i2, j1]) / det;
            }
        return inv;
    }

    private static double[,] Times(double[,] a, double[,] l)
    {
        int n = a.GetLength(0), r = l.GetLength(1);
        var c = new double[n, r];
        for (int i = 0; i < n; i++)
            for (int k = 0; k < r; k++)
                for (int q = 0; q < l.GetLength(0); q++) c[i, k] += a[i, q] * l[q, k];
        return c;
    }

    // Lᵀ S L.
    private static double[,] Congruent(double[,] s, double[,] l)
    {
        int r = l.GetLength(0);
        var c = new double[r, r];
        for (int a = 0; a < r; a++)
            for (int b = 0; b < r; b++)
                for (int i = 0; i < r; i++)
                    for (int j = 0; j < r; j++) c[a, b] += l[i, a] * s[i, j] * l[j, b];
        return c;
    }

    private static double Relative(double a, double b) => Math.Abs(a - b) / Math.Max(1e-300, Math.Max(Math.Abs(a), Math.Abs(b)));

    private static double MaxAbs(double[,] a)
    {
        double worst = 0.0;
        foreach (double v in a) worst = Math.Max(worst, Math.Abs(v));
        return worst;
    }

    private static double MaxDifference(double[,] a, double[,] b)
    {
        Assert.Equal(a.GetLength(0), b.GetLength(0));
        Assert.Equal(a.GetLength(1), b.GetLength(1));
        double worst = 0.0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) worst = Math.Max(worst, Math.Abs(a[i, j] - b[i, j]));
        return worst;
    }

    private static double MaxDifference(double[,,] a, double[,,] b)
    {
        double worst = 0.0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                for (int k = 0; k < a.GetLength(2); k++) worst = Math.Max(worst, Math.Abs(a[i, j, k] - b[i, j, k]));
        return worst;
    }
}
