using System;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// A basis fit's pointwise log density (stage 5, O4; OQ17; D07): in coefficient space, the sum over components of each
/// component's own retained-parameters density at its paired draw, by case or under any partition, with the joint chains;
/// under the declared diagonal Gaussian, the y-space density against a dense Gaussian in physical units; WAIC and
/// PSIS-LOO reading it; and the refusals.
/// </summary>
public sealed class BassOutputDensityTests
{
    private readonly ITestOutputHelper _out;
    public BassOutputDensityTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Coefficients_SumEachComponentsOwnDensity_ByCase()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        BassDecomposition d = fit.Decomposition;
        int n = d.Cases, r = d.Rank, s = fit.JointDraws.Count;
        PointwiseLogDensity density = fit.PointwiseLogDensity();
        Assert.Equal((s, n), (density.Draws, density.Units.Units));
        Assert.Equal(fit.JointDraws.ChainLengths, density.ChainLengths);
        Assert.False(density.LeavesUnitsOut);
        PointwiseLogDensity[] own = [.. Enumerable.Range(0, r).Select(k => fit.Components[k].PointwiseLogDensity(fit.ComponentData(k)))];
        for (int j = 0; j < s; j++)
            for (int i = 0; i < n; i++)
            {
                double sum = Enumerable.Range(0, r).Sum(k => own[k].Values[fit.JointDraws.Draw(k, j), i]);
                Assert.Equal(sum, density.Values[j, i], 12);
            }

        // Any partition of the n·r coefficients: here each (case, component) its own unit, then pairs of cases.
        PointwiseLogDensity each = fit.PointwiseLogDensity(UnitPartition.EachObservation(n * r));
        for (int j = 0; j < s; j += 7)
            for (int i = 0; i < n; i++)
                for (int k = 0; k < r; k++) Assert.Equal(own[k].Values[fit.JointDraws.Draw(k, j), i], each.Values[j, i * r + k], 12);
        PointwiseLogDensity pairs = fit.PointwiseLogDensity(new UnitPartition([.. Enumerable.Range(0, n * r).Select(o => o / (2 * r))]));
        for (int j = 0; j < s; j += 7)
            for (int u = 0; u < n / 2; u++) Assert.Equal(density.Values[j, 2 * u] + density.Values[j, 2 * u + 1], pairs.Values[j, u], 11);

        // The model-free consumers read it.
        PredictiveCheck waic = PredictiveChecks.Waic(density);
        PredictiveCheck loo = PredictiveChecks.PsisLoo(density);
        _out.WriteLine($"WAIC elpd {waic.Elpd:F2}, PSIS-LOO elpd {loo.Elpd:F2}");
        Assert.True(double.IsFinite(waic.Elpd) && double.IsFinite(loo.Elpd));
    }

    [Fact]
    public void DiagonalTruncation_IsADenseGaussianInPhysicalUnits()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        BassDecomposition d = fit.Decomposition;
        (double[,] y, _) = BassOutputComponentsTests.Problem();
        int n = d.Cases, r = d.Rank, m = d.Locations, s = fit.JointDraws.Count;
        PointwiseLogDensity density = fit.PointwiseLogDensity(space: BasisDensitySpace.ResponseWithDiagonalTruncation);
        Assert.Equal((s, n), (density.Draws, density.Units.Units));
        Assert.Contains(density.Latents, l => l.Name.Contains("truncation") && l.Treatment == LatentTreatment.Integrated);

        double[] eBar = d.Truncation.Mean();
        double[,] te = d.Truncation.Covariance();
        var fitted = new double[r][];
        for (int j = 0; j < s; j += 11)
        {
            var sigma = new double[m, m];
            for (int k = 0; k < r; k++)
            {
                BassFit component = fit.Components[k];
                BassPrediction own = component.Predict(fit.ComponentData(k), PredictionCases.Of(fit.Predictors),
                                                       new BassPredictOptions { Draws = [fit.JointDraws.Draw(k, j)], Predictive = false });
                fitted[k] = [.. Enumerable.Range(0, n).Select(i => own.Latent[0, i, 0])];
                double s2 = own.NoiseVariances[0];
                for (int l = 0; l < m; l++)
                    for (int q = 0; q < m; q++) sigma[l, q] += d.Map.Basis(l, k) * d.Map.Basis(q, k) * s2;
            }
            // Dense: Σ_phys = diag(y_s)(K diag(σ²)Kᵀ + diag(v))diag(y_s), the mean y_m + y_s ∘ (K w(x) + ē).
            for (int l = 0; l < m; l++)
            {
                sigma[l, l] += te[l, l];
                for (int q = 0; q < m; q++) sigma[l, q] *= d.Map.Scales[l] * d.Map.Scales[q];
            }
            var factor = new DenseCholesky(m);
            factor.Decompose(sigma);
            for (int i = 0; i < n; i++)
            {
                var x = new double[m];
                for (int l = 0; l < m; l++)
                    x[l] = y[i, l] - d.Map.Offsets[l] - d.Map.Scales[l] * (Enumerable.Range(0, r).Sum(k => d.Map.Basis(l, k) * fitted[k][i]) + eBar[l]);
                double[] solved = factor.Solve(x);
                double quadratic = Enumerable.Range(0, m).Sum(l => x[l] * solved[l]);
                double dense = -0.5 * (m * Math.Log(2 * Math.PI) + factor.LogDet + quadratic);
                Assert.Equal(dense, density.Values[j, i], 1e-9 * Math.Abs(dense));
            }
        }
    }

    [Fact]
    public void Refusals_SayWhy()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        int n = fit.Decomposition.Cases, m = fit.Decomposition.Locations;
        Assert.Contains("cover", Assert.Throws<ArgumentException>(() => fit.PointwiseLogDensity(UnitPartition.EachObservation(n))).Message);
        Assert.Contains("whole cases", Assert.Throws<ArgumentException>(() =>
            fit.PointwiseLogDensity(UnitPartition.EachObservation(n * m), BasisDensitySpace.ResponseWithDiagonalTruncation)).Message);

        // A constant column has zero truncation variance; supplied residuals carry no responses.
        (double[,] y, BassPredictor[] predictors) = BassOutputComponentsTests.Problem();
        for (int i = 0; i < n; i++) y[i, 3] = 2.5;
        BassDecomposition constant = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 0.95 });
        BassBasisFit flat = BassComponents.Fit(constant, predictors, BassOutputComponentsTests.Config(samples: 20));
        Assert.Contains("location 3", Assert.Throws<ArgumentException>(() =>
            flat.PointwiseLogDensity(space: BasisDensitySpace.ResponseWithDiagonalTruncation)).Message);
        BassDecomposition d = fit.Decomposition;
        BassDecomposition supplied = BassDecomposition.FromBasis(d.Map.BasisMatrix(), Coefficients(d), d.Map.Offsets.ToArray(), d.Map.Scales.ToArray(),
                                                                 residuals: d.Truncation.Values);
        BassBasisFit given = BassComponents.Assemble(supplied, fit.Predictors, fit.Components, fit.MasterSeed);
        Assert.Contains("training responses", Assert.Throws<ArgumentException>(() =>
            given.PointwiseLogDensity(space: BasisDensitySpace.ResponseWithDiagonalTruncation)).Message);
        Assert.Equal(fit.PointwiseLogDensity().Values.Cast<double>(), given.PointwiseLogDensity().Values.Cast<double>());
    }

    private static double[,] Coefficients(BassDecomposition d)
    {
        var w = new double[d.Rank, d.Cases];
        for (int k = 0; k < d.Rank; k++)
            for (int i = 0; i < d.Cases; i++) w[k, i] = d.Coefficient(k, i);
        return w;
    }
}
