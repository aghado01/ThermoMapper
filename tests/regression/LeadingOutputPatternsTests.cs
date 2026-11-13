using System;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The leading-patterns output readout (spectral-solver brief, E): E[G]'s leading eigenpairs by LOBPCG beside the full readout's
/// stacked SVD, with each draw's spectrum and patterns read alike.
/// </summary>
public sealed class LeadingOutputPatternsTests
{
    /// <summary>
    /// On smooth output factors (twelve shared cosine modes, variances falling as 1/j², per-draw noise) the leading readout's draw
    /// quantities are the full readout's to the bit; its eigenvalues match to 10⁻¹⁰ and its trace to 10⁻¹²; and its patterns lie in
    /// the full readout's leading eigenspace to the Davis–Kahan bound √k·r/(λ_k − λ_{k+1}), r the solve's reported residual.
    /// </summary>
    [Fact]
    public void Leading_AgreesWithTheFullReadout()
    {
        const int m = 200, k = 3;
        double[][,] factors = SmoothFactors(m, draws: 30, width: 5, seed: 17);
        double[] weights = [.. Enumerable.Range(1, 30).Select(s => 1.0 + 0.1 * s)];
        OutputPosterior full = OutputPosterior.Of(factors, weights, k);
        OutputPosterior leading = OutputPosterior.Leading(factors, weights, k);

        Assert.Equal(OutputSpectrum.Full, full.Spectrum);
        Assert.Null(full.MeanSolve);
        Assert.Equal(OutputSpectrum.Leading, leading.Spectrum);
        Assert.NotNull(leading.MeanSolve);
        Assert.True(leading.MeanSolve!.Converged);
        Assert.Equal(full.Weights, leading.Weights);
        Assert.Equal(full.Traces, leading.Traces);
        for (int s = 0; s < factors.Length; s++)
        {
            Assert.Equal(full.DrawEigenvalues[s], leading.DrawEigenvalues[s]);
            Assert.Equal(full.DrawPatterns[s], leading.DrawPatterns[s]);
        }
        Assert.Equal(k, leading.MeanEigenvalues.Length);
        for (int j = 0; j < k; j++)
            Assert.InRange(Math.Abs(leading.MeanEigenvalues[j] / full.MeanEigenvalues[j] - 1), 0.0, 1e-10);
        Assert.InRange(Math.Abs(leading.MeanTrace / full.MeanTrace - 1), 0.0, 1e-12);
        double gap = full.MeanEigenvalues[k - 1] - full.MeanEigenvalues[k];
        Assert.InRange(SinTheta(full.MeanPatterns, leading.MeanPatterns), 0.0, Math.Sqrt(k) * leading.MeanSolve.ResidualNorm / gap + 1e-12);
    }

    /// <summary>
    /// The brief's witness: two equally weighted draws with G = diag(4, 3, 0) and diag(0, 3, 4) on the first three of 40 locations.
    /// E[G] = diag(2, 3, 2), so its leading pattern is the second location's, with eigenvalue 3, though each draw leads with another.
    /// </summary>
    [Fact]
    public void Leading_ReadsTheMeansPattern_NotTheDraws()
    {
        double[][,] factors = Witness();
        OutputPosterior leading = OutputPosterior.Leading(factors, [1.0, 1.0], 1);

        Assert.NotNull(leading.MeanSolve);
        Assert.True(leading.MeanSolve!.Converged);
        Assert.InRange(Math.Abs(leading.MeanEigenvalues.Single() - 3.0), 0.0, 1e-9);
        for (int l = 0; l < 40; l++) Assert.InRange(Math.Abs(leading.MeanPatterns[l, 0] - (l == 1 ? 1.0 : 0.0)), 0.0, 1e-9);
        Assert.InRange(Math.Abs(leading.DrawEigenvalues[0][0] - 4.0), 0.0, 1e-12);
        Assert.InRange(Math.Abs(leading.DrawPatterns[0][0, 0] - 1.0), 0.0, 1e-12);
        Assert.InRange(Math.Abs(leading.DrawPatterns[1][2, 0] - 1.0), 0.0, 1e-12);
        Assert.InRange(Math.Abs(leading.MeanTrace - 7.0), 0.0, 1e-12);
    }

    /// <summary>
    /// Five patterns asked of the witness, whose E[G] has rank 3 over four compressed columns: too small to iterate, so the stacked
    /// SVD reads it, cut to the three nonzero eigenvalues 3, 2 and 2, and no solve is reported.
    /// </summary>
    [Fact]
    public void Leading_TooSmallToIterate_ReadsTheStackedSvd()
    {
        OutputPosterior leading = OutputPosterior.Leading(Witness(), [1.0, 1.0], 5);

        Assert.Null(leading.MeanSolve);
        Assert.Equal(OutputSpectrum.Leading, leading.Spectrum);
        Assert.Equal(3, leading.MeanEigenvalues.Length);
        Assert.Equal(3, leading.MeanPatterns.GetLength(1));
        double[] expected = [3.0, 2.0, 2.0];
        for (int j = 0; j < 3; j++) Assert.InRange(Math.Abs(leading.MeanEigenvalues[j] - expected[j]), 0.0, 1e-12);
        Assert.InRange(Math.Abs(leading.MeanTrace - 7.0), 0.0, 1e-12);
    }

    /// <summary>
    /// The spec's switch: on the functional fixture, asking for the leading spectrum leaves every draw-level output and the input-side
    /// readout as the full spectrum reads them, and E[G]'s two leading eigenvalues and patterns agree.
    /// </summary>
    [Fact]
    public void Spec_SwitchesTheMeanStageOnly()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        var spec = new ActiveSubspaceSpec { Draws = [.. Enumerable.Range(0, 8)], OutputPatterns = 2 };
        ActiveSubspaceAnalysis full = fit.ActiveSubspace(data, spec);
        ActiveSubspaceAnalysis leading = fit.ActiveSubspace(data, spec with { OutputSpectrum = OutputSpectrum.Leading });

        Assert.Equal(OutputSpectrum.Full, full.Outputs!.Spectrum);
        Assert.Equal(OutputSpectrum.Leading, leading.Outputs!.Spectrum);
        Assert.Equal(full.Outputs.Traces, leading.Outputs.Traces);
        for (int d = 0; d < full.Posterior.Draws.Count; d++)
        {
            Assert.Equal(full.Posterior.Draws[d].C, leading.Posterior.Draws[d].C);
            Assert.Equal(full.Outputs.DrawPatterns[d], leading.Outputs.DrawPatterns[d]);
        }
        for (int j = 0; j < 2; j++)
            Assert.InRange(Math.Abs(leading.Outputs.MeanEigenvalues[j] / full.Outputs.MeanEigenvalues[j] - 1), 0.0, 1e-8);
        Assert.InRange(SinTheta(full.Outputs.MeanPatterns, leading.Outputs.MeanPatterns), 0.0, 1e-5);
    }

    // ‖(I − UUᵀ)V‖_F for orthonormal columns U and V.
    private static double SinTheta(double[,] u, double[,] v)
    {
        int m = u.GetLength(0);
        double total = 0.0;
        for (int b = 0; b < v.GetLength(1); b++)
        {
            var rest = new double[m];
            for (int l = 0; l < m; l++) rest[l] = v[l, b];
            for (int a = 0; a < u.GetLength(1); a++)
            {
                double along = 0.0;
                for (int l = 0; l < m; l++) along += u[l, a] * v[l, b];
                for (int l = 0; l < m; l++) rest[l] -= along * u[l, a];
            }
            total += rest.Sum(x => x * x);
        }
        return Math.Sqrt(total);
    }

    private static double[][,] Witness()
    {
        var first = new double[40, 3];
        var second = new double[40, 3];
        first[0, 0] = 2.0;
        first[1, 1] = Math.Sqrt(3.0);
        second[1, 1] = Math.Sqrt(3.0);
        second[2, 2] = 2.0;
        return [first, second];
    }

    private static double[][,] SmoothFactors(int m, int draws, int width, int seed)
    {
        var random = new Random(seed);
        const int modes = 12;
        var basis = new double[modes, m];
        for (int j = 0; j < modes; j++)
        {
            double norm = 0.0;
            for (int l = 0; l < m; l++) norm += Math.Pow(basis[j, l] = Math.Cos(Math.PI * j * l / (m - 1)), 2);
            for (int l = 0; l < m; l++) basis[j, l] /= Math.Sqrt(norm);
        }
        var shared = new double[modes, width];
        for (int j = 0; j < modes; j++)
            for (int c = 0; c < width; c++) shared[j, c] = Gaussian(random) / (1 + j);
        var factors = new double[draws][,];
        for (int s = 0; s < draws; s++)
        {
            var f = new double[m, width];
            for (int c = 0; c < width; c++)
            {
                var a = new double[modes];
                for (int j = 0; j < modes; j++) a[j] = shared[j, c] + 0.3 * Gaussian(random) / (1 + j);
                for (int l = 0; l < m; l++)
                {
                    double v = 1e-3 * Gaussian(random);
                    for (int j = 0; j < modes; j++) v += basis[j, l] * a[j];
                    f[l, c] = v;
                }
            }
            factors[s] = f;
        }
        return factors;
    }

    private static double Gaussian(Random r) => Math.Sqrt(-2.0 * Math.Log(1.0 - r.NextDouble())) * Math.Cos(2.0 * Math.PI * r.NextDouble());
}
