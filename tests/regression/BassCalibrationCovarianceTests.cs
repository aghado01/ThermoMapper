using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// A synthetic emulator for calibration's tests: y_m + y_s ∘ (K w(θ, levels)) with w given, its draws' component noise
/// over y_s ∘ K (or σ²_e·I when it has no basis) and declared truncation residuals, all in physical units.
/// </summary>
internal sealed class TestEmulator : ICalibrationEmulator
{
    private readonly double[] _offset, _scale;
    private readonly double[,] _basis, _factor;
    private readonly Func<int, double[], int[], double[]> _w;
    private readonly double[][] _variances;
    private readonly double[,]? _residuals;
    private readonly bool _independent;

    /// <param name="w">Draw s's coefficients at θ (physical continuous values, level indices).</param>
    /// <param name="variances">Each draw's component variances, or its one σ²_e when <paramref name="independent"/>.</param>
    /// <param name="residuals">The truncation residuals in physical units, n × m; null for none.</param>
    public TestEmulator(IReadOnlyList<CalibrationParameter> parameters, double[] offset, double[] scale, double[,] basis,
                        Func<int, double[], int[], double[]> w, double[][] variances, double[,]? residuals = null, bool independent = false)
    {
        Parameters = parameters;
        _offset = offset;
        _scale = scale;
        _basis = basis;
        _w = w;
        _variances = variances;
        _residuals = residuals;
        _independent = independent;
        int m = offset.Length, r = basis.GetLength(1);
        _factor = new double[m, r];
        for (int l = 0; l < m; l++)
            for (int k = 0; k < r; k++) _factor[l, k] = scale[l] * basis[l, k];
    }

    public IReadOnlyList<CalibrationParameter> Parameters { get; }
    public int DrawCount => _variances.Length;
    public IReadOnlyList<double>? DrawWeights => null;
    public int Outputs => _offset.Length;
    public int Predictions;   // calls, for the caches' checks

    public double[,] Factor => _factor;

    public EmulatorNoise Noise(int draw)
        => _independent ? new IndependentNoise(_variances[draw][0]) : new ComponentNoise(_factor, (double[])_variances[draw].Clone());

    public double[,]? TruncationResiduals() => (double[,]?)_residuals?.Clone();

    public void Predict(int draw, ReadOnlySpan<double> continuous, ReadOnlySpan<int> levels, Span<double> output)
    {
        System.Threading.Interlocked.Increment(ref Predictions);
        double[] w = _w(draw, continuous.ToArray(), levels.ToArray());
        for (int l = 0; l < Outputs; l++)
        {
            double s = 0.0;
            for (int k = 0; k < w.Length; k++) s += _basis[l, k] * w[k];
            output[l] = _offset[l] + _scale[l] * s;
        }
    }
}

/// <summary>
/// Calibration's covariance formulations (stage 7, C1; D09) against dense forms built independently: the functional
/// formulation's residual-set term, low rank through <see cref="DiagonalPlusLowRank"/>, against (1/n)·Σ rᵢrᵢᵀ formed
/// whole; its diagonal option against R's centred sample variance; σ²_e·I for an emulator without a basis; the
/// multivariate formulation's dense sum; units (every output's scale multiplied by c moves the log-likelihood by
/// −m ln c and no further); and the refusals at construction.
/// </summary>
public sealed class BassCalibrationCovarianceTests
{
    private readonly ITestOutputHelper _out;
    public BassCalibrationCovarianceTests(ITestOutputHelper output) => _out = output;

    private const int M = 14, N = 6;

    private static (double[] Offset, double[] Scale, double[,] Basis, double[,] Residuals) Pieces(int seed = 3)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var offset = new double[M];
        var scale = new double[M];
        var basis = new double[M, 2];
        var residuals = new double[N, M];
        for (int l = 0; l < M; l++)
        {
            double t = l / (M - 1.0);
            offset[l] = 4 + t;
            scale[l] = 0.5 + 2 * t;
            basis[l, 0] = Math.Sin(Math.PI * t);
            basis[l, 1] = Math.Cos(Math.PI * t);
            for (int i = 0; i < N; i++) residuals[i, l] = scale[l] * (0.03 * Normal.SampleStandard(rng) + 0.01 * Math.Sin(3 * t + i));
        }
        return (offset, scale, basis, residuals);
    }

    internal static TestEmulator Emulator(double[,]? residuals, bool independent = false, double factorScale = 1.0)
    {
        (double[] offset, double[] scale, double[,] basis, double[,] r) = Pieces();
        double[] s = [.. scale.Select(v => v * factorScale)];
        double[] o = [.. offset.Select(v => v * factorScale)];
        double[,]? res = residuals is null ? null : (double[,])residuals.Clone();
        if (res is not null)
            for (int i = 0; i < res.GetLength(0); i++)
                for (int l = 0; l < M; l++) res[i, l] *= factorScale;
        return new TestEmulator([new ContinuousParameter("theta", 0, 1)], o, s, basis,
                                (_, theta, _) => [Math.Sin(3 * theta[0]), Math.Cos(2 * theta[0])],
                                [[0.3, 0.7], [0.05, 0.2]], res, independent);
    }

    private static MeasurementModel Groups(double ig = 0.1, double hc = 0.2)
        => new(new GaussianMeasurement([new VarianceGroup([.. Enumerable.Range(0, M / 2)], new InverseGammaVariance(10, ig)),
                                        new VarianceGroup([.. Enumerable.Range(M / 2, M / 2)], new HalfCauchyDeviation(hc))]), M);

    private static (double LogDet, double[] Solve) Dense(double[,] sigma, double[] b)
    {
        var f = new DenseCholesky(M);
        f.Decompose(sigma);
        Assert.True(f.Status.PositiveDefinite);
        return (f.LogDet, f.Solve(b));
    }

    private static double RelativeError(double[] a, double[] reference)
        => a.Zip(reference, (x, y) => Math.Abs(x - y)).Max() / reference.Max(Math.Abs);

    private static double[] Rhs() => [.. Enumerable.Range(0, M).Select(j => Math.Cos(0.7 * j) + 0.1 * j)];

    // diag(perOutput) + U diag(s) Uᵀ (+ the residual set's (1/n)·Σ rᵢrᵢᵀ), formed whole.
    private static double[,] DenseFunctional(double[] diagonal, double[,] u, double[] s, double[,]? residuals)
    {
        var sigma = new double[M, M];
        for (int j = 0; j < M; j++)
        {
            sigma[j, j] += diagonal[j];
            for (int l = 0; l < M; l++)
            {
                for (int k = 0; k < s.Length; k++) sigma[j, l] += u[j, k] * s[k] * u[l, k];
                if (residuals is null) continue;
                int n = residuals.GetLength(0);
                for (int i = 0; i < n; i++) sigma[j, l] += residuals[i, j] * residuals[i, l] / n;
            }
        }
        return sigma;
    }

    [Fact]
    public void Functional_ResidualSet_IsTheLowRankFormOfTheDenseSecondMoment()
    {
        (_, _, _, double[,] residuals) = Pieces();
        TestEmulator emulator = Emulator(residuals);
        MeasurementModel measurement = Groups();
        CalibrationCovariance covariance = CalibrationCovariance.Create(new FunctionalFormulation(), measurement, emulator);
        double[] variances = [0.02, 0.08];
        var diagonal = new double[M];
        measurement.OutputVariances(variances, diagonal);
        double[] b = Rhs();
        for (int s = 0; s < emulator.DrawCount; s++)
        {
            var noise = (ComponentNoise)emulator.Noise(s);
            CovarianceFactor factor = covariance.Factor(variances, noise);
            (double logDet, double[] solve) = Dense(DenseFunctional(diagonal, noise.Factor, noise.Variances, residuals), b);
            _out.WriteLine($"draw {s}: low rank vs dense ldet {Math.Abs(factor.LogDet - logDet):E1}, solve {RelativeError(factor.Solve(b), solve):E1}");
            Assert.Equal(logDet, factor.LogDet, 1e-11);
            Assert.InRange(RelativeError(factor.Solve(b), solve), 0, 1e-11);
            Assert.Equal(b.Zip(solve).Sum(t => t.First * t.Second), factor.InverseQuadraticForm(b), 1e-10);
        }
    }

    [Fact]
    public void Functional_ResidualSet_WiderThanTheOutputs_IsStillItsDenseForm()
    {
        // More residual vectors than outputs (n > m): the low-rank factor is wide, and the thin SVD takes its transpose.
        var rng = new Xoshiro256PlusPlus(19);
        var residuals = new double[3 * M, M];
        for (int i = 0; i < 3 * M; i++)
            for (int l = 0; l < M; l++) residuals[i, l] = 0.05 * Normal.SampleStandard(rng);
        TestEmulator emulator = Emulator(residuals);
        MeasurementModel measurement = Groups();
        double[] variances = [1e-4, 3e-4];
        var diagonal = new double[M];
        measurement.OutputVariances(variances, diagonal);
        var noise = (ComponentNoise)emulator.Noise(0);
        CovarianceFactor factor = CalibrationCovariance.Create(new FunctionalFormulation(), measurement, emulator).Factor(variances, noise);
        double[] b = Rhs();
        (double logDet, double[] solve) = Dense(DenseFunctional(diagonal, noise.Factor, noise.Variances, residuals), b);
        Assert.Equal(logDet, factor.LogDet, 1e-10);
        Assert.InRange(RelativeError(factor.Solve(b), solve), 0, 1e-10);
    }

    [Fact]
    public void Functional_Diagonal_IsRsCentredSampleVariance()
    {
        (_, _, _, double[,] residuals) = Pieces();
        TestEmulator emulator = Emulator(residuals);
        MeasurementModel measurement = Groups();
        CalibrationCovariance covariance = CalibrationCovariance.Create(new FunctionalFormulation { Truncation = TruncationTerm.Diagonal }, measurement, emulator);
        double[] variances = [0.02, 0.08];
        var diagonal = new double[M];
        measurement.OutputVariances(variances, diagonal);
        for (int j = 0; j < M; j++)
        {
            double[] column = [.. Enumerable.Range(0, N).Select(i => residuals[i, j])];
            double mean = column.Average();
            diagonal[j] += column.Sum(v => (v - mean) * (v - mean)) / (N - 1);   // R's var()
        }
        var noise = (ComponentNoise)emulator.Noise(0);
        CovarianceFactor factor = covariance.Factor(variances, noise);
        double[] b = Rhs();
        (double logDet, double[] solve) = Dense(DenseFunctional(diagonal, noise.Factor, noise.Variances, null), b);
        Assert.Equal(logDet, factor.LogDet, 1e-11);
        Assert.InRange(RelativeError(factor.Solve(b), solve), 0, 1e-11);
    }

    [Fact]
    public void WithoutABasis_TheNoiseIsOnTheDiagonal()
    {
        (_, _, _, double[,] residuals) = Pieces();
        MeasurementModel measurement = Groups();
        double[] variances = [0.02, 0.08];
        double[] b = Rhs();
        foreach (double[,]? res in new[] { null, residuals })
        {
            TestEmulator emulator = Emulator(res, independent: true);
            var diagonal = new double[M];
            measurement.OutputVariances(variances, diagonal);
            for (int j = 0; j < M; j++) diagonal[j] += 0.3;   // σ²_e of draw 0
            CovarianceFactor functional = CalibrationCovariance.Create(new FunctionalFormulation(), measurement, emulator).Factor(variances, emulator.Noise(0));
            (double logDet, double[] solve) = Dense(DenseFunctional(diagonal, new double[M, 0], [], res), b);
            Assert.Equal(logDet, functional.LogDet, 1e-11);
            Assert.InRange(RelativeError(functional.Solve(b), solve), 0, 1e-11);

            double[][] identity = [.. Enumerable.Range(0, M).Select(i => Enumerable.Range(0, M).Select(j => i == j ? 1.0 : 0.0).ToArray())];
            CovarianceFactor multivariate = CalibrationCovariance.Create(new MultivariateFormulation(identity), measurement, emulator).Factor(variances, emulator.Noise(0));
            Assert.Equal(logDet, multivariate.LogDet, 1e-11);
            Assert.InRange(RelativeError(multivariate.Solve(b), solve), 0, 1e-11);
        }
    }

    [Fact]
    public void Multivariate_IsTheDenseSum()
    {
        (_, _, _, double[,] residuals) = Pieces();
        TestEmulator emulator = Emulator(residuals);
        MeasurementModel measurement = Groups();
        // R: an AR(1) correlation within each group; Σ_δ: a squared-exponential covariance with a nugget.
        var correlation = new double[M][];
        var delta = new double[M][];
        for (int i = 0; i < M; i++)
        {
            correlation[i] = new double[M];
            delta[i] = new double[M];
            for (int j = 0; j < M; j++)
            {
                correlation[i][j] = measurement.GroupOf(i) == measurement.GroupOf(j) ? Math.Pow(0.6, Math.Abs(i - j)) : 0.0;
                double d = (i - j) / (M - 1.0);
                delta[i][j] = 0.01 * Math.Exp(-d * d / (2 * 0.2 * 0.2)) + (i == j ? 1e-8 : 0.0);
            }
        }
        double[] mean = [.. Enumerable.Range(0, M).Select(j => 0.1 * j)];
        CalibrationCovariance covariance = CalibrationCovariance.Create(new MultivariateFormulation(correlation, delta, mean), measurement, emulator);
        Assert.Equal(mean, covariance.Mean);
        Assert.Null(covariance.DiscrepancyBasis);

        double[] variances = [0.02, 0.08];
        var noise = (ComponentNoise)emulator.Noise(1);
        var sigma = DenseFunctional(new double[M], noise.Factor, noise.Variances, residuals);
        var sd = new double[M];
        measurement.OutputVariances(variances, sd);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < M; j++) sigma[i, j] += Math.Sqrt(sd[i]) * correlation[i][j] * Math.Sqrt(sd[j]) + delta[i][j];
        double[] b = Rhs();
        (double logDet, double[] solve) = Dense(sigma, b);
        CovarianceFactor factor = covariance.Factor(variances, noise);
        Assert.Equal(logDet, factor.LogDet, 1e-11);
        Assert.InRange(RelativeError(factor.Solve(b), solve), 0, 1e-11);

        // The diagonal option in place of the residual set: R's centred sample variance on the diagonal.
        CovarianceFactor diagonalOption = CalibrationCovariance.Create(new MultivariateFormulation(correlation, delta, mean) { Truncation = TruncationTerm.Diagonal },
                                                                       measurement, emulator).Factor(variances, noise);
        for (int j = 0; j < M; j++)
        {
            double[] column = [.. Enumerable.Range(0, N).Select(i => residuals[i, j])];
            double average = column.Average();
            for (int l = 0; l < M; l++)
                sigma[j, l] += (j == l ? column.Sum(v => (v - average) * (v - average)) / (N - 1) : 0.0)
                             - Enumerable.Range(0, N).Sum(i => residuals[i, j] * residuals[i, l]) / N;
        }
        (logDet, solve) = Dense(sigma, b);
        Assert.Equal(logDet, diagonalOption.LogDet, 1e-11);
        Assert.InRange(RelativeError(diagonalOption.Solve(b), solve), 0, 1e-11);
    }

    [Fact]
    public void Units_ScaledOutputs_MoveTheLogLikelihoodByAConstant()
    {
        // Every output, residual and noise factor in units c times larger, the group variances c² times: the Gaussian
        // log-likelihood moves by −m ln c at every θ, so the θ posterior is unchanged.
        (_, _, _, double[,] residuals) = Pieces();
        double[] y = [.. Enumerable.Range(0, M).Select(j => 4.5 + 0.1 * Math.Sin(j))];
        foreach (TruncationTerm term in new[] { TruncationTerm.ResidualSet, TruncationTerm.Diagonal })
            foreach (double c in new[] { 0.01, 37.0 })
            {
                double[] shifts = [.. new[] { 0.1, 0.45, 0.8 }.Select(theta =>
                {
                    double LogLikelihood(double factorScale)
                    {
                        TestEmulator emulator = Emulator(residuals, factorScale: factorScale);
                        CalibrationCovariance covariance = CalibrationCovariance.Create(new FunctionalFormulation { Truncation = term }, Groups(), emulator);
                        CovarianceFactor f = covariance.Factor([0.02 * factorScale * factorScale, 0.08 * factorScale * factorScale], emulator.Noise(0));
                        var pred = new double[M];
                        emulator.Predict(0, [theta], [], pred);
                        double[] r = [.. Enumerable.Range(0, M).Select(j => y[j] * factorScale - pred[j])];
                        return -0.5 * (f.LogDet + f.InverseQuadraticForm(r));
                    }
                    return LogLikelihood(c) - LogLikelihood(1.0);
                })];
                foreach (double shift in shifts) Assert.Equal(-M * Math.Log(c), shift, 1e-9);
            }
    }

    [Fact]
    public void MismatchedOrInvalidDeclarations_AreRefusedAtConstruction()
    {
        (_, _, _, double[,] residuals) = Pieces();
        TestEmulator emulator = Emulator(residuals);
        MeasurementModel measurement = Groups();
        double[][] Rows(int rows, int columns, Func<int, int, double> f)
            => [.. Enumerable.Range(0, rows).Select(i => Enumerable.Range(0, columns).Select(j => f(i, j)).ToArray())];
        double[][] identity = Rows(M, M, (i, j) => i == j ? 1 : 0);
        string Refusal(CovarianceFormulation f, ICalibrationEmulator? e = null)
            => Assert.Throws<ArgumentException>(() => CalibrationCovariance.Create(f, measurement, e ?? emulator)).Message;

        Assert.Contains("predicts 6 outputs", Refusal(new FunctionalFormulation(), new TestEmulator([new ContinuousParameter("theta", 0, 1)],
            new double[6], [.. Enumerable.Repeat(1.0, 6)], new double[6, 1], (_, _, _) => [0.0], [[1.0]])));
        Assert.Contains("rows", Refusal(new FunctionalFormulation { Discrepancy = new FunctionalDiscrepancy(Rows(M - 1, 2, (i, j) => 1), 1.0) }));
        Assert.Contains("no default", Refusal(new FunctionalFormulation { Discrepancy = new FunctionalDiscrepancy(Rows(M, 2, (i, j) => 1), null) }));
        Assert.Contains("positive", Refusal(new FunctionalFormulation { Discrepancy = new FunctionalDiscrepancy(Rows(M, 2, (i, j) => 1), -1.0) }));
        Assert.Contains("mean", Refusal(new FunctionalFormulation { Discrepancy = new FunctionalDiscrepancy(Rows(M, 2, (i, j) => 1), 1.0, new double[3]) }));
        Assert.Contains("not finite", Refusal(new FunctionalFormulation { Discrepancy = new FunctionalDiscrepancy(Rows(M, 2, (i, j) => double.NaN), 1.0) }));
        Assert.Contains("it has none", Assert.Throws<ArgumentException>(() =>
            CalibrationCovariance.Create(new FunctionalFormulation { Truncation = TruncationTerm.Diagonal }, measurement, Emulator(null))).Message);
        Assert.Contains("symmetric", Refusal(new MultivariateFormulation(Rows(M, M, (i, j) => i == j ? 1 : i < j ? 0.1 : 0))));
        Assert.Contains("diagonal", Refusal(new MultivariateFormulation(Rows(M, M, (i, j) => i == j ? 2 : 0))));
        Assert.Contains("positive definite", Refusal(new MultivariateFormulation(Rows(M, M, (i, j) => i == j ? 1 : -0.5))));
        // A discrepancy covariance of another shape is refused, not read as a functional basis (R's branch by shape).
        Assert.Contains($"hold {M} values", Refusal(new MultivariateFormulation(identity, Rows(M, 2, (i, j) => 1))));
        Assert.Contains("positive definite", Refusal(new MultivariateFormulation(identity, Rows(M, M, (i, j) => 1))));

        // A square discrepancy basis is a functional formulation's, with m coefficients: the declared type holds.
        CalibrationCovariance square = CalibrationCovariance.Create(new FunctionalFormulation { Discrepancy = new FunctionalDiscrepancy(identity, 1.0) }, measurement, emulator);
        Assert.Equal(M, square.DiscrepancyBasis!.GetLength(1));
        Assert.Equal(1.0, square.DiscrepancyScale);
    }
}
