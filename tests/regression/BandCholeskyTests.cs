using System;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Band Cholesky is the engine the banded spline marginal consumes: on a symmetric positive-definite matrix
/// with half-bandwidth p it must reproduce the dense <see cref="CholeskyDecomposition"/> exactly — same solve
/// A⁻¹b, same log-determinant — at O(n·p²) instead of O(n³). Both factorizations (L·Lᵀ and root-free L·D·Lᵀ)
/// must agree with the dense reference.
/// </summary>
public sealed class BandCholeskyTests
{
    private readonly ITestOutputHelper _out;
    public BandCholeskyTests(ITestOutputHelper output) => _out = output;

    // Random SPD matrix banded at half-bandwidth p: symmetric band entries, diagonally dominant diagonal.
    private static double[,] BandedSpd(int n, int p, Xoshiro256PlusPlus rng)
    {
        var a = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = Math.Max(0, i - p); j < i; j++)
            {
                double v = rng.NextDouble() - 0.5;
                a[i, j] = v;
                a[j, i] = v;
            }
        for (int i = 0; i < n; i++)
        {
            double off = 0.0;
            for (int j = Math.Max(0, i - p); j <= Math.Min(n - 1, i + p); j++)
                if (j != i) off += Math.Abs(a[i, j]);
            a[i, i] = off + 1.0;   // strict diagonal dominance ⇒ SPD
        }
        return a;
    }

    [Theory]
    [InlineData(30, 3, BandFactorization.Cholesky)]   // cubic-spline band (heptadiagonal)
    [InlineData(30, 3, BandFactorization.Ldlt)]
    [InlineData(50, 1, BandFactorization.Cholesky)]   // tridiagonal
    [InlineData(50, 1, BandFactorization.Ldlt)]
    [InlineData(40, 5, BandFactorization.Cholesky)]
    [InlineData(40, 5, BandFactorization.Ldlt)]
    [InlineData(8, 7, BandFactorization.Cholesky)]    // p ≥ n−1: degenerates to dense
    [InlineData(8, 7, BandFactorization.Ldlt)]
    public void BandCholesky_MatchesDenseSolveAndLogDet(int n, int p, BandFactorization mode)
    {
        var rng = new Xoshiro256PlusPlus(seed: 20260612);
        double[,] a = BandedSpd(n, p, rng);
        var b = new double[n];
        for (int i = 0; i < n; i++) b[i] = rng.NextDouble() - 0.5;

        var dense = new CholeskyDecomposition(n);
        dense.Decompose(a);
        var aInv = new double[n, n];
        dense.WriteInverseTo(aInv);
        var xDense = new double[n];
        for (int i = 0; i < n; i++)
        {
            double s = 0.0;
            for (int k = 0; k < n; k++) s += aInv[i, k] * b[k];
            xDense[i] = s;
        }

        var band = new BandCholesky(n, p, mode);
        band.Decompose(a);
        double[] xBand = band.Solve(b);

        double maxErr = 0.0;
        for (int i = 0; i < n; i++) maxErr = Math.Max(maxErr, Math.Abs(xBand[i] - xDense[i]));
        _out.WriteLine($"n={n} p={p} {mode} maxSolveErr={maxErr:E3} logDet band={band.LogDet:F6} dense={dense.LogDet:F6}");

        Assert.True(maxErr < 1e-9, $"{mode} band solve deviates from dense by {maxErr:E3}");
        Assert.Equal(dense.LogDet, band.LogDet, 9);
    }

    [Fact]
    public void CholeskyDecomposition_Solve_MatchesInverseAndResidual()
    {
        const int n = 24;
        var rng = new Xoshiro256PlusPlus(seed: 12321);
        double[,] a = BandedSpd(n, n - 1, rng);   // p = n−1 ⇒ a fully dense SPD matrix
        var b = new double[n];
        for (int i = 0; i < n; i++) b[i] = rng.NextDouble() - 0.5;

        var chol = new CholeskyDecomposition(n);
        chol.Decompose(a);
        double[] x = chol.Solve(b);

        // Agree with A⁻¹b from the full inverse, and have a tiny residual A·x − b.
        var aInv = new double[n, n];
        chol.WriteInverseTo(aInv);
        double maxErr = 0.0, maxResid = 0.0;
        for (int i = 0; i < n; i++)
        {
            double inv = 0.0, ax = 0.0;
            for (int k = 0; k < n; k++) { inv += aInv[i, k] * b[k]; ax += a[i, k] * x[k]; }
            maxErr = Math.Max(maxErr, Math.Abs(x[i] - inv));
            maxResid = Math.Max(maxResid, Math.Abs(ax - b[i]));
        }
        _out.WriteLine($"[dense-solve] maxErr vs inverse·b={maxErr:E3}, max residual |Ax−b|={maxResid:E3}");

        Assert.True(maxErr < 1e-9, $"Solve disagrees with inverse·b by {maxErr:E3}");
        Assert.True(maxResid < 1e-9, $"Solve residual {maxResid:E3}");
    }

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void SampleInnovation_HasInverseCovariance()
    {
        const int n = 12, p = 3;
        var rng = new Xoshiro256PlusPlus(seed: 4242);
        double[,] a = BandedSpd(n, p, rng);

        var dense = new CholeskyDecomposition(n);
        dense.Decompose(a);
        var aInv = new double[n, n];
        dense.WriteInverseTo(aInv);

        var band = new BandCholesky(n, p, BandFactorization.Ldlt);
        band.Decompose(a);

        // Empirical covariance of v = L⁻ᵀD^{-½}z must match A⁻¹.
        const int m = 100_000;
        var cov = new double[n, n];
        var z = new double[n];
        var v = new double[n];
        for (int s = 0; s < m; s++)
        {
            for (int i = 0; i < n; i++) z[i] = Gaussian(rng);
            band.SampleInnovation(z, v);
            for (int i = 0; i < n; i++)
                for (int j = 0; j <= i; j++)
                    cov[i, j] += v[i] * v[j];
        }

        double maxErr = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
                maxErr = Math.Max(maxErr, Math.Abs(cov[i, j] / m - aInv[i, j]));
        _out.WriteLine($"SampleInnovation cov maxErr={maxErr:E3}");
        Assert.True(maxErr < 0.02, $"empirical innovation covariance deviates from A⁻¹ by {maxErr:E3}");
    }
}
