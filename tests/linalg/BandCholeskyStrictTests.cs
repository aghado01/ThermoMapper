using System;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;

namespace Maths.LinAlg.Tests;

/// <summary>
/// The strict band factor (step 3.7): no pivot floor, a stop at the first pivot that is not finite and positive, and
/// solves that refuse after it, as <see cref="DenseCholesky"/>. On a matrix whose pivots all clear the floor it computes
/// the floored factor's numbers exactly.
/// </summary>
public sealed class BandCholeskyStrictTests
{
    // A = L·Lᵀ for a random lower-banded L with a positive diagonal: symmetric positive definite, half-bandwidth p.
    private static double[,] BandSpd(int n, int p, Xoshiro256PlusPlus rng)
    {
        var l = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = Math.Max(0, i - p); j <= i; j++)
                l[i, j] = i == j ? 1.0 + rng.NextDouble() : rng.NextDouble() - 0.5;
        var a = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                double s = 0.0;
                for (int k = 0; k < n; k++) s += l[i, k] * l[j, k];
                a[i, j] = s;
            }
        return a;
    }

    [Theory]
    [InlineData(BandFactorization.Cholesky)]
    [InlineData(BandFactorization.Ldlt)]
    public void PositiveDefinite_MatchesTheFlooredFactorExactly(BandFactorization mode)
    {
        var rng = new Xoshiro256PlusPlus(seed: 3);
        double[,] a = BandSpd(12, 3, rng);
        var floored = new BandCholesky(12, 3, mode);
        var strict = new BandCholesky(12, 3, mode, strict: true);
        floored.Decompose(a);
        strict.Decompose(a);

        Assert.True(strict.Strict);
        Assert.Equal(-1, strict.FailedAt);
        Assert.Equal(floored.Status, strict.Status);
        Assert.Equal(floored.LogDet, strict.LogDet);
        double[] b = new double[12];
        for (int i = 0; i < 12; i++) b[i] = rng.NextDouble();
        Assert.Equal(floored.Solve(b), strict.Solve(b));
        if (mode == BandFactorization.Ldlt)
        {
            var u = new double[12];
            var v = new double[12];
            floored.SampleInnovation(b, u);
            strict.SampleInnovation(b, v);
            Assert.Equal(u, v);
        }
    }

    [Theory]
    [InlineData(BandFactorization.Cholesky, 1.0, 1.0, 1.0, 1, 0.0)]    // singular: the second pivot is exactly 0
    [InlineData(BandFactorization.Ldlt, 1.0, 2.0, 1.0, 1, -3.0)]       // indefinite: the second pivot is 1 − 4
    [InlineData(BandFactorization.Ldlt, -1.0, 0.0, 1.0, 0, -1.0)]      // negative first diagonal
    public void NonPositiveDefinite_StopsAtTheFailingPivot_AndRefusesToSolve(
        BandFactorization mode, double a00, double a10, double a11, int failedAt, double pivot)
    {
        var f = new BandCholesky(2, 1, mode, strict: true);
        f.Decompose(new double[,] { { a00, a10 }, { a10, a11 } });
        Assert.True(f.Status.Finite);
        Assert.False(f.Status.PositiveDefinite);
        Assert.Equal(0, f.Status.FlooredPivots);
        Assert.Equal(failedAt, f.FailedAt);
        Assert.Equal(pivot, f.Status.MinPivot);
        Assert.True(double.IsNaN(f.LogDet));
        Assert.Throws<InvalidOperationException>(() => f.Solve([1.0, 1.0]));
        if (mode == BandFactorization.Ldlt)
            Assert.Throws<InvalidOperationException>(() => f.SampleInnovation([1.0, 1.0], new double[2]));

        f.Decompose(new double[,] { { 2.0, 0.5 }, { 0.5, 2.0 } });   // a later factorization starts afresh
        Assert.Equal(-1, f.FailedAt);
        Assert.True(f.Status.PositiveDefinite);
        Assert.Equal(2, f.Solve([1.0, 1.0]).Length);
    }

    [Fact]
    public void NonFiniteEntry_FailsWhereItEntersAPivot()
    {
        var f = new BandCholesky(3, 1, BandFactorization.Ldlt, strict: true);
        f.Decompose(new double[,] { { 2, 0, 0 }, { 0, 2, double.NaN }, { 0, double.NaN, 2 } });
        Assert.False(f.Status.Finite);
        Assert.Equal(2, f.FailedAt);
        Assert.Throws<InvalidOperationException>(() => f.Solve([1.0, 1.0, 1.0]));
    }

    [Theory]
    [InlineData(BandFactorization.Cholesky)]
    [InlineData(BandFactorization.Ldlt)]
    public void TinyScaleMatrix_IsFactoredWithoutFlooring(BandFactorization mode)
    {
        // 1e-20·(tridiagonal 2, −1): perfectly conditioned, but every pivot is far below the 1e-12 floor.
        int n = 5;
        var a = new double[n, n];
        for (int i = 0; i < n; i++)
        {
            a[i, i] = 2e-20;
            if (i > 0) { a[i, i - 1] = -1e-20; a[i - 1, i] = -1e-20; }
        }
        var floored = new BandCholesky(n, 1, mode);
        var strict = new BandCholesky(n, 1, mode, strict: true);
        floored.Decompose(a);
        strict.Decompose(a);

        Assert.Equal(n, floored.Status.FlooredPivots);
        Assert.Equal(0, strict.Status.FlooredPivots);
        Assert.True(strict.Status.PositiveDefinite);
        Assert.True(strict.Status.MinRelativePivot > 0.5);   // the equilibrated pivots are (k+1)/(2k) ≥ 0.6
        double[] x = strict.Solve([1e-20, 0.0, 0.0, 0.0, 1e-20]);   // A·[1,1,1,1,1] = [1e-20, 0, 0, 0, 1e-20]
        foreach (double v in x) Assert.Equal(1.0, v, 10);
    }
}
