using System;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;
using static Maths.LinAlg.Tests.FactorReferenceTests;

namespace Maths.LinAlg.Tests;

/// <summary>
/// The thin SVD by its defining properties: A = U·diag(σ)·Vᵀ to working precision, orthonormal U and V, σ descending
/// and nonnegative. Against constructed spectra it recovers each σ to an absolute ε·σ₁, including values eleven
/// orders below σ₁ that a Gram eigendecomposition would lose, and it reads the numerical rank of products of known
/// rank, of centred curves (the BASS PCA shape, one rank lost to centring) and of zero matrices, completing U where
/// σ is exactly zero.
/// </summary>
public sealed class ThinSvdTests
{
    private const double Epsilon = 2.220446049250313e-16;

    [Theory]
    [InlineData(1, 1)]
    [InlineData(5, 1)]
    [InlineData(1, 5)]
    [InlineData(8, 8)]
    [InlineData(40, 7)]
    [InlineData(7, 40)]
    [InlineData(200, 5)]
    public void Decomposition_ReconstructsTheMatrix_WithOrthonormalFactors(int m, int n)
    {
        var rng = new Xoshiro256PlusPlus(seed: 800 + 10 * m + n);
        double[,] a = RandomMatrix(m, n, rng);
        ThinSvdResult svd = ThinSvd.Decompose(a);
        AssertValid(a, svd);
        Assert.Equal(Math.Min(m, n), svd.NumericalRank);
    }

    [Fact]
    public void KnownSpectrum_IsRecoveredToEpsilonTimesTheLargest()
    {
        const int m = 30, n = 12;
        var rng = new Xoshiro256PlusPlus(seed: 811);
        var sigma = new double[n];
        for (int k = 0; k < n; k++) sigma[k] = Math.Pow(10.0, -k);       // 1 down to 10⁻¹¹
        double[,] u = Multiply(Householder(Vector(m, rng)), Householder(Vector(m, rng)));
        double[,] v = Multiply(Householder(Vector(n, rng)), Householder(Vector(n, rng)));
        var a = new double[m, n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                for (int k = 0; k < n; k++) a[i, j] += u[i, k] * sigma[k] * v[j, k];

        ThinSvdResult svd = ThinSvd.Decompose(a);
        AssertValid(a, svd);
        for (int k = 0; k < n; k++)
            Assert.True(Math.Abs(svd.SingularValues[k] - sigma[k]) <= 10.0 * m * Epsilon, $"σ[{k}] = {svd.SingularValues[k]:R}, expected {sigma[k]:R}");
        Assert.Equal(n, svd.NumericalRank);
    }

    [Fact]
    public void ProductOfKnownRank_HasThatNumericalRank()
    {
        const int m = 50, n = 20, r = 7;
        var rng = new Xoshiro256PlusPlus(seed: 812);
        double[,] a = Multiply(RandomMatrix(m, r, rng), RandomMatrix(r, n, rng));
        ThinSvdResult svd = ThinSvd.Decompose(a);
        AssertValid(a, svd);
        Assert.Equal(r, svd.NumericalRank);
        Assert.True(svd.SingularValues[r - 1] > 1e6 * svd.RankTolerance);
    }

    [Fact]
    public void CentredCurves_LoseOneRankToCentring()
    {
        // Five curves on 200 locations, centred across cases: Y_cᵀ is 200 × 5 of rank 4, as BASS's PCA sees it.
        const int cases = 5, locations = 200;
        var rng = new Xoshiro256PlusPlus(seed: 813);
        var y = new double[cases, locations];
        for (int c = 0; c < cases; c++)
        {
            double phase = rng.NextDouble(), amplitude = 1.0 + rng.NextDouble();
            for (int t = 0; t < locations; t++)
            {
                double s = (double)t / (locations - 1);
                y[c, t] = amplitude * Math.Sin(2.0 * Math.PI * (s - 0.3 * phase)) + 0.01 * (rng.NextDouble() - 0.5);
            }
        }
        var yct = new double[locations, cases];
        for (int t = 0; t < locations; t++)
        {
            double mean = 0.0;
            for (int c = 0; c < cases; c++) mean += y[c, t];
            mean /= cases;
            for (int c = 0; c < cases; c++) yct[t, c] = y[c, t] - mean;
        }
        ThinSvdResult svd = ThinSvd.Decompose(yct);
        AssertValid(yct, svd);
        Assert.Equal(cases - 1, svd.NumericalRank);
        Assert.True(svd.SingularValues[cases - 1] <= svd.RankTolerance);
    }

    [Theory]
    [InlineData(6, 4)]
    [InlineData(4, 6)]
    public void ZeroMatrix_HasRankZero_AndOrthonormalFactors(int m, int n)
    {
        var a = new double[m, n];
        ThinSvdResult svd = ThinSvd.Decompose(a);
        AssertValid(a, svd);
        Assert.Equal(0, svd.NumericalRank);
        Assert.Equal(0.0, svd.RankTolerance);
    }

    [Fact]
    public void ZeroColumn_IsCompletedInTheMiddle()
    {
        var rng = new Xoshiro256PlusPlus(seed: 814);
        double[,] a = RandomMatrix(9, 4, rng);
        for (int i = 0; i < 9; i++) a[i, 2] = 0.0;
        ThinSvdResult svd = ThinSvd.Decompose(a);
        AssertValid(a, svd);
        Assert.Equal(3, svd.NumericalRank);
        Assert.Equal(0.0, svd.SingularValues[3]);
    }

    private static void AssertValid(double[,] a, ThinSvdResult svd)
    {
        int m = a.GetLength(0), n = a.GetLength(1), k = Math.Min(m, n);
        Assert.Equal(k, svd.SingularValues.Length);
        Assert.Equal(k, svd.LeftVectors.Length);
        Assert.Equal(k, svd.RightVectors.Length);
        double s1 = k > 0 ? svd.SingularValues[0] : 0.0;
        for (int t = 0; t < k; t++)
        {
            Assert.True(svd.SingularValues[t] >= 0.0);
            if (t > 0) Assert.True(svd.SingularValues[t] <= svd.SingularValues[t - 1]);
            Assert.Equal(m, svd.LeftVectors[t].Length);
            Assert.Equal(n, svd.RightVectors[t].Length);
        }
        AssertOrthonormal(svd.LeftVectors);
        AssertOrthonormal(svd.RightVectors);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
            {
                double s = 0.0;
                for (int t = 0; t < k; t++) s += svd.LeftVectors[t][i] * svd.SingularValues[t] * svd.RightVectors[t][j];
                Assert.True(Math.Abs(s - a[i, j]) <= 10.0 * Math.Max(m, n) * Epsilon * s1, $"(UΣVᵀ)[{i},{j}] = {s:R} vs {a[i, j]:R}");
            }
    }

    private static void AssertOrthonormal(double[][] q)
    {
        for (int i = 0; i < q.Length; i++)
            for (int j = 0; j <= i; j++)
            {
                double d = PairwiseSum.Dot(q[i], q[j]);
                Assert.True(Math.Abs(d - (i == j ? 1.0 : 0.0)) <= 1e-13, $"q[{i}]·q[{j}] = {d:R}");
            }
    }

    private static double[,] RandomMatrix(int m, int n, Xoshiro256PlusPlus rng)
    {
        var a = new double[m, n];
        for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) a[i, j] = rng.NextDouble() - 0.5;
        return a;
    }
}
