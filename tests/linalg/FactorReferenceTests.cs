using System;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;

namespace Maths.LinAlg.Tests;

/// <summary>
/// B06 (R11): the dense and banded factorizations against references that share no code with them — an LU
/// solve with partial pivoting, the residual A·x − b, A·A⁻¹ = I, L·Lᵀ = A, and matrices with known eigenvalues —
/// plus the <see cref="FactorStatus"/> each reports on singular, tiny-scale, ill-conditioned and non-finite input.
/// </summary>
public sealed class FactorReferenceTests
{
    // ── Dense ──────────────────────────────────────────────────────────────────

    [Fact]
    public void Dense_SolveLogDetInverseAndFactor_MatchIndependentReferences()
    {
        const int n = 12;
        double[,] a = DenseSpd(n, new Xoshiro256PlusPlus(seed: 101));
        double[] b = Vector(n, new Xoshiro256PlusPlus(seed: 102));

        var chol = new CholeskyDecomposition(n);
        chol.Decompose(a);

        (double[] xRef, double logDetRef) = LuSolve(a, b);
        AssertClose(xRef, chol.Solve(b), 1e-10);
        Assert.Equal(logDetRef, chol.LogDet, 10);

        var inv = new double[n, n];
        chol.WriteInverseTo(inv);
        AssertIdentity(Multiply(a, inv), 1e-10);

        var lFlat = new double[n * n];
        chol.WriteLTo(lFlat);
        var l = new double[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) l[i, j] = lFlat[i * n + j];
        AssertClose(Flatten(a), Flatten(Multiply(l, Transpose(l))), 1e-10);

        Assert.True(chol.Status.PositiveDefinite);
        Assert.Equal(0, chol.Status.FlooredPivots);
    }

    // ── Banded ─────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(30, 3, BandFactorization.Cholesky)]
    [InlineData(30, 3, BandFactorization.Ldlt)]
    [InlineData(25, 1, BandFactorization.Ldlt)]
    [InlineData(8, 7, BandFactorization.Cholesky)]
    public void Band_SolveAndLogDet_MatchLu(int n, int p, BandFactorization mode)
    {
        double[,] a = BandedSpd(n, p, new Xoshiro256PlusPlus(seed: 201));
        double[] b = Vector(n, new Xoshiro256PlusPlus(seed: 202));

        var band = new BandCholesky(n, p, mode);
        band.Decompose(a);

        (double[] xRef, double logDetRef) = LuSolve(a, b);
        AssertClose(xRef, band.Solve(b), 1e-10);
        Assert.Equal(logDetRef, band.LogDet, 10);
        Assert.True(band.Status.PositiveDefinite);
        Assert.Equal(0, band.Status.FlooredPivots);
    }

    [Fact]
    public void Band_InnovationMap_HasExactlyTheInverseCovariance()
    {
        // v = M·z with M = L⁻ᵀD^{−½}; applying the map to each unit vector gives M's columns, and A·M·Mᵀ = I.
        const int n = 14, p = 3;
        double[,] a = BandedSpd(n, p, new Xoshiro256PlusPlus(seed: 301));
        var band = new BandCholesky(n, p, BandFactorization.Ldlt);
        band.Decompose(a);

        var m = new double[n, n];
        var e = new double[n];
        var col = new double[n];
        for (int k = 0; k < n; k++)
        {
            Array.Clear(e);
            e[k] = 1.0;
            band.SampleInnovation(e, col);
            for (int i = 0; i < n; i++) m[i, k] = col[i];
        }
        AssertIdentity(Multiply(a, Multiply(m, Transpose(m))), 1e-10);
    }

    // ── Status ─────────────────────────────────────────────────────────────────

    public static TheoryData<string> Factorizations => new() { "dense", "band-chol", "band-ldlt" };

    [Theory]
    [MemberData(nameof(Factorizations))]
    public void Status_SingularMatrix_IsNotPositiveDefinite(string kind)
    {
        FactorStatus s = Factor(kind, new double[,] { { 1, 1 }, { 1, 1 } });
        Assert.True(s.Finite);
        Assert.False(s.PositiveDefinite);
        Assert.Equal(1, s.FlooredPivots);
    }

    [Theory]
    [MemberData(nameof(Factorizations))]
    public void Status_TinyScaleIdentity_IsFlooredButPerfectlyConditioned(string kind)
    {
        // The absolute floor clamps every pivot of 1e-20·I; the pivot ratio shows it is a scale effect, not rank.
        FactorStatus s = Factor(kind, new double[,] { { 1e-20, 0 }, { 0, 1e-20 } });
        Assert.True(s.PositiveDefinite);
        Assert.Equal(2, s.FlooredPivots);
        Assert.Equal(1.0, s.PivotRatio);
    }

    [Theory]
    [MemberData(nameof(Factorizations))]
    public void Status_IllConditionedMatrix_ShowsASmallPivotRatioWithoutFlooring(string kind)
    {
        FactorStatus s = Factor(kind, new double[,] { { 1, 0 }, { 0, 1e-9 } });
        Assert.True(s.PositiveDefinite);
        Assert.Equal(0, s.FlooredPivots);
        Assert.Equal(1e-9, s.PivotRatio, 15);
    }

    [Theory]
    [MemberData(nameof(Factorizations))]
    public void Status_NonFiniteEntry_IsNotFinite(string kind)
    {
        FactorStatus s = Factor(kind, new double[,] { { 2, 0 }, { double.NaN, 2 } });
        Assert.False(s.Finite);
        Assert.False(s.PositiveDefinite);
    }

    [Fact]
    public void PivotRatio_IsBoundedBelowByTheReciprocalConditionNumber()
    {
        // A = Q·diag(λ)·Qᵀ with a Householder Q: κ(A) = λ_max/λ_min exactly, and every pivot lies in [λ_min, λ_max].
        double[] lambda = [1e-6, 0.01, 0.3, 1.0, 2.5, 40.0];
        int n = lambda.Length;
        double[,] q = Householder(Vector(n, new Xoshiro256PlusPlus(seed: 401)));
        var d = new double[n, n];
        for (int i = 0; i < n; i++) d[i, i] = lambda[i];
        double[,] a = Multiply(q, Multiply(d, Transpose(q)));

        var chol = new CholeskyDecomposition(n);
        chol.Decompose(a);
        double kappa = 40.0 / 1e-6;
        Assert.True(chol.Status.PivotRatio >= (1.0 / kappa) * (1.0 - 1e-6), $"ratio {chol.Status.PivotRatio:E3} < 1/κ");
        Assert.True(chol.Status.MinPivot >= 1e-6 * (1.0 - 1e-6) && chol.Status.MaxPivot <= 40.0 * (1.0 + 1e-9));
    }

    // ── Helpers ────────────────────────────────────────────────────────────────

    private static FactorStatus Factor(string kind, double[,] a)
    {
        int n = a.GetLength(0);
        switch (kind)
        {
            case "dense":
                var chol = new CholeskyDecomposition(n);
                chol.Decompose(a);
                return chol.Status;
            default:
                var band = new BandCholesky(n, n - 1, kind == "band-ldlt" ? BandFactorization.Ldlt : BandFactorization.Cholesky);
                band.Decompose(a);
                return band.Status;
        }
    }

    internal static double[,] DenseSpd(int n, Xoshiro256PlusPlus rng)
    {
        var b = new double[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) b[i, j] = rng.NextDouble() - 0.5;
        double[,] a = Multiply(b, Transpose(b));
        for (int i = 0; i < n; i++) a[i, i] += 0.5;
        return a;
    }

    private static double[,] BandedSpd(int n, int p, Xoshiro256PlusPlus rng)
    {
        var a = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = Math.Max(0, i - p); j < i; j++) { double v = rng.NextDouble() - 0.5; a[i, j] = v; a[j, i] = v; }
        for (int i = 0; i < n; i++)
        {
            double off = 0.0;
            for (int j = 0; j < n; j++) if (j != i) off += Math.Abs(a[i, j]);
            a[i, i] = off + 1.0;
        }
        return a;
    }

    internal static double[] Vector(int n, Xoshiro256PlusPlus rng)
    {
        var v = new double[n];
        for (int i = 0; i < n; i++) v[i] = rng.NextDouble() - 0.5;
        return v;
    }

    internal static double[,] Householder(double[] v)
    {
        int n = v.Length;
        double vv = 0.0;
        foreach (double x in v) vv += x * x;
        var q = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) q[i, j] = (i == j ? 1.0 : 0.0) - 2.0 * v[i] * v[j] / vv;
        return q;
    }

    // Gaussian elimination with partial pivoting: x = A⁻¹b and ln|det A|.
    internal static (double[] X, double LogAbsDet) LuSolve(double[,] a0, double[] b0)
    {
        int n = b0.Length;
        var a = (double[,])a0.Clone();
        var b = (double[])b0.Clone();
        double logDet = 0.0;
        for (int c = 0; c < n; c++)
        {
            int piv = c;
            for (int r = c + 1; r < n; r++) if (Math.Abs(a[r, c]) > Math.Abs(a[piv, c])) piv = r;
            if (piv != c)
            {
                for (int k = 0; k < n; k++) (a[c, k], a[piv, k]) = (a[piv, k], a[c, k]);
                (b[c], b[piv]) = (b[piv], b[c]);
            }
            logDet += Math.Log(Math.Abs(a[c, c]));
            for (int r = c + 1; r < n; r++)
            {
                double f = a[r, c] / a[c, c];
                for (int k = c; k < n; k++) a[r, k] -= f * a[c, k];
                b[r] -= f * b[c];
            }
        }
        var x = new double[n];
        for (int r = n - 1; r >= 0; r--)
        {
            double s = b[r];
            for (int k = r + 1; k < n; k++) s -= a[r, k] * x[k];
            x[r] = s / a[r, r];
        }
        return (x, logDet);
    }

    internal static double[,] Multiply(double[,] x, double[,] y)
    {
        int n = x.GetLength(0), m = y.GetLength(1), k = x.GetLength(1);
        var r = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
            {
                double s = 0.0;
                for (int t = 0; t < k; t++) s += x[i, t] * y[t, j];
                r[i, j] = s;
            }
        return r;
    }

    internal static double[,] Transpose(double[,] x)
    {
        int n = x.GetLength(0), m = x.GetLength(1);
        var r = new double[m, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) r[j, i] = x[i, j];
        return r;
    }

    private static double[] Flatten(double[,] x)
    {
        var r = new double[x.Length];
        int i = 0;
        foreach (double v in x) r[i++] = v;
        return r;
    }

    internal static void AssertClose(double[] expected, double[] actual, double tol)
    {
        Assert.Equal(expected.Length, actual.Length);
        for (int i = 0; i < expected.Length; i++)
            Assert.True(Math.Abs(expected[i] - actual[i]) <= tol * (1.0 + Math.Abs(expected[i])),
                        $"index {i}: expected {expected[i]:R}, got {actual[i]:R}");
    }

    internal static void AssertIdentity(double[,] x, double tol)
    {
        int n = x.GetLength(0);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                Assert.True(Math.Abs(x[i, j] - (i == j ? 1.0 : 0.0)) <= tol, $"[{i},{j}] = {x[i, j]:R}");
    }
}
