using System;
using Maths.LinAlg;
using Xunit;

namespace Maths.LinAlg.Tests;

/// <summary>
/// Functions of a symmetric matrix from its eigensystem (spectral-solver brief, B6): f(A) = V·f(Λ)·Vᵀ against a matrix of known
/// spectrum, the square root squaring back, the inverse inverting, and the factor forms reproducing A as F·Fᵀ and Rᵀ·R over the kept
/// eigenpairs. The selection is the caller's: a dropped eigenpair leaves its share out.
/// </summary>
public sealed class MatrixFunctionsTests
{
    // Q·diag(4, 1, 0.25)·Qᵀ for a rotation Q, so √A has eigenvalues 2, 1, 0.5 and A⁻¹ has 0.25, 1, 4.
    private static readonly double[,] Q = Rotation(0.4, 1.1);
    private static readonly double[,] A = Compose(Q, [4.0, 1.0, 0.25]);

    [Fact]
    public void Symmetric_IsTheFunctionOfTheSpectrum()
    {
        EigenResult e = Eigen.DecomposeSymmetric(A);
        double[,] root = MatrixFunctions.Symmetric(e, Math.Sqrt);
        double[,] inverse = MatrixFunctions.Symmetric(e, x => 1.0 / x);

        AssertClose(Compose(Q, [2.0, 1.0, 0.5]), root, 1e-13);
        AssertClose(A, Multiply(root, root), 1e-13);
        AssertClose(Identity(3), Multiply(A, inverse), 1e-13);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) Assert.Equal(root[i, j], root[j, i]);
    }

    /// <summary>A zero weight drops its term; a diagonal matrix's eigenvectors are coordinate axes and its root is exact.</summary>
    [Fact]
    public void Symmetric_SkipsZeroTerms_ExactOnADiagonal()
    {
        EigenResult diagonal = Eigen.DecomposeSymmetric(new[,] { { 9.0, 0.0 }, { 0.0, 0.0 } });
        double[,] root = MatrixFunctions.Symmetric(diagonal, x => Math.Sqrt(Math.Max(0.0, x)));

        Assert.Equal(new[,] { { 3.0, 0.0 }, { 0.0, 0.0 } }, root);
    }

    [Fact]
    public void Factors_ReproduceTheMatrixOverTheKeptEigenpairs()
    {
        EigenResult e = Eigen.DecomposeSymmetric(A);
        double[,] columns = MatrixFunctions.FactorColumns(e, [0, 1, 2], Math.Sqrt);
        double[,] rows = MatrixFunctions.FactorRows(e, [0, 1, 2], Math.Sqrt);
        double[,] leading = MatrixFunctions.FactorColumns(e, [0], Math.Sqrt);

        Assert.Equal(3, columns.GetLength(1));
        Assert.Equal(3, rows.GetLength(0));
        AssertClose(A, Multiply(columns, Transpose(columns)), 1e-13);
        AssertClose(A, Multiply(Transpose(rows), rows), 1e-13);
        AssertClose(Compose(Q, [4.0, 0.0, 0.0]), Multiply(leading, Transpose(leading)), 1e-13);
    }

    private static double[,] Rotation(double a, double b)
    {
        double ca = Math.Cos(a), sa = Math.Sin(a), cb = Math.Cos(b), sb = Math.Sin(b);
        double[,] x = { { ca, -sa, 0.0 }, { sa, ca, 0.0 }, { 0.0, 0.0, 1.0 } };
        double[,] y = { { 1.0, 0.0, 0.0 }, { 0.0, cb, -sb }, { 0.0, sb, cb } };
        return Multiply(x, y);
    }

    private static double[,] Compose(double[,] q, double[] spectrum)
    {
        int n = spectrum.Length;
        var m = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int k = 0; k < n; k++) m[i, j] += q[i, k] * spectrum[k] * q[j, k];
        return m;
    }

    private static double[,] Multiply(double[,] a, double[,] b)
    {
        int n = a.GetLength(0), inner = a.GetLength(1), m = b.GetLength(1);
        var c = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int k = 0; k < inner; k++)
                for (int j = 0; j < m; j++) c[i, j] += a[i, k] * b[k, j];
        return c;
    }

    private static double[,] Transpose(double[,] a)
    {
        var t = new double[a.GetLength(1), a.GetLength(0)];
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) t[j, i] = a[i, j];
        return t;
    }

    private static double[,] Identity(int n)
    {
        var m = new double[n, n];
        for (int i = 0; i < n; i++) m[i, i] = 1.0;
        return m;
    }

    private static void AssertClose(double[,] expected, double[,] actual, double tolerance)
    {
        for (int i = 0; i < expected.GetLength(0); i++)
            for (int j = 0; j < expected.GetLength(1); j++)
                Assert.InRange(Math.Abs(expected[i, j] - actual[i, j]), 0.0, tolerance);
    }
}
