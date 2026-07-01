using System;

namespace Maths.LinAlg;

public enum DenseEigenFastVariant
{
    Default = 0,
    Fma = 1,
}

public readonly record struct DenseEigenOptions(
    DenseEigenFastVariant FastVariant = DenseEigenFastVariant.Default);

/// <summary>
/// Compile-time dispatch point for dense symmetric eigendecomposition.
/// </summary>
public static class DenseEigen
{
    public static EigenResult DecomposeSymmetric(
        double[,] matrix,
        int maxSweeps = 256,
        double tol = 1e-12,
        DenseEigenOptions options = default)
    {
#if EIGEN_REFERENCE
        EnsureReferenceCompatible(options);
        return Eigen.DecomposeSymmetric(matrix, maxSweeps, tol);
#else
        return EigenFast.DecomposeSymmetric(matrix, maxSweeps, tol, options.FastVariant);
#endif
    }

    public static EigenResult DecomposeSymmetric(
        ReadOnlySpan<double> flatColumnMajorMatrix,
        int n,
        int maxSweeps = 256,
        double tol = 1e-12,
        DenseEigenOptions options = default)
    {
        if (flatColumnMajorMatrix.Length != n * n)
            throw new ArgumentException("Matrix dimensions must match N x N.", nameof(flatColumnMajorMatrix));

#if EIGEN_REFERENCE
        EnsureReferenceCompatible(options);

        var matrix = new double[n, n];
        for (int c = 0; c < n; c++)
        {
            for (int r = 0; r < n; r++)
                matrix[r, c] = flatColumnMajorMatrix[c * n + r];
        }

        return Eigen.DecomposeSymmetric(matrix, maxSweeps, tol);
#else
        return EigenFast.DecomposeSymmetric(flatColumnMajorMatrix, n, maxSweeps, tol, options.FastVariant);
#endif
    }

    private static void EnsureReferenceCompatible(DenseEigenOptions options)
    {
        if (options.FastVariant != DenseEigenFastVariant.Default)
        {
            throw new NotSupportedException(
                "Fast-family variants are unavailable when EIGEN_REFERENCE is active.");
        }
    }
}
