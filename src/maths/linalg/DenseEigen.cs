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
/// How a dense Jacobi solve ended, read from the diagonalization it returns (spectral-solver brief, A2). A solve stops after a sweep
/// that rotates nothing. <see cref="Converged"/> holds when every off-diagonal entry passes the solver's pivot test, so it can hold
/// while <see cref="OffDiagonalNorm"/> exceeds the tolerance: many entries, each below it.
/// </summary>
/// <param name="Sweeps">Sweeps run, the last one included.</param>
/// <param name="Rotations">Rotations applied over all sweeps.</param>
/// <param name="OffDiagonalNorm">off(A) = (Σ_{i≠j} a_ij²)^½ of the returned diagonalization VᵀAV, as the solver holds it.</param>
/// <param name="LargestOffDiagonal">max_{i≠j} |a_ij| of the returned diagonalization.</param>
/// <param name="Converged">Whether every off-diagonal entry passes the pivot test: below the tolerance for <see cref="EigenFast"/>, at
/// or below it for the reference <see cref="Eigen"/>.</param>
public readonly record struct DenseEigenConvergence(int Sweeps, long Rotations, double OffDiagonalNorm, double LargestOffDiagonal, bool Converged)
{
    // off(A) and the largest entry from the strict upper triangle of a symmetric A, read by entry(row, column); scaled by the
    // largest entry so the sum of squares stays in range.
    internal static (double Norm, double Largest) OffDiagonal(int n, Func<int, int, double> entry)
    {
        double largest = 0.0;
        for (int q = 1; q < n; q++)
            for (int p = 0; p < q; p++) largest = Math.Max(largest, Math.Abs(entry(p, q)));
        if (!(largest > 0.0) || double.IsInfinity(largest)) return (largest, largest);
        double squares = 0.0;
        for (int q = 1; q < n; q++)
            for (int p = 0; p < q; p++)
            {
                double scaled = entry(p, q) / largest;
                squares += scaled * scaled;
            }
        return (largest * Math.Sqrt(2.0 * squares), largest);
    }
}

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
