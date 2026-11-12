using System;
using System.Buffers;

namespace Maths.LinAlg;

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
    // off(A) and the largest entry from the strict upper triangle of a symmetric A held flat column-major, (p, q) at q·n + p; scaled
    // by the largest entry so the sum of squares stays in range.
    internal static (double Norm, double Largest) OffDiagonal(ReadOnlySpan<double> columnMajor, int n)
    {
        double largest = 0.0;
        for (int q = 1; q < n; q++)
            for (int p = 0; p < q; p++) largest = Math.Max(largest, Math.Abs(columnMajor[q * n + p]));
        if (!(largest > 0.0) || double.IsInfinity(largest)) return (largest, largest);
        double squares = 0.0;
        for (int q = 1; q < n; q++)
            for (int p = 0; p < q; p++)
            {
                double scaled = columnMajor[q * n + p] / largest;
                squares += scaled * scaled;
            }
        return (largest * Math.Sqrt(2.0 * squares), largest);
    }

    // The same from the strict upper triangle of a symmetric double[,].
    internal static (double Norm, double Largest) OffDiagonal(double[,] a)
    {
        int n = a.GetLength(0);
        double largest = 0.0;
        for (int q = 1; q < n; q++)
            for (int p = 0; p < q; p++) largest = Math.Max(largest, Math.Abs(a[p, q]));
        if (!(largest > 0.0) || double.IsInfinity(largest)) return (largest, largest);
        double squares = 0.0;
        for (int q = 1; q < n; q++)
            for (int p = 0; p < q; p++)
            {
                double scaled = a[p, q] / largest;
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
        double tol = 1e-12)
    {
#if EIGEN_REFERENCE
        return Eigen.DecomposeSymmetric(matrix, maxSweeps, tol);
#else
        return EigenFast.DecomposeSymmetric(matrix, maxSweeps, tol);
#endif
    }

    public static EigenResult DecomposeSymmetric(
        ReadOnlySpan<double> flatColumnMajorMatrix,
        int n,
        int maxSweeps = 256,
        double tol = 1e-12)
    {
        if (flatColumnMajorMatrix.Length != n * n)
            throw new ArgumentException("Matrix dimensions must match N x N.", nameof(flatColumnMajorMatrix));

#if EIGEN_REFERENCE

        var matrix = new double[n, n];
        for (int c = 0; c < n; c++)
        {
            for (int r = 0; r < n; r++)
                matrix[r, c] = flatColumnMajorMatrix[c * n + r];
        }

        return Eigen.DecomposeSymmetric(matrix, maxSweeps, tol);
#else
        return EigenFast.DecomposeSymmetric(flatColumnMajorMatrix, n, maxSweeps, tol);
#endif
    }

    /// <summary>
    /// Decomposes the symmetric n × n matrix held flat column-major into caller-owned buffers, with no allocation once the pools are
    /// warm: eigenvalues descending into <paramref name="eigenvalues"/> (length n), and eigenvectors as the columns of the flat
    /// column-major <paramref name="eigenvectors"/> (length n²), or none when it is empty, which skips their rotations and leaves the
    /// eigenvalues the same bits (spectral-solver brief, B2).
    /// </summary>
    public static DenseEigenConvergence DecomposeSymmetric(
        ReadOnlySpan<double> flatColumnMajorMatrix,
        int n,
        Span<double> eigenvalues,
        Span<double> eigenvectors,
        int maxSweeps = 256,
        double tol = 1e-12)
    {
        if (flatColumnMajorMatrix.Length != n * n)
            throw new ArgumentException("Matrix dimensions must match N x N.", nameof(flatColumnMajorMatrix));
#if EIGEN_REFERENCE
        if (eigenvalues.Length != n) throw new ArgumentException("Eigenvalues buffer size must equal N.", nameof(eigenvalues));
        if (eigenvectors.Length != 0 && eigenvectors.Length != n * n)
            throw new ArgumentException("Eigenvectors buffer size must equal N * N, or 0 for eigenvalues alone.", nameof(eigenvectors));
        var matrix = new double[n, n];
        for (int c = 0; c < n; c++)
            for (int r = 0; r < n; r++) matrix[r, c] = flatColumnMajorMatrix[c * n + r];
        EigenResult result = Eigen.DecomposeSymmetric(matrix, maxSweeps, tol);
        result.Eigenvalues.CopyTo(eigenvalues);
        if (eigenvectors.Length != 0)
            for (int j = 0; j < n; j++) result.Eigenvectors[j].CopyTo(eigenvectors.Slice(j * n, n));
        return result.Convergence;
#else
        return EigenFast.DecomposeSymmetricInPlace(flatColumnMajorMatrix, n, eigenvalues, eigenvectors, maxSweeps, tol);
#endif
    }

    /// <summary>
    /// Decomposes the symmetric <paramref name="matrix"/> into caller-owned outputs: eigenvalues descending into
    /// <paramref name="eigenvalues"/>, and eigenvectors as the columns of <paramref name="eigenvectorColumns"/> (n × n), the layout
    /// callers read, or none when it is null. The input is copied to a pooled column-major buffer and the eigenvectors written once,
    /// with no rows allocated in between (spectral-solver brief, B2).
    /// </summary>
    public static DenseEigenConvergence DecomposeSymmetricInto(
        double[,] matrix,
        Span<double> eigenvalues,
        double[,]? eigenvectorColumns,
        int maxSweeps = 256,
        double tol = 1e-12)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        int n = matrix.GetLength(0);
        if (matrix.GetLength(1) != n) throw new ArgumentException("Matrix must be square.", nameof(matrix));
        double[] flat = ArrayPool<double>.Shared.Rent(n * n);
        try
        {
            for (int c = 0; c < n; c++)
                for (int r = 0; r < n; r++) flat[c * n + r] = matrix[r, c];
            return DecomposeSymmetricInto(flat.AsSpan(0, n * n), n, eigenvalues, eigenvectorColumns, maxSweeps, tol);
        }
        finally
        {
            ArrayPool<double>.Shared.Return(flat);
        }
    }

    /// <summary>
    /// The flat column-major matrix into caller-owned outputs, eigenvectors written as the columns of
    /// <paramref name="eigenvectorColumns"/> (n × n), or none when it is null: a pooled column-major buffer receives them and is
    /// transposed into place once.
    /// </summary>
    public static DenseEigenConvergence DecomposeSymmetricInto(
        ReadOnlySpan<double> flatColumnMajorMatrix,
        int n,
        Span<double> eigenvalues,
        double[,]? eigenvectorColumns,
        int maxSweeps = 256,
        double tol = 1e-12)
    {
        if (eigenvectorColumns is null) return DecomposeSymmetric(flatColumnMajorMatrix, n, eigenvalues, Span<double>.Empty, maxSweeps, tol);
        if (eigenvectorColumns.GetLength(0) != n || eigenvectorColumns.GetLength(1) != n)
            throw new ArgumentException("The eigenvector columns are n × n.", nameof(eigenvectorColumns));
        double[] vectors = ArrayPool<double>.Shared.Rent(n * n);
        try
        {
            DenseEigenConvergence convergence = DecomposeSymmetric(flatColumnMajorMatrix, n, eigenvalues, vectors.AsSpan(0, n * n), maxSweeps, tol);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < n; i++) eigenvectorColumns[i, j] = vectors[j * n + i];
            return convergence;
        }
        finally
        {
            ArrayPool<double>.Shared.Return(vectors);
        }
    }
}
