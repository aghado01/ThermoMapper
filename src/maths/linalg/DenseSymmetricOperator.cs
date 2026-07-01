// src/maths/linalg/DenseSymmetricOperator.cs
#nullable enable
using System;

namespace Maths.LinAlg;

/// <summary>
/// Wraps a materialised symmetric matrix as an <see cref="ILinearOperator"/> so a
/// dense <c>double[,]</c> (or flat column-major span) can be fed to
/// <see cref="LOBPCG"/>. This is the dense peer of the matrix-free graph operator
/// in <c>Graphs.Spectral</c>: it lets LOBPCG stand alongside
/// <see cref="Eigen"/>/<see cref="EigenFast"/> on the same input for
/// apples-to-apples comparison.
/// </summary>
/// <remarks>
/// The matrix is assumed symmetric, so the column-major store is indexed as if
/// row-major in the matvec (<c>A[i,j] == A[j,i]</c>) to keep the inner loop
/// contiguous. For small dense problems a full Jacobi decomposition
/// (<see cref="Eigen"/>/<see cref="EigenFast"/>) is usually faster and returns the
/// whole spectrum; LOBPCG over a dense operator earns its keep when <c>n</c> is
/// large and only a few extremal pairs are wanted — and as the common surface for
/// benchmarking the solvers against one another.
/// </remarks>
public sealed class DenseSymmetricOperator : ILinearOperator
{
    private readonly double[] _matrix; // flat column-major, length n*n
    private readonly int _n;

    /// <summary>
    /// Wraps an idiomatic <c>double[,]</c> symmetric matrix. The values are copied
    /// into an internal flat buffer (one-time O(n²), negligible against the solve).
    /// </summary>
    public DenseSymmetricOperator(double[,] matrix)
    {
        if (matrix is null) throw new ArgumentNullException(nameof(matrix));
        _n = matrix.GetLength(0);
        if (matrix.GetLength(1) != _n)
            throw new ArgumentException("Matrix must be square.", nameof(matrix));

        _matrix = new double[_n * _n];
        for (int j = 0; j < _n; j++)
            for (int i = 0; i < _n; i++)
                _matrix[i + j * _n] = matrix[i, j];
    }

    /// <summary>
    /// Wraps a flat column-major symmetric matrix of dimension <paramref name="n"/>.
    /// The array is used by reference (not copied); do not mutate it for the
    /// lifetime of the operator.
    /// </summary>
    public DenseSymmetricOperator(double[] flatColumnMajor, int n)
    {
        if (flatColumnMajor is null) throw new ArgumentNullException(nameof(flatColumnMajor));
        if (n < 0 || flatColumnMajor.Length != n * n)
            throw new ArgumentException("Array length must equal n × n.", nameof(flatColumnMajor));

        _matrix = flatColumnMajor;
        _n = n;
    }

    public int Dimension => _n;

    public void Apply(ReadOnlySpan<double> block, Span<double> result, int columns)
    {
        int n = _n;
        for (int c = 0; c < columns; c++)
        {
            int colOffset = c * n;
            for (int i = 0; i < n; i++)
            {
                // (Ax)_i = Σ_j A[i,j] x_j. A is symmetric, so the column-major store
                // A[j + i*n] = A[i,j]; row i is then contiguous (stride 1) in j.
                int rowOffset = i * n;
                double sum = 0.0;
                for (int j = 0; j < n; j++)
                    sum += _matrix[rowOffset + j] * block[colOffset + j];
                result[colOffset + i] = sum;
            }
        }
    }
}
