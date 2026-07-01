using System;
using System.Collections.Generic;
using Maths.LinAlg;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// Parity guards for LOBPCG's dense entry points (<see cref="LOBPCG.BottomK(double[,], int, LOBPCG.Options)"/> /
/// <see cref="LOBPCG.TopK(double[,], int, LOBPCG.Options)"/>): the iterative solver must
/// recover the same extremal spectrum as the direct dense decomposition on the
/// same matrix, so the three solvers can be reasoned about — and benchmarked —
/// apples-to-apples.
/// </summary>
public sealed class LobpcgDenseTests
{
    [Fact]
    public void BottomK_MatchesDenseReference()
    {
        double[,] matrix = BuildRandomSymmetricMatrix(size: 64, seed: 7);
        const int k = 4;
        var options = new LOBPCG.Options { MaxIterations = 2000, Tolerance = 1e-11 };

        IReadOnlyList<EigenPair> direct = SpectralMath.BottomK(matrix, k);
        IReadOnlyList<EigenPair> iterative = LOBPCG.BottomK(matrix, k, options);

        Assert.Equal(direct.Count, iterative.Count);
        for (int i = 0; i < k; i++)
        {
            Assert.InRange(Math.Abs(direct[i].Lambda - iterative[i].Lambda), 0.0, 1e-4);
            Assert.InRange(ResidualNorm(matrix, iterative[i].Lambda, iterative[i].Vector), 0.0, 1e-4);
        }
    }

    [Fact]
    public void TopK_MatchesDenseReference()
    {
        double[,] matrix = BuildRandomSymmetricMatrix(size: 64, seed: 11);
        const int k = 4;
        var options = new LOBPCG.Options { MaxIterations = 2000, Tolerance = 1e-11 };

        // Reference: full decomposition is sorted descending, so the head is the top-k.
        EigenResult full = Eigen.DecomposeSymmetric(matrix);
        IReadOnlyList<EigenPair> iterative = LOBPCG.TopK(matrix, k, options);

        Assert.Equal(k, iterative.Count);
        for (int i = 0; i < k; i++)
        {
            Assert.InRange(Math.Abs(full.Eigenvalues[i] - iterative[i].Lambda), 0.0, 1e-4);
            Assert.InRange(ResidualNorm(matrix, iterative[i].Lambda, iterative[i].Vector), 0.0, 1e-4);
        }
    }

    private static double[,] BuildRandomSymmetricMatrix(int size, int seed)
    {
        var random = new Random(seed);
        var matrix = new double[size, size];
        for (int i = 0; i < size; i++)
        {
            for (int j = i; j < size; j++)
            {
                double value = random.NextDouble() * 2.0 - 1.0;
                matrix[i, j] = value;
                matrix[j, i] = value;
            }
        }

        return matrix;
    }

    private static double ResidualNorm(double[,] matrix, double eigenvalue, double[] eigenvector)
    {
        int n = eigenvector.Length;
        double sumSquares = 0.0;
        for (int row = 0; row < n; row++)
        {
            double projected = 0.0;
            for (int col = 0; col < n; col++)
                projected += matrix[row, col] * eigenvector[col];

            double residual = projected - eigenvalue * eigenvector[row];
            sumSquares += residual * residual;
        }

        return Math.Sqrt(sumSquares);
    }
}
