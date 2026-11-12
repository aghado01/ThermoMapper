using System;
using System.Collections.Generic;
using System.Linq;
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

    /// <summary>
    /// The reported residual and converged flag describe the pairs returned, at the iteration cap and on convergence (spectral-solver
    /// brief, A1). The witness is the tridiagonal matrix with diagonal [4, 3, 2] and unit couplings, its largest pair from seed 173:
    /// capped at one iteration, the block was updated after its residual was measured, and 0.76200174 was reported for a pair whose
    /// residual is 0.34933228. The operator is applied to the returned pairs independently here.
    /// </summary>
    [Fact]
    public void Solve_ResidualDescribesReturnedPairs()
    {
        var tridiagonal = new double[,] { { 4.0, 1.0, 0.0 }, { 1.0, 3.0, 1.0 }, { 0.0, 1.0, 2.0 } };
        var largest = new LOBPCG.Options { WantLargest = true, Seed = 173, Tolerance = 1e-12 };
        LOBPCG.Result capped = LOBPCG.Solve(new DenseSymmetricOperator(tridiagonal), 1, largest with { MaxIterations = 1 });
        LOBPCG.Result converged = LOBPCG.Solve(new DenseSymmetricOperator(tridiagonal), 1, largest);
        double cappedResidual = MaxResidual(tridiagonal, capped);
        double convergedResidual = MaxResidual(tridiagonal, converged);

        Assert.Equal(1, capped.Iterations);
        Assert.InRange(Math.Abs(capped.ResidualNorm - cappedResidual), 0.0, 1e-12);
        Assert.InRange(cappedResidual, 0.3493, 0.3494);
        Assert.False(capped.Converged);
        Assert.True(converged.Converged);
        Assert.InRange(converged.ResidualNorm, 0.0, 1e-12);
        Assert.InRange(convergedResidual, 0.0, 1e-12);

        // A block of four capped at three iterations: the largest residual over the returned pairs.
        double[,] matrix = BuildRandomSymmetricMatrix(size: 64, seed: 7);
        LOBPCG.Result block = LOBPCG.Solve(new DenseSymmetricOperator(matrix), 4, new LOBPCG.Options { MaxIterations = 3 });
        Assert.InRange(Math.Abs(block.ResidualNorm - MaxResidual(matrix, block)), 0.0, 1e-12);
        Assert.False(block.Converged);
    }

    private static double MaxResidual(double[,] matrix, LOBPCG.Result result) =>
        result.Eigenpairs.Max(pair => ResidualNorm(matrix, pair.Lambda, pair.Vector));

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
