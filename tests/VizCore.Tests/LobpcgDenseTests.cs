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

    /// <summary>
    /// The Rayleigh–Ritz step runs on the production dense kernel (spectral-solver brief, B3), and after the first iteration only on
    /// the trial subspace, the operator images carried (B5): a solve that converges at iteration I decomposes I projections, the first
    /// iteration's and one trial projection per iteration after it, none of them by the reference Jacobi; a solve capped at three
    /// decomposes four. The build's dense dispatch decides the kernel, the reference one only under EIGEN_REFERENCE.
    /// </summary>
    [Fact]
    public void Solve_RayleighRitzRunsOnTheDenseDispatch()
    {
#if EIGEN_REFERENCE
        const DecompositionKernel dense = DecompositionKernel.ReferenceJacobi;
#else
        const DecompositionKernel dense = DecompositionKernel.FastJacobi;
#endif
        double[,] matrix = BuildRandomSymmetricMatrix(size: 64, seed: 7);
        var op = new DenseSymmetricOperator(matrix);
        LOBPCG.Result converged, capped;
        IReadOnlyList<Decomposition> convergedWork, cappedWork;
        using (DecompositionLog log = DecompositionLog.Begin())
        {
            converged = LOBPCG.Solve(op, 4, new LOBPCG.Options { MaxIterations = 2000, Tolerance = 1e-9 });
            convergedWork = log.Entries;
        }
        using (DecompositionLog log = DecompositionLog.Begin())
        {
            capped = LOBPCG.Solve(op, 4, new LOBPCG.Options { MaxIterations = 3 });
            cappedWork = log.Entries;
        }

        Assert.True(converged.Converged);
        Assert.All(convergedWork.Concat(cappedWork), d => Assert.Equal(dense, d.Kernel));
        Assert.Equal(converged.Iterations, convergedWork.Count);
        Assert.False(capped.Converged);
        Assert.Equal(4, cappedWork.Count);
        Assert.InRange(Math.Abs(converged.ResidualNorm - MaxResidual(matrix, converged)), 0.0, 1e-12);
    }

    /// <summary>
    /// With the operator images carried across iterations (spectral-solver brief, B5), the residual reported is still the returned
    /// pairs': a solve below the tolerance on its carried images is confirmed on a fresh application. Checked on a constrained and
    /// preconditioned run, where the carried images pass through the deflation and the preconditioner. The constraint is the lowest
    /// eigenvector, an invariant subspace as a Laplacian's null vector is, so the solve finds the next five eigenpairs: they match the
    /// dense spectrum, stay orthogonal to the constraint, and the reported residual matches the operator applied to them independently.
    /// The constraint is itself an approximate eigenvector, which leaves a residual floor near 10⁻⁸ that fresh and carried images reach
    /// alike (7.0·10⁻⁹ carried, 8.7·10⁻⁹ fresh, after 2000 iterations at a tolerance of 10⁻¹⁰), so the tolerance here is 10⁻⁷; and the
    /// pairs lean on it alike, by 3.8·10⁻⁹ carried and 4.6·10⁻⁹ fresh, so orthogonality is read at 10⁻⁸.
    /// </summary>
    [Fact]
    public void Solve_CarriedImages_ReportTheReturnedPairsResidual()
    {
        const int n = 80;
        double[,] matrix = BuildRandomSymmetricMatrix(size: n, seed: 23);
        for (int i = 0; i < n; i++) matrix[i, i] += 4.0 + i * 0.05;   // a varied diagonal for the preconditioner to read
        EigenResult full = DenseEigen.DecomposeSymmetric(matrix);
        double[] constraint = full.Eigenvectors[n - 1];
        var options = new LOBPCG.Options
        {
            MaxIterations = 2000, Tolerance = 1e-7, Constraints = constraint, ConstraintColumns = 1,
            Preconditioner = new DiagonalPreconditioner(matrix),
        };
        LOBPCG.Result result = LOBPCG.Solve(new DenseSymmetricOperator(matrix), 5, options);

        Assert.True(result.Converged, $"iterations {result.Iterations}, residual {result.ResidualNorm:E3}, independent {MaxResidual(matrix, result):E3}");
        Assert.InRange(result.ResidualNorm, 0.0, 1e-7);
        Assert.InRange(Math.Abs(result.ResidualNorm - MaxResidual(matrix, result)), 0.0, 1e-12);
        for (int j = 0; j < 5; j++)
        {
            Assert.InRange(Math.Abs(result.Eigenpairs[j].Lambda - full.Eigenvalues[n - 2 - j]), 0.0, 1e-9);
            Assert.InRange(Math.Abs(result.Eigenpairs[j].Vector.Zip(constraint, (a, b) => a * b).Sum()), 0.0, 1e-8);
        }
    }

    // D⁻¹ for the matrix's diagonal D.
    private sealed class DiagonalPreconditioner(double[,] matrix) : IPreconditioner
    {
        public void Apply(ReadOnlySpan<double> residualBlock, Span<double> resultBlock, int columns)
        {
            int n = matrix.GetLength(0);
            for (int c = 0; c < columns; c++)
                for (int i = 0; i < n; i++) resultBlock[c * n + i] = residualBlock[c * n + i] / matrix[i, i];
        }
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
