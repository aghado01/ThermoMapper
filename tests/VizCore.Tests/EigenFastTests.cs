using System;
using Graphs.Primitives;
using Maths.LinAlg;
using Xunit;

namespace VizCore.Tests;

public sealed class EigenFastTests
{
    [Fact]
    public void DecomposeSymmetric_DefaultParameters_MatchReference()
    {
        double[,] matrix = new double[,]
        {
            { 2.0, 1.0, 0.0 },
            { 1.0, 2.0, 0.0 },
            { 0.0, 0.0, 5.0 }
        };

        AssertMatchesReference(matrix, eigenvalueTolerance: 1e-9, residualTolerance: 1e-9);
    }

    [Fact]
    public void DecomposeSymmetric_ExceptionContract_MatchesReference()
    {
        double[,] matrix = new double[2, 3];

        Assert.Throws<ArgumentException>(() => Eigen.DecomposeSymmetric(matrix));
        Assert.Throws<ArgumentException>(() => EigenFast.DecomposeSymmetric(matrix));
    }

    [Fact]
    public void DecomposeSymmetric_ThetaZeroCase_MatchesReference()
    {
        double[,] matrix = new double[,]
        {
            { 2.0, 1.0, 0.0 },
            { 1.0, 2.0, 0.0 },
            { 0.0, 0.0, 5.0 }
        };

        AssertMatchesReference(matrix, eigenvalueTolerance: 1e-9, residualTolerance: 1e-9);
    }

    [Fact]
    public void DecomposeSymmetric_KnownSpectrumFixtures_MatchReference()
    {
        AssertMatchesReference(
            BuildSymmetricMatrix(
                new[] { 5.0, 3.0, 1.0 },
                BuildPlanarRotationBasis3(Math.PI / 6.0)),
            eigenvalueTolerance: 1e-9,
            residualTolerance: 1e-9);

        AssertMatchesReference(
            BuildSymmetricMatrix(
                new[] { 7.0, 4.0, 2.0, -1.0 },
                BuildNormalizedHadamard(4)),
            eigenvalueTolerance: 1e-9,
            residualTolerance: 1e-9);

        AssertMatchesReference(
            BuildSymmetricMatrix(
                new[] { 12.0, 9.0, 7.0, 5.0, 4.0, 3.0, 2.0, 0.5 },
                BuildNormalizedHadamard(8)),
            eigenvalueTolerance: 1e-9,
            residualTolerance: 1e-9);
    }

    [Fact]
    public void DecomposeSymmetric_DisconnectedCombinatorialLaplacian_MatchesReference()
    {
        DisconnectedControl.Fixture control = DisconnectedControl.Generate();
        CsrGraph graph = FixtureGoldenHelpers.BuildDraftGraph(control.Points, k: 4, ensureConnected: false);
        double[,] laplacian = BuildCombinatorialLaplacian(graph);

        EigenResult reference = Eigen.DecomposeSymmetric(laplacian);
        EigenResult fast = EigenFast.DecomposeSymmetric(laplacian);

        AssertMatches(reference, fast, laplacian, eigenvalueTolerance: 1e-9, residualTolerance: 1e-9);
        Assert.InRange(Math.Abs(fast.Eigenvalues[^1]), 0.0, 1e-10);
        Assert.InRange(Math.Abs(fast.Eigenvalues[^2]), 0.0, 1e-10);
    }

    [Fact]
    public void DecomposeSymmetric_DeterministicRandomSymmetric_MatchesReference()
    {
        double[,] matrix = BuildRandomSymmetricMatrix(size: 50, seed: 1729);

        AssertMatchesReference(matrix, eigenvalueTolerance: 1e-8, residualTolerance: 1e-8);
    }

    [Fact]
    public void DecomposeSymmetric_NearDegenerateSpectrum_MatchesReference()
    {
        double[,] matrix = BuildSymmetricMatrix(
            new[] { 6.0, 6.0, 2.0, 1.0 },
            BuildNormalizedHadamard(4));

        AssertMatchesReference(matrix, eigenvalueTolerance: 1e-8, residualTolerance: 1e-8);
    }

    private static void AssertMatchesReference(double[,] matrix, double eigenvalueTolerance, double residualTolerance)
    {
        EigenResult reference = Eigen.DecomposeSymmetric(matrix);
        EigenResult fast = EigenFast.DecomposeSymmetric(matrix);

        AssertMatches(reference, fast, matrix, eigenvalueTolerance, residualTolerance);
    }

    private static void AssertMatches(
        EigenResult reference,
        EigenResult fast,
        double[,] matrix,
        double eigenvalueTolerance,
        double residualTolerance)
    {
        Assert.Equal(reference.Eigenvalues.Length, fast.Eigenvalues.Length);
        Assert.Equal(reference.Eigenvectors.Length, fast.Eigenvectors.Length);

        AssertDescending(reference.Eigenvalues);
        AssertDescending(fast.Eigenvalues);

        for (int i = 0; i < reference.Eigenvalues.Length; i++)
        {
            Assert.InRange(Math.Abs(reference.Eigenvalues[i] - fast.Eigenvalues[i]), 0.0, eigenvalueTolerance);
            Assert.InRange(ComputeResidualNorm(matrix, reference.Eigenvalues[i], reference.Eigenvectors[i]), 0.0, residualTolerance);
            Assert.InRange(ComputeResidualNorm(matrix, fast.Eigenvalues[i], fast.Eigenvectors[i]), 0.0, residualTolerance);
        }
    }

    private static void AssertDescending(double[] eigenvalues)
    {
        for (int i = 1; i < eigenvalues.Length; i++)
            Assert.True(eigenvalues[i - 1] >= eigenvalues[i], "Eigenvalues must be sorted in descending order.");
    }

    private static double ComputeResidualNorm(double[,] matrix, double eigenvalue, double[] eigenvector)
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

    private static double[,] BuildSymmetricMatrix(double[] eigenvalues, double[,] basis)
    {
        int n = eigenvalues.Length;
        var matrix = new double[n, n];

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                double sum = 0.0;
                for (int k = 0; k < n; k++)
                    sum += basis[i, k] * eigenvalues[k] * basis[j, k];

                matrix[i, j] = sum;
            }
        }

        return matrix;
    }

    private static double[,] BuildPlanarRotationBasis3(double theta)
    {
        double c = Math.Cos(theta);
        double s = Math.Sin(theta);

        return new[,]
        {
            { c, -s, 0.0 },
            { s,  c, 0.0 },
            { 0.0, 0.0, 1.0 }
        };
    }

    private static double[,] BuildNormalizedHadamard(int size)
    {
        if (size <= 0 || (size & (size - 1)) != 0)
            throw new ArgumentOutOfRangeException(nameof(size), "Hadamard size must be a positive power of two.");

        int[,] hadamard = BuildHadamard(size);
        var normalized = new double[size, size];
        double scale = 1.0 / Math.Sqrt(size);

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
                normalized[i, j] = hadamard[i, j] * scale;
        }

        return normalized;
    }

    private static int[,] BuildHadamard(int size)
    {
        if (size == 1)
            return new[,] { { 1 } };

        int half = size / 2;
        int[,] smaller = BuildHadamard(half);
        var hadamard = new int[size, size];

        for (int i = 0; i < half; i++)
        {
            for (int j = 0; j < half; j++)
            {
                int value = smaller[i, j];
                hadamard[i, j] = value;
                hadamard[i, j + half] = value;
                hadamard[i + half, j] = value;
                hadamard[i + half, j + half] = -value;
            }
        }

        return hadamard;
    }

    private static double[,] BuildCombinatorialLaplacian(CsrGraph graph)
    {
        int n = graph.NodeCount;
        var laplacian = new double[n, n];

        for (int i = 0; i < n; i++)
        {
            int start = graph.RowPointers[i];
            int end = graph.RowPointers[i + 1];
            double degree = 0.0;

            for (int edge = start; edge < end; edge++)
            {
                degree += graph.Weights[edge];
                laplacian[i, graph.Targets[edge]] -= graph.Weights[edge];
            }

            laplacian[i, i] = degree;
        }

        return laplacian;
    }
}
