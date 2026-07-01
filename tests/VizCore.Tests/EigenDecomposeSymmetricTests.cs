using System;
using Graphs.Primitives;
using Maths.LinAlg;
using Xunit;

namespace VizCore.Tests;

public sealed class EigenDecomposeSymmetricTests
{
    [Fact]
    public void DecomposeSymmetric_MatchesKnownEigenvalues()
    {
        AssertSpectrumMatches(
            BuildSymmetricMatrix(
                new[] { 5.0, 3.0, 1.0 },
                BuildPlanarRotationBasis3(Math.PI / 6.0)),
            new[] { 5.0, 3.0, 1.0 });

        AssertSpectrumMatches(
            BuildSymmetricMatrix(
                new[] { 7.0, 4.0, 2.0, -1.0 },
                BuildNormalizedHadamard(4)),
            new[] { 7.0, 4.0, 2.0, -1.0 });

        AssertSpectrumMatches(
            BuildSymmetricMatrix(
                new[] { 12.0, 9.0, 7.0, 5.0, 4.0, 3.0, 2.0, 0.5 },
                BuildNormalizedHadamard(8)),
            new[] { 12.0, 9.0, 7.0, 5.0, 4.0, 3.0, 2.0, 0.5 });
    }

    [Fact]
    public void DecomposeSymmetric_DisconnectedCombinatorialLaplacian_HasTwoZeroModes()
    {
        DisconnectedControl.Fixture control = DisconnectedControl.Generate();
        CsrGraph graph = FixtureGoldenHelpers.BuildDraftGraph(control.Points, k: 4, ensureConnected: false);
        double[,] laplacian = BuildCombinatorialLaplacian(graph);

        EigenResult result = Eigen.DecomposeSymmetric(laplacian);

        Assert.InRange(Math.Abs(result.Eigenvalues[^1]), 0.0, 1e-9);
        Assert.InRange(Math.Abs(result.Eigenvalues[^2]), 0.0, 1e-9);
    }

    private static void AssertSpectrumMatches(double[,] matrix, double[] expectedEigenvalues)
    {
        EigenResult result = Eigen.DecomposeSymmetric(matrix);

        Assert.Equal(expectedEigenvalues.Length, result.Eigenvalues.Length);
        for (int i = 0; i < expectedEigenvalues.Length; i++)
            Assert.InRange(Math.Abs(result.Eigenvalues[i] - expectedEigenvalues[i]), 0.0, 1e-9);
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
