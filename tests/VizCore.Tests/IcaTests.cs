using System;
using System.Linq;
using Maths.LinAlg;
using Xunit;

namespace VizCore.Tests;

public sealed class IcaTests
{
    [Fact]
    public void Compute_SymmetricMode_ProducesOrthogonalUnmixingAndRecoversSourcesUpToSignedPermutation()
    {
        (double[][] mixed, double[][] sources) = BuildMixedSignals(sampleCount: 2048);

        IcaResult result = FastIca.Compute(
            mixed,
            numComponents: 2,
            mode: FastIca.Mode.Symmetric,
            maxIter: 400,
            tol: 1e-7,
            seed: 7);

        double[,] gram = Multiply(result.UnmixingMatrix, Transpose(result.UnmixingMatrix));
        for (int i = 0; i < 2; i++)
        {
            Assert.InRange(Math.Abs(gram[i, i] - 1.0), 0.0, 5e-3);
            for (int j = 0; j < 2; j++)
            {
                if (i == j)
                    continue;

                Assert.InRange(Math.Abs(gram[i, j]), 0.0, 5e-3);
            }
        }

        double[][] recovered = RecoverSources(mixed, result);
        double[,] correlations = AbsoluteCorrelationMatrix(recovered, sources);

        Assert.True(Math.Max(correlations[0, 0], correlations[1, 0]) >= 0.95);
        Assert.True(Math.Max(correlations[0, 1], correlations[1, 1]) >= 0.95);
    }

    private static (double[][] Mixed, double[][] Sources) BuildMixedSignals(int sampleCount)
    {
        var sources = new double[sampleCount][];
        var mixed = new double[sampleCount][];

        for (int i = 0; i < sampleCount; i++)
        {
            double t = 2.0 * Math.PI * i / sampleCount;
            double s1 = Math.Tanh(2.5 * Math.Sin(t)) + 0.15 * Math.Sin(3.0 * t);
            double s2 = ((i / 32) % 2 == 0 ? 1.0 : -1.0) + 0.2 * Math.Sin(5.0 * t);

            sources[i] = [s1, s2];
            mixed[i] =
            [
                1.0 * s1 + 0.5 * s2,
                0.3 * s1 + 1.2 * s2
            ];
        }

        return (mixed, sources);
    }

    private static double[][] RecoverSources(double[][] mixed, IcaResult result)
    {
        Func<double[], double[]> projector = Pca.MakeProjector(result.WhiteningPca);
        double[][] whitened = mixed.Select(projector).ToArray();
        var recovered = new double[whitened.Length][];

        for (int i = 0; i < whitened.Length; i++)
        {
            recovered[i] = new double[result.UnmixingMatrix.Length];
            for (int j = 0; j < result.UnmixingMatrix.Length; j++)
                recovered[i][j] = Dot(whitened[i], result.UnmixingMatrix[j]);
        }

        return recovered;
    }

    private static double[,] AbsoluteCorrelationMatrix(double[][] left, double[][] right)
    {
        int components = left[0].Length;
        var matrix = new double[components, components];

        for (int i = 0; i < components; i++)
        {
            double[] leftColumn = left.Select(row => row[i]).ToArray();
            for (int j = 0; j < components; j++)
            {
                double[] rightColumn = right.Select(row => row[j]).ToArray();
                matrix[i, j] = Math.Abs(Correlation(leftColumn, rightColumn));
            }
        }

        return matrix;
    }

    private static double Correlation(double[] left, double[] right)
    {
        double leftMean = left.Average();
        double rightMean = right.Average();
        double numerator = 0.0;
        double leftNorm = 0.0;
        double rightNorm = 0.0;

        for (int i = 0; i < left.Length; i++)
        {
            double dx = left[i] - leftMean;
            double dy = right[i] - rightMean;
            numerator += dx * dy;
            leftNorm += dx * dx;
            rightNorm += dy * dy;
        }

        return numerator / Math.Sqrt(leftNorm * rightNorm);
    }

    private static double[,] Multiply(double[][] left, double[][] right)
    {
        int rows = left.Length;
        int inner = right.Length;
        int cols = right[0].Length;
        var product = new double[rows, cols];

        for (int i = 0; i < rows; i++)
        {
            for (int k = 0; k < inner; k++)
            {
                for (int j = 0; j < cols; j++)
                    product[i, j] += left[i][k] * right[k][j];
            }
        }

        return product;
    }

    private static double[][] Transpose(double[][] matrix)
    {
        int rows = matrix.Length;
        int cols = matrix[0].Length;
        var transpose = new double[cols][];

        for (int i = 0; i < cols; i++)
        {
            transpose[i] = new double[rows];
            for (int j = 0; j < rows; j++)
                transpose[i][j] = matrix[j][i];
        }

        return transpose;
    }

    private static double Dot(double[] left, double[] right)
    {
        double sum = 0.0;
        for (int i = 0; i < left.Length; i++)
            sum += left[i] * right[i];

        return sum;
    }
}
