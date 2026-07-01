using System;
using System.Linq;
using Maths.LinAlg;
using Xunit;

namespace VizCore.Tests;

public sealed class PcaTests
{
    [Fact]
    public void Compute_DeterministicData_ProducesOrderedComponentsAndWhitenedProjection()
    {
        double[][] data = BuildDeterministicData(sampleCount: 256, seed: 1234);

        PcaResult pca = Pca.Compute(data, numComponents: 3, center: true, whiten: true);

        Assert.Equal(3, pca.Components.Length);
        Assert.All(pca.Components, component => Assert.Equal(3, component.Length));

        for (int i = 1; i < pca.ExplainedVarianceRatio.Length; i++)
            Assert.True(pca.ExplainedVarianceRatio[i - 1] >= pca.ExplainedVarianceRatio[i]);

        double[][] whitened = data.Select(Pca.MakeProjector(pca)).ToArray();
        double[,] covariance = ComputeCovariance(whitened);

        for (int i = 0; i < 3; i++)
        {
            Assert.InRange(Math.Abs(covariance[i, i] - 1.0), 0.0, 5e-3);
            for (int j = 0; j < 3; j++)
            {
                if (i == j)
                    continue;

                Assert.InRange(Math.Abs(covariance[i, j]), 0.0, 5e-3);
            }
        }
    }

    private static double[][] BuildDeterministicData(int sampleCount, int seed)
    {
        var random = new Random(seed);
        var data = new double[sampleCount][];

        for (int i = 0; i < sampleCount; i++)
        {
            double z1 = NextGaussian(random);
            double z2 = NextGaussian(random);
            double z3 = NextGaussian(random);

            data[i] =
            [
                2.0 * z1 + 0.5 * z2,
                -1.0 * z1 + 1.5 * z2 + 0.25 * z3,
                0.75 * z1 - 0.5 * z2 + 1.25 * z3
            ];
        }

        return data;
    }

    private static double NextGaussian(Random random)
    {
        double u1 = 1.0 - random.NextDouble();
        double u2 = 1.0 - random.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double[,] ComputeCovariance(double[][] data)
    {
        int n = data.Length;
        int d = data[0].Length;
        double[] mean = new double[d];
        var covariance = new double[d, d];

        for (int i = 0; i < n; i++)
            for (int j = 0; j < d; j++)
                mean[j] += data[i][j];

        for (int j = 0; j < d; j++)
            mean[j] /= n;

        for (int i = 0; i < d; i++)
        {
            for (int j = i; j < d; j++)
            {
                double sum = 0.0;
                for (int k = 0; k < n; k++)
                    sum += (data[k][i] - mean[i]) * (data[k][j] - mean[j]);

                covariance[i, j] = sum / n;
                covariance[j, i] = covariance[i, j];
            }
        }

        return covariance;
    }
}
