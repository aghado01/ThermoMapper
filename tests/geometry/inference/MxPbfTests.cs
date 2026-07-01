using System;
using Maths.Geometry.Inference;
using Xunit;

namespace Maths.Geometry.Tests;

public sealed class MxPbfTests
{
    [Fact]
    public void TwoSampleMean_SeparatedMeans_GiveStrongerEvidenceThanIdentical()
    {
        double[][] x = Gaussian(rows: 30, cols: 5, seed: 1, shift: 0.0);
        double[][] ySame = Gaussian(rows: 30, cols: 5, seed: 2, shift: 0.0);
        double[][] ySep = Gaussian(rows: 30, cols: 5, seed: 2, shift: 5.0);

        double bfSame = MxPbf.TwoSampleMean(x, ySame);
        double bfSep = MxPbf.TwoSampleMean(x, ySep);

        Assert.True(double.IsFinite(bfSame) && double.IsFinite(bfSep));
        Assert.True(
            bfSep > bfSame,
            $"Separated means should raise the log mxPBF: sep={bfSep:F2} vs same={bfSame:F2}.");
    }

    [Fact]
    public void TwoSampleCovariance_DifferentScale_GivesStrongerEvidenceThanIdentical()
    {
        double[][] x = Gaussian(rows: 40, cols: 4, seed: 3, shift: 0.0);
        double[][] ySame = Gaussian(rows: 40, cols: 4, seed: 4, shift: 0.0);
        double[][] yScaled = Scale(Gaussian(rows: 40, cols: 4, seed: 4, shift: 0.0), 3.0);

        double bfSame = MxPbf.TwoSampleCovariance(x, ySame);
        double bfDiff = MxPbf.TwoSampleCovariance(x, yScaled);

        Assert.True(
            bfDiff > bfSame,
            $"A covariance-scale gap should raise the log mxPBF: diff={bfDiff:F2} vs same={bfSame:F2}.");
    }

    [Fact]
    public void TwoSampleMean_DimensionMismatch_Throws()
    {
        double[][] x = Gaussian(rows: 10, cols: 5, seed: 1, shift: 0.0);
        double[][] y = Gaussian(rows: 10, cols: 4, seed: 2, shift: 0.0);
        Assert.Throws<ArgumentException>(() => MxPbf.TwoSampleMean(x, y));
    }

    private static double[][] Gaussian(int rows, int cols, int seed, double shift)
    {
        var rng = new Random(seed);
        double[][] data = new double[rows][];
        for (int i = 0; i < rows; i++)
        {
            data[i] = new double[cols];
            for (int j = 0; j < cols; j++)
                data[i][j] = NextGaussian(rng) + shift;
        }
        return data;
    }

    private static double[][] Scale(double[][] data, double factor)
    {
        for (int i = 0; i < data.Length; i++)
            for (int j = 0; j < data[i].Length; j++)
                data[i][j] *= factor;
        return data;
    }

    private static double NextGaussian(Random rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }
}
