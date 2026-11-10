using System;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Regression.Tests;

public sealed class MultiOutputPredictivePrimitiveTests
{
    [Theory]
    [InlineData(1e-4)]
    [InlineData(1e-8)]
    public void DensityTransformsInFarUnitsAndWithNearlySingularOutputScale(double epsilon)
    {
        var mean = new double[2, 2]; var rows = new double[,] { { 1, 0.8 }, { 0.8, 2 } };
        var psi = new double[,] { { 1, 1 - epsilon }, { 1 - epsilon, 1 } };
        var y = new double[,] { { 0.2, 0.2 }, { 0.5, 0.5 } };
        var original = new MatrixTPredictive(mean, rows, psi, 9);
        double[] units = [1e-6, -1e6]; var transformed = new double[2, 2]; var scale = new double[2, 2];
        for (int i = 0; i < 2; i++) for (int k = 0; k < 2; k++) { transformed[i, k] = y[i, k] * units[k]; scale[i, k] = psi[i, k] * units[i] * units[k]; }
        var target = new MatrixTPredictive(mean, rows, scale, 9);
        Assert.InRange(Math.Abs(target.LogDensity(transformed) - original.LogDensity(y) + 2 * Math.Log(Math.Abs(units[0] * units[1]))), 0, 1e-5);
    }

    [Fact]
    public void MatrixTJointDrawsHaveBothRowAndOutputCovariance()
    {
        var mean = new double[,] { { 1.2, -0.3 }, { 0.7, 2.0 } };
        var rows = new double[,] { { 1.0, 0.6 }, { 0.6, 2.0 } };
        var scale = new double[,] { { 2.0, -0.7 }, { -0.7, 1.0 } };
        var law = new MatrixTPredictive(mean, rows, scale, 9);
        var rng = new Xoshiro256PlusPlus(20261222); var moments = new RunningMoments[5];
        for (int i = 0; i < 80_000; i++)
        {
            var value = law.Sample(rng);
            moments[0].Add(value[0, 0]); moments[1].Add(value[1, 1]);
            moments[2].Add((value[0, 0] - mean[0, 0]) * (value[1, 0] - mean[1, 0]));
            moments[3].Add((value[0, 0] - mean[0, 0]) * (value[1, 1] - mean[1, 1]));
            moments[4].Add((value[1, 1] - mean[1, 1]) * (value[1, 1] - mean[1, 1]));
        }
        double[] expected = [1.2, 2, 0.6 * 2 / 6, 0.6 * -0.7 / 6, 2.0 / 6];
        for (int j = 0; j < moments.Length; j++) Assert.InRange(Math.Abs(moments[j].Mean - expected[j]) / Math.Sqrt(moments[j].Variance / moments[j].Count), 0, 4.5);
        Assert.Equal(0.6 * -0.7 / 6, law.Covariance(0, 0, 1, 1));
    }

    [Fact]
    public void ScalarPredictiveIsStudentT_AndMomentExistenceIsExplicit()
    {
        foreach (double nu in new[] { 0.5, 1.0, 2, 7, 1000 })
        {
            var law = new MatrixTPredictive(new double[,] { { 0.7 } }, new double[,] { { 1.3 } }, new double[,] { { 2.1 } }, nu);
            var t = StudentT.LocationScale(nu, 0.7, Math.Sqrt(1.3 * 2.1 / nu));
            foreach (double value in new[] { -12, 0.7, 5.0 }) Assert.InRange(Math.Abs(law.LogDensity(new double[,] { { value } }) - t.LogPdf(value)), 0, 2e-11);
            Assert.Equal(nu > 1, law.HasMean); Assert.Equal(nu > 2, law.HasCovariance);
            Assert.Equal(nu > 2, law.Covariance(0, 0, 0, 0) is not null);
            Assert.Equal(t.Quantile(0.95), law.Marginal(0, 0).Quantile(0.95));
        }
        Assert.Throws<ArgumentException>(() => new MatrixTPredictive(new double[1, 2], new double[,] { { 0 } }, new double[,] { { 1, 0 }, { 0, 1 } }, 3));
    }
}
