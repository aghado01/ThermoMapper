using System;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class RelativeJacobiStageTwoTests(ITestOutputHelper output)
{
    [Fact]
    public void Ridge_NormalizationAndWeakDirectionSurviveColumnRescaling()
    {
        const double mass = 16;
        double rho = 1 - Math.ScaleB(1, -30);
        double[] expected = [mass * (1 + rho), mass, mass * (1 - rho)];
        double worst = 0;
        foreach (int exponent in new[] { -300, 0, 300 })
        {
            double[] scale = [Math.ScaleB(1, exponent), 1, Math.ScaleB(1, -exponent)];
            var gram = new double[3, 3];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    gram[i, j] = mass * scale[i] * scale[j] * (i == j ? 1 : i + j == 1 ? rho : 0);
            RidgeSpectrum spectrum = RidgeSpectrum.Of(gram, [scale[0], -scale[1], 0], rows: 19, mass: mass);
            Assert.Equal(19, spectrum.Rows);
            Assert.Equal(mass, spectrum.Mass);
            for (int j = 0; j < 3; j++)
            {
                double gap = Math.Abs(spectrum.Eigenvalues[j] / expected[j] - 1);
                worst = Math.Max(worst, gap);
                Assert.InRange(gap, 0, 2e-13);
            }
            foreach (double tau in new[] { 1e-10, 1.0, 1e10 })
                foreach (double b in new[] { 0.0, 0.3, 1.0 })
                {
                    RidgeTerms terms = spectrum.At(tau, b);
                    double h = tau / (tau + b * expected[2]), q = 2 / expected[2];
                    double logDet = -0.5 * expected.Sum(lambda => double.LogP1(b * lambda / tau));
                    Close(logDet, terms.LogDeterminant);
                    Close(q * h, terms.Shrinkage);
                    Close(-0.5 * expected.Sum(lambda => lambda / (b * lambda + tau)), terms.LogDeterminantSlope);
                    Close(q * h * h, terms.ShrinkageSlope);
                }
        }
        output.WriteLine($"Column scales 2^-300 through 2^300; weak eigenvalue {expected[2]:R}; worst eigenvalue relative gap {worst:E3}.");
    }

    [Fact]
    public void UnbalancedIntercepts_WhitenTheIndependentDenseCovariance()
    {
        foreach (int multiple in new[] { 1, 3 })
            foreach (double scale in new[] { 1e-3, 1.0, 1e3 })
            {
                int[] counts = [multiple, 2 * multiple, 5 * multiple, 20 * multiple];
                int[] codes = [.. counts.SelectMany((count, group) => Enumerable.Repeat(group, count))];
                int n = codes.Length;
                var grouping = new CategoricalPredictor("g", codes, ["a", "b", "c", "d"]);
                RandomInterceptBlock block = RandomInterceptBlock.Of(grouping, new BassRandomIntercept("g", scale));
                double norm = counts.Sum(c => (double)c * c);
                var covariance = new double[n, n];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                    {
                        double pi = (codes[i] == codes[j] ? 1 : 0) - (double)counts[codes[i]] * counts[codes[j]] / norm;
                        covariance[i, j] = (i == j ? 1 : 0) + scale * scale * pi;
                    }
                var cholesky = new DenseCholesky(n);
                cholesky.Decompose(covariance);
                Assert.True(cholesky.Status.PositiveDefinite);
                var w = new double[n, n];
                for (int j = 0; j < n; j++)
                {
                    var e = new double[n];
                    e[j] = 1;
                    block.Whiten(e);
                    for (int i = 0; i < n; i++) w[i, j] = e[i];
                }
                double[,] whitened = Multiply(Multiply(w, covariance), w);
                double worst = 0;
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        worst = Math.Max(worst, Math.Abs(whitened[i, j] - (i == j ? 1 : 0)));
                double tolerance = scale > 1 ? 1e-7 : 2e-12;
                Assert.InRange(worst, 0, tolerance);
                Assert.InRange(Math.Abs(cholesky.LogDet - block.LogDeterminant), 0, tolerance);
                double[] ones = [.. Enumerable.Repeat(1.0, n)];
                block.Whiten(ones);
                Assert.All(ones, value => Assert.InRange(Math.Abs(value - 1), 0, 2e-12));
                output.WriteLine($"n={n}, s={scale:G}: W V W - I {worst:E3}, log determinant gap {Math.Abs(cholesky.LogDet - block.LogDeterminant):E3}.");
            }
    }

    private static void Close(double expected, double actual) =>
        Assert.InRange(Math.Abs(actual - expected), 0, 2e-12 * Math.Max(1, Math.Abs(expected)));

    private static double[,] Multiply(double[,] a, double[,] b)
    {
        int n = a.GetLength(0);
        var c = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int k = 0; k < n; k++) c[i, j] += a[i, k] * b[k, j];
        return c;
    }
}
