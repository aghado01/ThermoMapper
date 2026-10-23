using System;
using Maths.Regression.Bgp;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The GP's conditional variance with the constant mean integrated out (<see cref="GpRegression.PredictVariance"/>),
/// checked against an independent route: a proper N(0, v) prior on the mean with v = 10⁸, the flat prior's limit,
/// conditioned as one joint Gaussian, Var = (k** + v) − (k* + v1)ᵀ(K + v11ᵀ)⁻¹(k* + v1), and solved by Gaussian
/// elimination in <see cref="decimal"/> so the cancellation against v costs nothing.
/// </summary>
public sealed class BgpVarianceTests
{
    [Fact]
    public void PredictVariance_MatchesTheProperPriorLimit()
    {
        var rng = new Xoshiro256PlusPlus(seed: 12);
        int n = 15;
        var x = new double[n, 1];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = (i + 0.5) / n;
            y[i] = Math.Sin(2.0 * Math.PI * x[i, 0]) + 0.1 * (rng.NextDouble() - 0.5);
        }
        const double sigma2 = 0.01, t = 0.02;
        var kernel = new SquaredExponentialKernel();
        var gp = new GpRegression(x, y, sigma2, kernel, GpInputScaling.None);
        double[] test = [0.0, 0.31, 0.5, 0.97, 1.4, 2.5];   // inside, at the edges, and well outside the data
        var xTest = new double[test.Length, 1];
        for (int a = 0; a < test.Length; a++) xTest[a, 0] = test[a];

        double[] variance = gp.PredictVariance(gp.Fit(t), xTest);

        const decimal v = 100_000_000m;
        for (int a = 0; a < test.Length; a++)
        {
            var big = new decimal[n, n];
            var kStar = new decimal[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    double s = (x[i, 0] - x[j, 0]) * (x[i, 0] - x[j, 0]);
                    big[i, j] = (decimal)kernel.Evaluate(s, t) + v + (i == j ? (decimal)sigma2 : 0m);
                }
                double sa = (test[a] - x[i, 0]) * (test[a] - x[i, 0]);
                kStar[i] = (decimal)kernel.Evaluate(sa, t) + v;
            }
            decimal reference = (decimal)kernel.Evaluate(0.0, t) + v - Quadratic(big, kStar);
            double expected = (double)reference;
            Assert.True(Math.Abs(variance[a] - expected) <= 1e-9 * expected + 1e-12,
                        $"x* = {test[a]}: {variance[a]:R} against the proper-prior limit {expected:R}");
        }
    }

    // bᵀA⁻¹b by Gaussian elimination with partial pivoting, in decimal.
    private static decimal Quadratic(decimal[,] a, decimal[] b)
    {
        int n = b.Length;
        var m = (decimal[,])a.Clone();
        var r = (decimal[])b.Clone();
        for (int c = 0; c < n; c++)
        {
            int p = c;
            for (int i = c + 1; i < n; i++) if (Math.Abs(m[i, c]) > Math.Abs(m[p, c])) p = i;
            if (p != c)
            {
                for (int j = 0; j < n; j++) (m[c, j], m[p, j]) = (m[p, j], m[c, j]);
                (r[c], r[p]) = (r[p], r[c]);
            }
            for (int i = c + 1; i < n; i++)
            {
                decimal f = m[i, c] / m[c, c];
                for (int j = c; j < n; j++) m[i, j] -= f * m[c, j];
                r[i] -= f * r[c];
            }
        }
        var sol = new decimal[n];
        for (int i = n - 1; i >= 0; i--)
        {
            decimal s = r[i];
            for (int j = i + 1; j < n; j++) s -= m[i, j] * sol[j];
            sol[i] = s / m[i, i];
        }
        decimal q = 0m;
        for (int i = 0; i < n; i++) q += b[i] * sol[i];
        return q;
    }
}
