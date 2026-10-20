using System;
using Maths.Regression.Bgp;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The GP's constant mean is integrated out under a flat prior. The evidence is
/// log ∫ N(y; m1, K) dm = −½[yᵀK⁻¹y − (1ᵀK⁻¹y)²/(1ᵀK⁻¹1)] − ½log|K| − ½log(1ᵀK⁻¹1) − ((n−1)/2)log 2π, and the
/// predictive mean is m̂ + k*ᵀK⁻¹(y − m̂1) with the GLS mean m̂ = 1ᵀK⁻¹y / 1ᵀK⁻¹1. The reference is independent of
/// that algebra: a dense GP with a wide Normal prior N(0, v) on m, whose evidence plus ½log(2πv) and whose
/// predictive mean (k* + v1)ᵀ(K + v11ᵀ)⁻¹y converge to the flat-prior values as v grows.
/// </summary>
public sealed class BgpConstantMeanTests
{
    private const double Sigma2 = 0.01, V = 1e6;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void EvidenceAndPrediction_MatchTheWidePriorLimit()
    {
        // Inputs already span [0, 1] from 0, so the fitted similarity map is the identity and the reference can use
        // raw coordinates. The response has a mean of about 3 that the GP must integrate out.
        const int n = 30;
        var rng = new Xoshiro256PlusPlus(seed: 77);
        var x = new double[n, 1];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = i / (n - 1.0);
            y[i] = 3.0 + Math.Sin(4.0 * x[i, 0]) + 0.1 * Gaussian(rng);
        }
        var xTest = new double[,] { { 0.1 }, { 0.55 }, { 1.4 } };   // the last point extrapolates
        var kernel = new SquaredExponentialKernel();
        var gp = new GpRegression(x, y, Sigma2, kernel);
        Assert.Equal(1.0, gp.Transform.InputScale);

        foreach (double t in new[] { 0.005, 0.05, 0.4 })
        {
            GpFit fit = gp.Fit(t);
            (double logEvidence, double mean, double[] pred) = WidePrior(x, y, xTest, kernel, t);
            Assert.True(Math.Abs(fit.LogMarginal - logEvidence) < 1e-5, $"evidence at t={t}: {fit.LogMarginal} vs wide-prior {logEvidence}");
            Assert.True(Math.Abs(fit.Mean - mean) < 1e-5, $"GLS mean at t={t}: {fit.Mean} vs wide-prior {mean}");
            double[] p = gp.PredictMean(fit, xTest);
            for (int i = 0; i < p.Length; i++)
                Assert.True(Math.Abs(p[i] - pred[i]) < 1e-5, $"prediction {i} at t={t}: {p[i]} vs wide-prior {pred[i]}");
        }
    }

    [Fact]
    public void FarFromTheData_PredictionReturnsToTheGlsMean()
    {
        const int n = 25;
        var x = new double[n, 1];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = i / (n - 1.0);
            y[i] = 7.0 + 0.5 * Math.Cos(6.0 * x[i, 0]);
        }
        var gp = new GpRegression(x, y, Sigma2, new SquaredExponentialKernel());
        GpFit fit = gp.Fit(0.01);
        double far = gp.PredictMean(fit, new double[,] { { 10.0 } })[0];
        Assert.Equal(fit.Mean, far, 9);
        Assert.InRange(fit.Mean, 6.5, 7.5);
    }

    // Dense wide-prior reference: y ~ N(0, K + v11ᵀ), and f* = m + g(x*) has cov(f*, y) = k* + v1.
    private static (double LogEvidence, double Mean, double[] Pred) WidePrior(
        double[,] x, double[] y, double[,] xTest, IGpKernel kernel, double t)
    {
        int n = y.Length;
        var c = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                double d = x[i, 0] - x[j, 0];
                c[i, j] = kernel.Evaluate(d * d, t) + V + (i == j ? Sigma2 : 0.0);
            }
        double[,] l = Cholesky(c);
        double logDet = 0.0;
        for (int i = 0; i < n; i++) logDet += 2.0 * Math.Log(l[i, i]);
        double[] w = Solve(l, y);
        double quad = 0.0;
        for (int i = 0; i < n; i++) quad += y[i] * w[i];
        double logEvidence = -0.5 * quad - 0.5 * logDet - 0.5 * n * Math.Log(2.0 * Math.PI) + 0.5 * Math.Log(2.0 * Math.PI * V);

        double mean = 0.0;   // E[m | y] = v·1ᵀ(K + v11ᵀ)⁻¹y
        for (int i = 0; i < n; i++) mean += V * w[i];

        var pred = new double[xTest.GetLength(0)];
        for (int a = 0; a < pred.Length; a++)
            for (int j = 0; j < n; j++)
            {
                double d = xTest[a, 0] - x[j, 0];
                pred[a] += (kernel.Evaluate(d * d, t) + V) * w[j];
            }
        return (logEvidence, mean, pred);
    }

    private static double[,] Cholesky(double[,] a)
    {
        int n = a.GetLength(0);
        var l = new double[n, n];
        for (int j = 0; j < n; j++)
        {
            double s = a[j, j];
            for (int k = 0; k < j; k++) s -= l[j, k] * l[j, k];
            if (!(s > 0.0)) throw new InvalidOperationException("matrix is not positive definite");
            l[j, j] = Math.Sqrt(s);
            for (int i = j + 1; i < n; i++)
            {
                double r = a[i, j];
                for (int k = 0; k < j; k++) r -= l[i, k] * l[j, k];
                l[i, j] = r / l[j, j];
            }
        }
        return l;
    }

    private static double[] Solve(double[,] l, double[] b)
    {
        int n = b.Length;
        var z = new double[n];
        for (int i = 0; i < n; i++)
        {
            double s = b[i];
            for (int k = 0; k < i; k++) s -= l[i, k] * z[k];
            z[i] = s / l[i, i];
        }
        for (int i = n - 1; i >= 0; i--)
        {
            double s = z[i];
            for (int k = i + 1; k < n; k++) s -= l[k, i] * z[k];
            z[i] = s / l[i, i];
        }
        return z;
    }
}
