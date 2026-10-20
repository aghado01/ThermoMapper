using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// R04 / Q2: under robust weights the coefficient prior is the <i>unweighted</i> g-prior β ~ N(0, mσ²(ZᵀZ)⁻¹), so the
/// prior on the curve does not depend on the weights. Checked against a dense reference:
/// log p(y | Z, w) = ½log|G/m| − ½log|Zᵀ(W + I/m)Z| + ½Σ log wᵢ − (m/2)·log S, S = yᵀWy − hᵀ(Zᵀ(W + I/m)Z)⁻¹h, h = ZᵀWy.
/// </summary>
public sealed class RobustWeightTests
{
    [Fact]
    public void WeightedMarginal_IsTheUnweightedGPriorMarginal()
    {
        (double[,] z, double[] y, double[] w) = Problem();
        Dense d = Dense.Of(z, y, w);

        double expected = 0.5 * (d.LogDetG - z.GetLength(1) * Math.Log(y.Length)) - 0.5 * d.LogDetA
                          + 0.5 * d.SumLogW - 0.5 * y.Length * Math.Log(d.S);

        Assert.Equal(expected, new WeightedNormalModel().LogMarginalLikelihood(z, y, w), 9);
    }

    [Fact]
    public void WeightedDraws_HaveTheUnweightedGPriorConditionalMeanAndCovariance()
    {
        (double[,] z, double[] y, double[] w) = Problem();
        Dense d = Dense.Of(z, y, w);
        int m = y.Length, nu = z.GetLength(1);
        double scale = d.S / (m - 2);   // E[σ²] under IG(m/2, S/2)

        IObservationModel model = new WeightedNormalModel();
        var rng = new Xoshiro256PlusPlus(seed: 5);
        const int n = 40_000;
        var mean = new double[nu];
        var second = new double[nu, nu];
        for (int s = 0; s < n; s++)
        {
            double[] beta = model.SampleCoefficients(z, y, w, rng);
            for (int i = 0; i < nu; i++)
            {
                mean[i] += beta[i];
                for (int j = 0; j < nu; j++) second[i, j] += beta[i] * beta[j];
            }
        }
        for (int i = 0; i < nu; i++) mean[i] /= n;

        for (int i = 0; i < nu; i++)
        {
            double sd = Math.Sqrt(scale * d.InvA[i, i]);
            Assert.True(Math.Abs(mean[i] - d.Mean[i]) < 5.0 * sd / Math.Sqrt(n), $"mean[{i}] {mean[i]:R} vs {d.Mean[i]:R}");
            for (int j = 0; j < nu; j++)
            {
                double c = second[i, j] / n - mean[i] * mean[j];
                double expected = scale * d.InvA[i, j];
                double tol = 0.05 * scale * Math.Sqrt(d.InvA[i, i] * d.InvA[j, j]);
                Assert.True(Math.Abs(c - expected) < tol, $"cov[{i},{j}] {c:R} vs {expected:R}");
            }
        }
    }

    [Fact]
    public void Resampler_WithANonNormalModel_IsRejected()
    {
        Assert.Throws<ArgumentException>(() =>
            new BarsEnsemble(new SplineBasis(3), new PoissonModel(), new PoissonPrior(4.0), new LocalBetaKernel(50.0),
                             resampler: new StudentTWeights(4.0)));
    }

    [Fact]
    public void StudentTWeights_FollowTheirFullConditional()
    {
        // wᵢ | rᵢ, σ² ~ Gamma((ν+1)/2, scale 2/(ν + rᵢ²/σ²)), so E[wᵢ] = (ν+1)/(ν + rᵢ²/σ²).
        double[] r = [0.0, 1.0, 3.0, 10.0];
        const double sigma2 = 1.0, nu = 4.0;
        var resampler = new StudentTWeights(nu);
        var rng = new Xoshiro256PlusPlus(seed: 21);
        const int n = 40_000;
        var mean = new double[r.Length];
        for (int s = 0; s < n; s++)
        {
            double[] w = resampler.Resample(r, sigma2, rng);
            for (int i = 0; i < r.Length; i++) mean[i] += w[i] / n;
        }
        for (int i = 0; i < r.Length; i++)
        {
            double expected = (nu + 1.0) / (nu + r[i] * r[i] / sigma2);
            Assert.True(Math.Abs(mean[i] - expected) < 0.02 * expected, $"E[w{i}] {mean[i]:R} vs {expected:R}");
        }
    }

    [Fact]
    public void ContaminatedNormalWeights_FollowTheirPosteriorResponsibility()
    {
        // P(outlier | r, σ²) = π·N(r; 0, κ²σ²) / (π·N(r; 0, κ²σ²) + (1−π)·N(r; 0, σ²)); outliers get weight 1/κ².
        double[] r = [0.0, 2.0, 3.5, 8.0];
        const double sigma2 = 1.0, pi = 0.1, kappa = 12.0;
        var resampler = new ContaminatedNormalWeights(pi, kappa);
        var rng = new Xoshiro256PlusPlus(seed: 22);
        const int n = 40_000;
        var outlierFrac = new double[r.Length];
        for (int s = 0; s < n; s++)
        {
            double[] w = resampler.Resample(r, sigma2, rng);
            for (int i = 0; i < r.Length; i++) if (w[i] < 1.0) outlierFrac[i] += 1.0 / n;
        }
        for (int i = 0; i < r.Length; i++)
        {
            double good = (1.0 - pi) * Math.Exp(-0.5 * r[i] * r[i]);
            double bad = pi * Math.Exp(-0.5 * r[i] * r[i] / (kappa * kappa)) / kappa;
            double expected = bad / (good + bad);
            Assert.True(Math.Abs(outlierFrac[i] - expected) < 0.01, $"P(outlier {i}) {outlierFrac[i]:R} vs {expected:R}");
        }
    }

    // The CoefficientDrawTests design with fixed, uneven positive weights (a few points downweighted hard).
    private static (double[,] Z, double[] Y, double[] W) Problem()
    {
        (double[,] z, double[] y) = CoefficientDrawTests.Problem();
        var w = new double[y.Length];
        for (int i = 0; i < w.Length; i++) w[i] = i % 7 == 3 ? 0.02 : 0.6 + 0.8 * ((i * 37) % 11) / 10.0;
        return (z, y, w);
    }

    private sealed record Dense(double LogDetG, double LogDetA, double SumLogW, double S, double[] Mean, double[,] InvA)
    {
        public static Dense Of(double[,] z, double[] y, double[] w)
        {
            int m = y.Length, nu = z.GetLength(1);
            var g = new double[nu, nu];
            var a = new double[nu, nu];
            var h = new double[nu];
            double yWy = 0.0, sumLogW = 0.0;
            for (int i = 0; i < m; i++)
            {
                yWy += w[i] * y[i] * y[i];
                sumLogW += Math.Log(w[i]);
                for (int p = 0; p < nu; p++)
                {
                    h[p] += z[i, p] * w[i] * y[i];
                    for (int q = 0; q < nu; q++)
                    {
                        g[p, q] += z[i, p] * z[i, q];
                        a[p, q] += (w[i] + 1.0 / m) * z[i, p] * z[i, q];
                    }
                }
            }
            double[,] invA = CoefficientDrawTests.Invert(a);
            var mean = new double[nu];
            double quad = 0.0;
            for (int p = 0; p < nu; p++)
            {
                double s = 0.0;
                for (int q = 0; q < nu; q++) s += invA[p, q] * h[q];
                mean[p] = s;
                quad += h[p] * s;
            }
            return new Dense(LogDet(g), LogDet(a), sumLogW, yWy - quad, mean, invA);
        }

        private static double LogDet(double[,] a0)
        {
            int n = a0.GetLength(0);
            var a = (double[,])a0.Clone();
            double logDet = 0.0;
            for (int c = 0; c < n; c++)
            {
                int piv = c;
                for (int r = c + 1; r < n; r++) if (Math.Abs(a[r, c]) > Math.Abs(a[piv, c])) piv = r;
                for (int k = 0; k < n; k++) (a[c, k], a[piv, k]) = (a[piv, k], a[c, k]);
                logDet += Math.Log(Math.Abs(a[c, c]));
                for (int r = c + 1; r < n; r++)
                {
                    double f = a[r, c] / a[c, c];
                    for (int k = c; k < n; k++) a[r, k] -= f * a[c, k];
                }
            }
            return logDet;
        }
    }
}
