using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// R12 / Q1: the Normal model's coefficient draws come from the exact conditional posterior given the knots —
/// σ² ~ IG(m/2, a/2), then β ~ N(μ, σ²·(m/(m+1))·A⁻¹) — checked against a dense reference that shares no code with
/// the band factor.
/// </summary>
public sealed class CoefficientDrawTests
{
    [Fact]
    public void NormalDraws_HaveTheConditionalPosteriorMeanAndCovariance()
    {
        (double[,] z, double[] y) = Problem();
        int m = y.Length, nu = z.GetLength(1);
        (double[] mu, double[,] cov) = Reference(z, y);

        IObservationModel model = new WeightedNormalModel();
        var rng = new Xoshiro256PlusPlus(seed: 77);
        const int n = 40_000;
        var mean = new double[nu];
        var second = new double[nu, nu];
        for (int s = 0; s < n; s++)
        {
            double[] beta = model.SampleCoefficients(z, y, null, rng);
            for (int i = 0; i < nu; i++)
            {
                mean[i] += beta[i];
                for (int j = 0; j < nu; j++) second[i, j] += beta[i] * beta[j];
            }
        }
        for (int i = 0; i < nu; i++) mean[i] /= n;

        for (int i = 0; i < nu; i++)
        {
            double sd = Math.Sqrt(cov[i, i]);
            Assert.True(Math.Abs(mean[i] - mu[i]) < 5.0 * sd / Math.Sqrt(n), $"mean[{i}] {mean[i]:R} vs {mu[i]:R}");
            for (int j = 0; j < nu; j++)
            {
                double c = second[i, j] / n - mean[i] * mean[j];
                double scale = Math.Sqrt(cov[i, i] * cov[j, j]);
                Assert.True(Math.Abs(c - cov[i, j]) < 0.05 * scale, $"cov[{i},{j}] {c:R} vs {cov[i, j]:R}");
            }
        }
    }

    [Fact]
    public void NormalConditionalVariance_MatchesTheDenseReference()
    {
        (double[,] z, double[] y) = Problem();
        (_, double[,] cov) = Reference(z, y);
        double[,] readout = new SplineBasis(3).Design(new KnotConfig(new[] { 0.5 }), new[] { 0.05, 0.37, 0.5, 0.81, 0.99 });

        double[]? v = new WeightedNormalModel().ConditionalVariance(z, y, null, readout);

        Assert.NotNull(v);
        for (int t = 0; t < v!.Length; t++)
            Assert.Equal(QuadraticForm(readout, t, cov), v[t], 12);
    }

    [Fact]
    public void WhittleAndPoissonConditionalVariance_MatchTheirLaplaceReferences()
    {
        (double[,] z, double[] s) = Problem();
        int m = s.Length, nu = z.GetLength(1);
        double[,] readout = new SplineBasis(3).Design(new KnotConfig(new[] { 0.5 }), new[] { 0.1, 0.5, 0.9 });

        // Whittle: J = BᵀB.
        var periodogram = new double[m];
        for (int i = 0; i < m; i++) periodogram[i] = Math.Exp(s[i]);
        double[,] whittleCov = Invert(Gram(z, null));
        double[]? vw = new WhittleModel().ConditionalVariance(z, periodogram, null, readout);
        for (int t = 0; t < readout.GetLength(0); t++) Assert.Equal(QuadraticForm(readout, t, whittleCov), vw![t], 10);

        // Poisson: H = BᵀŴB + Λ at the mode, Ŵ = exp(η̂), with the g-prior's Λ = (ȳ/m)·BᵀB, or Λ = 0 at BIC's MLE.
        var counts = new double[m];
        for (int i = 0; i < m; i++) counts[i] = Math.Round(5.0 + 3.0 * s[i]);
        double priorScale = 0.0;
        foreach (double c in counts) priorScale += c / ((double)m * m);   // ȳ/m
        foreach (PoissonMarginal marginal in new[] { PoissonMarginal.GPrior, PoissonMarginal.Bic })
        {
            var poisson = new PoissonModel(marginal);
            double[] betaHat = poisson.PosteriorMeanCoefficients(z, counts, null);
            var w = new double[m];
            for (int i = 0; i < m; i++)
            {
                double eta = 0.0;
                for (int j = 0; j < nu; j++) eta += z[i, j] * betaHat[j];
                w[i] = Math.Exp(eta);
            }
            double[,] h = Gram(z, w);
            if (marginal == PoissonMarginal.GPrior)
            {
                double[,] gram = Gram(z, null);
                for (int p = 0; p < nu; p++) for (int q = 0; q < nu; q++) h[p, q] += priorScale * gram[p, q];
            }
            double[,] poissonCov = Invert(h);
            double[]? vp = poisson.ConditionalVariance(z, counts, null, readout);
            for (int t = 0; t < readout.GetLength(0); t++) Assert.Equal(QuadraticForm(readout, t, poissonCov), vp![t], 10);
        }
    }

    [Fact]
    public void UncertaintyLabels_MatchEachModel()
    {
        Assert.Equal(CoefficientUncertainty.Exact, ((IObservationModel)new WeightedNormalModel()).Uncertainty);
        Assert.Equal(CoefficientUncertainty.Approximate, ((IObservationModel)new WhittleModel()).Uncertainty);
        Assert.Equal(CoefficientUncertainty.Approximate, ((IObservationModel)new PoissonModel()).Uncertainty);
        Assert.Equal(CoefficientUncertainty.KnotOnly, ((IObservationModel)new GibbsLossModel()).Uncertainty);
    }

    [Fact]
    public void BarsFitSd_IsAFinitePosteriorSpread_AndCarriesTheModelsLabel()
    {
        const int n = 100;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.15 * Math.Sin(53.0 * x[i] + 0.7);
        }
        BarsResult r = new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(4.0),
                                        new LocalBetaKernel(50.0))
            .Run(x, y, grid: x, chains: 2, masterSeed: 13, burn: 800, samples: 600);

        Assert.Equal(CoefficientUncertainty.Exact, r.Uncertainty);
        double mean = 0.0;
        foreach (double sd in r.FitSd) { Assert.True(double.IsFinite(sd) && sd > 0.0); mean += sd; }
        mean /= r.FitSd.Length;
        // A smooth fit of noise amplitude 0.15 over 100 points: pointwise sd well below the noise, well above zero.
        Assert.InRange(mean, 0.01, 0.15);
    }

    private static double QuadraticForm(double[,] rows, int t, double[,] cov)
    {
        int nu = rows.GetLength(1);
        double q = 0.0;
        for (int i = 0; i < nu; i++)
            for (int j = 0; j < nu; j++) q += rows[t, i] * cov[i, j] * rows[t, j];
        return q;
    }

    private static double[,] Gram(double[,] z, double[]? w)
    {
        int m = z.GetLength(0), nu = z.GetLength(1);
        var a = new double[nu, nu];
        for (int i = 0; i < m; i++)
        {
            double wi = w is null ? 1.0 : w[i];
            for (int p = 0; p < nu; p++) for (int q = 0; q < nu; q++) a[p, q] += wi * z[i, p] * z[i, q];
        }
        return a;
    }

    // A cubic spline design with one interior knot on 30 points, and a fixed noisy response.
    internal static (double[,] Z, double[] Y) Problem()
    {
        const int m = 30;
        var x = new double[m];
        var y = new double[m];
        for (int i = 0; i < m; i++)
        {
            x[i] = (i + 0.5) / m;
            y[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.2 * Math.Sin(37.0 * x[i] + 1.3);
        }
        return (new SplineBasis(3).Design(new KnotConfig(new[] { 0.5 }), x), y);
    }

    // The unit-information conditional posterior from a dense Gauss–Jordan inverse of A = ZᵀZ:
    // μ = (m/(m+1))A⁻¹b and, with σ² integrated out, Cov β = (a/(m−2))·(m/(m+1))·A⁻¹.
    internal static (double[] Mu, double[,] Cov) Reference(double[,] z, double[] y)
    {
        int m = y.Length, nu = z.GetLength(1);
        var a = new double[nu, nu];
        var b = new double[nu];
        double yy = 0.0;
        for (int i = 0; i < m; i++)
        {
            yy += y[i] * y[i];
            for (int p = 0; p < nu; p++)
            {
                b[p] += z[i, p] * y[i];
                for (int q = 0; q < nu; q++) a[p, q] += z[i, p] * z[i, q];
            }
        }
        double[,] inv = Invert(a);
        double shrink = (double)m / (m + 1);
        var mu = new double[nu];
        double quad = 0.0;
        for (int p = 0; p < nu; p++)
        {
            double s = 0.0;
            for (int q = 0; q < nu; q++) s += inv[p, q] * b[q];
            mu[p] = shrink * s;
            quad += b[p] * s;
        }
        double aScalar = yy - shrink * quad;
        double scale = aScalar / (m - 2) * shrink;
        var cov = new double[nu, nu];
        for (int p = 0; p < nu; p++) for (int q = 0; q < nu; q++) cov[p, q] = scale * inv[p, q];
        return (mu, cov);
    }

    internal static double[,] Invert(double[,] a0)
    {
        int n = a0.GetLength(0);
        var a = (double[,])a0.Clone();
        var inv = new double[n, n];
        for (int i = 0; i < n; i++) inv[i, i] = 1.0;
        for (int c = 0; c < n; c++)
        {
            int piv = c;
            for (int r = c + 1; r < n; r++) if (Math.Abs(a[r, c]) > Math.Abs(a[piv, c])) piv = r;
            for (int k = 0; k < n; k++) { (a[c, k], a[piv, k]) = (a[piv, k], a[c, k]); (inv[c, k], inv[piv, k]) = (inv[piv, k], inv[c, k]); }
            double d = a[c, c];
            for (int k = 0; k < n; k++) { a[c, k] /= d; inv[c, k] /= d; }
            for (int r = 0; r < n; r++)
            {
                if (r == c) continue;
                double f = a[r, c];
                if (f == 0.0) continue;
                for (int k = 0; k < n; k++) { a[r, k] -= f * a[c, k]; inv[r, k] -= f * inv[c, k]; }
            }
        }
        return inv;
    }
}
