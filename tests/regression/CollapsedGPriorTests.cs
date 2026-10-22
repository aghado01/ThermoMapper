using System;
using Maths.LinAlg;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The collapsed g-prior block (step 2.3). With β | σ² ~ N(0, σ²(τG)⁻¹), G = ZᵀZ, and the likelihood raised to the
/// inverse temperature b, the block's closed forms are checked against references computed a different way:
/// <list type="bullet">
/// <item>the scale and mean from dense Gauss–Jordan algebra, S = b·yᵀy − b²hᵀ((b + τ)G)⁻¹h, rather than from the
/// block's explicit residual;</item>
/// <item>the complete log marginal against numerical quadrature over σ²;</item>
/// <item>the draws against their conditional moments;</item>
/// <item>the cancellation-free scale against a constructed near-interpolating fit whose residual is known.</item>
/// </list>
/// </summary>
public sealed class CollapsedGPriorTests
{
    public static TheoryData<double, double> Settings => new()
    {
        { 1.0 / 30.0, 1.0 },   // BARS: τ = 1/m, untempered
        { 0.01, 1.0 },
        { 0.3, 0.25 },
        { 2.0, 0.7 },
    };

    [Theory]
    [MemberData(nameof(Settings))]
    public void Unweighted_MatchesDenseAlgebra_AtAnyTauAndTemperature(double tau, double b)
    {
        (double[,] z, double[] y) = CoefficientDrawTests.Problem();
        CollapsedGPrior block = BandedGPrior.Fit(new BandedDesign(z), y, PositiveDefinite)!;
        GPriorConditional c = block.At(tau, b);
        Dense d = Dense.Of(z, y, tau, b);

        Close(d.S, c.S, 1e-12);
        Close(d.LogDeterminantTerm, c.LogDeterminantTerm, 1e-12);
        Close(1.0 / (tau + b), c.CovarianceScale, 1e-15);
        for (int p = 0; p < d.Mean.Length; p++) Close(d.Mean[p], c.Mean[p], 1e-11);
    }

    public static TheoryData<double, double, double, double> MarginalSettings => new()
    {
        { 1.0 / 30.0, 1.0, 0.0, 0.0 },   // improper σ² prior (g1 = g2 = 0)
        { 1.0 / 30.0, 1.0, 2.0, 0.5 },
        { 0.3, 0.25, 3.0, 0.4 },
        { 2.0, 0.7, 0.5, 2.0 },
    };

    [Theory]
    [MemberData(nameof(MarginalSettings))]
    public void LogMarginal_MatchesQuadratureOverSigma2(double tau, double b, double g1, double g2)
    {
        // ∫ p_b(y | σ², τ) π(σ²) dσ², with the coefficient integral from dense algebra and σ² integrated numerically
        // on u = ln σ² (trapezoid, log-sum-exp). The improper prior is π(σ²) ∝ 1/σ², whose constant the block drops.
        (double[,] z, double[] y) = CoefficientDrawTests.Problem();
        int m = y.Length;
        Dense d = Dense.Of(z, y, tau, b);
        InverseGamma? prior = g1 > 0.0 ? InverseGamma.ShapeScale(g1, g2) : null;

        double LogIntegrand(double u)
        {
            double s2 = Math.Exp(u);
            double logLik = -0.5 * m * b * Math.Log(2.0 * Math.PI * s2) + d.LogDeterminantTerm - d.S / (2.0 * s2);
            double logPrior = prior is InverseGamma ig ? ig.LogPdf(s2) : -Math.Log(s2);
            return logLik + logPrior + u;   // + u: dσ² = σ² du
        }

        const int n = 400_000;
        const double lo = -40.0, hi = 40.0;
        double h = (hi - lo) / n, max = double.NegativeInfinity;
        var vals = new double[n + 1];
        for (int i = 0; i <= n; i++) { vals[i] = LogIntegrand(lo + i * h); max = Math.Max(max, vals[i]); }
        double sum = 0.0;
        for (int i = 0; i <= n; i++) sum += (i == 0 || i == n ? 0.5 : 1.0) * Math.Exp(vals[i] - max);
        double quadrature = max + Math.Log(sum * h);

        CollapsedGPrior block = BandedGPrior.Fit(new BandedDesign(z), y, PositiveDefinite)!;
        Assert.True(Math.Abs(block.At(tau, b).LogMarginal(prior) - quadrature) < 1e-8,
                    $"log marginal {block.At(tau, b).LogMarginal(prior):R} vs quadrature {quadrature:R}");
    }

    [Fact]
    public void Weighted_MatchesDenseAlgebra_AndIsFixedAtItsTauAndUntempered()
    {
        (double[,] z, double[] y) = CoefficientDrawTests.Problem();
        int m = y.Length, nu = z.GetLength(1);
        var w = new double[m];
        for (int i = 0; i < m; i++) w[i] = i % 7 == 3 ? 0.02 : 0.6 + 0.8 * ((i * 37) % 11) / 10.0;
        double tau = 1.0 / m;

        // Pa = Zᵀ(W + τI)Z, h = ZᵀWy; S = yᵀWy − hᵀPa⁻¹h; ½log|τG| − ½log|Pa| + ½Σ ln wᵢ.
        var g = new double[nu, nu];
        var pa = new double[nu, nu];
        var hv = new double[nu];
        double yWy = 0.0, sumLogW = 0.0;
        for (int i = 0; i < m; i++)
        {
            yWy += w[i] * y[i] * y[i];
            sumLogW += Math.Log(w[i]);
            for (int p = 0; p < nu; p++)
            {
                hv[p] += z[i, p] * w[i] * y[i];
                for (int q = 0; q < nu; q++) { g[p, q] += z[i, p] * z[i, q]; pa[p, q] += (w[i] + tau) * z[i, p] * z[i, q]; }
            }
        }
        double[,] paInv = CoefficientDrawTests.Invert(pa);
        var mean = new double[nu];
        double quad = 0.0;
        for (int p = 0; p < nu; p++) { for (int q = 0; q < nu; q++) mean[p] += paInv[p, q] * hv[q]; quad += hv[p] * mean[p]; }
        double logDet = 0.5 * (LogDet(g) + nu * Math.Log(tau)) - 0.5 * LogDet(pa) + 0.5 * sumLogW;

        CollapsedGPrior block = BandedGPrior.Fit(new BandedDesign(z), y, w, tau, PositiveDefinite)!;
        GPriorConditional c = block.At(tau);
        Close(yWy - quad, c.S, 1e-11);
        Close(logDet, c.LogDeterminantTerm, 1e-11);
        Assert.Equal(1.0, c.CovarianceScale);
        for (int p = 0; p < nu; p++) Close(mean[p], c.Mean[p], 1e-11);

        Assert.Throws<ArgumentException>(() => block.At(2.0 * tau));
        Assert.Throws<ArgumentException>(() => block.At(tau, 0.5));
    }

    [Fact]
    public void Draws_HaveTheirConditionalMoments_AtATemperedRungWithAProperPrior()
    {
        // σ² | y ~ IG(g1 + mb/2, g2 + S/2); β | σ², y ~ N(b/(b+τ)·G⁻¹h, σ²G⁻¹/(b+τ)), so Cov β = E[σ²]·G⁻¹/(b+τ).
        (double[,] z, double[] y) = CoefficientDrawTests.Problem();
        int m = y.Length, nu = z.GetLength(1);
        const double tau = 0.05, b = 0.5;
        var prior = InverseGamma.ShapeScale(3.0, 0.4);
        Dense d = Dense.Of(z, y, tau, b);
        InverseGamma post = InverseGamma.ShapeScale(3.0 + 0.5 * m * b, 0.4 + 0.5 * d.S);

        GPriorConditional c = BandedGPrior.Fit(new BandedDesign(z), y, PositiveDefinite)!.At(tau, b);
        Close(post.Mean, c.Sigma2Posterior(prior).Mean, 1e-12);

        var rng = new Xoshiro256PlusPlus(seed: 23);
        const int n = 40_000;
        var beta = new double[nu];
        var mean = new double[nu];
        var second = new double[nu, nu];
        double s2Sum = 0.0;
        for (int s = 0; s < n; s++)
        {
            s2Sum += c.Sample(prior, rng, beta);
            for (int i = 0; i < nu; i++) { mean[i] += beta[i] / n; for (int j = 0; j < nu; j++) second[i, j] += beta[i] * beta[j] / n; }
        }
        Assert.True(Math.Abs(s2Sum / n - post.Mean) < 5.0 * Math.Sqrt(post.Variance / n), $"E[σ²] {s2Sum / n:R} vs {post.Mean:R}");
        double covScale = post.Mean / (tau + b);
        for (int i = 0; i < nu; i++)
        {
            double sd = Math.Sqrt(covScale * d.GramInverse[i, i]);
            Assert.True(Math.Abs(mean[i] - d.Mean[i]) < 5.0 * sd / Math.Sqrt(n), $"mean[{i}] {mean[i]:R} vs {d.Mean[i]:R}");
            for (int j = 0; j < nu; j++)
            {
                double cov = second[i, j] - mean[i] * mean[j];
                double expected = covScale * d.GramInverse[i, j];
                double tol = 0.05 * covScale * Math.Sqrt(d.GramInverse[i, i] * d.GramInverse[j, j]);
                Assert.True(Math.Abs(cov - expected) < tol, $"cov[{i},{j}] {cov:R} vs {expected:R}");
            }
        }
    }

    [Fact]
    public void ConditionalVariance_IsThePosteriorMeanVarianceTimesTheQuadratic()
    {
        (double[,] z, double[] y) = CoefficientDrawTests.Problem();
        int m = y.Length, nu = z.GetLength(1);
        double tau = 1.0 / m;
        Dense d = Dense.Of(z, y, tau, 1.0);
        double eSigma2 = 0.5 * d.S / (0.5 * m - 1.0);

        GPriorConditional c = BandedGPrior.Fit(new BandedDesign(z), y, PositiveDefinite)!.At(tau);
        double[] v = c.ConditionalVariance(null, z);
        for (int t = 0; t < m; t++)
        {
            double q = 0.0;
            for (int i = 0; i < nu; i++) for (int j = 0; j < nu; j++) q += z[t, i] * d.GramInverse[i, j] * z[t, j];
            Close(eSigma2 * q / (1.0 + tau), v[t], 1e-10);
        }
    }

    [Fact]
    public void Scale_IsCancellationFree_ForANearInterpolatingFitAtSmallTau()
    {
        // y = Zβ₀ + e with e orthogonal to the column space, so ‖r_OLS‖² = ‖e‖² and Q = ‖Zβ₀‖² are known, and
        // S = ‖e‖² + τ/(1+τ)·Q. At τ = 1e-18 and ‖e‖ = 1e-8 the reference is about 1e-16 while yᵀy is about 50, so
        // the subtracted form yᵀy − Q/(1+τ) is a multiple of ulp(yᵀy) ≈ 7e-15: its relative error is 1 or above 40.
        // The explicit residual keeps the digits.
        (double[,] z, double[] _) = CoefficientDrawTests.Problem();
        int m = z.GetLength(0), nu = z.GetLength(1);
        var rng = new Xoshiro256PlusPlus(seed: 24);
        double[] beta0 = [1.3, -0.4, 2.1, 0.7, -1.6];
        var fit = new double[m];
        for (int i = 0; i < m; i++) for (int p = 0; p < nu; p++) fit[i] += z[i, p] * beta0[p];
        var e = new double[m];
        for (int i = 0; i < m; i++) e[i] = rng.NextDouble() - 0.5;
        for (int pass = 0; pass < 3; pass++) OrthogonalizeAgainstColumns(e, z);
        double eNorm = 0.0;
        foreach (double v in e) eNorm += v * v;
        double scale = 1e-8 / Math.Sqrt(eNorm);
        for (int i = 0; i < m; i++) e[i] *= scale;

        var y = new double[m];
        double q = 0.0, ee = 0.0, yy = 0.0;
        for (int i = 0; i < m; i++) { y[i] = fit[i] + e[i]; q += fit[i] * fit[i]; ee += e[i] * e[i]; yy += y[i] * y[i]; }
        const double tau = 1e-18;
        double reference = ee + tau / (1.0 + tau) * q;
        double subtracted = yy - q / (1.0 + tau);

        GPriorConditional c = BandedGPrior.Fit(new BandedDesign(z), y, PositiveDefinite)!.At(tau);
        Assert.True(Math.Abs(c.S - reference) <= 1e-3 * reference, $"S {c.S:R} vs {reference:R}");
        Assert.True(Math.Abs(subtracted - reference) >= 0.5 * reference, $"the subtracted form {subtracted:R} should have lost its digits"); 
    }

    [Fact]
    public void Settings_AreValidated()
    {
        (double[,] z, double[] y) = CoefficientDrawTests.Problem();
        CollapsedGPrior block = BandedGPrior.Fit(new BandedDesign(z), y, PositiveDefinite)!;
        Assert.Throws<ArgumentOutOfRangeException>(() => block.At(0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => block.At(0.1, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => block.At(double.NaN));
    }

    private static bool PositiveDefinite(ISpdFactor factor) => factor.Status.PositiveDefinite;

    private static void Close(double expected, double actual, double rel)
        => Assert.True(Math.Abs(actual - expected) <= rel * Math.Max(1.0, Math.Abs(expected)), $"expected {expected:R}, got {actual:R}");

    // Dense reference for the unweighted block: G = ZᵀZ, h = Zᵀy, posterior precision (b + τ)G.
    private sealed record Dense(double S, double LogDeterminantTerm, double[] Mean, double[,] GramInverse)
    {
        public static Dense Of(double[,] z, double[] y, double tau, double b)
        {
            int m = y.Length, nu = z.GetLength(1);
            var g = new double[nu, nu];
            var h = new double[nu];
            double yy = 0.0;
            for (int i = 0; i < m; i++)
            {
                yy += y[i] * y[i];
                for (int p = 0; p < nu; p++) { h[p] += z[i, p] * y[i]; for (int q = 0; q < nu; q++) g[p, q] += z[i, p] * z[i, q]; }
            }
            double[,] gInv = CoefficientDrawTests.Invert(g);
            var gh = new double[nu];
            double hGh = 0.0;
            for (int p = 0; p < nu; p++) { for (int q = 0; q < nu; q++) gh[p] += gInv[p, q] * h[q]; hGh += h[p] * gh[p]; }
            var mean = new double[nu];
            for (int p = 0; p < nu; p++) mean[p] = b / (b + tau) * gh[p];
            return new Dense(b * yy - b * b / (b + tau) * hGh, 0.5 * nu * Math.Log(tau / (tau + b)), mean, gInv);
        }
    }

    private static void OrthogonalizeAgainstColumns(double[] e, double[,] z)
    {
        int m = z.GetLength(0), nu = z.GetLength(1);
        var basis = new double[nu][];
        for (int p = 0; p < nu; p++)
        {
            var col = new double[m];
            for (int i = 0; i < m; i++) col[i] = z[i, p];
            for (int k = 0; k < p; k++) Project(col, basis[k]);
            double n2 = 0.0;
            foreach (double v in col) n2 += v * v;
            double inv = 1.0 / Math.Sqrt(n2);
            for (int i = 0; i < m; i++) col[i] *= inv;
            basis[p] = col;
        }
        foreach (double[] q in basis) Project(e, q);

        static void Project(double[] v, double[]? unit)
        {
            if (unit is null) return;
            double d = 0.0;
            for (int i = 0; i < v.Length; i++) d += v[i] * unit[i];
            for (int i = 0; i < v.Length; i++) v[i] -= d * unit[i];
        }
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
