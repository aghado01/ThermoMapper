using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The Poisson g-prior marginal (step 3.10, OQ12). The prior is N(β₀, Λ⁻¹), with Λ = (ȳ/m)·BᵀB and β₀ = log ȳ·1, the
/// null fit in a B-spline basis, whose functions sum to one. It keeps the mode finite where the counts under a basis
/// function are all zero, its Laplace marginal matches importance sampling of the same integral, and a zero response is
/// refused. Every reference here is computed densely from the prior's definition, not from the model's banded code.
/// </summary>
public sealed class PoissonGPriorTests
{
    private readonly ITestOutputHelper _out;
    public PoissonGPriorTests(ITestOutputHelper output) => _out = output;

    private static double[] Even(int n) => Enumerable.Range(0, n).Select(i => (i + 0.5) / n).ToArray();

    // Knuth's Poisson sampler (adequate for the moderate λ here).
    private static int PoissonSample(Xoshiro256PlusPlus rng, double lambda)
    {
        double l = Math.Exp(-lambda), p = 1.0;
        int k = 0;
        do { k++; p *= rng.NextDouble(); } while (p > l);
        return k - 1;
    }

    // OQ12's case: about three counts per bin below 0.797 and none beyond, with a knot at 0.8, so the last basis
    // function sees only zeros.
    private static (double[,] Design, double[] Y) Separated()
    {
        double[] x = Even(100);
        double[] y = x.Select((t, i) => t < 0.797 ? 2.0 + i % 3 : 0.0).ToArray();
        return (new SplineBasis(3).Design(new KnotConfig([0.8]), x), y);
    }

    [Fact]
    public void Separation_TheGPriorModeStaysBounded_WhereBicRunsToTheClamp()
    {
        (double[,] design, double[] y) = Separated();
        double[] gPrior = new PoissonModel().PosteriorMeanCoefficients(design, y, null);
        double[] bic = new PoissonModel(PoissonMarginal.Bic).PosteriorMeanCoefficients(design, y, null);
        _out.WriteLine($"[separation] g-prior β = [{string.Join(", ", gPrior.Select(b => b.ToString("F3")))}]");
        _out.WriteLine($"[separation] BIC β = [{string.Join(", ", bic.Select(b => b.ToString("F1")))}]");

        Assert.All(gPrior, b => Assert.InRange(b, -15.0, 5.0));
        Assert.True(double.IsFinite(new PoissonModel().LogMarginalLikelihood(design, y, null)));
        // BIC's documented caveat: the MLE lies at −∞ in the separated direction, and IRLS stops at the η clamp.
        Assert.True(bic.Min() < -100.0, $"the BIC pathway's separated coefficient was {bic.Min():F1}");

        // The mode is the stationary point of the log posterior: Bᵀ(y − μ̂) − Λ(β̂ − β₀) = 0.
        double[] gradient = Gradient(design, y, gPrior);
        Assert.All(gradient, g => Assert.InRange(g, -1e-6, 1e-6));
    }

    [Fact]
    public void ZeroResponse_IsRefused()
    {
        (double[,] design, double[] y) = Separated();
        var zeros = new double[y.Length];
        var model = new PoissonModel();
        Assert.Equal(DesignRejection.ZeroResponse, model.TemperedMarginal(design, zeros, null).Rejection);
        Assert.Equal(double.NegativeInfinity, model.LogMarginalLikelihood(design, zeros, null));
        // Without a g-prior the draw falls back to the fit, as when H is refused, rather than to BIC's ridge target.
        Assert.Equal(model.PosteriorMeanCoefficients(design, zeros, null),
                     model.SampleCoefficients(design, zeros, new Xoshiro256PlusPlus(seed: 1)));
        Assert.Equal(DesignRejection.None, model.TemperedMarginal(design, y, null).Rejection);
    }

    [Fact]
    public void LaplaceMarginal_MatchesImportanceSampling()
    {
        // A small design, ν = 5 and 6 at m = 40, where the Laplace approximation is at its least accurate.
        double[] x = Even(40);
        var rng = new Xoshiro256PlusPlus(seed: 310);
        double[] y = x.Select(t => (double)PoissonSample(rng, Math.Exp(1.2 + 0.8 * Math.Sin(2.0 * Math.PI * t)))).ToArray();
        var model = new PoissonModel();
        var estimates = new (double Laplace, double Sampled, double StdError)[2];
        double[][] configurations = [[0.5], [0.3, 0.7]];
        for (int c = 0; c < configurations.Length; c++)
        {
            double[,] design = new SplineBasis(3).Design(new KnotConfig(configurations[c]), x);
            double laplace = model.LogMarginalLikelihood(design, y, null);
            (double sampled, double se) = ImportanceSampled(design, y, model.PosteriorMeanCoefficients(design, y, null),
                                                            new Xoshiro256PlusPlus(seed: 3100 + c));
            estimates[c] = (laplace, sampled, se);
            _out.WriteLine($"[laplace-is] ν={design.GetLength(1)}: Laplace {laplace:F4}, sampled {sampled:F4} ± {se:F4}");
        }
        // The Laplace error at this size, measured: −0.017 and −0.024 for the two marginals and −0.008 for their ratio,
        // against standard errors of about 0.0015. The bounds hold it a hundred times below the Occam factor per
        // coefficient, ½·log(1 + m) ≈ 1.86.
        foreach ((double laplace, double sampled, double se) in estimates)
        {
            Assert.True(se < 0.003, $"importance-sampling standard error {se:F4}");
            Assert.InRange(laplace - sampled, -0.05, 0.05);
        }
        double ratioError = (estimates[1].Laplace - estimates[0].Laplace) - (estimates[1].Sampled - estimates[0].Sampled);
        Assert.InRange(ratioError, -0.02, 0.02);
    }

    // log ∫ L(β)·N(β; β₀, Λ⁻¹) dβ, dropping Σ log yᵢ! as the model does, by importance sampling from a multivariate
    // Student t (4 degrees of freedom) centred on the mode with scale H⁻¹, H = BᵀŴB + Λ. Returns the estimate and its
    // delta-method standard error.
    private static (double LogZ, double StdError) ImportanceSampled(double[,] design, double[] y, double[] mode,
                                                                    Xoshiro256PlusPlus rng)
    {
        const int draws = 200_000;
        const int df = 4;
        int m = y.Length, nu = mode.Length;
        double mean = y.Average(), scale = mean / m, level = Math.Log(mean);
        double[,] precision = Scaled(DenseGram(design, null), scale);
        double[,] lambdaFactor = Cholesky(precision);
        double logDetLambda = LogDet(lambdaFactor);
        var w = new double[m];
        for (int i = 0; i < m; i++) w[i] = Math.Exp(Row(design, mode, i));
        double[,] h = DenseGram(design, w);
        for (int p = 0; p < nu; p++) for (int q = 0; q < nu; q++) h[p, q] += precision[p, q];
        double[,] hFactor = Cholesky(h);
        double logDetH = LogDet(hFactor);
        double logT = SpecialFunctions.LogGamma(0.5 * (df + nu)) - SpecialFunctions.LogGamma(0.5 * df)
                      - 0.5 * nu * Math.Log(df * Math.PI) + 0.5 * logDetH;

        var logWeights = new double[draws];
        var z = new double[nu];
        var beta = new double[nu];
        for (int d = 0; d < draws; d++)
        {
            double zz = 0.0, chi = 0.0;
            for (int j = 0; j < nu; j++) { z[j] = Normal.SampleStandard(rng); zz += z[j] * z[j]; }
            for (int j = 0; j < df; j++) { double e = Normal.SampleStandard(rng); chi += e * e; }
            double[] offset = SolveUpper(hFactor, z);              // Lᵀ·offset = z, so offset ~ N(0, H⁻¹)
            double stretch = Math.Sqrt(df / chi);
            for (int j = 0; j < nu; j++) beta[j] = mode[j] + stretch * offset[j];
            double logProposal = logT - 0.5 * (df + nu) * Math.Log(1.0 + zz / chi);

            double logLik = 0.0, squares = 0.0;
            for (int i = 0; i < m; i++)
            {
                double eta = Row(design, beta, i);
                logLik += y[i] * eta - Math.Exp(eta);
                squares += (eta - level) * (eta - level);   // ‖B(β − β₀)‖², with Bβ₀ = log ȳ
            }
            double logPrior = -0.5 * nu * Math.Log(2.0 * Math.PI) + 0.5 * logDetLambda - 0.5 * scale * squares;
            logWeights[d] = logLik + logPrior - logProposal;
        }
        double max = logWeights.Max();
        double sum = 0.0, sumSq = 0.0;
        foreach (double lw in logWeights) { double r = Math.Exp(lw - max); sum += r; sumSq += r * r; }
        double meanRatio = sum / draws;
        double variance = sumSq / draws - meanRatio * meanRatio;
        return (max + Math.Log(meanRatio), Math.Sqrt(variance / draws) / meanRatio);
    }

    // Bᵀ(y − μ) − Λ(β − β₀) at β, with Λ = (ȳ/m)·BᵀB and Bβ₀ = log ȳ.
    private static double[] Gradient(double[,] design, double[] y, double[] beta)
    {
        int m = y.Length, nu = beta.Length;
        double mean = y.Average(), scale = mean / m, level = Math.Log(mean);
        var g = new double[nu];
        for (int i = 0; i < m; i++)
        {
            double eta = Row(design, beta, i);
            double r = y[i] - Math.Exp(eta) - scale * (eta - level);
            for (int j = 0; j < nu; j++) g[j] += design[i, j] * r;
        }
        return g;
    }

    private static double Row(double[,] design, double[] beta, int i)
    {
        double eta = 0.0;
        for (int j = 0; j < beta.Length; j++) eta += design[i, j] * beta[j];
        return eta;
    }

    private static double[,] DenseGram(double[,] z, double[]? w)
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

    private static double[,] Scaled(double[,] a, double s)
    {
        var r = (double[,])a.Clone();
        for (int p = 0; p < r.GetLength(0); p++) for (int q = 0; q < r.GetLength(1); q++) r[p, q] *= s;
        return r;
    }

    // The lower Cholesky factor L of a = LLᵀ.
    private static double[,] Cholesky(double[,] a)
    {
        int n = a.GetLength(0);
        var l = new double[n, n];
        for (int j = 0; j < n; j++)
        {
            double d = a[j, j];
            for (int k = 0; k < j; k++) d -= l[j, k] * l[j, k];
            l[j, j] = Math.Sqrt(d);
            for (int i = j + 1; i < n; i++)
            {
                double s = a[i, j];
                for (int k = 0; k < j; k++) s -= l[i, k] * l[j, k];
                l[i, j] = s / l[j, j];
            }
        }
        return l;
    }

    private static double LogDet(double[,] l)
    {
        double s = 0.0;
        for (int j = 0; j < l.GetLength(0); j++) s += 2.0 * Math.Log(l[j, j]);
        return s;
    }

    // Solves Lᵀx = b.
    private static double[] SolveUpper(double[,] l, double[] b)
    {
        int n = b.Length;
        var x = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            double s = b[i];
            for (int k = i + 1; k < n; k++) s -= l[k, i] * x[k];
            x[i] = s / l[i, i];
        }
        return x;
    }
}
