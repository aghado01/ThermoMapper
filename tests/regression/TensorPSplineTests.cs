using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The tensor-product P-spline (He, Yang &amp; Kang) smooths a bivariate surface on the same banded machinery as
/// the 1-D case: the flattened tensor design plus the tensor difference penalty, with the smoothing chosen by
/// REML. On a noisy 2-D surface the REML evidence must peak at a finite λ and the fit there must recover the
/// truth better than both the unpenalized overfit and the over-smoothed limit.
/// </summary>
public sealed class TensorPSplineTests
{
    private readonly ITestOutputHelper _out;
    public TensorPSplineTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void TensorPSpline_RecoversSurface_ViaReml()
    {
        var rng = new Xoshiro256PlusPlus(seed: 202);
        const int gx = 18, gy = 18;
        int n = gx * gy;
        var xs = new double[n];
        var ys = new double[n];
        var f = new double[n];
        var y = new double[n];
        for (int a = 0; a < gx; a++)
            for (int b = 0; b < gy; b++)
            {
                int i = a * gy + b;
                xs[i] = (a + 0.5) / gx;
                ys[i] = (b + 0.5) / gy;
                f[i] = Math.Sin(2.0 * Math.PI * xs[i]) * Math.Sin(2.0 * Math.PI * ys[i]);
                y[i] = f[i] + 0.10 * Gaussian(rng);
            }

        const int nKnots = 8;
        var knots = new double[nKnots];
        for (int k = 0; k < nKnots; k++) knots[k] = (k + 1.0) / (nKnots + 1);
        var basis = new SplineBasis(3);
        TensorDesign td = TensorProductDesign.Build(basis, new KnotConfig(knots), xs,
                                                    basis, new KnotConfig(knots), ys);
        int nu = td.NuX * td.NuY;
        var penalty = new TensorPenalty(td.NuX, td.NuY, 2, 2);
        var ps = new PenalizedSpline(td.Design, y, penalty);

        double Mse(double lambda)
        {
            double[] beta = ps.Coefficients(lambda);
            double s = 0.0;
            for (int i = 0; i < n; i++)
            {
                double fit = 0.0;
                for (int p = 0; p < nu; p++) fit += td.Design[i, p] * beta[p];
                double d = fit - f[i];
                s += d * d;
            }
            return s / n;
        }

        const int g = 25;
        var lambdas = new double[g];
        int best = 0;
        var ev = new double[g];
        for (int j = 0; j < g; j++)
        {
            lambdas[j] = Math.Pow(10.0, -3.0 + 7.0 * j / (g - 1));   // 1e-3 … 1e4
            ev[j] = ps.RemlLogEvidence(lambdas[j]);
            if (ev[j] > ev[best]) best = j;
        }

        double mseStar = Mse(lambdas[best]);
        double mseUnder = Mse(lambdas[0]);
        double mseOver = Mse(lambdas[g - 1]);
        _out.WriteLine($"[tensor] ν={nu} (={td.NuX}×{td.NuY}) λ*={lambdas[best]:E2} " +
                       $"mse λ*={mseStar:F5} under={mseUnder:F5} over={mseOver:F5}");

        Assert.True(best > 0 && best < g - 1, $"REML optimum at grid edge (index {best}).");
        Assert.True(mseStar < mseUnder, $"REML fit {mseStar:F5} should beat the overfit {mseUnder:F5}");
        Assert.True(mseStar < mseOver, $"REML fit {mseStar:F5} should beat the over-smoothed {mseOver:F5}");
        Assert.True(mseStar < 0.01, $"REML fit {mseStar:F5} should recover the surface");
    }
}
