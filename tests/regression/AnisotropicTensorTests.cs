using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Anisotropic tensor smoothing: on a surface that is gentle in x but oscillatory in y, a single (isotropic)
/// smoothing parameter must compromise, while the 2-D REML over (λ_x, λ_y) picks heavier x-smoothing than y and
/// recovers the surface better. Validates both the anisotropic REML (built on the 1-D penalty eigenvalues) and
/// that it detects the asymmetry (λ_x* &gt; λ_y*).
/// </summary>
public sealed class AnisotropicTensorTests
{
    private readonly ITestOutputHelper _out;
    public AnisotropicTensorTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Fact]
    public void Anisotropic_BeatsIsotropic_OnAsymmetricSurface()
    {
        var rng = new Xoshiro256PlusPlus(seed: 808);
        const int gx = 22, gy = 22;
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
                f[i] = Math.Sin(1.5 * Math.PI * xs[i]) * Math.Sin(6.0 * Math.PI * ys[i]);   // gentle x, fast y
                y[i] = f[i] + 0.10 * Gaussian(rng);
            }

        const int nKnots = 9;
        var knots = new double[nKnots];
        for (int k = 0; k < nKnots; k++) knots[k] = (k + 1.0) / (nKnots + 1);
        var basis = new SplineBasis(3);
        TensorDesign td = TensorProductDesign.Build(basis, new KnotConfig(knots), xs,
                                                    basis, new KnotConfig(knots), ys);
        int nu = td.NuX * td.NuY;

        double Mse(double[] beta)
        {
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

        static double Grid(int j, int g) => Math.Pow(10.0, -2.0 + 6.0 * j / (g - 1));   // 1e-2 … 1e4

        // Isotropic baseline: single-λ tensor penalty, 1-D REML.
        var iso = new PenalizedSpline(td.Design, y, new TensorPenalty(td.NuX, td.NuY, 2, 2));
        double isoLam = 1.0, isoBest = double.NegativeInfinity;
        for (int j = 0; j < 25; j++)
        {
            double lam = Grid(j, 25);
            double ev = iso.RemlLogEvidence(lam);
            if (ev > isoBest) { isoBest = ev; isoLam = lam; }
        }
        double isoMse = Mse(iso.Coefficients(isoLam));

        // Anisotropic: 2-D REML over (λ_x, λ_y).
        var aniso = new AnisotropicTensorSpline(td, y);
        double lamX = 1.0, lamY = 1.0, anisoBest = double.NegativeInfinity;
        const int g = 13;
        for (int ix = 0; ix < g; ix++)
            for (int iy = 0; iy < g; iy++)
            {
                double lx = Grid(ix, g), ly = Grid(iy, g);
                double ev = aniso.RemlLogEvidence(lx, ly);
                if (ev > anisoBest) { anisoBest = ev; lamX = lx; lamY = ly; }
            }
        double anisoMse = Mse(aniso.Coefficients(lamX, lamY));

        _out.WriteLine($"[aniso] iso: λ={isoLam:E1} mse={isoMse:F5} | aniso: λx={lamX:E1} λy={lamY:E1} " +
                       $"mse={anisoMse:F5} (ratio {anisoMse / isoMse:F2})");

        Assert.True(anisoMse < isoMse, $"anisotropic MSE {anisoMse:F5} should beat isotropic {isoMse:F5}");
        Assert.True(lamX > lamY, $"x is smoother ⇒ λx {lamX:E1} should exceed λy {lamY:E1}");
    }
}
