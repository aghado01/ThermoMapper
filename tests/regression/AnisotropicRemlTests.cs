using System;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The anisotropic penalty seam and the REML optimizer. The tensor penalty's components assemble to the explicit
/// Kronecker form and its isotropic case; the joint spectrum's log pseudo-determinant matches a dense eigensolve of
/// λ_x P_x + λ_y P_y, and its derivatives match finite differences; Nelder–Mead finds a known minimum and stops at a box
/// edge; and the REML optimum beats a grid around it, detects the asymmetry of a surface gentle in x and fast in y, and
/// flags an axis along which the surface has no curvature.
/// </summary>
public sealed class AnisotropicRemlTests
{
    private readonly ITestOutputHelper _out;
    public AnisotropicRemlTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    // P_x ⊗ I and I ⊗ P_y densely, coefficients flattened as i·ν_y + k.
    private static (double[,] Px, double[,] Py) Kronecker(int nuX, int nuY, int orderX, int orderY)
    {
        double[,] dx = new DifferencePenalty(orderX).ToDense(nuX), dy = new DifferencePenalty(orderY).ToDense(nuY);
        int nu = nuX * nuY;
        var px = new double[nu, nu];
        var py = new double[nu, nu];
        for (int i = 0; i < nuX; i++)
            for (int k = 0; k < nuY; k++)
                for (int i2 = 0; i2 < nuX; i2++)
                    for (int k2 = 0; k2 < nuY; k2++)
                    {
                        if (k == k2) px[i * nuY + k, i2 * nuY + k2] = dx[i, i2];
                        if (i == i2) py[i * nuY + k, i2 * nuY + k2] = dy[k, k2];
                    }
        return (px, py);
    }

    private static double[,] Dense(double[,] band, int nu)
    {
        var m = new double[nu, nu];
        for (int d = 0; d < band.GetLength(0); d++)
            for (int j = 0; j + d < nu; j++) { m[j + d, j] = band[d, j]; m[j, j + d] = band[d, j]; }
        return m;
    }

    [Fact]
    public void TensorComponents_AssembleTheKroneckerForm_AndItsIsotropicCase()
    {
        const int nuX = 5, nuY = 6;
        var penalty = new TensorPenalty(nuX, nuY, orderX: 2, orderY: 1);
        (double[,] px, double[,] py) = Kronecker(nuX, nuY, 2, 1);
        int nu = nuX * nuY;
        var band = new double[penalty.Bandwidth + 1, nu];
        penalty.AccumulateComponents(band, nu, [0.7, 3.0]);
        double[,] m = Dense(band, nu);
        for (int a = 0; a < nu; a++)
            for (int b = 0; b < nu; b++) Assert.Equal(0.7 * px[a, b] + 3.0 * py[a, b], m[a, b], 12);

        var iso = new double[penalty.Bandwidth + 1, nu];
        var both = new double[penalty.Bandwidth + 1, nu];
        penalty.AccumulateInto(iso, nu, 1.9);
        penalty.AccumulateComponents(both, nu, [1.9, 1.9]);
        Assert.Equal(iso.Cast<double>(), both.Cast<double>());

        var rng = new Xoshiro256PlusPlus(seed: 3);
        double[] beta = Enumerable.Range(0, nu).Select(_ => Gaussian(rng)).ToArray();
        double qx = 0.0, qy = 0.0;
        for (int a = 0; a < nu; a++)
            for (int b = 0; b < nu; b++) { qx += beta[a] * px[a, b] * beta[b]; qy += beta[a] * py[a, b] * beta[b]; }
        Assert.Equal(qx, penalty.ComponentRoughness(beta, 0), 10);
        Assert.Equal(qy, penalty.ComponentRoughness(beta, 1), 10);
        Assert.Equal(penalty.Roughness(beta), qx + qy, 10);
    }

    [Fact]
    public void Spectrum_MatchesADenseEigensolve_AndItsDerivativesMatchFiniteDifferences()
    {
        const int nuX = 5, nuY = 6;
        var penalty = new TensorPenalty(nuX, nuY, orderX: 2, orderY: 1);
        (double[,] px, double[,] py) = Kronecker(nuX, nuY, 2, 1);
        int nu = nuX * nuY;
        PenaltySpectrum spectrum = penalty.Spectrum;
        Assert.Equal(2, spectrum.Nullity);   // degree < 2 in x times degree < 1 in y
        Assert.Equal(nu - 2, spectrum.Rank);

        foreach ((double lx, double ly) in new[] { (0.3, 2.0), (5.0, 0.01), (1.0, 1.0) })
        {
            var m = new double[nu, nu];
            for (int a = 0; a < nu; a++)
                for (int b = 0; b < nu; b++) m[a, b] = lx * px[a, b] + ly * py[a, b];
            double[] eig = DenseEigen.DecomposeSymmetric(m).Eigenvalues.OrderBy(e => e).ToArray();
            double logDet = eig.Skip(2).Sum(Math.Log);
            Assert.Equal(logDet, spectrum.LogPseudoDeterminant([lx, ly]), 9);

            var g = new double[2];
            var h = new double[4];
            double value = spectrum.LogPseudoDeterminant([lx, ly], g, h);
            Assert.Equal(spectrum.LogPseudoDeterminant([lx, ly]), value, 12);
            const double e = 1e-5;
            for (int j = 0; j < 2; j++)
            {
                double[] up = [lx, ly], down = [lx, ly];
                up[j] *= Math.Exp(e);
                down[j] *= Math.Exp(-e);
                Assert.Equal((spectrum.LogPseudoDeterminant(up) - spectrum.LogPseudoDeterminant(down)) / (2 * e), g[j], 6);
                var gu = new double[2];
                var gd = new double[2];
                spectrum.LogPseudoDeterminant(up, gu, new double[4]);
                spectrum.LogPseudoDeterminant(down, gd, new double[4]);
                for (int i = 0; i < 2; i++) Assert.Equal((gu[i] - gd[i]) / (2 * e), h[j * 2 + i], 6);
            }
            Assert.True(h[0] >= 0.0 && h[3] >= 0.0 && h[0] * h[3] - h[1] * h[2] >= -1e-12, "convex in ln λ");
        }
    }

    [Fact]
    public void NelderMead_FindsAMinimum_AndStopsAtABoxEdge()
    {
        static double Rosenbrock(double[] u) => 100.0 * Math.Pow(u[1] - u[0] * u[0], 2) + Math.Pow(1.0 - u[0], 2);
        NelderMead.Result r = NelderMead.Minimize(Rosenbrock, [-1.2, 1.0], [0.5, 0.5], [-5.0, -5.0], [5.0, 5.0], 1e-10, 2000);
        Assert.True(r.Converged);
        Assert.Equal(1.0, r.Point[0], 4);
        Assert.Equal(1.0, r.Point[1], 4);

        // (u − 3)² + (v + 1)² with u ≤ 2: the minimum within the box is (2, −1), on the edge.
        NelderMead.Result edge = NelderMead.Minimize(u => Math.Pow(u[0] - 3.0, 2) + Math.Pow(u[1] + 1.0, 2),
                                                     [0.0, 0.0], [1.0, 1.0], [-5.0, -5.0], [2.0, 5.0], 1e-10, 2000);
        Assert.True(edge.Converged);
        Assert.Equal(2.0, edge.Point[0], 12);
        Assert.Equal(-1.0, edge.Point[1], 4);
    }

    private static (TensorDesign Td, double[] Y) Surface(Func<double, double, double> f, int seed, double noise)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        const int g = 20;
        var xs = new double[g * g];
        var ys = new double[g * g];
        var y = new double[g * g];
        for (int a = 0; a < g; a++)
            for (int b = 0; b < g; b++)
            {
                int i = a * g + b;
                xs[i] = (a + 0.5) / g;
                ys[i] = (b + 0.5) / g;
                y[i] = f(xs[i], ys[i]) + noise * Gaussian(rng);
            }
        var knots = new KnotConfig(Enumerable.Range(1, 7).Select(k => k / 8.0).ToArray());
        TensorDesign td = TensorProductDesign.Build(new SplineBasis(3), knots, xs, new SplineBasis(3), knots, ys);
        return (td, y);
    }

    [Fact]
    public void RemlOptimum_BeatsAGridAroundIt_AndSeesTheAsymmetry()
    {
        (TensorDesign td, double[] y) = Surface((x, v) => Math.Sin(1.5 * Math.PI * x) * Math.Sin(5.0 * Math.PI * v), seed: 11, noise: 0.1);
        var spline = new AnisotropicTensorSpline(td, y);
        AnisotropicReml fit = spline.RemlOptimum();
        _out.WriteLine($"[aniso-reml] λx={fit.LambdaX:E3} λy={fit.LambdaY:E3} ℓ={fit.LogEvidence:F6} after {fit.Evaluations} evaluations");

        Assert.True(fit.Converged);
        Assert.False(fit.LambdaXAtBound || fit.LambdaYAtBound);
        Assert.True(fit.LambdaX > fit.LambdaY, "x is smoother, so λx should exceed λy");
        Assert.Equal(spline.RemlLogEvidence(fit.LambdaX, fit.LambdaY), fit.LogEvidence, 12);

        double best = double.NegativeInfinity;
        for (int a = -12; a <= 12; a++)
            for (int b = -12; b <= 12; b++)
                best = Math.Max(best, spline.RemlLogEvidence(fit.LambdaX * Math.Exp(0.25 * a), fit.LambdaY * Math.Exp(0.25 * b)));
        Assert.True(fit.LogEvidence >= best - 1e-9 * Math.Abs(best), $"optimum {fit.LogEvidence} against the grid's {best}");
    }

    [Fact]
    public void RemlOptimum_FlagsAnAxisWithoutCurvature()
    {
        // Constant in x: the x-penalty's null space holds the truth, so λx runs to the box edge or far above λy.
        (TensorDesign td, double[] y) = Surface((x, v) => Math.Sin(4.0 * Math.PI * v), seed: 12, noise: 0.1);
        AnisotropicReml fit = new AnisotropicTensorSpline(td, y).RemlOptimum();
        _out.WriteLine($"[aniso-reml] flat x: λx={fit.LambdaX:E3} (at bound {fit.LambdaXAtBound}) λy={fit.LambdaY:E3}");
        Assert.True(fit.LambdaXAtBound || fit.LambdaX > 1e3 * fit.LambdaY, $"λx {fit.LambdaX:E2} against λy {fit.LambdaY:E2}");
        Assert.False(fit.LambdaYAtBound);
    }
}
