using System;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The smoothing-variance priors and prior scaling. Each prior's log density in ρ = ln τ² is its distribution's log density
/// at e^ρ plus the Jacobian, with the derivatives it reports; the tensor penalty's prior draws have covariance K⁺ and no
/// null-space component; a scaled Weibull prior's function draws exceed the bound at the stated rate on fresh draws, and
/// its scale is quadratic in the bound; and BAPS's anisotropic scheme under a scaled Weibull prior separates the axes.
/// </summary>
public sealed class SmoothingVariancePriorTests
{
    private readonly ITestOutputHelper _out;
    public SmoothingVariancePriorTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    public static TheoryData<int> Priors => new() { 0, 1, 2 };

    private static (ISmoothingVariancePrior Prior, Func<double, double> LogPdf) Prior(int which) => which switch
    {
        0 => (new InverseGammaVariancePrior(0.5, 0.2), InverseGamma.ShapeScale(0.5, 0.2).LogPdf),
        1 => (new WeibullVariancePrior(3.0), Weibull.ShapeScale(0.5, 3.0).LogPdf),
        _ => (new WeibullVariancePrior(0.4, shape: 2.0), Weibull.ShapeScale(2.0, 0.4).LogPdf),
    };

    [Theory]
    [MemberData(nameof(Priors))]
    public void LogDensityOfLog_IsTheDistributionPlusTheJacobian_WithItsDerivatives(int which)
    {
        (ISmoothingVariancePrior prior, Func<double, double> logPdf) = Prior(which);
        double reference = prior.LogDensityOfLog(0.0, out _, out _) - logPdf(1.0);
        foreach (double rho in new[] { -3.0, -0.5, 1.2, 2.5 })
        {
            double value = prior.LogDensityOfLog(rho, out double first, out double second);
            Assert.Equal(reference, value - (logPdf(Math.Exp(rho)) + rho), 9);
            const double e = 1e-5;
            double up = prior.LogDensityOfLog(rho + e, out double firstUp, out _), down = prior.LogDensityOfLog(rho - e, out double firstDown, out _);
            Assert.Equal((up - down) / (2 * e), first, 6);
            Assert.Equal((firstUp - firstDown) / (2 * e), second, 6);
        }
    }

    [Fact]
    public void TensorPriorDraws_HaveCovarianceKPlus_AndNoNullSpaceComponent()
    {
        const int nuX = 4, nuY = 5;
        var penalty = new TensorPenalty(nuX, nuY);
        int nu = nuX * nuY;
        double[] w = [0.7, 2.5];
        var band = new double[penalty.Bandwidth + 1, nu];
        penalty.AccumulateComponents(band, nu, w);
        var k = new double[nu, nu];
        for (int d = 0; d < band.GetLength(0); d++)
            for (int j = 0; j + d < nu; j++) { k[j + d, j] = band[d, j]; k[j, j + d] = band[d, j]; }
        EigenResult eig = DenseEigen.DecomposeSymmetric(k);
        double top = eig.Eigenvalues.Max();
        var plus = new double[nu, nu];   // K⁺
        for (int m = 0; m < nu; m++)
        {
            if (eig.Eigenvalues[m] < 1e-9 * top) continue;
            double[] q = eig.Eigenvectors[m];
            for (int a = 0; a < nu; a++)
                for (int b = 0; b < nu; b++) plus[a, b] += q[a] * q[b] / eig.Eigenvalues[m];
        }

        var rng = new Xoshiro256PlusPlus(seed: 5);
        var beta = new double[nu];
        var cov = new double[nu, nu];
        const int n = 200_000;
        double worstNull = 0.0;
        for (int s = 0; s < n; s++)
        {
            penalty.SamplePrior(w, rng, beta);
            for (int a = 0; a < nu; a++)
                for (int b = 0; b < nu; b++) cov[a, b] += beta[a] * beta[b] / n;
            for (int px = 0; px < 2; px++)   // the null space: (1, i) ⊗ (1, k) in the coefficient indices
                for (int py = 0; py < 2; py++)
                {
                    double dot = 0.0;
                    for (int i = 0; i < nuX; i++)
                        for (int kk = 0; kk < nuY; kk++) dot += Math.Pow(i, px) * Math.Pow(kk, py) * beta[i * nuY + kk];
                    worstNull = Math.Max(worstNull, Math.Abs(dot));
                }
        }
        Assert.True(worstNull < 1e-9, $"null-space component {worstNull}");
        double scale = Enumerable.Range(0, nu).Max(a => plus[a, a]);
        for (int a = 0; a < nu; a++)
            for (int b = 0; b < nu; b++) Assert.True(Math.Abs(cov[a, b] - plus[a, b]) < 0.02 * scale, $"Cov[{a},{b}] {cov[a, b]} against K⁺ {plus[a, b]}");
    }

    private static (double[,] Z, double[] Y, double[] F, TensorPenalty Penalty) Surface(int seed, double noise = 0.1)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        const int g = 18;
        var xs = new double[g * g];
        var ys = new double[g * g];
        for (int a = 0; a < g; a++)
            for (int b = 0; b < g; b++) { xs[a * g + b] = (a + 0.5) / g; ys[a * g + b] = (b + 0.5) / g; }
        double[] f = xs.Select((x, i) => Math.Sin(1.5 * Math.PI * x) * Math.Sin(5.0 * Math.PI * ys[i])).ToArray();
        double[] y = f.Select(v => v + noise * Gaussian(rng)).ToArray();
        var knots = new KnotConfig(Enumerable.Range(1, 6).Select(k => k / 7.0).ToArray());
        TensorDesign td = TensorProductDesign.Build(new SplineBasis(3), knots, xs, new SplineBasis(3), knots, ys);
        return (td.Design, y, f, new TensorPenalty(td.NuX, td.NuY));
    }

    [Fact]
    public void ScaledPrior_HoldsItsBoundAtTheStatedRate_AndScalesWithTheBoundSquared()
    {
        (double[,] z, _, _, TensorPenalty penalty) = Surface(seed: 1);
        WeibullVariancePrior prior = WeibullVariancePrior.ScaledTo(z, penalty, bound: 2.0, tailProbability: 0.05, draws: 4000, seed: 3);
        WeibullVariancePrior doubled = WeibullVariancePrior.ScaledTo(z, penalty, bound: 4.0, tailProbability: 0.05, draws: 4000, seed: 3);
        Assert.Equal(4.0 * prior.Scale, doubled.Scale, 9);

        // Fresh draws under the scaled prior: τⱼ² ~ W(½, λ), β ~ N(0, K(τ)⁺), f = Zβ.
        Weibull tau2 = Weibull.ShapeScale(0.5, prior.Scale);
        var rng = new Xoshiro256PlusPlus(seed: 99);
        int nu = z.GetLength(1), n = z.GetLength(0), within = 0;
        const int draws = 4000;
        var beta = new double[nu];
        for (int d = 0; d < draws; d++)
        {
            penalty.SamplePrior([1.0 / tau2.Sample(rng), 1.0 / tau2.Sample(rng)], rng, beta);
            double max = 0.0;
            for (int i = 0; i < n; i++)
            {
                double f = 0.0;
                for (int k = 0; k < nu; k++) f += z[i, k] * beta[k];
                max = Math.Max(max, Math.Abs(f));
            }
            if (max <= 2.0) within++;
        }
        double rate = within / (double)draws;
        _out.WriteLine($"[prior-scaling] λ = {prior.Scale:E3}; P(max|f| ≤ 2) on fresh draws {rate:F4}");
        Assert.InRange(rate, 0.95 - 0.015, 0.95 + 0.015);
    }

    [Fact]
    public void AnisotropicBaps_UnderAScaledWeibullPrior_SeparatesTheAxes()
    {
        (double[,] z, double[] y, double[] f, TensorPenalty penalty) = Surface(seed: 31);
        WeibullVariancePrior prior = WeibullVariancePrior.ScaledTo(z, penalty, bound: 3.0);
        BapsResult r = new Baps(z, y, penalty, componentPrior: prior)
            .Run(BapsScheme.AnisotropicGibbs, chains: 4, masterSeed: 2, burn: 1000, samples: 1500);
        PosteriorInterval lx = r.ComponentLambdaInterval(0), ly = r.ComponentLambdaInterval(1);
        double mse = 0.0;
        for (int i = 0; i < f.Length; i++)
        {
            double fit = 0.0;
            for (int p = 0; p < r.Coefficients.Length; p++) fit += z[i, p] * r.Coefficients[p];
            mse += (fit - f[i]) * (fit - f[i]) / f.Length;
        }
        _out.WriteLine($"[baps-weibull] prior λ {prior.Scale:E3}; λx {lx.Median:E2} [{lx.Lo:E2}, {lx.Hi:E2}] λy {ly.Median:E2} [{ly.Lo:E2}, {ly.Hi:E2}]; " +
                       $"MSE {mse:F5}; R̂ {string.Join(", ", r.ComponentLogLambdaConvergence!.Select(c => c.RHat.Value.ToString("F3")))}");
        Assert.True(lx.Lo > ly.Hi, "the 95% intervals separate, λx above λy");
        Assert.All(r.ComponentLogLambdaConvergence!, c => Assert.True(c.RHat.Value < 1.05, $"R̂ {c.RHat}"));
        Assert.True(mse < 0.0045, $"MSE {mse:F5}");   // the isotropic scheme's on this surface is 0.00448
    }
}
