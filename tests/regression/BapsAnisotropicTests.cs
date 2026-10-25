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
/// BAPS's anisotropic scheme: one smoothing variance per penalty component, the log variances moved by a Newton-proposal
/// Metropolis step (Bach &amp; Klein 2025). The conditional reduces to the isotropic scheme's inverse gamma with one
/// component and has the derivatives it reports; the Metropolis block leaves that conditional invariant, against
/// numerical integration; on a surface gentle in x and fast in y the posterior puts λ_x above λ_y and fits better than
/// the isotropic scheme; a penalty without components is refused; and a resumed run continues exactly.
/// </summary>
public sealed class BapsAnisotropicTests
{
    private readonly ITestOutputHelper _out;
    public BapsAnisotropicTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    // A 1-D difference penalty as a one-component anisotropic penalty, for the reduction to the isotropic conditional.
    private sealed class OneComponent(int nu) : IAnisotropicPenalty
    {
        private readonly DifferencePenalty _d = new(2);
        public int Bandwidth => _d.Bandwidth;
        public int Nullity => _d.Nullity;
        public int Components => 1;
        public void AccumulateInto(double[,] band, int dim, double lambda) => _d.AccumulateInto(band, dim, lambda);
        public void AccumulateComponents(double[,] band, int dim, ReadOnlySpan<double> lambdas) => _d.AccumulateInto(band, dim, lambdas[0]);
        public double Roughness(double[] beta) => _d.Roughness(beta);
        public double ComponentRoughness(double[] beta, int component) => _d.Roughness(beta);
        public PenaltySpectrum Spectrum { get; } = new([DenseEigen.DecomposeSymmetric(new DifferencePenalty(2).ToDense(nu)).Eigenvalues
                                                         .OrderBy(e => e).Select((e, i) => i < 2 ? 0.0 : e).ToArray()]);
        public void SamplePrior(ReadOnlySpan<double> weights, Xoshiro256PlusPlus rng, Span<double> beta) => throw new NotSupportedException();
    }

    [Fact]
    public void Conditional_IsTheInverseGammaWithOneComponent_AndHasItsDerivatives()
    {
        const int nu = 12;
        var rng = new Xoshiro256PlusPlus(seed: 2);
        double[] beta = Enumerable.Range(0, nu).Select(_ => Gaussian(rng)).ToArray();
        const double a = 0.5, b = 0.2;

        var one = new ComponentSmoothing(new OneComponent(nu), new InverseGammaVariancePrior(a, b));
        one.Update(new double[1], beta, new Xoshiro256PlusPlus(1));   // reads the roughness
        double shape = a + (nu - 2) / 2.0, scale = b + new DifferencePenalty(2).Roughness(beta) / 2.0;
        double LogIg(double rho) => -shape * rho - scale * Math.Exp(-rho);   // IG(shape, scale) on τ², in ρ = ln τ²
        double reference = one.LogDensity([0.0], new double[1], new double[1]) - LogIg(0.0);
        foreach (double rho in new[] { -2.0, -0.3, 1.1, 3.0 })
            Assert.Equal(reference, one.LogDensity([rho], new double[1], new double[1]) - LogIg(rho), 9);

        var tensor = new TensorPenalty(4, 5);
        double[] surface = Enumerable.Range(0, 20).Select(_ => Gaussian(rng)).ToArray();
        var two = new ComponentSmoothing(tensor, new InverseGammaVariancePrior(a, b));
        two.Update(new double[2], surface, new Xoshiro256PlusPlus(1));
        double[] at = [0.4, -0.7];
        var g = new double[2];
        var h = new double[4];
        two.LogDensity(at, g, h);
        const double e = 1e-5;
        for (int j = 0; j < 2; j++)
        {
            double[] up = (double[])at.Clone(), down = (double[])at.Clone();
            up[j] += e;
            down[j] -= e;
            var gu = new double[2];
            var gd = new double[2];
            double fu = two.LogDensity(up, gu, new double[4]), fd = two.LogDensity(down, gd, new double[4]);
            Assert.Equal((fu - fd) / (2 * e), g[j], 6);
            for (int k = 0; k < 2; k++) Assert.Equal((gu[k] - gd[k]) / (2 * e), h[j * 2 + k], 6);
        }
    }

    // The conditional of a fixed β on a 0.01 grid: each cell's centre and unnormalized mass, with the running sum, so
    // exact draws can be taken from it.
    private static (double[] R0, double[] R1, double[] Cdf) Grid(ComponentSmoothing smoothing)
    {
        const double step = 0.01, lo0 = -4.0, lo1 = -5.0;
        const int m = 1400;
        var r0 = new double[m * m];
        var r1 = new double[m * m];
        var f = new double[m * m];
        double peak = double.NegativeInfinity;
        for (int a = 0; a < m; a++)
            for (int b = 0; b < m; b++)
            {
                int c = a * m + b;
                r0[c] = lo0 + (a + 0.5) * step;
                r1[c] = lo1 + (b + 0.5) * step;
                f[c] = smoothing.LogDensity([r0[c], r1[c]]);
                peak = Math.Max(peak, f[c]);
            }
        var cdf = new double[m * m];
        double total = 0.0;
        for (int c = 0; c < m * m; c++) cdf[c] = total += Math.Exp(f[c] - peak);
        return (r0, r1, cdf);
    }

    private static (ComponentSmoothing Smoothing, double[] Beta) FixedConditional(
        ISmoothingVariancePrior? prior = null, SmoothingVarianceUpdate update = SmoothingVarianceUpdate.Newton)
    {
        var rng = new Xoshiro256PlusPlus(seed: 7);
        double[] beta = Enumerable.Range(0, 30).Select(i => Gaussian(rng) * (i % 2 == 0 ? 0.3 : 1.0)).ToArray();
        var smoothing = new ComponentSmoothing(new TensorPenalty(5, 6), prior ?? new InverseGammaVariancePrior(), update);
        smoothing.Update(new double[2], beta, new Xoshiro256PlusPlus(1));   // reads the roughness
        return (smoothing, beta);
    }

    [Theory]
    [InlineData(false, SmoothingVarianceUpdate.Newton)]
    [InlineData(true, SmoothingVarianceUpdate.Newton)]
    [InlineData(false, SmoothingVarianceUpdate.Slice)]
    [InlineData(true, SmoothingVarianceUpdate.Slice)]
    public void Update_FromExactDraws_StaysOnTheConditional(bool weibull, SmoothingVarianceUpdate update)
    {
        (ComponentSmoothing smoothing, double[] beta) = FixedConditional(weibull ? new WeibullVariancePrior(2.0) : null, update);
        (double[] r0, double[] r1, double[] cdf) = Grid(smoothing);
        var rng = new Xoshiro256PlusPlus(seed: 9);
        const int n = 500_000;
        double b0 = 0, b00 = 0, a0 = 0, a00 = 0, a1 = 0;
        for (int i = 0; i < n; i++)
        {
            int c = Array.BinarySearch(cdf, rng.NextDouble() * cdf[^1]);
            if (c < 0) c = ~c;
            double[] rho = [r0[c] + 0.01 * (rng.NextDouble() - 0.5), r1[c] + 0.01 * (rng.NextDouble() - 0.5)];
            b0 += rho[0]; b00 += rho[0] * rho[0];
            smoothing.Update(rho, beta, rng);
            a0 += rho[0]; a00 += rho[0] * rho[0]; a1 += rho[1];
        }
        double mb = b0 / n, ma = a0 / n, vb = b00 / n - mb * mb, va = a00 / n - ma * ma;
        _out.WriteLine($"[aniso-block] {(weibull ? "Weibull" : "IG")} {update} exact draws: mean {mb:F5} var {vb:F5}; after one update: mean {ma:F5} var {va:F5}");
        double se = Math.Sqrt(vb / n);
        Assert.InRange(ma - mb, -5.0 * se, 5.0 * se);
        Assert.InRange(va / vb, 0.98, 1.02);   // the upper tail is exponential, so the variance estimate is noisy: ±2% is about 4σ here
    }

    [Theory]
    [InlineData(SmoothingVarianceUpdate.Newton)]
    [InlineData(SmoothingVarianceUpdate.Slice)]
    public void Chain_VisitsTheConditionalsUpperTail(SmoothingVarianceUpdate update)
    {
        (ComponentSmoothing smoothing, double[] beta) = FixedConditional(update: update);
        (double[] r0, _, double[] cdf) = Grid(smoothing);
        double tail = 0.0;
        for (int c = 0; c < cdf.Length; c++)
            if (r0[c] > 3.5) tail += cdf[c] - (c == 0 ? 0.0 : cdf[c - 1]);
        tail /= cdf[^1];

        long n = 0, inTail = 0;
        int longest = 0;
        for (int chain = 0; chain < 4; chain++)
        {
            var rng = new Xoshiro256PlusPlus(seed: 300 + chain);
            double[] rho = [1.95, 0.8];
            int run = 0;
            for (int i = 0; i < 250_000; i++)
            {
                double before = rho[0] + rho[1];
                smoothing.Update(rho, beta, rng);
                run = rho[0] + rho[1] == before ? run + 1 : 0;
                longest = Math.Max(longest, run);
                if (rho[0] > 3.5) inTail++;
                n++;
            }
        }
        double visited = inTail / (double)n;
        _out.WriteLine($"[aniso-block] {update} P(ρ0 > 3.5): chain {visited:E3}, conditional {tail:E3}; longest stay {longest}");
        Assert.InRange(visited / tail, 0.85, 1.15);
        Assert.True(longest < 100, $"the chain stayed {longest} updates at one point");
    }

    private static (TensorDesign Td, double[] Y, double[] F) Surface(int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        const int g = 18;
        var xs = new double[g * g];
        var ys = new double[g * g];
        var f = new double[g * g];
        var y = new double[g * g];
        for (int a = 0; a < g; a++)
            for (int b = 0; b < g; b++)
            {
                int i = a * g + b;
                xs[i] = (a + 0.5) / g;
                ys[i] = (b + 0.5) / g;
                f[i] = Math.Sin(1.5 * Math.PI * xs[i]) * Math.Sin(5.0 * Math.PI * ys[i]);   // gentle x, fast y
                y[i] = f[i] + 0.1 * Gaussian(rng);
            }
        var knots = new KnotConfig(Enumerable.Range(1, 6).Select(k => k / 7.0).ToArray());
        return (TensorProductDesign.Build(new SplineBasis(3), knots, xs, new SplineBasis(3), knots, ys), y, f);
    }

    private static double Mse(double[,] z, double[] beta, double[] f)
    {
        double s = 0.0;
        for (int i = 0; i < f.Length; i++)
        {
            double fit = 0.0;
            for (int p = 0; p < beta.Length; p++) fit += z[i, p] * beta[p];
            s += (fit - f[i]) * (fit - f[i]);
        }
        return s / f.Length;
    }

    [Fact]
    public void Anisotropic_SmoothsXMoreThanY_AndFitsBetterThanIsotropic()
    {
        (TensorDesign td, double[] y, double[] f) = Surface(seed: 31);
        var penalty = new TensorPenalty(td.NuX, td.NuY);
        var baps = new Baps(td.Design, y, penalty);
        BapsResult aniso = baps.Run(BapsScheme.AnisotropicGibbs, chains: 4, masterSeed: 2, burn: 1000, samples: 1500);
        BapsResult iso = baps.Run(BapsScheme.Gibbs, chains: 4, masterSeed: 2, burn: 1000, samples: 1500);

        PosteriorInterval lx = aniso.ComponentLambdaInterval(0), ly = aniso.ComponentLambdaInterval(1);
        AnisotropicReml reml = new AnisotropicTensorSpline(td, y).RemlOptimum();
        double anisoMse = Mse(td.Design, aniso.Coefficients, f), isoMse = Mse(td.Design, iso.Coefficients, f);
        _out.WriteLine($"[baps-aniso] λx {lx.Median:E2} [{lx.Lo:E2}, {lx.Hi:E2}] λy {ly.Median:E2} [{ly.Lo:E2}, {ly.Hi:E2}]; " +
                       $"REML ({reml.LambdaX:E2}, {reml.LambdaY:E2}); MSE aniso {anisoMse:F5} iso {isoMse:F5}; " +
                       $"acceptance {aniso.AcceptanceRate:F3}; R̂ {string.Join(", ", aniso.ComponentLogLambdaConvergence!.Select(c => c.RHat.Value.ToString("F3")))}");

        Assert.True(lx.Lo > ly.Hi, "the 95% intervals separate, λx above λy");
        Assert.True(anisoMse < isoMse, $"anisotropic MSE {anisoMse:F5} against isotropic {isoMse:F5}");
        Assert.Equal(2, aniso.ComponentLogLambdaConvergence!.Count);
        Assert.All(aniso.ComponentLogLambdaConvergence, c => Assert.True(c.RHat.Value < 1.05, $"R̂ {c.RHat}"));
        foreach (BapsDraw d in aniso.PosteriorDraws[0])
            Assert.Equal(d.Sigma2 / Math.Sqrt(d.ComponentTau2![0] * d.ComponentTau2[1]), d.Lambda, 10);   // the geometric mean
        Assert.Null(iso.ComponentLogLambdaConvergence);
        Assert.Throws<InvalidOperationException>(() => iso.ComponentLambdaInterval(0));
    }

    [Fact]
    public void SliceUpdate_AgreesWithTheNewtonUpdate()
    {
        (TensorDesign td, double[] y, _) = Surface(seed: 31);
        var penalty = new TensorPenalty(td.NuX, td.NuY);
        BapsResult newton = new Baps(td.Design, y, penalty).Run(BapsScheme.AnisotropicGibbs, chains: 4, masterSeed: 2, burn: 1000, samples: 1500);
        BapsResult slice = new Baps(td.Design, y, penalty, componentUpdate: SmoothingVarianceUpdate.Slice)
            .Run(BapsScheme.AnisotropicGibbs, chains: 4, masterSeed: 2, burn: 1000, samples: 1500);
        for (int j = 0; j < 2; j++)
        {
            double a = newton.ComponentLambdaInterval(j).Median, b = slice.ComponentLambdaInterval(j).Median;
            _out.WriteLine($"[baps-aniso] λ{j}: Newton {a:E3}, slice {b:E3}");
            Assert.True(Math.Abs(Math.Log(a / b)) < 0.1, $"component {j}: medians {a:E3} and {b:E3}");
        }
        Assert.Equal(1.0, slice.AcceptanceRate);
    }

    [Fact]
    public void Anisotropic_RefusesAPenaltyWithoutComponents_AndResumesExactly()
    {
        var rng = new Xoshiro256PlusPlus(seed: 3);
        double[] x = Enumerable.Range(0, 40).Select(i => (i + 0.5) / 40).ToArray();
        double[] y1 = x.Select(v => Math.Sin(6 * v) + 0.1 * Gaussian(rng)).ToArray();
        double[,] z1 = new SplineBasis(3).Design(new KnotConfig([0.25, 0.5, 0.75]), x);
        Assert.Throws<ArgumentException>(() => new Baps(z1, y1, new DifferencePenalty(2))
            .Run(BapsScheme.AnisotropicGibbs, chains: 1, burn: 1, samples: 2));

        (TensorDesign td, double[] y, _) = Surface(seed: 32);
        var baps = new Baps(td.Design, y, new TensorPenalty(td.NuX, td.NuY));
        BapsResult straight = baps.Run(BapsScheme.AnisotropicGibbs, chains: 2, masterSeed: 5, burn: 100, samples: 160);
        BapsResult first = baps.Run(BapsScheme.AnisotropicGibbs, chains: 2, masterSeed: 5, burn: 100, samples: 60, checkpoint: true);
        BapsResult resumed = baps.Run(BapsScheme.AnisotropicGibbs, chains: 2, masterSeed: 5, burn: 100, samples: 100,
                                      resumeFrom: first.Checkpoint);
        Assert.Equal(straight.Coefficients, resumed.Coefficients);
        for (int c = 0; c < 2; c++)
            for (int i = 0; i < straight.PosteriorDraws[c].Count; i++)
            {
                Assert.Equal(straight.PosteriorDraws[c][i].ComponentTau2, resumed.PosteriorDraws[c][i].ComponentTau2);
                Assert.Equal(straight.PosteriorDraws[c][i].Sigma2, resumed.PosteriorDraws[c][i].Sigma2);
            }
    }
}
