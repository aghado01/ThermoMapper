using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The hierarchical P-spline's Gibbs scan, checked against closed forms on a tiny model (ν = 5, two replicates, eight
/// points). The full conditionals are read through the public API: a one-sweep run from the known initial state (the
/// fixed start, <see cref="HierarchicalStart.Fixed"/>)
/// exposes the β_r and β_0 draws, and the end state of any sweep exposes everything the three variance draws condition
/// on. The stationary posterior means are checked against quadrature over the variance components with every
/// coefficient integrated out exactly, a reference that shares no code with the sampler.
/// </summary>
public sealed class HierarchicalPSplineConditionalTests
{
    private const int N = 8, Reps = 2, Order = 2;
    private const double PriorShape = 2.0, PriorScale = 0.2;
    private const int Draws = 4000;

    private readonly ITestOutputHelper _out;
    public HierarchicalPSplineConditionalTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    // Two replicates of an offset sine (offset keeps β far from the zero initial state), noise sd 0.3.
    private static (double[,] Z, double[][] Y) TinyModel()
    {
        var x = new double[N];
        for (int i = 0; i < N; i++) x[i] = (i + 0.5) / N;
        double[,] z = new SplineBasis(3).Design(new KnotConfig(new[] { 0.5 }), x);
        var rng = new Xoshiro256PlusPlus(seed: 303);
        var y = new double[Reps][];
        for (int r = 0; r < Reps; r++)
        {
            y[r] = new double[N];
            for (int i = 0; i < N; i++)
                y[r][i] = 1.5 + Math.Sin(2.0 * Math.PI * x[i]) + 0.3 * r * Math.Cos(Math.PI * x[i]) + 0.3 * Gaussian(rng);
        }
        return (z, y);
    }

    private static HierarchicalPSpline Sampler(double[,] z, double[][] y)
        => new(z, y, new DifferencePenalty(Order), PriorShape, PriorScale);

    [Fact]
    public void FirstSweep_CoefficientDraws_MatchTheirFullConditionals()
    {
        (double[,] z, double[][] y) = TinyModel();
        int nu = z.GetLength(1);
        double[,] ztz = Gram(z);
        double[,] p = PenaltyMatrix(nu);
        HierarchicalPSpline sampler = Sampler(z, y);

        // Initial state: β_r = β_0 = 0 and σ² = τ_u² = τ_0² = 1.
        // β_r | · ~ N(Q_r⁻¹ Zᵀy_r, Q_r⁻¹) with Q_r = ZᵀZ + I;  β_0 | β ~ N(Q_0⁻¹ Σ_r β_r, Q_0⁻¹) with Q_0 = R·I + P.
        var qr = new double[nu, nu];
        var q0 = new double[nu, nu];
        for (int i = 0; i < nu; i++)
            for (int j = 0; j < nu; j++)
            {
                qr[i, j] = ztz[i, j] + (i == j ? 1.0 : 0.0);
                q0[i, j] = p[i, j] + (i == j ? Reps : 0.0);
            }
        double[,] lr = Cholesky(qr), l0 = Cholesky(q0);
        var mr = new double[Reps][];
        for (int r = 0; r < Reps; r++) mr[r] = Solve(lr, MulTranspose(z, y[r]));

        int[] seeds = SeedTree.Derive(31, Draws);
        var wr = new double[Reps][][];
        for (int r = 0; r < Reps; r++) wr[r] = new double[Draws][];
        var w0 = new double[Draws][];
        for (int s = 0; s < Draws; s++)
        {
            HierarchicalResult h = sampler.Run(chains: 1, masterSeed: seeds[s], burn: 0, samples: 1, start: HierarchicalStart.Fixed);
            var sum = new double[nu];
            for (int r = 0; r < Reps; r++)
            {
                double[] br = h.ReplicateCoefficients[r];
                var dev = new double[nu];
                for (int j = 0; j < nu; j++) { dev[j] = br[j] - mr[r][j]; sum[j] += br[j]; }
                wr[r][s] = LowerTransposeTimes(lr, dev);   // Lᵀ(β − m) ~ N(0, I) when Cov = (LLᵀ)⁻¹
            }
            double[] m0 = Solve(l0, sum);
            var dev0 = new double[nu];
            for (int j = 0; j < nu; j++) dev0[j] = h.PopulationCoefficients[j] - m0[j];
            w0[s] = LowerTransposeTimes(l0, dev0);
        }

        for (int r = 0; r < Reps; r++) AssertStandardNormal(wr[r], $"β_r (replicate {r + 1}) | β_0, σ², τ_u²");
        AssertStandardNormal(w0, "β_0 | β_r, τ_u², τ_0²");
    }

    [Theory]
    [InlineData(0)]   // the first sweep, from the known initial state
    [InlineData(4)]   // a later sweep: the draws condition only on the end state, which the run reports
    public void VarianceDraws_MatchTheirFullConditionals(int burn)
    {
        (double[,] z, double[][] y) = TinyModel();
        int nu = z.GetLength(1);
        double[,] p = PenaltyMatrix(nu);
        HierarchicalPSpline sampler = Sampler(z, y);

        // Each variance v ~ IG(α, b + S/2) given the current coefficients, so (b + S/2)/v ~ Gamma(α, 1).
        int[] seeds = SeedTree.Derive(47 + burn, Draws);
        var gSigma = new double[Draws];
        var gDev = new double[Draws];
        var gPop = new double[Draws];
        for (int s = 0; s < Draws; s++)
        {
            HierarchicalResult h = sampler.Run(chains: 1, masterSeed: seeds[s], burn: burn, samples: 1, start: HierarchicalStart.Fixed);
            double[] b0 = h.PopulationCoefficients;
            double resid = 0.0, dev = 0.0;
            for (int r = 0; r < Reps; r++)
            {
                double[] br = h.ReplicateCoefficients[r];
                for (int i = 0; i < N; i++)
                {
                    double fit = 0.0;
                    for (int j = 0; j < nu; j++) fit += z[i, j] * br[j];
                    resid += (y[r][i] - fit) * (y[r][i] - fit);
                }
                for (int j = 0; j < nu; j++) dev += (br[j] - b0[j]) * (br[j] - b0[j]);
            }
            double rough = Quadratic(p, b0);
            gSigma[s] = (PriorScale + 0.5 * resid) / (h.NoiseSd * h.NoiseSd);
            gDev[s] = (PriorScale + 0.5 * dev) / (h.DeviationSd * h.DeviationSd);
            gPop[s] = (PriorScale + 0.5 * rough) / (h.PopulationSmoothingSd * h.PopulationSmoothingSd);
        }

        AssertUnitGamma(gSigma, PriorShape + 0.5 * Reps * N, $"σ² | β_r (burn {burn})");
        AssertUnitGamma(gDev, PriorShape + 0.5 * Reps * nu, $"τ_u² | β_r, β_0 (burn {burn})");
        AssertUnitGamma(gPop, PriorShape + 0.5 * (nu - Order), $"τ_0² | β_0 (burn {burn})");
    }

    [Fact]
    public void StationaryMeans_MatchQuadratureReference()
    {
        (double[,] z, double[][] y) = TinyModel();
        int nu = z.GetLength(1);
        (double sigma, double tauU, double tau0, double[] beta0) reference = QuadratureMeans(z, y);
        HierarchicalPSpline sampler = Sampler(z, y);

        // Independent chains give an honest Monte Carlo standard error for each reported mean.
        const int chains = 16;
        int[] seeds = SeedTree.Derive(59, chains);
        int dim = 3 + nu;
        var est = new double[chains][];
        for (int c = 0; c < chains; c++)
        {
            HierarchicalResult h = sampler.Run(chains: 1, masterSeed: seeds[c], burn: 2000, samples: 20000);
            est[c] = new double[dim];
            est[c][0] = h.NoiseSd;
            est[c][1] = h.DeviationSd;
            est[c][2] = h.PopulationSmoothingSd;
            for (int j = 0; j < nu; j++) est[c][3 + j] = h.PopulationCoefficients[j];
        }

        var refs = new double[dim];
        refs[0] = reference.sigma;
        refs[1] = reference.tauU;
        refs[2] = reference.tau0;
        for (int j = 0; j < nu; j++) refs[3 + j] = reference.beta0[j];
        string[] names = new string[dim];
        names[0] = "E[σ]"; names[1] = "E[τ_u]"; names[2] = "E[τ_0]";
        for (int j = 0; j < nu; j++) names[3 + j] = $"E[β_0,{j}]";

        bool ok = true;
        for (int k = 0; k < dim; k++)
        {
            double mean = 0.0;
            for (int c = 0; c < chains; c++) mean += est[c][k];
            mean /= chains;
            double ss = 0.0;
            for (int c = 0; c < chains; c++) ss += (est[c][k] - mean) * (est[c][k] - mean);
            double se = Math.Sqrt(ss / (chains - 1) / chains);
            double zScore = (mean - refs[k]) / se;
            bool pass = Math.Abs(mean - refs[k]) <= 4.5 * se + 1e-4 * Math.Abs(refs[k]);
            ok &= pass;
            _out.WriteLine($"[hier-stationary] {names[k]}: MCMC {mean:F5} ± {se:F5}  quadrature {refs[k]:F5}  z = {zScore:F2}{(pass ? "" : "  FAIL")}");
        }
        Assert.True(ok, "stationary means disagree with the quadrature reference (see output)");
    }

    // ── Assertions ────────────────────────────────────────────────────────────

    private void AssertStandardNormal(double[][] w, string label)
    {
        int count = w.Length, d = w[0].Length;
        var mean = new double[d];
        foreach (double[] v in w) for (int i = 0; i < d; i++) mean[i] += v[i] / count;
        double worstMean = 0.0, worstCov = 0.0;
        bool ok = true;
        for (int i = 0; i < d; i++)
        {
            double zm = mean[i] * Math.Sqrt(count);
            worstMean = Math.Max(worstMean, Math.Abs(zm));
            ok &= Math.Abs(zm) < 5.0;
            for (int j = 0; j <= i; j++)
            {
                double c = 0.0;
                foreach (double[] v in w) c += (v[i] - mean[i]) * (v[j] - mean[j]);
                c /= count - 1;
                double target = i == j ? 1.0 : 0.0;
                double zc = (c - target) / Math.Sqrt((1.0 + target) / count);
                worstCov = Math.Max(worstCov, Math.Abs(zc));
                ok &= Math.Abs(zc) < 5.0;
            }
        }
        _out.WriteLine($"[hier-conditional] {label}: worst |z| mean {worstMean:F2}, covariance {worstCov:F2}");
        Assert.True(ok, $"{label}: whitened draws are not N(0, I) (worst |z| mean {worstMean:F2}, covariance {worstCov:F2})");
    }

    private void AssertUnitGamma(double[] g, double shape, string label)
    {
        int count = g.Length;
        double mean = 0.0;
        foreach (double v in g) mean += v / count;
        double var = 0.0;
        foreach (double v in g) var += (v - mean) * (v - mean);
        var /= count - 1;
        double zMean = (mean - shape) / Math.Sqrt(shape / count);
        double zVar = (var - shape) / Math.Sqrt((2.0 * shape * shape + 6.0 * shape) / count);
        _out.WriteLine($"[hier-conditional] {label}: Gamma({shape}, 1) mean {mean:F4} (z {zMean:F2}), variance {var:F4} (z {zVar:F2})");
        Assert.True(Math.Abs(zMean) < 5.0 && Math.Abs(zVar) < 5.0,
            $"{label}: scaled draws are not Gamma({shape}, 1): mean {mean:F4} (z {zMean:F2}), variance {var:F4} (z {zVar:F2})");
    }

    // ── Quadrature reference ──────────────────────────────────────────────────

    // Posterior means of σ, τ_u, τ_0 and β_0 by quadrature over (ln σ², ln τ_u², ln τ_0²). Given the variances, all
    // coefficients x = (β_1, …, β_R, β_0) are jointly Gaussian with precision Λ and linear term h, so they integrate out
    // exactly: log p(θ | y) = ½hᵀΛ⁻¹h − ½ log|Λ| + the variance-only terms, and E[β_0 | θ, y] is the β_0 block of Λ⁻¹h.
    private static (double, double, double, double[]) QuadratureMeans(double[,] z, double[][] y)
    {
        int nu = z.GetLength(1);
        double[,] ztz = Gram(z);
        double[,] p = PenaltyMatrix(nu);
        var zty = new double[Reps][];
        double yy = 0.0;
        for (int r = 0; r < Reps; r++)
        {
            zty[r] = MulTranspose(z, y[r]);
            foreach (double v in y[r]) yy += v * v;
        }
        int dim = (Reps + 1) * nu;
        var lam = new double[dim, dim];
        var h = new double[dim];

        double LogPost(double u1, double u2, double u3, double[]? beta0)
        {
            double s2 = Math.Exp(u1), tu = Math.Exp(u2), t0 = Math.Exp(u3);
            Array.Clear(lam);
            Array.Clear(h);
            int o = Reps * nu;
            for (int r = 0; r < Reps; r++)
            {
                int b = r * nu;
                for (int i = 0; i < nu; i++)
                {
                    for (int j = 0; j < nu; j++) lam[b + i, b + j] = ztz[i, j] / s2;
                    lam[b + i, b + i] += 1.0 / tu;
                    lam[b + i, o + i] = -1.0 / tu;
                    lam[o + i, b + i] = -1.0 / tu;
                    h[b + i] = zty[r][i] / s2;
                }
            }
            for (int i = 0; i < nu; i++)
            {
                for (int j = 0; j < nu; j++) lam[o + i, o + j] = p[i, j] / t0;
                lam[o + i, o + i] += Reps / tu;
            }
            double[,] l = Cholesky(lam);
            double[] m = Solve(l, h);
            double hm = 0.0;
            for (int k = 0; k < dim; k++) hm += h[k] * m[k];
            double logDet = 0.0;
            for (int k = 0; k < dim; k++) logDet += 2.0 * Math.Log(l[k, k]);

            double lp = -0.5 * Reps * N * u1 - 0.5 * yy / s2
                        - 0.5 * Reps * nu * u2 - 0.5 * (nu - Order) * u3
                        + 0.5 * hm - 0.5 * logDet;
            // IG(a, b) priors in log coordinates (Jacobian included): −a·u − b·e^{−u}.
            lp += -PriorShape * (u1 + u2 + u3) - PriorScale * (1.0 / s2 + 1.0 / tu + 1.0 / t0);
            if (beta0 is not null) Array.Copy(m, o, beta0, 0, nu);
            return lp;
        }

        // Pass 1: coarse grid to locate the posterior mass.
        const double lo = -16.0, hi = 12.0, coarse = 0.5;
        int nc = (int)Math.Round((hi - lo) / coarse) + 1;
        var lpc = new double[nc, nc, nc];
        double maxC = double.NegativeInfinity;
        for (int a = 0; a < nc; a++)
            for (int b = 0; b < nc; b++)
                for (int c = 0; c < nc; c++)
                {
                    lpc[a, b, c] = LogPost(lo + a * coarse, lo + b * coarse, lo + c * coarse, null);
                    maxC = Math.Max(maxC, lpc[a, b, c]);
                }
        var boxLo = new[] { nc, nc, nc };
        var boxHi = new[] { -1, -1, -1 };
        for (int a = 0; a < nc; a++)
            for (int b = 0; b < nc; b++)
                for (int c = 0; c < nc; c++)
                    if (lpc[a, b, c] > maxC - 35.0)
                    {
                        int[] idx = { a, b, c };
                        for (int k = 0; k < 3; k++) { boxLo[k] = Math.Min(boxLo[k], idx[k]); boxHi[k] = Math.Max(boxHi[k], idx[k]); }
                    }
        for (int k = 0; k < 3; k++)
            if (boxLo[k] == 0 || boxHi[k] == nc - 1) throw new InvalidOperationException("quadrature box touches the grid edge");

        // Pass 2: fine trapezoid grid over the box widened by one coarse cell (the integrand is negligible at its edge).
        const int nf = 64;
        var start = new double[3];
        var step = new double[3];
        for (int k = 0; k < 3; k++)
        {
            start[k] = lo + (boxLo[k] - 1) * coarse;
            step[k] = (boxHi[k] - boxLo[k] + 2) * coarse / (nf - 1);
        }
        var bufBeta = new double[nu];
        double maxF = maxC, wSum = 0.0, eSigma = 0.0, eTauU = 0.0, eTau0 = 0.0;
        var eBeta = new double[nu];
        for (int a = 0; a < nf; a++)
            for (int b = 0; b < nf; b++)
                for (int c = 0; c < nf; c++)
                {
                    double u1 = start[0] + a * step[0], u2 = start[1] + b * step[1], u3 = start[2] + c * step[2];
                    double w = Math.Exp(LogPost(u1, u2, u3, bufBeta) - maxF);
                    wSum += w;
                    eSigma += w * Math.Exp(0.5 * u1);
                    eTauU += w * Math.Exp(0.5 * u2);
                    eTau0 += w * Math.Exp(0.5 * u3);
                    for (int j = 0; j < nu; j++) eBeta[j] += w * bufBeta[j];
                }
        for (int j = 0; j < nu; j++) eBeta[j] /= wSum;
        return (eSigma / wSum, eTauU / wSum, eTau0 / wSum, eBeta);
    }

    // ── Dense helpers (independent of the banded machinery under test) ────────

    private static double[,] Gram(double[,] z)
    {
        int n = z.GetLength(0), nu = z.GetLength(1);
        var g = new double[nu, nu];
        for (int i = 0; i < nu; i++)
            for (int j = 0; j < nu; j++)
            {
                double s = 0.0;
                for (int k = 0; k < n; k++) s += z[k, i] * z[k, j];
                g[i, j] = s;
            }
        return g;
    }

    private static double[] MulTranspose(double[,] z, double[] y)
    {
        int n = z.GetLength(0), nu = z.GetLength(1);
        var v = new double[nu];
        for (int j = 0; j < nu; j++)
            for (int k = 0; k < n; k++) v[j] += z[k, j] * y[k];
        return v;
    }

    // DᵀD for second differences (rows 1, −2, 1).
    private static double[,] PenaltyMatrix(int nu)
    {
        var p = new double[nu, nu];
        double[] c = { 1.0, -2.0, 1.0 };
        for (int row = 0; row + Order < nu; row++)
            for (int a = 0; a <= Order; a++)
                for (int b = 0; b <= Order; b++) p[row + a, row + b] += c[a] * c[b];
        return p;
    }

    private static double Quadratic(double[,] a, double[] v)
    {
        double s = 0.0;
        for (int i = 0; i < v.Length; i++)
            for (int j = 0; j < v.Length; j++) s += v[i] * a[i, j] * v[j];
        return s;
    }

    private static double[,] Cholesky(double[,] a)
    {
        int n = a.GetLength(0);
        var l = new double[n, n];
        for (int j = 0; j < n; j++)
        {
            double d = a[j, j];
            for (int k = 0; k < j; k++) d -= l[j, k] * l[j, k];
            if (!(d > 0.0)) throw new InvalidOperationException("matrix is not positive definite");
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

    private static double[] Solve(double[,] l, double[] b)
    {
        int n = b.Length;
        var x = new double[n];
        for (int i = 0; i < n; i++)
        {
            double s = b[i];
            for (int k = 0; k < i; k++) s -= l[i, k] * x[k];
            x[i] = s / l[i, i];
        }
        for (int i = n - 1; i >= 0; i--)
        {
            double s = x[i];
            for (int k = i + 1; k < n; k++) s -= l[k, i] * x[k];
            x[i] = s / l[i, i];
        }
        return x;
    }

    private static double[] LowerTransposeTimes(double[,] l, double[] v)
    {
        int n = v.Length;
        var r = new double[n];
        for (int i = 0; i < n; i++)
            for (int k = i; k < n; k++) r[i] += l[k, i] * v[k];
        return r;
    }
}
