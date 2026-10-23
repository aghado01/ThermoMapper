using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The <see cref="BapsLambdaUpdate.MarginalEvidence"/> target. The profiled REML evidence ℓ_p levels off as λ → ∞,
/// so a flat prior on u = ln λ cannot be normalized. The mode therefore puts a uniform prior on the effective degrees
/// of freedom edf(λ) = tr[(ZᵀZ + λP)⁻¹ZᵀZ], which induces the density |d edf/du| on u and decays like 1/λ at the top.
/// The target is one-dimensional, so the sampler is checked against quadrature of a dense reference that shares no
/// code with the banded path: ℓ_p from an explicit residual, and the edf slope as λ·‖Z A⁻¹ Dᵀ‖²_F.
/// </summary>
public sealed class BapsMarginalEvidenceTests
{
    private const int N = 30, Knots = 8, Order = 2;

    private readonly ITestOutputHelper _out;
    public BapsMarginalEvidenceTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    // A cubic P-spline (ν = 12) on 30 points; "linear" data lies near the penalty null space, "sine" does not.
    private static (double[,] Z, double[] Y) Data(string shape)
    {
        var x = new double[N];
        for (int i = 0; i < N; i++) x[i] = (i + 0.5) / N;
        var knots = new double[Knots];
        for (int k = 0; k < Knots; k++) knots[k] = (k + 1.0) / (Knots + 1);
        double[,] z = new SplineBasis(3).Design(new KnotConfig(knots), x);
        var rng = new Xoshiro256PlusPlus(seed: 11);
        var y = new double[N];
        for (int i = 0; i < N; i++)
            y[i] = shape == "linear" ? 1.0 + 2.0 * x[i] + 0.1 * Gaussian(rng)
                                     : Math.Sin(2.0 * Math.PI * x[i]) + 0.15 * Gaussian(rng);
        return (z, y);
    }

    [Theory]
    [InlineData("linear")]
    [InlineData("sine")]
    public void ProfiledEvidence_Plateaus_SoAFlatPriorIsImproper(string shape)
    {
        (double[,] z, double[] y) = Data(shape);
        var ps = new PenalizedSpline(z, y, new DifferencePenalty(Order));
        double e6 = ps.RemlLogEvidence(1e6), e8 = ps.RemlLogEvidence(1e8);
        _out.WriteLine($"[reml-plateau {shape}] ℓ_p(1e6)={e6:F6} ℓ_p(1e8)={e8:F6}");
        Assert.True(Math.Abs(e8 - e6) < 1e-3, $"ℓ_p should have levelled off: {e6:F6} vs {e8:F6}");
    }

    [Fact]
    public void EdfAndItsSlope_MatchDenseReference_AndThePriorIsProper()
    {
        (double[,] z, double[] y) = Data("sine");
        int nu = z.GetLength(1);
        var ps = new PenalizedSpline(z, y, new DifferencePenalty(Order));
        var reference = new DenseReference(z, y);

        foreach (double lambda in new[] { 1e-3, 1e-1, 1.0, 10.0, 1e3, 1e6, 1e9 })
        {
            (_, double edf, double slope) = reference.At(lambda);
            double edfCode = ps.EffectiveDegreesOfFreedom(lambda);
            double slopeCode = ps.EdfSlope(lambda);
            _out.WriteLine($"[edf] λ={lambda:E0}: edf {edfCode:F8} (dense {edf:F8}), |d edf/d ln λ| {slopeCode:E6} (dense {slope:E6})");
            Assert.True(Math.Abs(edfCode - edf) <= 1e-7 * nu, $"edf at λ={lambda}: {edfCode} vs dense {edf}");
            Assert.True(Math.Abs(slopeCode - slope) <= 1e-6 * slope, $"edf slope at λ={lambda}: {slopeCode} vs dense {slope}");
        }

        // Uniform on edf ∈ (r, ν): the induced density integrates to ν − r over u = ln λ.
        double integral = 0.0;
        const double lo = -40.0, hi = 50.0, h = 0.01;
        for (double u = lo; u <= hi; u += h) integral += ps.EdfSlope(Math.Exp(u)) * h;
        _out.WriteLine($"[edf] ∫ |d edf/du| du = {integral:F6} (ν − r = {nu - Order})");
        Assert.True(Math.Abs(integral - (nu - Order)) < 1e-4 * (nu - Order), $"edf prior mass {integral:F6} should be ν − r");
    }

    [Theory]
    [InlineData("linear")]
    [InlineData("sine")]
    public void MarginalEvidence_MatchesQuadratureOfTheEdfPosterior(string shape)
    {
        (double[,] z, double[] y) = Data(shape);
        var reference = new DenseReference(z, y);

        // Quadrature of p(u | y) ∝ exp(ℓ_p(e^u)) · |d edf/du| on a fine grid. Near the null space the upper tail decays
        // only like e^(−u), and the dense reference loses digits beyond λ ≈ 10¹², so the grid stops at u = 28 and the
        // ends must sit 15 nats below the mode (the omitted tail mass is then about e^(−15) of the total).
        const double lo = -25.0, hi = 28.0, h = 0.005;
        int g = (int)Math.Round((hi - lo) / h) + 1;
        var logp = new double[g];
        double max = double.NegativeInfinity;
        for (int i = 0; i < g; i++)
        {
            (double ell, _, double slope) = reference.At(Math.Exp(lo + i * h));
            logp[i] = ell + Math.Log(slope);
            max = Math.Max(max, logp[i]);
        }
        _out.WriteLine($"[edf-posterior {shape}] log density at the grid ends: {logp[0] - max:F1}, {logp[g - 1] - max:F1} (relative to the mode)");
        Assert.True(logp[0] < max - 15.0 && logp[g - 1] < max - 15.0, "quadrature range must hold the mass");
        var cdf = new double[g];
        for (int i = 1; i < g; i++)
            cdf[i] = cdf[i - 1] + 0.5 * h * (Math.Exp(logp[i - 1] - max) + Math.Exp(logp[i] - max));
        for (int i = 0; i < g; i++) cdf[i] /= cdf[g - 1];
        double Inverse(double p)
        {
            int i = 1;
            while (cdf[i] < p) i++;
            double f = (p - cdf[i - 1]) / (cdf[i] - cdf[i - 1]);
            return lo + (i - 1 + f) * h;
        }
        double[] probs = { 0.025, 0.5, 0.975 };
        var refQ = new double[3];
        for (int k = 0; k < 3; k++) refQ[k] = Inverse(probs[k]);

        // Independent runs give a Monte Carlo standard error for each pooled quantile of ln λ.
        const int runs = 10;
        int[] seeds = SeedTree.Derive(101, runs);
        var sampler = new BapsSampler(z, y, new DifferencePenalty(Order));
        var est = new double[3][];
        for (int k = 0; k < 3; k++) est[k] = new double[runs];
        for (int r = 0; r < runs; r++)
        {
            BapsResult res = sampler.Run(BapsLambdaUpdate.MarginalEvidence, chains: 4, masterSeed: seeds[r], burn: 500, samples: 4000);
            est[0][r] = Math.Log(res.LambdaLo);
            est[1][r] = Math.Log(res.LambdaMedian);
            est[2][r] = Math.Log(res.LambdaHi);
        }

        bool ok = true;
        string[] names = { "2.5%", "50%", "97.5%" };
        for (int k = 0; k < 3; k++)
        {
            double mean = 0.0;
            foreach (double v in est[k]) mean += v / runs;
            double ss = 0.0;
            foreach (double v in est[k]) ss += (v - mean) * (v - mean);
            double se = Math.Sqrt(ss / (runs - 1) / runs);
            bool pass = Math.Abs(mean - refQ[k]) <= 4.5 * se + 0.02;
            ok &= pass;
            _out.WriteLine($"[edf-posterior {shape}] ln λ {names[k]}: sampler {mean:F3} ± {se:F3}  quadrature {refQ[k]:F3}{(pass ? "" : "  FAIL")}");
        }
        Assert.True(ok, "MarginalEvidence quantiles disagree with quadrature (see output)");
    }

    [Theory]
    [InlineData("linear")]
    [InlineData("sine")]
    public void MarginalEvidence_ConditionalDraws_MatchQuadrature(string shape)
    {
        // Each retained λ draws σ² | λ, y ~ IG((n − r)/2, RSS_pen/2) and then β | σ², λ, y ~ N(β̂, σ²A⁻¹). Their marginal
        // moments are mixtures over p(u | y): E[σ² | λ] = RSS_pen/(n − r − 2), Var[σ² | λ] = E[σ² | λ]²/((n − r)/2 − 2),
        // E[f | λ] = zᵀβ̂ and Var[f | λ] = E[σ² | λ]·zᵀA⁻¹z for f = zᵀβ at a mid-design row, all from the dense reference.
        (double[,] z, double[] y) = Data(shape);
        var reference = new DenseReference(z, y);
        int nu = z.GetLength(1);
        var zRow = new double[nu];
        for (int j = 0; j < nu; j++) zRow[j] = z[N / 2, j];

        const double lo = -25.0, hi = 28.0, h = 0.005;   // the grid of the λ quadrature above
        int g = (int)Math.Round((hi - lo) / h) + 1;
        var logp = new double[g];
        double max = double.NegativeInfinity;
        for (int i = 0; i < g; i++)
        {
            (double ell, _, double slope) = reference.At(Math.Exp(lo + i * h));
            logp[i] = ell + Math.Log(slope);
            max = Math.Max(max, logp[i]);
        }
        double shapeIg = 0.5 * (N - Order);
        double mass = 0.0, s1 = 0.0, s2 = 0.0, f1 = 0.0, f2 = 0.0;
        for (int i = 0; i < g; i++)
        {
            double w = Math.Exp(logp[i] - max) * (i == 0 || i == g - 1 ? 0.5 : 1.0);   // trapezoid
            (double rss, double fitMean, double fitScale) = reference.Conditional(Math.Exp(lo + i * h), zRow);
            double e = rss / (N - Order - 2);
            mass += w;
            s1 += w * e;
            s2 += w * (e * e / (shapeIg - 2.0) + e * e);
            f1 += w * fitMean;
            f2 += w * (e * fitScale + fitMean * fitMean);
        }
        double refSigma2 = s1 / mass, refSigma2Var = s2 / mass - refSigma2 * refSigma2;
        double refFit = f1 / mass, refFitVar = f2 / mass - refFit * refFit;

        const int runs = 10;
        int[] seeds = SeedTree.Derive(202, runs);
        var sampler = new BapsSampler(z, y, new DifferencePenalty(Order));
        var est = new double[4][];
        for (int k = 0; k < 4; k++) est[k] = new double[runs];
        for (int r = 0; r < runs; r++)
        {
            BapsResult res = sampler.Run(BapsLambdaUpdate.MarginalEvidence, chains: 4, masterSeed: seeds[r], burn: 500, samples: 4000);
            double count = 0.0, m1 = 0.0, m2 = 0.0, q1 = 0.0, q2 = 0.0;
            foreach (var chain in res.PosteriorDraws)
                foreach (BapsDraw d in chain)
                {
                    double fit = 0.0;
                    for (int j = 0; j < nu; j++) fit += zRow[j] * d.Coefficients[j];
                    count++;
                    m1 += d.Sigma2; m2 += d.Sigma2 * d.Sigma2;
                    q1 += fit; q2 += fit * fit;
                }
            est[0][r] = m1 / count;
            est[1][r] = m2 / count - (m1 / count) * (m1 / count);
            est[2][r] = q1 / count;
            est[3][r] = q2 / count - (q1 / count) * (q1 / count);
            if (r == 0) AssertDrawsMatchTheirConditionals(res, reference, zRow, shape);
        }

        double[] refs = { refSigma2, refSigma2Var, refFit, refFitVar };
        string[] names = { "E[σ²]", "Var[σ²]", "E[f]", "Var[f]" };
        bool ok = true;
        for (int k = 0; k < 4; k++)
        {
            double mean = 0.0;
            foreach (double v in est[k]) mean += v / runs;
            double ss = 0.0;
            foreach (double v in est[k]) ss += (v - mean) * (v - mean);
            double se = Math.Sqrt(ss / (runs - 1) / runs);
            bool pass = Math.Abs(mean - refs[k]) <= 4.5 * se + 1e-3 * Math.Abs(refs[k]);
            ok &= pass;
            _out.WriteLine($"[marginal-draws {shape}] {names[k]}: sampler {mean:E5} ± {se:E2}  quadrature {refs[k]:E5}  z = {(mean - refs[k]) / se:F2}{(pass ? "" : "  FAIL")}");
        }
        Assert.True(ok, "MarginalEvidence conditional draws disagree with quadrature (see output)");
    }

    // The conditional draws alone, separated from the λ chain: on the run's own λ sample, the draws' mean σ² and f
    // must match the mean of their exact conditional expectations. Given the λs the draws are independent (their own
    // stream), so the standard error is the root of the summed conditional variances over the draw count.
    private void AssertDrawsMatchTheirConditionals(BapsResult res, DenseReference reference, double[] zRow, string shape)
    {
        double shapeIg = 0.5 * (N - Order);
        double count = 0.0, drawSigma2 = 0.0, condSigma2 = 0.0, varSigma2 = 0.0, drawFit = 0.0, condFit = 0.0, varFit = 0.0;
        foreach (var chain in res.PosteriorDraws)
            foreach (BapsDraw d in chain)
            {
                (double rss, double fitMean, double fitScale) = reference.Conditional(d.Lambda, zRow);
                double e = rss / (N - Order - 2);
                double fit = 0.0;
                for (int j = 0; j < zRow.Length; j++) fit += zRow[j] * d.Coefficients[j];
                count++;
                drawSigma2 += d.Sigma2; condSigma2 += e; varSigma2 += e * e / (shapeIg - 2.0);
                drawFit += fit; condFit += fitMean; varFit += e * fitScale;
            }
        double zSigma2 = (drawSigma2 - condSigma2) / Math.Sqrt(varSigma2);
        double zFit = (drawFit - condFit) / Math.Sqrt(varFit);
        _out.WriteLine($"[marginal-draws {shape}] given the run's λs: σ² z = {zSigma2:F2}, f z = {zFit:F2} over {count} draws");
        Assert.InRange(zSigma2, -4.0, 4.0);
        Assert.InRange(zFit, -4.0, 4.0);
    }

    // Dense reference: ℓ_p from an explicit penalized residual, edf = tr(A⁻¹G), and |d edf/d ln λ| = λ‖Z A⁻¹ Dᵀ‖²_F.
    private sealed class DenseReference
    {
        private readonly double[,] _z;
        private readonly double[] _y;
        private readonly double[,] _g;
        private readonly double[,] _d;
        private readonly double[] _zty;
        private readonly int _n, _nu;

        public DenseReference(double[,] z, double[] y)
        {
            _z = z; _y = y;
            _n = z.GetLength(0);
            _nu = z.GetLength(1);
            _g = new double[_nu, _nu];
            for (int i = 0; i < _nu; i++)
                for (int j = 0; j < _nu; j++)
                    for (int k = 0; k < _n; k++) _g[i, j] += z[k, i] * z[k, j];
            _zty = new double[_nu];
            for (int j = 0; j < _nu; j++)
                for (int k = 0; k < _n; k++) _zty[j] += z[k, j] * y[k];
            _d = new double[_nu - Order, _nu];   // second differences: rows (1, −2, 1)
            for (int r = 0; r < _nu - Order; r++) { _d[r, r] = 1.0; _d[r, r + 1] = -2.0; _d[r, r + 2] = 1.0; }
        }

        // The Cholesky factor of A = ZᵀZ + λDᵀD.
        private double[,] FactorA(double lambda)
        {
            var a = new double[_nu, _nu];
            for (int i = 0; i < _nu; i++)
                for (int j = 0; j < _nu; j++)
                {
                    double p = 0.0;
                    for (int r = 0; r < _nu - Order; r++) p += _d[r, i] * _d[r, j];
                    a[i, j] = _g[i, j] + lambda * p;
                }
            return Cholesky(a);
        }

        // Given λ: RSS_pen from an explicit penalized residual, and for f = zᵀβ its conditional mean zᵀβ̂ and zᵀA⁻¹z,
        // so that Var[f | σ², λ, y] = σ²·zᵀA⁻¹z.
        public (double Rss, double FitMean, double FitScale) Conditional(double lambda, double[] zRow)
        {
            double[,] l = FactorA(lambda);
            double[] beta = Solve(l, _zty);
            double rss = 0.0;
            for (int k = 0; k < _n; k++)
            {
                double fit = 0.0;
                for (int j = 0; j < _nu; j++) fit += _z[k, j] * beta[j];
                rss += (_y[k] - fit) * (_y[k] - fit);
            }
            for (int r = 0; r < _nu - Order; r++)
            {
                double db = 0.0;
                for (int j = 0; j < _nu; j++) db += _d[r, j] * beta[j];
                rss += lambda * db * db;
            }
            double mean = 0.0, scale = 0.0;
            double[] aInvZ = Solve(l, zRow);
            for (int j = 0; j < _nu; j++) { mean += zRow[j] * beta[j]; scale += zRow[j] * aInvZ[j]; }
            return (rss, mean, scale);
        }

        public (double Ell, double Edf, double Slope) At(double lambda)
        {
            double[,] l = FactorA(lambda);
            double logDet = 0.0;
            for (int i = 0; i < _nu; i++) logDet += 2.0 * Math.Log(l[i, i]);

            double[] beta = Solve(l, _zty);
            double rss = 0.0;
            for (int k = 0; k < _n; k++)
            {
                double fit = 0.0;
                for (int j = 0; j < _nu; j++) fit += _z[k, j] * beta[j];
                rss += (_y[k] - fit) * (_y[k] - fit);
            }
            for (int r = 0; r < _nu - Order; r++)
            {
                double db = 0.0;
                for (int j = 0; j < _nu; j++) db += _d[r, j] * beta[j];
                rss += lambda * db * db;
            }
            double ell = -0.5 * ((_n - Order) * Math.Log(rss) - (_nu - Order) * Math.Log(lambda) + logDet);

            double edf = 0.0;
            var col = new double[_nu];
            for (int k = 0; k < _nu; k++)
            {
                for (int j = 0; j < _nu; j++) col[j] = _g[j, k];
                edf += Solve(l, col)[k];
            }

            double slope = 0.0;
            for (int r = 0; r < _nu - Order; r++)
            {
                for (int j = 0; j < _nu; j++) col[j] = _d[r, j];
                double[] w = Solve(l, col);
                for (int k = 0; k < _n; k++)
                {
                    double zw = 0.0;
                    for (int j = 0; j < _nu; j++) zw += _z[k, j] * w[j];
                    slope += zw * zw;
                }
            }
            return (ell, edf, lambda * slope);
        }

        private static double[,] Cholesky(double[,] a)
        {
            int n = a.GetLength(0);
            var l = new double[n, n];
            for (int j = 0; j < n; j++)
            {
                double s = a[j, j];
                for (int k = 0; k < j; k++) s -= l[j, k] * l[j, k];
                if (!(s > 0.0)) throw new InvalidOperationException("matrix is not positive definite");
                l[j, j] = Math.Sqrt(s);
                for (int i = j + 1; i < n; i++)
                {
                    double t = a[i, j];
                    for (int k = 0; k < j; k++) t -= l[i, k] * l[j, k];
                    l[i, j] = t / l[j, j];
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
    }
}
