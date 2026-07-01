using System;
using Maths.LinAlg;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The difference penalty is the P-spline (measure-side) smoother: on a deliberately over-rich fixed B-spline
/// basis the unpenalized fit chases noise, and adding λ·DᵀD into the (already banded) normal equations both
/// drives the coefficient roughness βᵀDᵀD β down monotonically in λ and recovers a fit closer to the truth.
/// </summary>
public sealed class PenaltyTests
{
    private readonly ITestOutputHelper _out;
    public PenaltyTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double Mse(double[] a, double[] b)
    {
        double s = 0.0;
        for (int i = 0; i < a.Length; i++) { double d = a[i] - b[i]; s += d * d; }
        return s / a.Length;
    }

    // Σ_i (Δ² β)_i² — the order-2 coefficient roughness the penalty shrinks.
    private static double Roughness(double[] beta)
    {
        double s = 0.0;
        for (int i = 0; i + 2 < beta.Length; i++)
        {
            double d = beta[i] - 2.0 * beta[i + 1] + beta[i + 2];
            s += d * d;
        }
        return s;
    }

    [Fact]
    public void DifferencePenalty_SmoothsOverRichBasis()
    {
        var rng = new Xoshiro256PlusPlus(seed: 31);
        int n = 120;
        var x = new double[n];
        var f = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            f[i] = Math.Sin(2.0 * Math.PI * x[i]);
            y[i] = f[i] + 0.15 * Gaussian(rng);
        }

        // Deliberately over-rich fixed basis (P-spline style): many equally-spaced interior knots.
        const int nKnots = 25;
        var knots = new double[nKnots];
        for (int k = 0; k < nKnots; k++) knots[k] = (k + 1.0) / (nKnots + 1);
        var basis = new SplineBasis(3);
        double[,] z = basis.Design(new KnotConfig(knots), x);
        int nu = z.GetLength(1);

        var bd = new BandedDesign(z);
        var penalty = new DifferencePenalty(2);
        int bw = Math.Max(bd.Bandwidth, penalty.Order);

        double[] Coefficients(double lambda)
        {
            var band = new double[bw + 1, nu];
            var b = new double[nu];
            bd.Accumulate(null, band, y, b);
            penalty.AccumulateInto(band, nu, lambda);
            var chol = new BandCholesky(nu, bw, BandFactorization.Ldlt);
            chol.DecomposeBanded(band);
            return chol.Solve(b);
        }

        double[] Fit(double[] beta)
        {
            var fit = new double[n];
            for (int i = 0; i < n; i++)
            {
                double s = 0.0;
                for (int j = 0; j < nu; j++) s += z[i, j] * beta[j];
                fit[i] = s;
            }
            return fit;
        }

        // Roughness strictly decreases as the penalty tightens — the operator does what it should.
        double r0 = Roughness(Coefficients(0.0));
        double r1 = Roughness(Coefficients(1.0));
        double r2 = Roughness(Coefficients(100.0));
        _out.WriteLine($"[penalty] roughness λ=0:{r0:E3} λ=1:{r1:E3} λ=100:{r2:E3}");
        Assert.True(r1 < r0, $"roughness should fall with λ: {r1:E3} !< {r0:E3}");
        Assert.True(r2 < r1, $"roughness should fall with λ: {r2:E3} !< {r1:E3}");

        // The best penalized fit beats the unpenalized rich-basis overfit on MSE-to-truth.
        double mseOver = Mse(Fit(Coefficients(0.0)), f);
        double mseBest = double.PositiveInfinity;
        double bestLambda = 0.0;
        foreach (double lambda in new[] { 0.3, 1.0, 3.0, 10.0, 30.0 })
        {
            double mse = Mse(Fit(Coefficients(lambda)), f);
            if (mse < mseBest) { mseBest = mse; bestLambda = lambda; }
        }
        _out.WriteLine($"[penalty] mseOverfit={mseOver:F5} mseBest={mseBest:F5} (λ*={bestLambda}) ratio={mseBest / mseOver:F3}");
        Assert.True(mseBest < mseOver, $"penalized MSE {mseBest:F5} should beat unpenalized {mseOver:F5}");
    }
}
