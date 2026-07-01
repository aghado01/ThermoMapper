using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The BAPS sampler infers the smoothing on a rich fixed P-spline basis. Both λ-update modes — conjugate Gibbs
/// over the variance components and Metropolis against the REML marginal evidence — must recover the smooth
/// truth, concentrate the λ posterior near the REML optimum, and mix (R̂ ≈ 1).
/// </summary>
public sealed class BapsTests
{
    private readonly ITestOutputHelper _out;
    public BapsTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Theory]
    [InlineData(BapsLambdaUpdate.Gibbs)]
    [InlineData(BapsLambdaUpdate.MarginalEvidence)]
    public void Baps_InfersSmoothing_BothModes(BapsLambdaUpdate mode)
    {
        var rng = new Xoshiro256PlusPlus(seed: 91);
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

        const int nKnots = 25;
        var knots = new double[nKnots];
        for (int k = 0; k < nKnots; k++) knots[k] = (k + 1.0) / (nKnots + 1);
        var basis = new SplineBasis(3);
        double[,] z = basis.Design(new KnotConfig(knots), x);
        int nu = z.GetLength(1);
        var penalty = new DifferencePenalty(2);

        // REML reference smoothing (grid argmax of the marginal evidence).
        var ps = new PenalizedSpline(z, y, penalty);
        double lamReml = 1.0, bestEv = double.NegativeInfinity;
        for (int j = 0; j < 31; j++)
        {
            double lam = Math.Pow(10.0, -3.0 + 8.0 * j / 30.0);
            double ev = ps.RemlLogEvidence(lam);
            if (ev > bestEv) { bestEv = ev; lamReml = lam; }
        }

        var sampler = new BapsSampler(z, y, penalty);
        BapsResult r = sampler.Run(mode, chains: 4, masterSeed: 7, burn: 500, samples: 1500);

        double mse = 0.0;
        for (int i = 0; i < n; i++)
        {
            double fit = 0.0;
            for (int p = 0; p < nu; p++) fit += z[i, p] * r.Coefficients[p];
            double d = fit - f[i];
            mse += d * d;
        }
        mse /= n;

        _out.WriteLine($"[baps {mode}] mse={mse:F5} λ̄={r.LambdaMean:F2} [{r.LambdaLo:F2},{r.LambdaHi:F2}] " +
                       $"λ_reml={lamReml:F2} σ̄={r.SigmaMean:F3} R̂={r.RHatLogLambda:F3}");

        // Recovers the smooth truth.
        Assert.True(mse < 0.01, $"BAPS fit MSE {mse:F5} should recover the smooth truth");
        // Posterior smoothing concentrates around the REML optimum (order-of-magnitude agreement).
        Assert.InRange(r.LambdaMean, lamReml / 5.0, lamReml * 5.0);
        // Noise level recovered near the truth (σ = 0.15).
        Assert.InRange(r.SigmaMean, 0.10, 0.22);
        // Chains mixed.
        Assert.True(r.RHatLogLambda < 1.2, $"R̂ {r.RHatLogLambda:F3} indicates poor mixing");
    }

    [Fact]
    public void AdaptiveBaps_BeatsGlobal_OnHeterogeneousSmoothness()
    {
        // A sharp bump on a flat background: a single global λ must either blur the bump or leave noise in the
        // flat region; locally-adaptive λ smooths the background hard and tracks the bump.
        var rng = new Xoshiro256PlusPlus(seed: 123);
        int n = 150;
        var x = new double[n];
        var f = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            f[i] = Math.Exp(-Math.Pow((x[i] - 0.5) / 0.04, 2.0));   // bump at 0.5, flat elsewhere
            y[i] = f[i] + 0.10 * Gaussian(rng);
        }

        const int nKnots = 30;
        var knots = new double[nKnots];
        for (int k = 0; k < nKnots; k++) knots[k] = (k + 1.0) / (nKnots + 1);
        var basis = new SplineBasis(3);
        double[,] z = basis.Design(new KnotConfig(knots), x);
        int nu = z.GetLength(1);
        var penalty = new DifferencePenalty(2);

        double Mse(double[] beta)
        {
            double s = 0.0;
            for (int i = 0; i < n; i++)
            {
                double fit = 0.0;
                for (int p = 0; p < nu; p++) fit += z[i, p] * beta[p];
                double d = fit - f[i];
                s += d * d;
            }
            return s / n;
        }

        var sampler = new BapsSampler(z, y, penalty);
        BapsResult global = sampler.Run(BapsLambdaUpdate.Gibbs, chains: 4, masterSeed: 5, burn: 500, samples: 1500);
        BapsResult adaptive = sampler.Run(BapsLambdaUpdate.AdaptiveGibbs, chains: 4, masterSeed: 5, burn: 500, samples: 1500);

        double mseGlobal = Mse(global.Coefficients);
        double mseAdaptive = Mse(adaptive.Coefficients);

        Assert.NotNull(adaptive.LocalSmoothing);
        double[] lam = adaptive.LocalSmoothing!;
        int len = lam.Length;
        // Smoothing in a window around the bump (centre) vs the flat ends.
        double central = 0.0; int cN = 0;
        double outer = 0.0; int oN = 0;
        for (int i = 0; i < len; i++)
        {
            double pos = (double)i / (len - 1);
            if (pos > 0.4 && pos < 0.6) { central += lam[i]; cN++; }
            else if (pos < 0.2 || pos > 0.8) { outer += lam[i]; oN++; }
        }
        central /= cN;
        outer /= oN;

        _out.WriteLine($"[adaptive] mseGlobal={mseGlobal:F5} mseAdaptive={mseAdaptive:F5} " +
                       $"λ̄(bump)={central:F2} λ̄(flat)={outer:F2} R̂={adaptive.RHatLogLambda:F3}");

        Assert.True(mseAdaptive < mseGlobal,
            $"adaptive MSE {mseAdaptive:F5} should beat global {mseGlobal:F5} on heterogeneous smoothness");
        Assert.True(central < outer,
            $"local smoothing at the bump ({central:F2}) should be lighter than the flat region ({outer:F2})");
        Assert.True(adaptive.RHatLogLambda < 1.3, $"R̂ {adaptive.RHatLogLambda:F3} indicates poor mixing");
    }
}
