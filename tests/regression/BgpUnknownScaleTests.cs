using System;
using System.Collections.Generic;
using Maths.Regression.Bgp;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BGP with an unknown scale: K = s²(Σ_t + ηI), s² integrated out under 1/s² and drawn exactly given (t, η), and a
/// joint walk on (log t, log η). The closed-form marginal is checked by quadrature over s² against a dense reference;
/// rescaling the response shifts it by a constant; the s² draws match their conditional; the noise variance of a
/// response far from unit scale is recovered; and checkpoints resume exactly.
/// </summary>
public sealed class BgpUnknownScaleTests
{
    private readonly ITestOutputHelper _out;
    public BgpUnknownScaleTests(ITestOutputHelper output) => _out = output;

    private static (double[,] X, double[] Y, double[] F) Data(int n, double amplitude, double noiseSd, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[n, 1];
        var y = new double[n];
        var f = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = (i + 0.5) / n;
            f[i] = amplitude * Math.Sin(2.0 * Math.PI * x[i, 0]);
            double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
            y[i] = f[i] + noiseSd * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }
        return (x, y, f);
    }

    [Fact]
    public void Marginal_MatchesQuadratureOverTheSignalVariance()
    {
        (double[,] x, double[] y, _) = Data(20, 3.0, 0.4, seed: 4);
        var kernel = new SquaredExponentialKernel();
        GpRegression gp = GpRegression.WithUnknownScale(x, y, kernel, GpInputScaling.None);
        int n = y.Length;
        foreach ((double t, double eta) in new[] { (0.01, 0.05), (0.03, 0.001), (0.005, 0.5) })
        {
            // Dense reference with the mean integrated out: log ∫ N(y; m1, s²R) dm
            //   = −((n−1)/2) log(2π s²) − ½ log|R| − ½ log 1ᵀR⁻¹1 − Q/(2s²),
            // then s² integrated against 1/s² by quadrature in u = ln s² (ds²/s² = du).
            var r = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    r[i, j] = kernel.Evaluate((x[i, 0] - x[j, 0]) * (x[i, 0] - x[j, 0]), t) + (i == j ? eta : 0.0);
            double[,] l = Cholesky(r);
            double logDet = 0.0;
            for (int i = 0; i < n; i++) logDet += 2.0 * Math.Log(l[i, i]);
            var ones = new double[n];
            Array.Fill(ones, 1.0);
            double[] rInvOne = Solve(l, ones), rInvY = Solve(l, y);
            double a = 0.0, b = 0.0, c = 0.0;
            for (int i = 0; i < n; i++) { a += rInvOne[i]; b += rInvY[i]; c += y[i] * rInvY[i]; }
            double q = c - b * b / a;

            double centre = Math.Log(q / (n - 1)), h = 1e-3, max = double.NegativeInfinity;
            var terms = new List<double>();
            for (double u = centre - 20.0; u <= centre + 20.0; u += h)
            {
                double s2 = Math.Exp(u);
                double term = -0.5 * (n - 1) * Math.Log(2.0 * Math.PI * s2) - 0.5 * logDet - 0.5 * Math.Log(a) - q / (2.0 * s2);
                terms.Add(term);
                max = Math.Max(max, term);
            }
            double sum = 0.0;
            foreach (double term in terms) sum += Math.Exp(term - max) * h;
            double reference = max + Math.Log(sum);

            double code = gp.Fit(t, eta).LogMarginal;
            _out.WriteLine($"t={t} η={eta}: closed form {code:F8}, quadrature {reference:F8}");
            Assert.Equal(reference, code, 6);
        }
    }

    [Fact]
    public void RescalingTheResponse_ShiftsTheMarginalByAConstant()
    {
        (double[,] x, double[] y, _) = Data(40, 1.0, 0.2, seed: 6);
        const double c = 37.5;
        var scaled = new double[y.Length];
        for (int i = 0; i < y.Length; i++) scaled[i] = c * y[i];
        GpRegression gp = GpRegression.WithUnknownScale(x, y, new SquaredExponentialKernel());
        GpRegression gpScaled = GpRegression.WithUnknownScale(x, scaled, new SquaredExponentialKernel());
        foreach ((double t, double eta) in new[] { (0.01, 0.05), (0.03, 0.001), (0.2, 2.0) })
        {
            GpFit fit = gp.Fit(t, eta), fitScaled = gpScaled.Fit(t, eta);
            Assert.Equal(-(y.Length - 1) * Math.Log(c), fitScaled.LogMarginal - fit.LogMarginal, 9);
            Assert.Equal(c * c * fit.Quadratic, fitScaled.Quadratic, 9 - (int)Math.Log10(c * c * fit.Quadratic + 1));
            Assert.Equal(c * fit.Mean, fitScaled.Mean, 9);
        }
    }

    [Fact]
    public void SignalVarianceDraws_MatchTheirConditional()
    {
        // Given each draw's own (t, η), s² ~ IG((n−1)/2, Q/2), so E[s² | ·] = Q/(n − 3) and Var[s² | ·] = E²/((n−1)/2 − 2).
        // The draws come from their own stream, independent given the chain, so their mean is tested against the mean of
        // the conditional means with the root of the summed conditional variances as the standard error.
        (double[,] x, double[] y, _) = Data(60, 4.0, 0.5, seed: 8);
        GpRegression gp = GpRegression.WithUnknownScale(x, y, new SquaredExponentialKernel());
        int n = y.Length;
        BgpResult r = new BgpSampler(gp, new EmpiricalBayesBandwidthPrior(gp)).Run(Grid(25), draws: 1500, burn: 300, seed: 3, chains: 2);
        double count = 0.0, drawn = 0.0, expected = 0.0, variance = 0.0;
        foreach (IReadOnlyList<BgpDraw> chain in r.BandwidthDraws)
            foreach (BgpDraw d in chain)
            {
                double e = d.Fit.Quadratic / (n - 3);
                count++;
                drawn += d.SignalVariance;
                expected += e;
                variance += e * e / (0.5 * (n - 1) - 2.0);
                Assert.Equal(d.SignalVariance * d.NoiseRatio, d.NoiseVariance, 12);
            }
        double z = (drawn - expected) / Math.Sqrt(variance);
        _out.WriteLine($"s² draws: mean {drawn / count:E4} against E[s² | ·] {expected / count:E4}, z = {z:F2}");
        Assert.InRange(z, -4.0, 4.0);
    }

    [Fact]
    public void NoiseVariance_IsRecovered_OnAResponseFarFromUnitScale()
    {
        // Signal amplitude 6, noise sd 0.5: a unit signal variance would be badly mis-scaled here.
        (double[,] x, double[] y, double[] f) = Data(80, 6.0, 0.5, seed: 12);
        GpRegression gp = GpRegression.WithUnknownScale(x, y, new SquaredExponentialKernel());
        BgpResult r = new BgpSampler(gp, new EmpiricalBayesBandwidthPrior(gp)).Run(x, draws: 1500, burn: 400, seed: 2, chains: 2);

        PosteriorInterval noise = r.NoiseVarianceInterval, signal = r.SignalVarianceInterval;
        double mse = 0.0;
        for (int i = 0; i < f.Length; i++) mse += (r.PosteriorMean[i] - f[i]) * (r.PosteriorMean[i] - f[i]);
        mse /= f.Length;
        _out.WriteLine($"σ² median {noise.Median:F4} [{noise.Lo:F4}, {noise.Hi:F4}] (true 0.25); s² median {signal.Median:F2} " +
                       $"[{signal.Lo:F2}, {signal.Hi:F2}]; η median {r.NoiseRatioInterval.Median:E3}; fit MSE {mse:F4}; " +
                       $"t R̂ {r.BandwidthConvergence.RHat}, η R̂ {r.NoiseRatioConvergence!.Value.RHat}, acceptance {r.AcceptanceRate:F2}");
        Assert.InRange(0.25, noise.Lo, noise.Hi);
        Assert.True(signal.Median > 1.0, "the signal variance should be far above the unit the known-σ² model assumes");
        Assert.True(mse < 0.05, $"posterior mean MSE {mse:F4}");
    }

    [Fact]
    public void ResumedRun_MatchesTheUninterruptedRun()
    {
        (double[,] x, double[] y, _) = Data(50, 3.0, 0.3, seed: 9);
        GpRegression gp = GpRegression.WithUnknownScale(x, y, new SquaredExponentialKernel());
        var prior = new EmpiricalBayesBandwidthPrior(gp);
        double[,] grid = Grid(15);
        BgpResult straight = new BgpSampler(gp, prior).Run(grid, draws: 400, burn: 100, seed: 5, chains: 2);
        EnsembleCheckpoint checkpoint = new BgpSampler(gp, prior).Run(grid, draws: 150, burn: 100, seed: 5, chains: 2, checkpoint: true).Checkpoint!;
        BgpResult resumed = new BgpSampler(gp, prior).Run(grid, draws: 250, burn: 100, seed: 5, chains: 2, resumeFrom: checkpoint);

        Assert.Equal(straight.PosteriorMean, resumed.PosteriorMean);
        Assert.Equal(straight.PosteriorSd, resumed.PosteriorSd);
        Assert.Equal(straight.NoiseRatioConvergence, resumed.NoiseRatioConvergence);
        for (int c = 0; c < straight.BandwidthDraws.Count; c++)
            for (int d = 0; d < straight.BandwidthDraws[c].Count; d++)
            {
                BgpDraw a = straight.BandwidthDraws[c][d], b = resumed.BandwidthDraws[c][d];
                Assert.Equal(a.Bandwidth, b.Bandwidth);
                Assert.Equal(a.NoiseRatio, b.NoiseRatio);
                Assert.Equal(a.SignalVariance, b.SignalVariance);
            }
        Assert.Throws<ArgumentException>(() => new BgpSampler(gp, prior, noiseRatioUpper: 10.0)
            .Run(grid, draws: 250, burn: 100, seed: 5, chains: 2, resumeFrom: checkpoint));
    }

    [Fact]
    public void KnownNoise_DrawsReportUnitSignalAndTheGivenNoise()
    {
        (double[,] x, double[] y, _) = Data(40, 1.0, 0.1, seed: 3);
        var gp = new GpRegression(x, y, 0.01, new SquaredExponentialKernel());
        BgpResult r = new BgpSampler(gp, new EmpiricalBayesBandwidthPrior(gp)).Run(Grid(10), draws: 200, burn: 50, seed: 1, chains: 2);
        Assert.Null(r.NoiseRatioConvergence);
        Assert.Equal(new PosteriorInterval(0.01, 0.01, 0.01), r.NoiseVarianceInterval);
        Assert.Equal(new PosteriorInterval(1.0, 1.0, 1.0), r.SignalVarianceInterval);
    }

    [Fact]
    public void TheFitMustMatchTheNoiseModel()
    {
        (double[,] x, double[] y, _) = Data(20, 1.0, 0.1, seed: 3);
        GpRegression unknown = GpRegression.WithUnknownScale(x, y, new SquaredExponentialKernel());
        var known = new GpRegression(x, y, 0.01, new SquaredExponentialKernel());
        Assert.Equal(GpNoiseModel.UnknownScale, unknown.NoiseModel);
        Assert.Throws<InvalidOperationException>(() => unknown.Fit(0.01));
        Assert.Throws<InvalidOperationException>(() => known.Fit(0.01, 0.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => unknown.Fit(0.01, 0.0));
        var constant = new double[y.Length];
        Array.Fill(constant, 2.5);
        Assert.Throws<ArgumentException>(() => GpRegression.WithUnknownScale(x, constant, new SquaredExponentialKernel()));
    }

    private static double[,] Grid(int points)
    {
        var grid = new double[points, 1];
        for (int i = 0; i < points; i++) grid[i, 0] = i / (points - 1.0);
        return grid;
    }

    private static double[,] Cholesky(double[,] a)
    {
        int n = a.GetLength(0);
        var l = new double[n, n];
        for (int j = 0; j < n; j++)
        {
            double s = a[j, j];
            for (int k = 0; k < j; k++) s -= l[j, k] * l[j, k];
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
