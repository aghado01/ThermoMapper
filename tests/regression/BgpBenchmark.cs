using System;
using Maths.Regression.Bgp;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Reproduction of the Swiss-Roll rate experiment (Tang, Wu, Cheng &amp; Dunson 2025, Fig. 1): the empirical-Bayes
/// GP — which never sees the intrinsic dimension d=2 — should track the evidence-optimal ("oracle") bandwidth GP
/// in out-of-sample error as the training size grows, both generalizing on the manifold. Guarded by the
/// <c>BGP_BENCH</c> environment variable so the default suite stays fast; run explicitly to produce the table.
/// </summary>
public sealed class BgpBenchmark
{
    private readonly ITestOutputHelper _out;
    public BgpBenchmark(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static (double[,] X, double[] F) SwissRoll(int n, Xoshiro256PlusPlus rng)
    {
        var x = new double[n, 3];
        var f = new double[n];
        for (int i = 0; i < n; i++)
        {
            double u = rng.NextDouble(), v = rng.NextDouble();
            double phi = 1.5 * Math.PI * (1.0 + 2.0 * u);
            x[i, 0] = phi * Math.Cos(phi);
            x[i, 1] = phi * Math.Sin(phi);
            x[i, 2] = 12.0 * v;
            f[i] = Math.Sin(2.0 * Math.PI * u) + 0.5 * Math.Cos(2.0 * Math.PI * v);
        }
        return (x, f);
    }

    private static void Normalize(double[,] x, double[] lo, double[] span)
    {
        for (int i = 0; i < x.GetLength(0); i++)
            for (int k = 0; k < x.GetLength(1); k++) x[i, k] = (x[i, k] - lo[k]) / span[k];
    }

    private static double Rmse(double[] pred, double[] truth)
    {
        double s = 0.0;
        for (int i = 0; i < pred.Length; i++) { double e = pred[i] - truth[i]; s += e * e; }
        return Math.Sqrt(s / pred.Length);
    }

    [Fact]
    public void RateRecovery_SwissRoll_EbTracksOracle()
    {
        if (Environment.GetEnvironmentVariable("BGP_BENCH") is null) return;   // heavy; set BGP_BENCH to run

        int[] ns = { 50, 100, 200, 400 };
        const int reps = 8, nTest = 1500;
        const double sigma = 0.1;
        var eb = new double[ns.Length];
        var oracle = new double[ns.Length];

        _out.WriteLine("Swiss Roll (d=2 manifold in R³), σ=0.1 — RMSE on a fixed 1500-point test set, mean over 8 reps:");
        for (int gi = 0; gi < ns.Length; gi++)
        {
            int n = ns[gi];
            double ebSum = 0.0, orSum = 0.0;
            for (int rep = 0; rep < reps; rep++)
            {
                var rng = new Xoshiro256PlusPlus(1000 + n * 7 + rep);
                (double[,] xTr, double[] fTr) = SwissRoll(n, rng);
                (double[,] xTe, double[] fTe) = SwissRoll(nTest, rng);

                var lo = new double[3]; var span = new double[3];
                for (int k = 0; k < 3; k++)
                {
                    double mn = double.MaxValue, mx = double.MinValue;
                    for (int i = 0; i < n; i++) { mn = Math.Min(mn, xTr[i, k]); mx = Math.Max(mx, xTr[i, k]); }
                    lo[k] = mn; span[k] = mx - mn;
                }
                Normalize(xTr, lo, span);
                Normalize(xTe, lo, span);

                var y = new double[n];
                for (int i = 0; i < n; i++) y[i] = fTr[i] + sigma * Gaussian(rng);

                var gp = new GpRegression(xTr, y, sigma * sigma, new SquaredExponentialKernel());

                // GP-EB (ours): the empirical-Bayes prior, d never supplied.
                var prior = new EmpiricalBayesBandwidthPrior(xTr, gp);
                BgpResult r = new BgpSampler(gp, prior).Run(xTe, draws: 120, burn: 100, proposalSd: 0.5, seed: rep + 1);
                ebSum += Rmse(r.PosteriorMean, fTe);

                // Oracle: the evidence-optimal bandwidth (grid-max marginal likelihood, no prior).
                GpFit best = gp.Fit(0.05);
                for (int s = 0; s < 18; s++)
                {
                    double t = 0.01 * Math.Pow(80.0, s / 17.0);
                    GpFit fit = gp.Fit(t);
                    if (fit.LogMarginal > best.LogMarginal) best = fit;
                }
                orSum += Rmse(gp.PredictMean(best, xTe), fTe);
            }
            eb[gi] = ebSum / reps;
            oracle[gi] = orSum / reps;
            _out.WriteLine($"  n={n,4}: GP-EB RMSE={eb[gi]:F4}   oracle(marg-max) RMSE={oracle[gi]:F4}   ratio={eb[gi] / oracle[gi]:F3}");
        }

        // Empirical convergence exponent of the EB error.
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (int gi = 0; gi < ns.Length; gi++)
        {
            double lx = Math.Log(ns[gi]), ly = Math.Log(eb[gi]);
            sx += lx; sy += ly; sxx += lx * lx; sxy += lx * ly;
        }
        double slope = (ns.Length * sxy - sx * sy) / (ns.Length * sxx - sx * sx);
        _out.WriteLine($"GP-EB convergence: RMSE ∼ n^{slope:F3} (intrinsic d=2 vs ambient D=3; for C^∞ f* both rates → −1/2).");

        // The EB prior tracks the evidence-optimal bandwidth without knowing d, and the error falls with n.
        for (int gi = 0; gi < ns.Length; gi++)
            Assert.True(eb[gi] < 1.4 * oracle[gi], $"EB should track the oracle bandwidth at n={ns[gi]} (EB {eb[gi]:F4} vs oracle {oracle[gi]:F4})");
        Assert.True(eb[^1] < 0.7 * eb[0], $"EB error should fall with n ({eb[0]:F4} → {eb[^1]:F4})");
        Assert.True(slope < -0.1, $"EB error should decay with n (exponent {slope:F3})");
    }
}
