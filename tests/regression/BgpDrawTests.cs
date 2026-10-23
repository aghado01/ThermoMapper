using System;
using System.Collections.Generic;
using Maths.Regression.Bgp;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The retained bandwidth draws (<see cref="BgpResult.BandwidthDraws"/>): every draw of every chain is kept with the GP
/// fit the chain held there, and the posterior mean and the bandwidth interval are exactly what those draws give.
/// </summary>
public sealed class BgpDrawTests
{
    [Fact]
    public void BandwidthDraws_KeepEveryDraw_AndGiveTheSummaries()
    {
        var rng = new Xoshiro256PlusPlus(seed: 31);
        int n = 60;
        var x = new double[n, 1];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = (i + 0.5) / n;
            double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
            y[i] = Math.Sin(2.0 * Math.PI * x[i, 0]) + 0.1 * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }
        var xTest = new double[25, 1];
        for (int i = 0; i < 25; i++) xTest[i, 0] = i / 24.0;
        var gp = new GpRegression(x, y, 0.01, new SquaredExponentialKernel());
        BgpResult r = new BgpSampler(gp, new EmpiricalBayesBandwidthPrior(gp))
            .Run(xTest, draws: 300, burn: 100, proposalSd: 0.5, seed: 5, chains: 3);

        Assert.Equal(3, r.BandwidthDraws.Count);
        var bandwidths = new List<double>();
        var mean = new double[25];
        int shared = 0;
        foreach (IReadOnlyList<BgpDraw> chain in r.BandwidthDraws)
        {
            Assert.Equal(r.Draws / 3, chain.Count);
            for (int d = 0; d < chain.Count; d++)
            {
                BgpDraw draw = chain[d];
                // The fit is the one at the draw's bandwidth (the start's t = exp(log t0) may differ from t0 in the last bit).
                Assert.Equal(draw.Bandwidth, draw.Fit.Bandwidth, 1e-14 * draw.Bandwidth);
                if (d > 0 && ReferenceEquals(draw.Fit, chain[d - 1].Fit)) shared++;
                bandwidths.Add(draw.Bandwidth);
                double[] pred = gp.PredictMean(draw.Fit, xTest);
                for (int i = 0; i < 25; i++) mean[i] += pred[i];
            }
        }
        Assert.True(shared > 0, "rejected proposals should share the previous fit");

        for (int i = 0; i < 25; i++) Assert.Equal(mean[i] / r.Draws, r.PosteriorMean[i], 12);
        bandwidths.Sort();
        Assert.Equal(Type7(bandwidths, 0.025), r.Bandwidth.Lo);
        Assert.Equal(Type7(bandwidths, 0.975), r.Bandwidth.Hi);
    }

    private static double Type7(List<double> sorted, double q)
    {
        double pos = q * (sorted.Count - 1);
        int lo = (int)Math.Floor(pos);
        double frac = pos - lo;
        return lo + 1 < sorted.Count ? sorted[lo] * (1.0 - frac) + sorted[lo + 1] * frac : sorted[lo];
    }
}
