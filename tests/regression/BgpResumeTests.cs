using System;
using System.Collections.Generic;
using Maths.Regression.Bgp;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// BGP checkpoints and the early stop. A run checkpointed partway and resumed gives exactly the uninterrupted run's
/// draws and summaries; a checkpoint is independent of the live run; another problem's checkpoint (data, noise variance,
/// prior, test points or proposal) is refused; and an R̂ target ends the run once the bandwidth has converged.
/// </summary>
public sealed class BgpResumeTests
{
    private static (double[,] X, double[] Y, double[,] XTest) Data(double shift = 0.0)
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
        y[5] += shift;
        var xTest = new double[25, 1];
        for (int i = 0; i < 25; i++) xTest[i, 0] = i / 24.0;
        return (x, y, xTest);
    }

    private static BgpSampler Sampler(double[,] x, double[] y, double sigma2 = 0.01, double a0 = 1.0)
    {
        var gp = new GpRegression(x, y, sigma2, new SquaredExponentialKernel());
        return new BgpSampler(gp, new EmpiricalBayesBandwidthPrior(gp, a0: a0));
    }

    [Fact]
    public void ResumedRun_MatchesTheUninterruptedRun()
    {
        (double[,] x, double[] y, double[,] xTest) = Data();
        BgpResult straight = Sampler(x, y).Run(xTest, draws: 500, burn: 100, proposalSd: 0.5, seed: 5, chains: 3);

        BgpResult first = Sampler(x, y).Run(xTest, draws: 200, burn: 100, proposalSd: 0.5, seed: 5, chains: 3, checkpoint: true);
        EnsembleCheckpoint checkpoint = first.Checkpoint!;
        Assert.Equal(200, checkpoint.DrawsPerChain);
        BgpResult resumed = Sampler(x, y).Run(xTest, draws: 300, burn: 100, proposalSd: 0.5, seed: 5, chains: 3, resumeFrom: checkpoint);
        AssertSame(straight, resumed);

        BgpResult again = Sampler(x, y).Run(xTest, draws: 300, burn: 100, proposalSd: 0.5, seed: 5, chains: 3, resumeFrom: checkpoint);
        AssertSame(straight, again);   // resuming did not change the checkpoint
    }

    [Fact]
    public void Checkpoint_OfAnotherProblem_IsRefused()
    {
        (double[,] x, double[] y, double[,] xTest) = Data();
        EnsembleCheckpoint checkpoint = Sampler(x, y)
            .Run(xTest, draws: 100, burn: 50, proposalSd: 0.5, seed: 5, chains: 2, checkpoint: true).Checkpoint!;

        (double[,] x2, double[] y2, _) = Data(shift: 1e-9);
        var otherTest = (double[,])xTest.Clone();
        otherTest[3, 0] += 1e-9;
        Assert.Throws<ArgumentException>(() => Sampler(x2, y2).Run(xTest, draws: 100, burn: 50, proposalSd: 0.5, seed: 5, chains: 2, resumeFrom: checkpoint));
        Assert.Throws<ArgumentException>(() => Sampler(x, y, sigma2: 0.02).Run(xTest, draws: 100, burn: 50, proposalSd: 0.5, seed: 5, chains: 2, resumeFrom: checkpoint));
        Assert.Throws<ArgumentException>(() => Sampler(x, y, a0: 2.0).Run(xTest, draws: 100, burn: 50, proposalSd: 0.5, seed: 5, chains: 2, resumeFrom: checkpoint));
        Assert.Throws<ArgumentException>(() => Sampler(x, y).Run(otherTest, draws: 100, burn: 50, proposalSd: 0.5, seed: 5, chains: 2, resumeFrom: checkpoint));
        Assert.Throws<ArgumentException>(() => Sampler(x, y).Run(xTest, draws: 100, burn: 50, proposalSd: 0.4, seed: 5, chains: 2, resumeFrom: checkpoint));
    }

    [Fact]
    public void EarlyStop_EndsOnceTheBandwidthConverges()
    {
        (double[,] x, double[] y, double[,] xTest) = Data();
        const int chains = 4, cap = 20000;
        BgpResult r = Sampler(x, y).Run(xTest, draws: cap, burn: 100, proposalSd: 0.5, seed: 9, chains: chains,
                                        rHatTarget: 1.01, batchSize: 250);

        int perChain = r.Draws / chains;
        Assert.True(perChain < cap, $"the run should stop before its cap ({perChain} of {cap} draws per chain)");
        foreach (IReadOnlyList<BgpDraw> chain in r.BandwidthDraws) Assert.Equal(perChain, chain.Count);
        Assert.True(r.BandwidthConvergence.RHat.MeetsTarget(1.01), $"R̂ {r.BandwidthConvergence.RHat}");
        Assert.True(r.BandwidthConvergence.Ess.MeetsMinimum(ChainEnsemble.MinEssPerChain * chains), $"ESS {r.BandwidthConvergence.Ess}");
    }

    private static void AssertSame(BgpResult a, BgpResult b)
    {
        Assert.Equal(a.Draws, b.Draws);
        Assert.Equal(a.PosteriorMean, b.PosteriorMean);
        Assert.Equal(a.PosteriorSd, b.PosteriorSd);
        Assert.Equal(a.Bandwidth, b.Bandwidth);
        Assert.Equal(a.BandwidthOriginal, b.BandwidthOriginal);
        Assert.Equal(a.BandwidthConvergence, b.BandwidthConvergence);
        VerdictAssert.Same(a.Verdict, b.Verdict);
        Assert.Equal(a.AcceptanceRate, b.AcceptanceRate);
        for (int c = 0; c < a.BandwidthDraws.Count; c++)
        {
            Assert.Equal(a.BandwidthDraws[c].Count, b.BandwidthDraws[c].Count);
            for (int d = 0; d < a.BandwidthDraws[c].Count; d++)
            {
                BgpDraw x = a.BandwidthDraws[c][d], w = b.BandwidthDraws[c][d];
                Assert.Equal(x.Bandwidth, w.Bandwidth);
                Assert.Equal(x.Fit.LogMarginal, w.Fit.LogMarginal);
                Assert.Equal(x.Fit.Alpha, w.Fit.Alpha);
                Assert.Equal(x.ConditionalVariance, w.ConditionalVariance);
            }
        }
    }
}
