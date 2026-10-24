using System;
using System.Collections.Generic;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// BAPS on the ensemble shell: rank-normalized diagnostics on log λ and log σ², an early stop gated on them, and
/// checkpoints. A run checkpointed partway and resumed gives exactly the uninterrupted run's draws and summaries in
/// every mode, a checkpoint is independent of the live run, and another problem's or mode's checkpoint is refused.
/// </summary>
public sealed class BapsEnsembleTests
{
    private static (double[,] Z, double[] Y) Data(double shift = 0.0)
    {
        var rng = new Xoshiro256PlusPlus(seed: 31);
        int n = 60;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
            y[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.1 * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }
        y[3] += shift;
        var knots = new double[12];
        for (int k = 0; k < knots.Length; k++) knots[k] = (k + 1.0) / (knots.Length + 1);
        return (new SplineBasis(3).Design(new KnotConfig(knots), x), y);
    }

    public static TheoryData<BapsScheme> Modes => new()
    {
        BapsScheme.Gibbs, BapsScheme.Collapsed, BapsScheme.LocallyAdaptiveGibbs,
    };

    [Theory]
    [MemberData(nameof(Modes))]
    public void ResumedRun_MatchesTheUninterruptedRun(BapsScheme mode)
    {
        (double[,] z, double[] y) = Data();
        var sampler = new Baps(z, y, new DifferencePenalty(2));
        BapsResult straight = sampler.Run(mode, chains: 3, masterSeed: 5, burn: 200, samples: 500);

        BapsResult first = sampler.Run(mode, chains: 3, masterSeed: 5, burn: 200, samples: 200, checkpoint: true);
        EnsembleCheckpoint checkpoint = first.Checkpoint!;
        Assert.Equal(200, checkpoint.DrawsPerChain);
        BapsResult resumed = sampler.Run(mode, chains: 3, masterSeed: 5, burn: 200, samples: 300, resumeFrom: checkpoint);
        AssertSame(straight, resumed);

        BapsResult again = sampler.Run(mode, chains: 3, masterSeed: 5, burn: 200, samples: 300, resumeFrom: checkpoint);
        AssertSame(straight, again);   // resuming did not change the checkpoint
    }

    [Fact]
    public void Checkpoint_OfAnotherProblemOrMode_IsRefused()
    {
        (double[,] z, double[] y) = Data();
        EnsembleCheckpoint checkpoint = new Baps(z, y, new DifferencePenalty(2))
            .Run(BapsScheme.Gibbs, chains: 2, masterSeed: 5, burn: 100, samples: 100, checkpoint: true).Checkpoint!;

        (double[,] z2, double[] y2) = Data(shift: 1e-9);
        Assert.Throws<ArgumentException>(() => new Baps(z2, y2, new DifferencePenalty(2))
            .Run(BapsScheme.Gibbs, chains: 2, masterSeed: 5, burn: 100, samples: 100, resumeFrom: checkpoint));
        Assert.Throws<ArgumentException>(() => new Baps(z, y, new DifferencePenalty(2))
            .Run(BapsScheme.LocallyAdaptiveGibbs, chains: 2, masterSeed: 5, burn: 100, samples: 100, resumeFrom: checkpoint));
        Assert.Throws<ArgumentException>(() => new Baps(z, y, new DifferencePenalty(2), phi: 2.0)
            .Run(BapsScheme.Gibbs, chains: 2, masterSeed: 5, burn: 100, samples: 100, resumeFrom: checkpoint));
    }

    [Theory]
    [MemberData(nameof(Modes))]
    public void EarlyStop_EndsOnceLogLambdaAndLogSigma2Converge(BapsScheme mode)
    {
        (double[,] z, double[] y) = Data();
        const int chains = 4, cap = 20000;
        BapsResult r = new Baps(z, y, new DifferencePenalty(2))
            .Run(mode, chains: chains, masterSeed: 9, burn: 300, samples: cap, rHatTarget: 1.01, batchSize: 250);

        int perChain = r.Draws / chains;
        Assert.True(perChain < cap, $"the run should stop before its cap ({perChain} of {cap} draws per chain)");
        foreach (IReadOnlyList<BapsDraw> chain in r.PosteriorDraws) Assert.Equal(perChain, chain.Count);
        foreach (ConvergenceDiagnostic d in new[] { r.LogLambdaConvergence, r.LogSigma2Convergence })
        {
            Assert.True(d.RHat.MeetsTarget(1.01), $"R̂ {d.RHat}");
            Assert.True(d.Ess.MeetsMinimum(ChainEnsemble.MinEssPerChain * chains), $"ESS {d.Ess}");
        }
    }

    [Theory]
    [MemberData(nameof(Modes))]
    public void Diagnostics_AndAcceptance_AreReported(BapsScheme mode)
    {
        (double[,] z, double[] y) = Data();
        BapsResult r = new Baps(z, y, new DifferencePenalty(2)).Run(mode, chains: 3, masterSeed: 2, burn: 200, samples: 600);

        Assert.Equal(3, r.ChainSeeds.Length);
        Assert.Null(r.Checkpoint);
        Assert.Equal(DiagnosticState.Estimated, r.LogLambdaConvergence.RHat.State);
        Assert.Equal(DiagnosticState.Estimated, r.LogSigma2Convergence.RHat.State);
        Assert.True(r.LogLambdaConvergence.BulkEss.Value > 0.0);
        if (mode == BapsScheme.Collapsed)
            Assert.InRange(r.AcceptanceRate, 0.05, 0.95);   // Metropolis on log λ
        else
            Assert.Equal(1.0, r.AcceptanceRate);              // every Gibbs sweep is accepted
    }

    private static void AssertSame(BapsResult a, BapsResult b)
    {
        Assert.Equal(a.Draws, b.Draws);
        Assert.Equal(a.Coefficients, b.Coefficients);
        Assert.Equal(a.LambdaMedian, b.LambdaMedian);
        Assert.Equal(a.LambdaLo, b.LambdaLo);
        Assert.Equal(a.LambdaHi, b.LambdaHi);
        Assert.Equal(a.SigmaMean, b.SigmaMean);
        Assert.Equal(a.RHatLogLambda, b.RHatLogLambda);
        Assert.Equal(a.LogLambdaConvergence, b.LogLambdaConvergence);
        Assert.Equal(a.LogSigma2Convergence, b.LogSigma2Convergence);
        Assert.Equal(a.AcceptanceRate, b.AcceptanceRate);
        Assert.Equal(a.LocalSmoothing, b.LocalSmoothing);
        for (int c = 0; c < a.PosteriorDraws.Count; c++)
        {
            Assert.Equal(a.PosteriorDraws[c].Count, b.PosteriorDraws[c].Count);
            for (int d = 0; d < a.PosteriorDraws[c].Count; d++)
            {
                BapsDraw x = a.PosteriorDraws[c][d], w = b.PosteriorDraws[c][d];
                Assert.Equal(x.Coefficients, w.Coefficients);
                Assert.Equal(x.Lambda, w.Lambda);
                Assert.Equal(x.Sigma2, w.Sigma2);
                Assert.Equal(x.Tau2, w.Tau2);
                Assert.Equal(x.Eta, w.Eta);
            }
        }
    }
}
