using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// Checkpoints (step 3.4): a BARS run checkpointed partway and resumed in a fresh ensemble gives exactly the result of
/// the run that never stopped, across the untempered chain, the robust weight block, the adaptive kernel and a tempered
/// ladder with burn-in adaptation. A checkpoint is independent of the live run, and resuming another problem's
/// checkpoint is refused.
/// </summary>
public sealed class BarsResumeTests
{
    private static (double[] X, double[] Y) Data()
    {
        var rng = new Xoshiro256PlusPlus(seed: 19);
        int n = 80;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
            y[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.15 * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2)
                 + (i % 13 == 0 ? 1.5 : 0.0);
        }
        return (x, y);
    }

    public static TheoryData<string> Setups => new() { "plain", "robust", "adaptive", "tempered" };

    private static BarsEnsemble Ensemble(string setup) => setup switch
    {
        "robust" => new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0),
                                     resampler: new StudentTWeights(4.0)),
        "adaptive" => new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0),
                                       adaptKernelTau: 20.0),
        "tempered" => new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0),
                                       temperLevels: 3),
        _ => new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0)),
    };

    [Theory]
    [MemberData(nameof(Setups))]
    public void ResumedRun_MatchesTheUninterruptedRun(string setup)
    {
        (double[] x, double[] y) = Data();
        BarsResult straight = Ensemble(setup).Run(x, y, grid: x, chains: 2, masterSeed: 6, burn: 300, samples: 500);

        BarsResult first = Ensemble(setup).Run(x, y, grid: x, chains: 2, masterSeed: 6, burn: 300, samples: 200, checkpoint: true);
        EnsembleCheckpoint checkpoint = first.Checkpoint!;
        Assert.Equal(200, checkpoint.DrawsPerChain);
        BarsResult resumed = Ensemble(setup).Run(x, y, grid: x, chains: 2, masterSeed: 6, burn: 300, samples: 300,
                                                 resumeFrom: checkpoint);
        AssertSame(straight, resumed);

        BarsResult again = Ensemble(setup).Run(x, y, grid: x, chains: 2, masterSeed: 6, burn: 300, samples: 300,
                                               resumeFrom: checkpoint);
        AssertSame(straight, again);   // resuming did not change the checkpoint
    }

    [Fact]
    public void Checkpoint_OfAnotherProblem_IsRefused()
    {
        (double[] x, double[] y) = Data();
        EnsembleCheckpoint checkpoint = Ensemble("plain").Run(x, y, grid: x, chains: 2, masterSeed: 6, burn: 100, samples: 100,
                                                              checkpoint: true).Checkpoint!;
        var other = (double[])y.Clone();
        other[3] += 1e-9;
        Assert.Throws<ArgumentException>(() => Ensemble("plain").Run(x, other, grid: x, chains: 2, masterSeed: 6, burn: 100,
                                                                     samples: 100, resumeFrom: checkpoint));
        Assert.Throws<ArgumentException>(() => Ensemble("tempered").Run(x, y, grid: x, chains: 2, masterSeed: 6, burn: 100,
                                                                        samples: 100, resumeFrom: checkpoint));
    }

    private static void AssertSame(BarsResult a, BarsResult b)
    {
        Assert.Equal(a.SamplesUsed, b.SamplesUsed);
        Assert.Equal(a.Fit, b.Fit);
        Assert.Equal(a.FitSd, b.FitSd);
        Assert.Equal(a.MeanKnots, b.MeanKnots);
        Assert.Equal(a.Knots, b.Knots);
        Assert.Equal(a.RHatFit, b.RHatFit);
        Assert.Equal(a.PeakCountMean, b.PeakCountMean);
        Assert.Equal(a.PeakIntensity, b.PeakIntensity);
        Assert.Equal(a.SpanCoverage, b.SpanCoverage);
        Assert.Equal(a.Peak, b.Peak);
        Assert.Equal(a.PeakModes, b.PeakModes);
        Assert.Equal(a.PeakDraws.Count, b.PeakDraws.Count);
        for (int c = 0; c < a.PeakDraws.Count; c++)
        {
            Assert.Equal(a.PeakDraws[c].Count, b.PeakDraws[c].Count);
            for (int d = 0; d < a.PeakDraws[c].Count; d++)
            {
                Assert.Equal(a.PeakDraws[c][d].ArgmaxLocation, b.PeakDraws[c][d].ArgmaxLocation);
                Assert.Equal(a.PeakDraws[c][d].ArgmaxHeight, b.PeakDraws[c][d].ArgmaxHeight);
                Assert.Equal(a.PeakDraws[c][d].Spans, b.PeakDraws[c][d].Spans);
            }
        }
        Assert.Equal(a.AcceptanceRate, b.AcceptanceRate);
        Assert.Equal(a.MoveStats, b.MoveStats);
        Assert.Equal(a.Evaluations.Entries, b.Evaluations.Entries);
        Assert.Equal(a.Tempering is null, b.Tempering is null);
        if (a.Tempering is TemperingReadout ta && b.Tempering is TemperingReadout tb)
        {
            Assert.Equal(ta.RoundTripsPerScan, tb.RoundTripsPerScan);
            Assert.Equal(ta.PairRejection, tb.PairRejection);
            for (int c = 0; c < ta.Ladders.Length; c++) Assert.Equal(ta.Ladders[c], tb.Ladders[c]);
        }
    }
}
