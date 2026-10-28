using System;
using System.Threading;
using System.Threading.Tasks;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Samplers.Tests;

/// <summary>
/// The ensemble shell's parallelism bound: a run's result does not depend on how many chains run at once, and the bound's
/// cancellation stops the run.
/// </summary>
public sealed class ChainEnsembleParallelTests
{
    [Fact]
    public void Run_DoesNotDependOnTheParallelismBound()
    {
        EnsembleRun Run(int degree) => new ChainEnsemble().Run(new WalkModel(), chains: 5, masterSeed: 7, burn: 50, samples: 400,
            parallel: new ParallelOptions { MaxDegreeOfParallelism = degree }).Summary;
        EnsembleRun reference = new ChainEnsemble().Run(new WalkModel(), chains: 5, masterSeed: 7, burn: 50, samples: 400).Summary;
        foreach (int degree in new[] { 1, 2, -1 })
        {
            EnsembleRun run = Run(degree);
            Assert.Equal(reference.FunctionalMean, run.FunctionalMean);
            Assert.Equal(reference.ChainSeeds, run.ChainSeeds);
            for (int c = 0; c < 5; c++) Assert.Equal(reference.Traces[c].Draws(0), run.Traces[c].Draws(0));
        }
    }

    [Fact]
    public void Run_StopsOnCancellation()
    {
        using var cancelled = new CancellationTokenSource();
        cancelled.Cancel();
        Assert.Throws<OperationCanceledException>(() => new ChainEnsemble().Run(new WalkModel(), chains: 3, masterSeed: 1, burn: 10,
            samples: 10, parallel: new ParallelOptions { CancellationToken = cancelled.Token }));
    }

    // A Gaussian random walk per chain, monitored on its position.
    private sealed class WalkModel : IEnsembleModel<double>
    {
        public int FunctionalDim => 1;
        public int MonitoredDim => 1;
        public IChainRun<double> StartChain(int chainSeed, int readoutSeed) => new Walk(new Xoshiro256PlusPlus(chainSeed));
    }

    private sealed class Walk(Xoshiro256PlusPlus rng) : IChainRun<double>, IChain<double>
    {
        private double _x;
        public IChain<double> Chain => this;
        public long Accepted { get; private set; }
        public long Attempts { get; private set; }
        public void Burn(int steps) { for (int i = 0; i < steps; i++) Step(); }
        public void Accumulate(in double draw, Span<double> functionals) => functionals[0] = draw;

        public double Step()
        {
            Attempts++;
            Accepted++;
            return _x = 0.9 * _x + rng.NextDouble() - 0.5;
        }
    }
}
