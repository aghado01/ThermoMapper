using System;
using System.Collections.Generic;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Samplers.Tests;

/// <summary>
/// Restart state (step 3.4): a chain or ladder restored from a snapshot continues exactly as the live one does, and a
/// snapshot is independent of the live sampler. The ensemble shell refuses to checkpoint a chain run that cannot save
/// its state, or to resume a model that cannot restore one.
/// </summary>
public sealed class SnapshotTests
{
    // −0.3·|state − 3| over the integers, explored by unit steps.
    private static readonly ISamplingTarget<int> Peaked = SamplingTarget.FromLogDensity<int>(s => -0.3 * Math.Abs(s - 3));
    private static IRjMove<int>[] Moves(int replica) => [new UnitStep("up", "down", +1), new UnitStep("down", "up", -1)];

    [Fact]
    public void Chain_RestoredFromASnapshot_ContinuesIdentically()
    {
        var live = new ReversibleJumpChain<int>(Moves(0), Peaked, 0, new Xoshiro256PlusPlus(5), new Temperature(0.7, TemperingPath.Likelihood));
        for (int i = 0; i < 300; i++) live.Step();
        ChainSnapshot<int> snapshot = live.Snapshot();
        var restored = ReversibleJumpChain<int>.Restore(Moves(0), Peaked, snapshot);

        for (int i = 0; i < 500; i++) Assert.Equal(live.Step(), restored.Step());
        Assert.Equal(live.Attempts, restored.Attempts);
        Assert.Equal(live.Accepted, restored.Accepted);
        Assert.Equal(live.MoveStats(), restored.MoveStats());
        Assert.Equal(live.Temperature, restored.Temperature);
        Assert.Equal(300, snapshot.MoveAttempts[0] + snapshot.MoveAttempts[1]);   // the snapshot did not move on
    }

    [Theory]
    [InlineData(SwapSchedule.DeterministicEvenOdd)]
    [InlineData(SwapSchedule.StochasticEvenOdd)]
    public void Ladder_RestoredFromASnapshot_ContinuesIdentically(SwapSchedule schedule)
    {
        var live = new ParallelTempering<int>(Moves, Peaked, 0, [1.0, 0.6, 0.3, 0.1], TemperingPath.Likelihood, masterSeed: 8, schedule);
        for (int round = 1; round <= 64; round *= 2)
        {
            for (int i = 0; i < round; i++) live.Step();
            live.AdaptLadder();
        }
        for (int i = 0; i < 37; i++) live.Step();   // mid-round, with statistics gathered since the last adaptation
        var restored = ParallelTempering<int>.Restore(Moves, Peaked, live.Snapshot());

        for (int i = 0; i < 400; i++) Assert.Equal(live.Step(), restored.Step());
        Assert.Equal(live.InverseTemperatures, restored.InverseTemperatures);
        Assert.Equal(live.RoundTrips, restored.RoundTrips);
        Assert.Equal(live.AnnealedRestarts, restored.AnnealedRestarts);
        Assert.Equal(live.SwapAccepts, restored.SwapAccepts);
        Assert.Equal(live.PairStats, restored.PairStats);
        Assert.Equal(live.AdaptLadder(), restored.AdaptLadder());   // the post-adaptation statistics were restored too
        Assert.Equal(live.InverseTemperatures, restored.InverseTemperatures);
        for (int r = 0; r < 4; r++) Assert.Equal(live.Replicas[r].Current, restored.Replicas[r].Current);
    }

    [Fact]
    public void Shell_RefusesWhatItCannotSaveOrRestore()
    {
        var model = new PlainModel();
        Assert.Throws<InvalidOperationException>(() =>
            new ChainEnsemble().Run(model, chains: 2, masterSeed: 1, burn: 0, samples: 10, checkpoint: true));
        var checkpoint = new EnsembleCheckpoint(new RunIdentity(1, 2, ""), [1, 2], 10,
            [new ChainCheckpoint(new NoState(), [], new ChainTrace(1), 0, 0), new ChainCheckpoint(new NoState(), [], new ChainTrace(1), 0, 0)]);
        Assert.Throws<NotSupportedException>(() =>
            new ChainEnsemble().Run(model, chains: 2, masterSeed: 1, burn: 0, samples: 10, resumeFrom: checkpoint));
        Assert.Throws<ArgumentException>(() =>
            new ChainEnsemble().Run(model, chains: 2, masterSeed: 2, burn: 0, samples: 10, resumeFrom: checkpoint));
    }

    private sealed class NoState : IChainState { }

    private sealed class PlainModel : IEnsembleModel<int>
    {
        public int FunctionalDim => 1;
        public int MonitoredDim => 1;
        public IChainRun<int> StartChain(int chainSeed, int readoutSeed)
            => new PlainRun(new ReversibleJumpChain<int>(Moves(0), Peaked, 0, new Xoshiro256PlusPlus(chainSeed)));
    }

    private sealed class PlainRun : IChainRun<int>, IChain<int>
    {
        private readonly ReversibleJumpChain<int> _chain;
        public PlainRun(ReversibleJumpChain<int> chain) => _chain = chain;
        public IChain<int> Chain => this;
        public int Step() => _chain.Step();
        public long Accepted => _chain.Accepted;
        public long Attempts => _chain.Attempts;
        public void Burn(int steps) { for (int i = 0; i < steps; i++) _chain.Step(); }
        public void Accumulate(in int draw, Span<double> functionals) => functionals[0] = draw;
    }

    private sealed class UnitStep : IRjMove<int>
    {
        private readonly int _delta;
        public UnitStep(string key, string reverseKey, int delta) { Key = key; ReverseKey = reverseKey; _delta = delta; }
        public string Key { get; }
        public string ReverseKey { get; }
        public double Weight(int state) => 1.0;
        public Proposal<int>? Propose(int current, Xoshiro256PlusPlus rng) => new(current + _delta, 0.0, 0.0);
    }
}
