using System;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// Reversible-jump chain contracts: move pairings must be reciprocal, and replacing the target keeps the chain's
/// state and cumulative counters while re-evaluating the cached log target.
/// </summary>
public sealed class ChainContractTests
{
    [Fact]
    public void Palette_WithNonReciprocalReverse_IsRejected()
    {
        // "up" names "down" as its reverse, but "down" names itself.
        IRjMove<int>[] moves = [new StepMove("up", "down", +1), new StepMove("down", "down", -1)];
        Assert.Throws<ArgumentException>(() =>
            new ReversibleJumpChain<int>(moves, new FlatTarget(), 0, new Xoshiro256PlusPlus(1)));
    }

    [Fact]
    public void Palette_WithReciprocalAndSelfReverseMoves_IsAccepted()
    {
        IRjMove<int>[] moves = [new StepMove("up", "down", +1), new StepMove("down", "up", -1), new StepMove("stay", "stay", 0)];
        var chain = new ReversibleJumpChain<int>(moves, new FlatTarget(), 0, new Xoshiro256PlusPlus(1));
        Assert.Equal(0, chain.Current);
    }

    [Fact]
    public void ReplaceTarget_KeepsStateAndCounters_AndReevaluatesTheCachedTarget()
    {
        IRjMove<int>[] moves = [new StepMove("up", "down", +1), new StepMove("down", "up", -1)];
        var chain = new ReversibleJumpChain<int>(moves, new FlatTarget(), 0, new Xoshiro256PlusPlus(3));
        for (int i = 0; i < 50; i++) chain.Step();

        int state = chain.Current;
        long attempts = chain.Attempts, accepted = chain.Accepted;
        MoveStat[] before = [.. chain.MoveStats()];

        chain.ReplaceTarget(new LinearTarget(0.5));

        Assert.Equal(state, chain.Current);
        Assert.Equal(0.5 * state, chain.CurrentLogTarget);
        Assert.Equal(attempts, chain.Attempts);
        Assert.Equal(accepted, chain.Accepted);
        Assert.Equal(before, chain.MoveStats());

        chain.Step();
        Assert.Equal(attempts + 1, chain.Attempts);
    }

    [Fact]
    public void ParallelTempering_BuildsEachReplicaItsOwnMoves()
    {
        var built = new System.Collections.Generic.List<(int Replica, IRjMove<int>[] Moves)>();
        IRjMove<int>[] MovesFor(int replica)
        {
            IRjMove<int>[] moves = [new StepMove("up", "down", +1), new StepMove("down", "up", -1)];
            built.Add((replica, moves));
            return moves;
        }

        double[] betas = ParallelTempering<int>.GeometricLadder(levels: 4, betaMin: 0.1);
        var pt = new ParallelTempering<int>(MovesFor, new FlatTarget(), 0, betas, masterSeed: 2);
        pt.Step();

        Assert.Equal(new[] { 0, 1, 2, 3 }, built.ConvertAll(b => b.Replica));
        for (int i = 0; i < built.Count; i++)
            for (int j = i + 1; j < built.Count; j++)
                Assert.NotSame(built[i].Moves[0], built[j].Moves[0]);
    }

    private sealed class FlatTarget : IRjTarget<int>
    {
        public double LogPosterior(int state) => 0.0;
    }

    private sealed class LinearTarget : IRjTarget<int>
    {
        private readonly double _slope;
        public LinearTarget(double slope) => _slope = slope;
        public double LogPosterior(int state) => _slope * state;
    }

    private sealed class StepMove : IRjMove<int>
    {
        private readonly int _delta;
        public StepMove(string key, string reverseKey, int delta) { Key = key; ReverseKey = reverseKey; _delta = delta; }
        public string Key { get; }
        public string ReverseKey { get; }
        public double Weight(int state) => 1.0;
        public Proposal<int>? Propose(int current, Xoshiro256PlusPlus rng) => new(current + _delta, 0.0, 0.0);
    }
}
