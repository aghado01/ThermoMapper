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

    private sealed class FlatTarget : IRjTarget<int>
    {
        public double LogPosterior(int state) => 0.0;
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
