using System;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Samplers.Tests;

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
            new ReversibleJumpChain<int>(moves, Flat, 0, new Xoshiro256PlusPlus(1)));
    }

    [Fact]
    public void Palette_WithReciprocalAndSelfReverseMoves_IsAccepted()
    {
        IRjMove<int>[] moves = [new StepMove("up", "down", +1), new StepMove("down", "up", -1), new StepMove("stay", "stay", 0)];
        var chain = new ReversibleJumpChain<int>(moves, Flat, 0, new Xoshiro256PlusPlus(1));
        Assert.Equal(0, chain.Current);
    }

    [Fact]
    public void ReplaceTarget_KeepsStateAndCounters_AndReevaluatesTheCachedTarget()
    {
        IRjMove<int>[] moves = [new StepMove("up", "down", +1), new StepMove("down", "up", -1)];
        var chain = new ReversibleJumpChain<int>(moves, Flat, 0, new Xoshiro256PlusPlus(3));
        for (int i = 0; i < 50; i++) chain.Step();

        int state = chain.Current;
        long attempts = chain.Attempts, accepted = chain.Accepted;
        MoveStat[] before = [.. chain.MoveStats()];

        chain.ReplaceTarget(Linear(0.5));

        Assert.Equal(state, chain.Current);
        Assert.Equal(0.5 * state, chain.CurrentLogDensity);
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
        var pt = new ParallelTempering<int>(MovesFor, Flat, 0, betas, masterSeed: 2);
        pt.Step();

        Assert.Equal(new[] { 0, 1, 2, 3 }, built.ConvertAll(b => b.Replica));
        for (int i = 0; i < built.Count; i++)
            for (int j = i + 1; j < built.Count; j++)
                Assert.NotSame(built[i].Moves[0], built[j].Moves[0]);
    }

    [Fact]
    public void ParallelTempering_SwapStream_DoesNotDependOnRoundingTies()
    {
        // A swap between states whose log targets tie up to rounding has a log A of rounding size, and rounding
        // decides its sign. Two ladders whose targets break the 2-vs-4 tie in opposite directions, by 10⁻¹⁵, must run
        // identically: the swap uniform is drawn on every attempt, so a tie's sign cannot shift the stream.
        IRjMove<int>[] MovesFor(int replica) => [new StepMove("up", "down", +1), new StepMove("down", "up", -1)];
        double[] betas = ParallelTempering<int>.GeometricLadder(levels: 4, betaMin: 0.05);
        var a = new ParallelTempering<int>(MovesFor, Peaked(tieBreak: 2), 0, betas, masterSeed: 9);
        var b = new ParallelTempering<int>(MovesFor, Peaked(tieBreak: 4), 0, betas, masterSeed: 9);
        for (int i = 0; i < 5000; i++) Assert.Equal(a.Step(), b.Step());
        Assert.Equal(a.SwapAccepts, b.SwapAccepts);
    }

    // −0.3·|state − 3|: a peaked target whose replicas meet at the mirror states 2 and 4, whose tie is broken by 10⁻¹⁵
    // in favour of tieBreak.
    private static ISamplingTarget<int> Peaked(int tieBreak)
        => SamplingTarget.FromLogDensity<int>(s => -0.3 * Math.Abs(s - 3) - 1.234567 + (s == tieBreak ? 1e-15 : 0.0));

    [Fact]
    public void Candidate_OffTheSupport_IsRejected_WhateverItsProposalRatio()
    {
        // A sampler that underflows onto the support's edge can report an infinite proposal ratio there; the
        // candidate's −∞ (or NaN) log target decides the move, not −∞ + ∞.
        IRjMove<int>[] moves = [new InfiniteRatioMove("up", "down", +1), new StepMove("down", "up", -1)];
        var chain = new ReversibleJumpChain<int>(moves, SupportOfZero, 0, new Xoshiro256PlusPlus(seed: 3));
        for (int i = 0; i < 1000; i++) Assert.Equal(0, chain.Step());
        Assert.Equal(0, chain.Accepted);
    }

    // Support {0}: −∞ above it, NaN below it.
    private static readonly ISamplingTarget<int> SupportOfZero
        = SamplingTarget.FromLogDensity<int>(s => s == 0 ? 0.0 : s > 0 ? double.NegativeInfinity : double.NaN);

    private sealed class InfiniteRatioMove : IRjMove<int>
    {
        private readonly int _delta;
        public InfiniteRatioMove(string key, string reverseKey, int delta) { Key = key; ReverseKey = reverseKey; _delta = delta; }
        public string Key { get; }
        public string ReverseKey { get; }
        public double Weight(int state) => 1.0;
        public Proposal<int>? Propose(int current, Xoshiro256PlusPlus rng) => new(current + _delta, double.PositiveInfinity, 0.0);
    }

    private static readonly ISamplingTarget<int> Flat = SamplingTarget.FromLogDensity<int>(_ => 0.0);

    private static ISamplingTarget<int> Linear(double slope) => SamplingTarget.FromLogDensity<int>(s => slope * s);

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
