using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;
using Xunit;
using static Maths.Samplers.Tests.ParallelTemperingBlockTests;

namespace Maths.Samplers.Tests;

/// <summary>
/// A ladder whose rungs step at once (<see cref="ParallelTempering{TState}.MaxDegreeOfParallelism"/>) runs as the
/// in-order ladder does, bit for bit: on the (k, θ) target of <see cref="ParallelTemperingBlockTests"/>, with each rung's
/// slice block and the reference rung's exact draw, on every swap schedule, through the burn-in's adaptations and a resize
/// and the traced phase, at many rungs and more workers than the machine has cores. The degree is not restart state, and
/// a rung's exception surfaces as the in-order loop throws it.
/// </summary>
public sealed class ParallelTemperingRungParallelismTests
{
    private const int Rungs = 32;

    public static TheoryData<SwapSchedule, int> Cases => new()
    {
        { SwapSchedule.DeterministicEvenOdd, 2 },
        { SwapSchedule.DeterministicEvenOdd, 7 },
        { SwapSchedule.DeterministicEvenOdd, 64 },
        { SwapSchedule.StochasticEvenOdd, 32 },
        { SwapSchedule.Sequential, 32 },
    };

    [Theory]
    [MemberData(nameof(Cases))]
    public void RungsAtOnce_RunAsInOrder_BitForBit(SwapSchedule schedule, int degree)
    {
        List<string> inOrder = Run(schedule, degree: 1), atOnce = Run(schedule, degree);
        AssertSameRecord(inOrder, atOnce);
    }

    [Fact]
    public void TheDegree_IsNotRestartState_AndAResumeAtAnotherDegreeContinuesBitForBit()
    {
        var target = new IndexAndCoordinate();
        ParallelTempering<State> whole = Ladder(target, SwapSchedule.DeterministicEvenOdd, degree: 8);
        for (int t = 0; t < 400; t++) whole.Step();
        LadderSnapshot<State> snapshot = whole.Snapshot();
        ParallelTempering<State> restored = ParallelTempering<State>.Restore(_ => Moves, target, snapshot, new PriorDraw(target),
                                                                             _ => [new CoordinateSlice(target)]);
        Assert.Equal(1, restored.MaxDegreeOfParallelism);
        restored.MaxDegreeOfParallelism = 3;
        for (int t = 0; t < 400; t++)
            Assert.Equal(Describe(whole.Step()), Describe(restored.Step()));
        AssertSameRecord(Final(whole), Final(restored));

        Assert.Throws<ArgumentOutOfRangeException>(() => whole.MaxDegreeOfParallelism = 0);
        whole.MaxDegreeOfParallelism = 1;
        Assert.Equal(1, whole.MaxDegreeOfParallelism);
    }

    [Fact]
    public void ARungsException_SurfacesAsTheInOrderLoopThrowsIt()
    {
        // Rungs 3 and 9 throw on their fifth step: in order, rung 3's exception ends the scan; at once, rung 3's still
        // surfaces, its own type and message, not the fan-out's aggregate.
        foreach (int degree in new[] { 1, 4, Rungs })
        {
            var target = new IndexAndCoordinate();
            var pt = new ParallelTempering<State>(i => i switch
            {
                3 => [new Throwing(() => new InvalidOperationException("rung 3"))],
                9 => [new Throwing(() => new ArgumentException("rung 9"))],
                _ => Moves,
            }, target, new State(0, 0.0), Ladder(), TemperingPath.Likelihood, masterSeed: 41) { MaxDegreeOfParallelism = degree };
            for (int t = 0; t < 4; t++) pt.Step();
            var thrown = Assert.Throws<InvalidOperationException>(() => pt.Step());
            Assert.Equal("rung 3", thrown.Message);
            Assert.Equal(4, pt.Scans);
        }
    }

    // 32 rungs from 1 down to b = 0 on the likelihood path, the last the reference rung.
    private static double[] Ladder() => ParallelTempering<State>.ReferenceLadder(Rungs, hotEnd: 0.01);

    private static ParallelTempering<State> Ladder(IndexAndCoordinate target, SwapSchedule schedule, int degree)
        => new(_ => Moves, target, new State(0, 0.0), Ladder(), TemperingPath.Likelihood, masterSeed: 23, schedule, new PriorDraw(target),
               _ => [new CoordinateSlice(target)])
        {
            MaxDegreeOfParallelism = degree,
        };

    // The ladder's run as a tempered fit drives it: burn-in rounds of 1, 2, 4, … scans, each followed by an adaptation, with
    // one resize to 40 rungs placed on the barrier; then the rungs traced. Every scan's cold state and swap count, then
    // every replica's state, density, moves and blocks, the pairs' statistics, the round trips and the traces.
    private static List<string> Run(SwapSchedule schedule, int degree)
    {
        var target = new IndexAndCoordinate();
        ParallelTempering<State> pt = Ladder(target, schedule, degree);
        var record = new List<string>();
        void Scan() => record.Add($"{Describe(pt.Step())} {pt.SwapAccepts}");
        for (int round = 1; round <= 256; round *= 2)
        {
            for (int t = 0; t < round; t++) Scan();
            if (round == 32)
            {
                double[] rejection = pt.RoundRejection, now = pt.InverseTemperatures;
                double[] next = TemperingLadder.Place(now, rejection, 40);
                pt.Resize(next, TemperingLadder.Sources(now, rejection, next), s => s, _ => Moves, _ => [new CoordinateSlice(target)]);
                Assert.Equal(degree, pt.MaxDegreeOfParallelism);
            }
            else pt.AdaptLadder();
        }
        pt.TraceRungs();
        for (int t = 0; t < 600; t++) Scan();
        record.AddRange(Final(pt));
        return record;
    }

    private static List<string> Final(ParallelTempering<State> pt)
    {
        var record = new List<string>();
        for (int i = 0; i < pt.Replicas.Count; i++)
        {
            ReversibleJumpChain<State> r = pt.Replicas[i];
            record.Add($"rung {i}: {Describe(r.Current)} {Bits(r.CurrentLogDensity)} {Bits(r.Temperature.InverseTemperature)} " +
                       string.Join(" ", r.MoveStats().Select(m => $"{m.Key} {m.Attempts}/{m.Accepted}")) + " " +
                       string.Join(" ", pt.BlockStats(i).Select(m => $"{m.Key} {m.Attempts}/{m.Accepted}")));
            if (pt.RungTraces is { } traces)
                record.Add($"trace {i}: " + string.Join(" ", traces[i].Draws(RungFunctional.Slope).Concat(traces[i].Draws(RungFunctional.LogRatio)).Select(Bits)));
        }
        foreach (SwapPairStat p in pt.PairStats) record.Add($"pair {p.Attempts}/{p.Accepts} {Bits(p.MeanRejection)}");
        record.Add($"scans {pt.Scans}, swaps {pt.SwapAttempts}/{pt.SwapAccepts}, round trips {pt.RoundTrips}, restarts {pt.AnnealedRestarts}");
        return record;
    }

    private static string Describe(State s) => $"{s.K} {Bits(s.Theta)}";

    private static string Bits(double x) => BitConverter.DoubleToInt64Bits(x).ToString("x16");

    private static void AssertSameRecord(List<string> expected, List<string> actual)
    {
        Assert.Equal(expected.Count, actual.Count);
        for (int i = 0; i < expected.Count; i++)
            if (expected[i] != actual[i]) Assert.Fail($"Entry {i} differs: in order '{expected[i]}', at once '{actual[i]}'.");
    }

    // A move that throws on its fifth proposal: one per replica, so its count is the rung's own.
    private sealed class Throwing(Func<Exception> failure) : IRjMove<State>
    {
        private int _proposals;
        public string Key => "stay";
        public string ReverseKey => "stay";
        public double Weight(State state) => 1.0;
        public Proposal<State>? Propose(State current, Xoshiro256PlusPlus rng) => ++_proposals == 5 ? throw failure() : new(current, 0.0, 0.0);
    }
}
