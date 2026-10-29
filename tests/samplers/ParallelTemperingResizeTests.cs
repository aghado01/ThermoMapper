using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;
using static Maths.Samplers.Tests.ParallelTemperingBlockTests;

namespace Maths.Samplers.Tests;

/// <summary>
/// A ladder's burn-in resize (<see cref="ParallelTempering{TState}.Resize"/>), on <see cref="ParallelTemperingBlockTests"/>'
/// target, whose replicas run moves on an index and a slice block on a coordinate: a ladder that never resizes is bit for
/// bit what it was before the resize existed; kept replicas keep their states and streams, and new ones take theirs from
/// the growth stream; after a resize every rung samples its own law, against quadrature; snapshots taken before and after
/// a resize continue bit for bit; and what cannot be resized is refused, the ladder unchanged.
/// </summary>
public sealed class ParallelTemperingResizeTests
{
    private readonly ITestOutputHelper _out;
    public ParallelTemperingResizeTests(ITestOutputHelper output) => _out = output;

    // Recorded at 2bdc22a, before the resize, the growth stream and the slope switch existed.
    [Theory]
    [InlineData(false, false, "6d304a81156cfe28")]
    [InlineData(true, false, "af901d9b84936bb7")]
    [InlineData(false, true, "56552bae58ef236d")]
    [InlineData(true, true, "1fecbcf6ecd328d8")]
    public void NeverResized_IsBitIdenticalToBefore(bool blocks, bool reference, string pinned)
    {
        ParallelTempering<State> ladder = Ladder(blocks, reference, seed: 31);
        for (int round = 1; round <= 32; round *= 2)
        {
            for (int t = 0; t < round; t++) ladder.Step();
            ladder.AdaptLadder();
        }
        ladder.TraceRungs();
        for (int t = 0; t < 300; t++) ladder.Step();
        Assert.Equal(pinned, Pin(ladder));
    }

    [Theory]
    [InlineData(false, false)]
    [InlineData(true, true)]
    public void Resize_KeepsKeptReplicas_AndNewOnesTakeTheGrowthStream(bool blocks, bool reference)
    {
        const int seed = 41;
        ParallelTempering<State> ladder = Ladder(blocks, reference, seed);
        int levels = ladder.Replicas.Count;
        // The growth stream is the SeedTree's last child, after every stream a ladder without it had.
        int[] seeds = SeedTree.Derive(seed, levels + 2 + (blocks ? levels : 0) + 1);
        Assert.Equal(new Xoshiro256PlusPlus(seeds[^1]).SaveState(), ladder.Snapshot().GrowthRng);

        for (int t = 0; t < 50; t++) ladder.Step();
        long scans = ladder.Scans, swaps = ladder.SwapAttempts, trips = ladder.RoundTrips;
        LadderSnapshot<State> before = ladder.Snapshot();
        ReversibleJumpChain<State>[] old = [.. ladder.Replicas];
        double[] rejection = ladder.RoundRejection;
        double[] inverseTemperatures = TemperingLadder.Place(ladder.InverseTemperatures, rejection, 2 * levels - 1);
        int[] sources = TemperingLadder.Sources(ladder.InverseTemperatures, rejection, inverseTemperatures);
        var forked = new List<State>();
        int[] kept = ladder.Resize(inverseTemperatures, sources, s => { forked.Add(s); return s; }, _ => Moves, blocks ? _ => [new CoordinateSlice(new IndexAndCoordinate())] : null);
        LadderSnapshot<State> after = ladder.Snapshot();

        Assert.Equal(inverseTemperatures, ladder.InverseTemperatures);
        Assert.Equal(inverseTemperatures.Length, ladder.Replicas.Count);
        Assert.Equal(0, kept[0]);
        if (reference) Assert.Equal(levels - 1, kept[^1]);
        Xoshiro256PlusPlus growth = Xoshiro256PlusPlus.Restore(before.GrowthRng.Item1, before.GrowthRng.Item2, before.GrowthRng.Item3, before.GrowthRng.Item4);
        int f = 0;
        for (int k = 0; k < inverseTemperatures.Length; k++)
        {
            ChainSnapshot<State> now = after.Replicas[k];
            if (kept[k] is int s and >= 0)
            {
                // A kept replica: the same chain, its state, stream, counters and blocks, at its new rung.
                Assert.Equal(s, sources[k]);
                Assert.Same(old[s], ladder.Replicas[k]);
                ChainSnapshot<State> was = before.Replicas[s];
                Assert.Equal((was.Current, was.Rng0, was.Rng1, was.Rng2, was.Rng3, was.Attempts, was.Accepted),
                             (now.Current, now.Rng0, now.Rng1, now.Rng2, now.Rng3, now.Attempts, now.Accepted));
                if (blocks)
                {
                    Assert.Equal(before.BlockRngs![s], after.BlockRngs![k]);
                    Assert.Equal(before.BlockAttempts![s], after.BlockAttempts![k]);
                }
                continue;
            }
            // A new replica: a fork of its source's state, its move stream and then its block stream from the growth stream.
            Assert.Equal(before.Replicas[sources[k]].Current, forked[f++]);
            Assert.Equal(before.Replicas[sources[k]].Current, now.Current);
            Assert.Equal(new Xoshiro256PlusPlus(Seed(growth)).SaveState(), (now.Rng0, now.Rng1, now.Rng2, now.Rng3));
            Assert.Equal((0L, 0L), (now.Attempts, now.Accepted));
            if (blocks)
            {
                Assert.Equal(new Xoshiro256PlusPlus(Seed(growth)).SaveState(), after.BlockRngs![k]);
                Assert.All(after.BlockAttempts![k], a => Assert.Equal(0, a));
            }
        }
        Assert.Equal(forked.Count, f);
        Assert.Equal(growth.SaveState(), after.GrowthRng);

        // The pairs' statistics, the round-trip labels and the round restart; the totals carry on.
        Assert.Equal((scans, swaps, trips), (ladder.Scans, ladder.SwapAttempts, ladder.RoundTrips));
        Assert.Equal(0, ladder.PairScans);
        Assert.All(ladder.PairStats, p => Assert.Equal((0L, 0L, double.NaN), (p.Attempts, p.Accepts, p.MeanRejection)));
        Assert.All(ladder.RoundRejection, r => Assert.True(double.IsNaN(r)));
        Assert.Equal(Enumerable.Range(0, inverseTemperatures.Length), after.LabelAt);
        Assert.Equal(new sbyte[inverseTemperatures.Length - 1].Append((sbyte)1), after.LastEnd);
        ladder.Step();
        Assert.Equal(1, ladder.PairScans);
        Assert.Equal(scans + 1, ladder.Scans);
    }

    [Fact]
    public void Resize_KeepsEachReplicaAtItsNearestNamingRung_AndDropsTheUnnamed()
    {
        ParallelTempering<State> ladder = Ladder(blocks: false, reference: false, seed: 5);   // 1, 0.5, 0.2, 0.05
        for (int t = 0; t < 10; t++) ladder.Step();
        ReversibleJumpChain<State>[] old = [.. ladder.Replicas];
        // Replica 1 (b = 0.5) is named by 0.625 and 0.4375, nearer the latter; replica 2 (0.2) by nothing; replica 3 by 0.1.
        int[] kept = ladder.Resize([1.0, 0.625, 0.4375, 0.1], [0, 1, 1, 3], s => s, _ => Moves);
        Assert.Equal([0, -1, 1, 3], kept);
        Assert.Same(old[0], ladder.Replicas[0]);
        Assert.Same(old[1], ladder.Replicas[2]);
        Assert.Same(old[3], ladder.Replicas[3]);
        Assert.DoesNotContain(old[2], ladder.Replicas);
        Assert.Equal(old[1].Current, ladder.Replicas[1].Current);
        // On a tie the colder naming rung keeps it: 0.4375 is 0.125 from both 0.5625 and 0.3125.
        Assert.Equal([0, 2, -1, 3], ladder.Resize([1.0, 0.5625, 0.3125, 0.1], [0, 2, 2, 3], s => s, _ => Moves));
    }

    // The split rule's resize (TemperingLadder.Split, F-M2's 4b): every rung stays where it is, so every replica is kept
    // at its own rung, untouched, and only the split rungs start new replicas, each from a neighbour of its pair.
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Resize_ToASplitLadder_KeepsEveryReplicaAtItsOwnRung(bool reference)
    {
        ParallelTempering<State> ladder = Ladder(blocks: true, reference, seed: 43);
        for (int t = 0; t < 40; t++) ladder.Step();
        ReversibleJumpChain<State>[] old = [.. ladder.Replicas];
        double[] before = ladder.InverseTemperatures;
        double[] rejection = [.. ladder.RoundRejection.Select((r, i) => i == 1 ? 0.97 : Math.Min(r, 0.5))];   // pair 1 read as steep
        double[] split = TemperingLadder.Split(before, rejection, 0.9, 64);
        Assert.Equal(before.Length + 1, split.Length);
        int[] sources = TemperingLadder.Sources(before, rejection, split);
        int[] kept = ladder.Resize(split, sources, s => s, _ => Moves, _ => [new CoordinateSlice(new IndexAndCoordinate())]);
        for (int k = 0, i = 0; k < split.Length; k++)
        {
            if (i < before.Length && split[k] == before[i])
            {
                Assert.Equal(i, kept[k]);
                Assert.Same(old[i], ladder.Replicas[k]);
                Assert.Equal(before[i], ladder.Replicas[k].Temperature.InverseTemperature);
                i++;
                continue;
            }
            Assert.Equal(-1, kept[k]);
            Assert.InRange(split[k], before[2], before[1]);
            Assert.Contains(sources[k], new[] { 1, 2 });
        }
    }

    [Theory]
    [InlineData(false, TemperingPath.Likelihood)]
    [InlineData(true, TemperingPath.Likelihood)]
    [InlineData(false, TemperingPath.Powered)]
    public void AfterAResize_EveryRungSamplesItsTemperedLaw(bool reference, TemperingPath path)
    {
        var target = new IndexAndCoordinate();
        ParallelTempering<State> ladder = Ladder(blocks: true, reference, seed: 23, path);
        for (int t = 0; t < 200; t++) ladder.Step();
        ladder.Resize(Grown(ladder, reference ? 8 : 7), Sources(ladder, reference ? 8 : 7), s => s, _ => Moves, _ => [new CoordinateSlice(target)]);
        double[] inverseTemperatures = ladder.InverseTemperatures;
        for (int t = 0; t < 5_000; t++) ladder.Step();

        const int n = 150_000;
        var counts = new double[inverseTemperatures.Length, IndexAndCoordinate.Indices];
        var sums = new double[inverseTemperatures.Length, IndexAndCoordinate.Indices];
        for (int i = 0; i < n; i++)
        {
            ladder.Step();
            for (int r = 0; r < inverseTemperatures.Length; r++)
            {
                State s = ladder.Replicas[r].Current;
                counts[r, s.K] += 1.0;
                sums[r, s.K] += s.Theta;
            }
        }
        for (int r = 0; r < inverseTemperatures.Length; r++)
        {
            (double[] law, double[] conditionalMean) = IndexAndCoordinate.Law(inverseTemperatures[r], path);
            for (int k = 0; k < IndexAndCoordinate.Indices; k++)
            {
                double observed = counts[r, k] / n, mean = sums[r, k] / counts[r, k];
                Assert.True(Math.Abs(observed - law[k]) < 0.012, $"b = {inverseTemperatures[r]:F4}, P(k = {k}): {observed:F4} vs {law[k]:F4}");
                Assert.True(Math.Abs(mean - conditionalMean[k]) < 0.04, $"b = {inverseTemperatures[r]:F4}, E[θ | k = {k}]: {mean:F4} vs {conditionalMean[k]:F4}");
            }
        }
        _out.WriteLine($"{path}, reference {reference}: {inverseTemperatures.Length} rungs [{string.Join(", ", inverseTemperatures.Select(b => b.ToString("F4")))}], " +
                       $"rejection [{string.Join(", ", ladder.PairStats.Select(p => p.MeanRejection.ToString("F2")))}]");
    }

    [Theory]
    [InlineData(false, false)]
    [InlineData(true, true)]
    public void SnapshotsAcrossAResize_ContinueBitForBit(bool blocks, bool reference)
    {
        var target = new IndexAndCoordinate();
        Func<int, IReadOnlyList<IBlock<ReversibleJumpChain<State>>>>? blockSet = blocks ? _ => [new CoordinateSlice(target)] : null;
        IReferenceSampler<State>? prior = reference ? new PriorDraw(target) : null;
        ParallelTempering<State> straight = Ladder(blocks, reference, seed: 61), restored = Ladder(blocks, reference, seed: 61);
        void Both(int scans)
        {
            for (int t = 0; t < scans; t++)
            {
                straight.Step();
                restored.Step();
            }
        }
        void Resize(ParallelTempering<State> l, int levels) => l.Resize(Grown(l, levels), Sources(l, levels), s => s, _ => Moves, blockSet);

        Both(30);
        straight.AdaptLadder();
        restored.AdaptLadder();
        Both(20);
        // Before a resize: the restored ladder resizes as the original does.
        restored = ParallelTempering<State>.Restore(_ => Moves, target, restored.Snapshot(), prior, blockSet);
        Resize(straight, 9);
        Resize(restored, 9);
        AssertSame(straight, restored);
        Both(100);
        // After it, and after a shrink, and once the rungs are traced.
        restored = ParallelTempering<State>.Restore(_ => Moves, target, restored.Snapshot(), prior, blockSet);
        Both(100);
        AssertSame(straight, restored);
        Resize(straight, 5);
        Resize(restored, 5);
        Both(50);
        straight.TraceRungs();
        restored.TraceRungs();
        Both(50);
        restored = ParallelTempering<State>.Restore(_ => Moves, target, restored.Snapshot(), prior, blockSet);
        Both(200);
        AssertSame(straight, restored);
        Assert.Equal(5, straight.Replicas.Count);
        for (int i = 0; i < 5; i++)
            for (int f = 0; f < 2; f++) Assert.Equal(straight.RungTraces![i].Draws(f), restored.RungTraces![i].Draws(f));
    }

    [Fact]
    public void Resize_RefusesWhatItCannotDo_AndLeavesTheLadderAsItWas()
    {
        var target = new IndexAndCoordinate();
        ParallelTempering<State> ladder = Ladder(blocks: true, reference: true, seed: 3);   // 1, 0.5, 0.2, 0.05, 0
        Func<int, IReadOnlyList<IBlock<ReversibleJumpChain<State>>>> blockSet = _ => [new CoordinateSlice(target)];
        for (int t = 0; t < 10; t++) ladder.Step();
        LadderSnapshot<State> before = ladder.Snapshot();
        double[] inverseTemperatures = [1.0, 0.7, 0.4, 0.1, 0.0];
        int[] sources = [0, 1, 2, 3, 4];

        Assert.Throws<ArgumentNullException>(() => ladder.Resize(null!, sources, s => s, _ => Moves, blockSet));
        Assert.Throws<ArgumentNullException>(() => ladder.Resize(inverseTemperatures, null!, s => s, _ => Moves, blockSet));
        Assert.Throws<ArgumentNullException>(() => ladder.Resize(inverseTemperatures, sources, null!, _ => Moves, blockSet));
        Assert.Throws<ArgumentNullException>(() => ladder.Resize(inverseTemperatures, sources, s => s, null!, blockSet));
        Assert.Equal("inverseTemperatures", Assert.Throws<ArgumentException>(() => ladder.Resize([1.0, 0.5, 0.1], [0, 1, 4], s => s, _ => Moves, blockSet)).ParamName);   // no b = 0
        Assert.Equal("inverseTemperatures", Assert.Throws<ArgumentException>(() => ladder.Resize([0.9, 0.5, 0.0], [0, 1, 4], s => s, _ => Moves, blockSet)).ParamName);
        Assert.Equal("blocksForReplica", Assert.Throws<ArgumentException>(() => ladder.Resize(inverseTemperatures, sources, s => s, _ => Moves)).ParamName);
        Assert.Equal("sources", Assert.Throws<ArgumentException>(() => ladder.Resize(inverseTemperatures, [0, 1, 4], s => s, _ => Moves, blockSet)).ParamName);
        Assert.Equal("sources", Assert.Throws<ArgumentOutOfRangeException>(() => ladder.Resize(inverseTemperatures, [0, 1, 2, 5, 4], s => s, _ => Moves, blockSet)).ParamName);
        Assert.Equal("sources", Assert.Throws<ArgumentOutOfRangeException>(() => ladder.Resize(inverseTemperatures, [0, -1, 2, 3, 4], s => s, _ => Moves, blockSet)).ParamName);
        Assert.Equal("sources", Assert.Throws<ArgumentException>(() => ladder.Resize(inverseTemperatures, [1, 1, 2, 3, 4], s => s, _ => Moves, blockSet)).ParamName);
        Assert.Equal("sources", Assert.Throws<ArgumentException>(() => ladder.Resize(inverseTemperatures, [0, 1, 2, 3, 3], s => s, _ => Moves, blockSet)).ParamName);
        // A factory that throws for a new replica leaves the ladder, its growth stream included, as it was.
        // The new rungs are 1 and 3, so the second's blocks fail after the first drew its streams.
        Assert.ThrowsAny<ArgumentException>(() => ladder.Resize([1.0, 0.7, 0.5, 0.3, 0.2, 0.05, 0.0], [0, 1, 1, 2, 2, 3, 4], s => s, _ => Moves,
                                                                i => i == 3 ? [new CoordinateSlice(target), null!] : blockSet(i)));
        Assert.Throws<ArgumentException>(() => ladder.Resize([1.0, 0.7, 0.5, 0.3, 0.2, 0.05, 0.0], [0, 1, 1, 2, 2, 3, 4], s => s, _ => [],
                                                             blockSet));
        AssertSame(ladder, ParallelTempering<State>.Restore(_ => Moves, target, before, new PriorDraw(target), blockSet));
        Assert.Equal(before.GrowthRng, ladder.Snapshot().GrowthRng);

        // A ladder without blocks takes none, and a traced ladder is fixed.
        ParallelTempering<State> plain = Ladder(blocks: false, reference: false, seed: 3);
        plain.Step();
        Assert.Equal("blocksForReplica", Assert.Throws<ArgumentException>(() =>
            plain.Resize([1.0, 0.5, 0.05], [0, 1, 3], s => s, _ => Moves, blockSet)).ParamName);
        plain.Resize([1.0, 0.5, 0.05], [0, 1, 3], s => s, _ => Moves);
        plain.TraceRungs();
        Assert.Throws<InvalidOperationException>(() => plain.Resize([1.0, 0.5, 0.05], [0, 1, 2], s => s, _ => Moves));
    }

    [Fact]
    public void TraceRungs_WithoutSlopes_TracesTheSameLogRatios()
    {
        // The powered path's slope is ln p + ℓ₁, which every target supplies.
        var target = new IndexAndCoordinate();
        ParallelTempering<State> Powered() => new(_ => Moves, target, new State(0, 0.0), [1.0, 0.5, 0.2, 0.05], TemperingPath.Powered, 13,
                                                  blocksForReplica: _ => [new CoordinateSlice(target)]);
        ParallelTempering<State> slopes = Powered(), none = Powered();
        for (int t = 0; t < 50; t++)
        {
            slopes.Step();
            none.Step();
        }
        slopes.TraceRungs();
        none.TraceRungs(traceSlopes: false);
        for (int t = 0; t < 200; t++)
        {
            slopes.Step();
            none.Step();
        }
        Assert.True(slopes.TracesSlopes);
        Assert.False(none.TracesSlopes);
        Assert.False(none.Snapshot().TracesSlopes);
        for (int i = 0; i < slopes.Replicas.Count; i++)
        {
            Assert.Equal(slopes.RungTraces![i].Draws(RungFunctional.LogRatio), none.RungTraces![i].Draws(RungFunctional.LogRatio));
            Assert.All(none.RungTraces[i].Draws(RungFunctional.Slope), v => Assert.True(double.IsNaN(v)));
            Assert.All(slopes.RungTraces[i].Draws(RungFunctional.Slope), v => Assert.True(double.IsFinite(v)));
            Assert.Equal(slopes.Replicas[i].Current, none.Replicas[i].Current);
        }
        // A restored ladder keeps the switch.
        Assert.False(ParallelTempering<State>.Restore(_ => Moves, target, none.Snapshot(), blocksForReplica: _ => [new CoordinateSlice(target)]).TracesSlopes);
    }

    // The seed a new replica's stream takes from the growth stream.
    private static int Seed(Xoshiro256PlusPlus growth) => unchecked((int)(growth.NextUInt64() >> 32));

    // `levels` rungs on the ladder's barrier from its round's rejections, and their sources.
    private static double[] Grown(ParallelTempering<State> ladder, int levels) => TemperingLadder.Place(ladder.InverseTemperatures, ladder.RoundRejection, levels);

    private static int[] Sources(ParallelTempering<State> ladder, int levels)
        => TemperingLadder.Sources(ladder.InverseTemperatures, ladder.RoundRejection, Grown(ladder, levels));

    private static ParallelTempering<State> Ladder(bool blocks, bool reference, int seed, TemperingPath path = TemperingPath.Likelihood)
    {
        var target = new IndexAndCoordinate();
        double[] inverseTemperatures = reference ? [1.0, 0.5, 0.2, 0.05, 0.0] : [1.0, 0.5, 0.2, 0.05];
        return new ParallelTempering<State>(_ => Moves, target, new State(0, 0.0), inverseTemperatures, path, seed,
                                            reference: reference ? new PriorDraw(target) : null,
                                            blocksForReplica: blocks ? _ => [new CoordinateSlice(target)] : null);
    }

    private static void AssertSame(ParallelTempering<State> a, ParallelTempering<State> b)
    {
        Assert.Equal(a.InverseTemperatures, b.InverseTemperatures);
        Assert.Equal((a.Scans, a.SwapAttempts, a.SwapAccepts, a.RoundTrips, a.AnnealedRestarts, a.PairScans),
                     (b.Scans, b.SwapAttempts, b.SwapAccepts, b.RoundTrips, b.AnnealedRestarts, b.PairScans));
        Assert.Equal(a.PairStats, b.PairStats);
        Assert.Equal(a.RoundRejection, b.RoundRejection);
        LadderSnapshot<State> x = a.Snapshot(), y = b.Snapshot();
        Assert.Equal(x.GrowthRng, y.GrowthRng);
        Assert.Equal((x.SwapRng0, x.SwapRng1, x.SwapRng2, x.SwapRng3), (y.SwapRng0, y.SwapRng1, y.SwapRng2, y.SwapRng3));
        Assert.Equal(x.LabelAt, y.LabelAt);
        Assert.Equal(x.LastEnd, y.LastEnd);
        for (int i = 0; i < x.Replicas.Count; i++)
        {
            ChainSnapshot<State> p = x.Replicas[i], q = y.Replicas[i];
            Assert.Equal((p.Current, p.Temperature, p.Rng0, p.Rng1, p.Rng2, p.Rng3, p.Attempts, p.Accepted),
                         (q.Current, q.Temperature, q.Rng0, q.Rng1, q.Rng2, q.Rng3, q.Attempts, q.Accepted));
            Assert.Equal(a.BlockStats(i), b.BlockStats(i));
        }
        Assert.Equal(x.BlockRngs, y.BlockRngs);
    }

    // A hash of a ladder's trajectory: every replica's state, rung and counters, the swap statistics, the round trips and
    // the traces.
    private static string Pin(ParallelTempering<State> ladder)
    {
        ulong h = 14695981039346656037UL;
        void Add(double v) => Word(BitConverter.DoubleToInt64Bits(v));
        void Word(long w)
        {
            unchecked
            {
                h ^= (ulong)w;
                h *= 1099511628211UL;
            }
        }
        for (int i = 0; i < ladder.Replicas.Count; i++)
        {
            ReversibleJumpChain<State> r = ladder.Replicas[i];
            Word(r.Current.K);
            Add(r.Current.Theta);
            Add(r.Temperature.InverseTemperature);
            Word(r.Attempts);
            Word(r.Accepted);
            foreach (MoveStat m in ladder.BlockStats(i)) Word(m.Accepted);
        }
        foreach (SwapPairStat s in ladder.PairStats)
        {
            Word(s.Attempts);
            Word(s.Accepts);
            Add(s.MeanRejection);
        }
        Word(ladder.RoundTrips);
        Word(ladder.AnnealedRestarts);
        foreach (ChainTrace trace in ladder.RungTraces!)
            for (int f = 0; f < 2; f++)
                foreach (double v in trace.Draws(f)) Add(v);
        return h.ToString("x16");
    }
}
