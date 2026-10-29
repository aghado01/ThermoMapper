using System;
using System.Collections.Generic;
using Maths.Rng;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Samplers.Tests;

/// <summary>
/// Swap schedules, round trips and ladder adaptation (step 3.3), against Syed, Bouchard-Côté, Deligiannidis and
/// Doucet (2022). Where every rung's local move is an exact independent draw (their efficient local exploration), the
/// round-trip rate must match their Theorem 1 from the observed rejection rates: 1/(2 + 2E) for deterministic even-odd
/// and 1/(2N + 2E) for stochastic even-odd. Adding rungs to an equalized ladder must not slow DEO, and must slow SEO.
/// Equalization must leave the barrier Σ r unchanged (their Corollary 2) and the rates equal.
/// <para>The sequential sweep's round trips are counted as mechanism, not merit. Within one swap phase a state can
/// cascade up several rungs with no local step between, so the sweep completes more round trips than the even-odd
/// schedules, each carrying less exploration; the count is not comparable with theirs.</para>
/// </summary>
public sealed class SwapScheduleTests
{
    private readonly ITestOutputHelper _out;
    public SwapScheduleTests(ITestOutputHelper output) => _out = output;

    // Prior N(0, 1); likelihood of y = 3 under N(x, 0.5²). Rung b's likelihood-tempered law is Gaussian with precision
    // 1 + b/0.25, so each rung can draw from it exactly.
    private const double Mu = 3.0, S2 = 0.25;

    private static readonly ISamplingTarget<double> Gaussian = SamplingTarget.Separable<double>(
        x => -0.5 * x * x, x => -(x - Mu) * (x - Mu) / (2.0 * S2));

    private static ParallelTempering<double> ExactLadder(double[] inverseTemperatures, SwapSchedule schedule, int seed)
    {
        ParallelTempering<double>? pt = null;
        pt = new ParallelTempering<double>(
            i => [new ExactDraw(() => pt!.Replicas[i].Temperature.InverseTemperature)],
            Gaussian, 0.0, inverseTemperatures, TemperingPath.Likelihood, seed, schedule);
        return pt;
    }

    private static (double Observed, double Ideal, double[] Rejection) Measure(ParallelTempering<double> pt, int scans)
    {
        for (int i = 0; i < 2_000; i++) pt.Step();
        long trips0 = pt.RoundTrips, scans0 = pt.Scans;
        IReadOnlyList<SwapPairStat> before = pt.PairStats;
        for (int i = 0; i < scans; i++) pt.Step();
        IReadOnlyList<SwapPairStat> after = pt.PairStats;
        var rejection = new double[after.Count];
        for (int i = 0; i < rejection.Length; i++)
            rejection[i] = (after[i].MeanRejection * pt.Scans - before[i].MeanRejection * scans0) / scans;
        double observed = (double)(pt.RoundTrips - trips0) / scans;
        return (observed, TemperingLadder.IdealRoundTripRate(rejection, pt.Schedule), rejection);
    }

    private static double MeasureTrips(ParallelTempering<double> pt, int scans)
    {
        for (int i = 0; i < 2_000; i++) pt.Step();
        long trips0 = pt.RoundTrips;
        for (int i = 0; i < scans; i++) pt.Step();
        return (double)(pt.RoundTrips - trips0) / scans;
    }

    [Theory]
    [InlineData(SwapSchedule.DeterministicEvenOdd)]
    [InlineData(SwapSchedule.StochasticEvenOdd)]
    public void RoundTripRate_MatchesTheorem1_UnderExactLocalExploration(SwapSchedule schedule)
    {
        double[] inverseTemperatures = [1.0, 0.55, 0.3, 0.16, 0.08, 0.035, 0.01, 0.0];
        (double observed, double ideal, double[] r) = Measure(ExactLadder(inverseTemperatures, schedule, seed: 21), scans: 200_000);
        _out.WriteLine($"[round-trips] {schedule}: observed {observed:F4}, Theorem 1 {ideal:F4}, Λ̂ = {TemperingLadder.Barrier(r):F3}");
        Assert.InRange(observed / ideal, 0.94, 1.06);
    }

    [Theory]
    [InlineData(SwapSchedule.DeterministicEvenOdd, 0.5)]
    [InlineData(SwapSchedule.StochasticEvenOdd, 0.1)]   // 1/(2N) with N = 5 pairs: a random walk
    [InlineData(SwapSchedule.Sequential, 1.0)]          // the sweep cascades a state up the whole ladder in one swap phase, then down a rung a scan
    public void RoundTripRate_WithoutRejections_IsTheSchedulesOwn(SwapSchedule schedule, double expected)
    {
        IRjMove<int>[] MovesFor(int replica) => [new Stay()];
        var pt = new ParallelTempering<int>(MovesFor, SamplingTarget.FromLogDensity<int>(_ => 0.0), 0,
                                            [1.0, 0.8, 0.6, 0.4, 0.2, 0.05], TemperingPath.Likelihood, masterSeed: 4, schedule);
        for (int i = 0; i < 40_000; i++) pt.Step();
        double rate = (double)pt.RoundTrips / pt.Scans;
        _out.WriteLine($"[round-trips] rejection-free {schedule}: {rate:F4}");
        Assert.InRange(rate, 0.95 * expected, 1.05 * expected);
    }

    [Fact]
    public void AddingRungs_KeepsDeoMoving_AndSlowsSeo()
    {
        var rate = new Dictionary<(SwapSchedule, int), double>();
        foreach (int levels in new[] { 5, 17 })
        {
            double[] ladder = Adapted(levels, seed: 30 + levels);
            // The sweep is logged for the record, not asserted: its cascades make its count incomparable.
            foreach (SwapSchedule schedule in new[] { SwapSchedule.DeterministicEvenOdd, SwapSchedule.StochasticEvenOdd, SwapSchedule.Sequential })
            {
                ParallelTempering<double> pt = ExactLadder(ladder, schedule, seed: 7);
                double observed = schedule == SwapSchedule.Sequential ? MeasureTrips(pt, 100_000) : Measure(pt, scans: 100_000).Observed;
                rate[(schedule, levels)] = observed;
                _out.WriteLine($"[round-trips] {schedule}, {levels} rungs: {observed:F4}");
            }
        }
        Assert.True(rate[(SwapSchedule.DeterministicEvenOdd, 17)] > 0.95 * rate[(SwapSchedule.DeterministicEvenOdd, 5)]);
        Assert.True(rate[(SwapSchedule.StochasticEvenOdd, 17)] < 0.6 * rate[(SwapSchedule.StochasticEvenOdd, 5)]);
        Assert.True(rate[(SwapSchedule.DeterministicEvenOdd, 17)] > 3.0 * rate[(SwapSchedule.StochasticEvenOdd, 17)]);
    }

    [Fact]
    public void AdaptLadder_EqualizesRejection_AndKeepsTheBarrier()
    {
        (double[] before, double[] adapted, double[] after) = AdaptAndMeasure(measureSeed: 3, adaptSeed: 5);

        double spreadBefore = Spread(before), spreadAfter = Spread(after);
        double barrierBefore = TemperingLadder.Barrier(before), barrierAfter = TemperingLadder.Barrier(after);
        _out.WriteLine($"[adapt] max/min rejection {spreadBefore:F2} → {spreadAfter:F2}; Λ̂ {barrierBefore:F3} → {barrierAfter:F3}; " +
                       $"ladder [{string.Join(", ", Array.ConvertAll(adapted, b => b.ToString("F4")))}]");

        Assert.Equal(1.0, adapted[0]);
        Assert.Equal(0.0, adapted[^1]);
        Assert.True(spreadAfter < 1.25, $"rates still unequal: max/min {spreadAfter:F2}");
        Assert.True(spreadAfter < spreadBefore);
        Assert.InRange(barrierAfter / barrierBefore, 0.93, 1.07);
    }

    // The test's computation at any seeds, which the seed sweep reruns (ComparativeSweep): the pair rejection rates of a
    // uniform ladder, its adapted rungs, and their rejection rates.
    internal static (double[] Before, double[] Adapted, double[] After) AdaptAndMeasure(int measureSeed, int adaptSeed)
    {
        double[] start = [1.0, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0.0];   // uniform: most barrier near b = 0
        (_, _, double[] before) = Measure(ExactLadder(start, SwapSchedule.DeterministicEvenOdd, seed: measureSeed), scans: 60_000);
        double[] adapted = Adapted(start, seed: adaptSeed);
        (_, _, double[] after) = Measure(ExactLadder(adapted, SwapSchedule.DeterministicEvenOdd, seed: measureSeed), scans: 60_000);
        return (before, adapted, after);
    }

    [Fact]
    public void AdaptLadder_RetempersWithoutReevaluating()
    {
        int evaluations = 0;
        ISamplingTarget<double> counting = SamplingTarget.Separable<double>(x => -0.5 * x * x,
            x => { evaluations++; return -(x - Mu) * (x - Mu) / (2.0 * S2); });
        ParallelTempering<double>? pt = null;
        pt = new ParallelTempering<double>(i => [new ExactDraw(() => pt!.Replicas[i].Temperature.InverseTemperature)],
                                           counting, 0.0, [1.0, 0.9, 0.5, 0.0], TemperingPath.Likelihood, masterSeed: 1);
        for (int i = 0; i < 500; i++) pt.Step();
        int before = evaluations;
        pt.AdaptLadder();
        Assert.Equal(before, evaluations);
        for (int r = 0; r < pt.Replicas.Count; r++)
        {
            ReversibleJumpChain<double> c = pt.Replicas[r];
            Assert.Equal(c.Temperature.LogDensity(c.CurrentEvaluation), c.CurrentLogDensity);
        }
        Assert.Throws<InvalidOperationException>(() => pt.AdaptLadder());   // no scans since
    }

    [Theory]
    [InlineData(SwapSchedule.DeterministicEvenOdd)]
    [InlineData(SwapSchedule.Sequential)]
    public void Schedules_AttemptTheirPairs(SwapSchedule schedule)
    {
        IRjMove<int>[] MovesFor(int replica) => [new Stay()];
        var pt = new ParallelTempering<int>(MovesFor, SamplingTarget.FromLogDensity<int>(_ => 0.0), 0,
                                            [1.0, 0.7, 0.4, 0.1], TemperingPath.Likelihood, masterSeed: 2, schedule);
        for (int i = 0; i < 100; i++) pt.Step();
        long[] attempts = Array.ConvertAll([.. pt.PairStats], s => s.Attempts);
        long[] expected = schedule == SwapSchedule.Sequential ? [100, 100, 100] : [50, 50, 50];   // even scans 0 and 2, odd 1
        Assert.Equal(expected, attempts);
        Assert.All(pt.PairStats, s => Assert.Equal(0.0, s.MeanRejection));
    }

    [Fact]
    public void Equalize_InvertsAKnownBarrier()
    {
        // Λ(b) = 2√b, observed on an uneven ladder: equal steps of Λ put rung k at b = (1 − k/8)².
        double[] inverseTemperatures = [1.0, 0.8, 0.6, 0.45, 0.3, 0.2, 0.1, 0.04, 0.0];
        var r = new double[inverseTemperatures.Length - 1];
        for (int i = 0; i < r.Length; i++) r[i] = (2.0 * Math.Sqrt(inverseTemperatures[i]) - 2.0 * Math.Sqrt(inverseTemperatures[i + 1])) / 4.0;   // scaled into [0, 1]
        double[] next = TemperingLadder.Equalize(inverseTemperatures, r);
        for (int k = 0; k < next.Length; k++)
            Assert.Equal(Math.Pow(1.0 - k / 8.0, 2), next[k], 0.02);
        Assert.Equal(1.0, next[0]);
        Assert.Equal(0.0, next[^1]);
    }

    [Fact]
    public void Equalize_LeavesWhatItCannotImprove()
    {
        double[] inverseTemperatures = [1.0, 0.5, 0.1];
        Assert.Equal(inverseTemperatures, TemperingLadder.Equalize(inverseTemperatures, [0.0, 0.0]));        // no rejections: no barrier to spread
        Assert.Equal([1.0, 0.2], TemperingLadder.Equalize([1.0, 0.2], [0.4]));   // no interior rung
        Assert.Throws<ArgumentException>(() => TemperingLadder.Equalize(inverseTemperatures, [0.1]));
        Assert.Throws<ArgumentOutOfRangeException>(() => TemperingLadder.Equalize(inverseTemperatures, [0.1, 1.5]));
    }

    [Fact]
    public void RecommendedLevels_IsTwiceTheBarrierPlusOne()
    {
        Assert.Equal(2, TemperingLadder.RecommendedLevels(0.0));
        Assert.Equal(6, TemperingLadder.RecommendedLevels(2.3));   // N ≈ 2Λ = 4.6 → 5 pairs
        Assert.Equal(1.6, TemperingLadder.Barrier([0.5, 0.25, 0.85]), 12);
    }

    // Algorithm 4's tuning: rounds of 1, 2, 4, … scans, adapting after each.
    private static double[] Adapted(int levels, int seed)
    {
        var start = new double[levels];
        for (int i = 0; i < levels; i++) start[i] = 1.0 - (double)i / (levels - 1);
        return Adapted(start, seed);
    }

    private static double[] Adapted(double[] start, int seed)
    {
        ParallelTempering<double> pt = ExactLadder(start, SwapSchedule.DeterministicEvenOdd, seed);
        for (int round = 1; round <= 8192; round *= 2)
        {
            for (int i = 0; i < round; i++) pt.Step();
            pt.AdaptLadder();
        }
        return pt.InverseTemperatures;
    }

    internal static double Spread(double[] r)
    {
        double lo = double.PositiveInfinity, hi = 0.0;
        foreach (double v in r) { lo = Math.Min(lo, v); hi = Math.Max(hi, v); }
        return hi / lo;
    }

    // An independence proposal from the replica's own tempered law, so every move is an exact draw (always accepted).
    private sealed class ExactDraw : IRjMove<double>
    {
        private readonly Func<double> _inverseTemperature;
        public ExactDraw(Func<double> inverseTemperature) => _inverseTemperature = inverseTemperature;
        public string Key => "draw";
        public string ReverseKey => "draw";
        public double Weight(double state) => 1.0;
        public Proposal<double>? Propose(double current, Xoshiro256PlusPlus rng)
        {
            double b = _inverseTemperature();
            double precision = 1.0 + b / S2, mean = b * Mu / S2 / precision;
            double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
            double z = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
            double next = mean + z / Math.Sqrt(precision);
            double logRatio = -0.5 * precision * ((current - mean) * (current - mean) - (next - mean) * (next - mean));
            return new Proposal<double>(next, logRatio, 0.0);
        }
    }

    private sealed class Stay : IRjMove<int>
    {
        public string Key => "stay";
        public string ReverseKey => "stay";
        public double Weight(int state) => 1.0;
        public Proposal<int>? Propose(int current, Xoshiro256PlusPlus rng) => new(current, 0.0, 0.0);
    }
}
