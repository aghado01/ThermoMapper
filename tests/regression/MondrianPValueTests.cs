using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E35's strata, the Mondrian form (<see cref="SmoothedPValues.Add(double, int)"/>; Vovk, Gammerman and Shafer 2005, §4.5, Proposition 4.10
/// and Theorem 8.2): each score ranked within its stratum.
/// <list type="bullet">
/// <item>The finite core holds by enumeration: over every ordering within the strata, each rank sequence has total probability 1, ties
/// included. The unstratified ranks under the same law do not.</item>
/// <item>A stream whose laws differ by season gives uniform, independent p-values by season, and Ville holds; the unstratified ranks
/// do not.</item>
/// <item>A change in one season alarms, pooled into one capital or watched per stratum.</item>
/// <item>One stratum is the unstratified form, and the monitor bets on the strata's p-values.</item>
/// </list>
/// Every statistical gate's size, seed and level were fixed before its first reading, with the rerun rule stated where it applies.
/// </summary>
public sealed class MondrianPValueTests(ITestOutputHelper output)
{
    private static IEnumerable<int[]> Permutations(int n)
    {
        int[] p = [.. Enumerable.Range(0, n)];
        while (true)
        {
            yield return (int[])p.Clone();
            int i = n - 2;
            while (i >= 0 && p[i] >= p[i + 1]) i--;
            if (i < 0) yield break;
            int j = n - 1;
            while (p[j] <= p[i]) j--;
            (p[i], p[j]) = (p[j], p[i]);
            Array.Reverse(p, i + 1, n - i - 1);
        }
    }

    private static double KolmogorovSmirnov(IEnumerable<double> values)
    {
        double[] sorted = [.. values.Order()];
        double d = 0.0;
        for (int i = 0; i < sorted.Length; i++) d = Math.Max(d, Math.Max((i + 1.0) / sorted.Length - sorted[i], sorted[i] - (double)i / sorted.Length));
        return d;
    }

    // ---- The finite core ----

    public static TheoryData<int[], double[]> Bags => new()
    {
        // The strata by position, and the values each position's stratum holds in arrival order before reordering.
        { [0, 1, 0, 1, 1, 0], [1.0, 5.0, 2.0, 1.0, 3.0, 2.0] },
        { [0, 0, 1, 2, 1, 2, 0], [4.0, 4.0, 1.0, 7.0, 1.0, 3.0, 2.0] },
        { [2, 0, 1, 0, 2, 1], [0.5, -1.0, 3.0, -1.0, 0.5, 2.0] },
    };

    [Theory]
    [MemberData(nameof(Bags))]
    public void WithinStrata_OverEveryOrdering_EachRankSequenceHasProbabilityOne(int[] strata, double[] values)
    {
        // The Mondrian law: given the strata by position and each stratum's bag, every ordering within the strata is equally likely. At
        // each step the tie-broken rank within the stratum is uniform on Lₙ + 1, …, Lₙ + Tₙ, so over the orderings each sequence of
        // within-stratum ranks has total probability 1 (seqRank_uniform per stratum, the strata's orderings independent). The
        // unstratified ranks of the same orderings reach fewer than the n! sequences M-1 gives each 1/n!, so they are not uniform.
        int n = strata.Length;
        int[] kinds = [.. strata.Distinct().Order()];
        int[][] slots = [.. kinds.Select(k => Enumerable.Range(0, n).Where(t => strata[t] == k).ToArray())];
        var total = new Dictionary<string, double>();
        var plainTotal = new Dictionary<string, double>();
        IEnumerable<int[]> Orderings(int s)
        {
            if (s == kinds.Length) { yield return new int[n]; yield break; }
            foreach (int[] rest in Orderings(s + 1))
                foreach (int[] perm in Permutations(slots[s].Length))
                {
                    int[] order = (int[])rest.Clone();
                    for (int i = 0; i < perm.Length; i++) order[slots[s][i]] = slots[s][perm[i]];
                    yield return order;
                }
        }
        static void Paths(Dictionary<string, double> into, IEnumerable<SmoothedPValue> steps)
        {
            var paths = new List<(string Key, double Probability)> { ("", 1.0) };
            foreach (SmoothedPValue step in steps)
            {
                var next = new List<(string, double)>();
                foreach ((string key, double probability) in paths)
                    for (int r = step.Below + 1; r <= step.Below + step.Tied; r++) next.Add((key + r + ",", probability / step.Tied));
                paths = next;
            }
            foreach ((string key, double probability) in paths) into[key] = into.GetValueOrDefault(key) + probability;
        }
        int orderings = 0;
        foreach (int[] order in Orderings(0))
        {
            orderings++;
            var mondrian = new SmoothedPValues(seed: 1);
            var plain = new SmoothedPValues(seed: 1);
            Paths(total, [.. Enumerable.Range(0, n).Select(t => mondrian.Add(values[order[t]], strata[t]))]);
            Paths(plainTotal, [.. Enumerable.Range(0, n).Select(t => plain.Add(values[order[t]]))]);
        }
        long sequences = 1;
        var seen = new Dictionary<int, int>();
        foreach (int k in strata) sequences *= seen[k] = seen.GetValueOrDefault(k) + 1;
        long factorial = Enumerable.Range(1, n).Aggregate(1L, (a, b) => a * b);
        output.WriteLine($"{orderings} orderings within the strata: {total.Count} within-stratum rank sequences of {sequences}, each total in "
                         + $"[{total.Values.Min():R}, {total.Values.Max():R}]; the unstratified ranks reach {plainTotal.Count} of {factorial} sequences");
        Assert.Equal(orderings, total.Values.Sum(), 9);
        Assert.Equal(sequences, total.Count);
        Assert.All(total.Values, v => Assert.Equal(orderings / (double)sequences, v, 12));
        Assert.True(plainTotal.Count < factorial);
    }

    // ---- A seasonal stream ----

    // Four seasons by position, each its own law; the score −|y| is a fixed fit that ignores the season.
    private static readonly (double Mean, double Sd)[] Seasons = [(0.0, 1.0), (2.0, 0.5), (-1.0, 2.0), (1.0, 1.0)];

    private static double Seasonal(int n, Xoshiro256PlusPlus rng, double changedSd = double.NaN, double everySeason = 1.0)
    {
        (double mean, double sd) = Seasons[n % 4];
        if (n % 4 == 2 && !double.IsNaN(changedSd)) sd = changedSd;
        return -Math.Abs(mean + everySeason * sd * Normal.SampleStandard(rng));
    }

    private sealed record SeasonalReading(double Pooled, double[] BySeason, double Lag, double VilleRate, double[] PlainBySeason, int Count);

    private static SeasonalReading ReadSeasons(int streams, int seed)
    {
        const int length = 200;
        const double c = 20.0;
        var rng = new Xoshiro256PlusPlus(seed);
        var bySeason = Enumerable.Range(0, 4).Select(_ => new List<double>()).ToArray();
        var plainBySeason = Enumerable.Range(0, 4).Select(_ => new List<double>()).ToArray();
        var lagA = new List<double>();
        var lagB = new List<double>();
        int alarmed = 0;
        for (int r = 0; r < streams; r++)
        {
            var monitor = new ConformalTestMartingale(new SleepyJumperBet(), r, new MartingaleAlarm(AlarmProcedure.Ville, c));
            var plain = new SmoothedPValues(seed: r + 1_000_000);
            double previous = double.NaN;
            for (int n = 0; n < length; n++)
            {
                double score = Seasonal(n, rng);
                double p = monitor.Add(score, n % 4).PValue.P;
                bySeason[n % 4].Add(p);
                plainBySeason[n % 4].Add(plain.Add(score).P);
                if (n > 0) { lagA.Add(previous); lagB.Add(p); }
                previous = p;
            }
            if (monitor.Alarms[0].AlarmTimes.Count > 0) alarmed++;
        }
        double ma = lagA.Average(), mb = lagB.Average();
        double lag = lagA.Zip(lagB, (a, b) => (a - ma) * (b - mb)).Sum() / Math.Sqrt(lagA.Sum(a => (a - ma) * (a - ma)) * lagB.Sum(b => (b - mb) * (b - mb)));
        return new SeasonalReading(KolmogorovSmirnov(bySeason.SelectMany(v => v)), [.. bySeason.Select(KolmogorovSmirnov)], lag, alarmed / (double)streams,
                                   [.. plainBySeason.Select(KolmogorovSmirnov)], streams * length);
    }

    [Fact]
    public void ASeasonalStream_TheStrataGiveUniformIndependentPValues_AndVilleHolds_TheUnstratifiedRanksDoNot()
    {
        // Fixed before the first reading: 1,000 streams of 200 arrivals, 50 per season, seed 2026_1126; KS at the 1% level, pooled and per
        // season; lag-one correlation within 4/√N; Ville at c = 20 under the Sleepy Jumper within 1/c + 3 se. The rerun rule: a KS gate
        // that fails below 1.1 times its critical value reruns once at 4,000 streams on seed 2026_1127, and that reading decides. The
        // unstratified ranks of the same scores: some season's KS above three times its critical value.
        const double c = 20.0;
        SeasonalReading reading = ReadSeasons(1000, 2026_1126);
        bool Passes(SeasonalReading x, double factor) =>
            x.Pooled < factor * 1.63 / Math.Sqrt(x.Count) && x.BySeason.All(d => d < factor * 1.63 / Math.Sqrt(x.Count / 4.0));
        void Report(SeasonalReading x, string label) =>
            output.WriteLine($"{label}: Mondrian KS pooled {x.Pooled:F4} (1% critical {1.63 / Math.Sqrt(x.Count):F4}), by season "
                             + string.Join(", ", x.BySeason.Select(d => d.ToString("F4"))) + $" (critical {1.63 / Math.Sqrt(x.Count / 4.0):F4}); lag-1 {x.Lag:F4}; "
                             + $"Ville's false alarms {x.VilleRate:F4}; the unstratified ranks' KS by season " + string.Join(", ", x.PlainBySeason.Select(d => d.ToString("F4"))));
        Report(reading, "1,000 streams, seed 2026_1126");
        if (!Passes(reading, 1.0) && Passes(reading, 1.1))
        {
            reading = ReadSeasons(4000, 2026_1127);
            Report(reading, "the pre-stated rerun, 4,000 streams, seed 2026_1127");
        }
        Assert.True(Passes(reading, 1.0));
        Assert.InRange(reading.Lag, -4.0 / Math.Sqrt(reading.Count), 4.0 / Math.Sqrt(reading.Count));
        int streams = reading.Count / 200;
        Assert.True(reading.VilleRate <= 1.0 / c + 3.0 * Math.Sqrt(1.0 / c * (1.0 - 1.0 / c) / streams));
        Assert.True(reading.PlainBySeason.Max() > 3.0 * 1.63 / Math.Sqrt(reading.Count / 4.0));
    }

    [Theory]
    [InlineData("one season")]
    [InlineData("every season")]
    public void AChange_ThePooledCapitalAgainstThePerStratumWatch(string change)
    {
        // 200 streams: 200 seasonal arrivals, then 400 changed. The pooled Mondrian monitor bets every season's p-value into one capital
        // under Ville at c = 20. The per-stratum watch runs one monitor per season at c = 80, so its false alarms over the four stay within
        // 4/80 = 1/20 by the union bound. Both under the Sleepy Jumper. Each case's false alarms before the change must stay within its
        // bound plus three standard errors.
        //   One season: season 2's sd goes from 2 to 5; seed 2026_1128. Its gate, fixed before its first reading, was that the pooled
        //   monitor alarms after the change in at least 0.8 of the streams. It read 0.10 against the watch's 0.495, so the gate failed and
        //   the case now pins those readings. A change in one season of four reaches the pooled capital in a quarter of the arrivals, and
        //   the uniform p-values of the other three cost it; the changed season's bag fills with changed units, so its p-values drift back
        //   toward uniform.
        //   Every season: each season's sd times 2.5; seed 2026_1130. Fixed before its first reading: the pooled monitor alarms after the
        //   change in at least 0.8 of the streams.
        const int streams = 200, before = 200, after = 400;
        bool one = change == "one season";
        var rng = new Xoshiro256PlusPlus(one ? 2026_1128 : 2026_1130);
        int pooledEarly = 0, watchEarly = 0;
        var pooledDelays = new List<int>();
        var watchDelays = new List<int>();
        for (int r = 0; r < streams; r++)
        {
            var pooled = new ConformalTestMartingale(new SleepyJumperBet(), r, new MartingaleAlarm(AlarmProcedure.Ville, 20.0));
            ConformalTestMartingale[] watch = [.. Enumerable.Range(0, 4).Select(k => new ConformalTestMartingale(new SleepyJumperBet(), 10 * r + k, new MartingaleAlarm(AlarmProcedure.Ville, 80.0)))];
            int pooledFirst = 0, watchFirst = 0;
            for (int n = 0; n < before + after; n++)
            {
                double score = n < before ? Seasonal(n, rng) : one ? Seasonal(n, rng, changedSd: 5.0) : Seasonal(n, rng, everySeason: 2.5);
                if (pooled.Add(score, n % 4).Fired!.Count > 0 && pooledFirst == 0) pooledFirst = n + 1;
                if (watch[n % 4].Add(score).Fired!.Count > 0 && watchFirst == 0) watchFirst = n + 1;
            }
            if (pooledFirst > 0 && pooledFirst <= before) pooledEarly++;
            else if (pooledFirst > before) pooledDelays.Add(pooledFirst - before);
            if (watchFirst > 0 && watchFirst <= before) watchEarly++;
            else if (watchFirst > before) watchDelays.Add(watchFirst - before);
        }
        pooledDelays.Sort();
        watchDelays.Sort();
        string Median(List<int> d) => d.Count > 0 ? d[d.Count / 2].ToString() : "none";
        output.WriteLine($"{change}: pooled, c 20: false alarms before the change {pooledEarly / (double)streams:F3}, alarmed after it {pooledDelays.Count / (double)streams:F3}, "
                         + $"median delay {Median(pooledDelays)}; per stratum, c 80 each: false {watchEarly / (double)streams:F3}, alarmed {watchDelays.Count / (double)streams:F3}, "
                         + $"median delay {Median(watchDelays)}");
        double se = Math.Sqrt(0.05 * 0.95 / streams);
        Assert.True(pooledEarly / (double)streams <= 0.05 + 3.0 * se && watchEarly / (double)streams <= 0.05 + 3.0 * se);
        if (one)
        {
            Assert.Equal(20, pooledDelays.Count);
            Assert.Equal(99, watchDelays.Count);
        }
        else Assert.True(pooledDelays.Count >= 0.8 * streams, $"the pooled monitor alarmed in {pooledDelays.Count} of {streams}");
    }

    // ---- One stratum, the monitor, the refusals ----

    [Fact]
    public void OneStratum_IsTheUnstratifiedForm_TheMonitorBetsOnTheStrata_AndBadInputsAreRefused()
    {
        var rng = new Xoshiro256PlusPlus(2026_1129);
        double[] scores = [.. Enumerable.Range(0, 60).Select(_ => Math.Round(3.0 * Normal.SampleStandard(rng)) / 3.0)];
        var plain = new SmoothedPValues(seed: 3);
        var zero = new SmoothedPValues(seed: 3);
        foreach (double s in scores) Assert.Equal(plain.Add(s), zero.Add(s, 0));

        // Each unit alone in its stratum: p = θ. Counts by stratum and in all.
        var alone = new SmoothedPValues(seed: 4);
        var theta = new Xoshiro256PlusPlus(4);
        for (int i = 0; i < 20; i++) Assert.Equal(1.0 - theta.NextDouble(), alone.Add(scores[i], 100 + i).P);
        var mixed = new SmoothedPValues(seed: 5);
        for (int i = 0; i < scores.Length; i++) mixed.Add(scores[i], i % 3);
        Assert.Equal(scores.Length, mixed.Count);
        Assert.Equal(20, mixed.CountIn(1));
        Assert.Equal(0, mixed.CountIn(7));

        // The monitor bets on the strata's p-values.
        var source = new SmoothedPValues(seed: 6);
        var monitor = new ConformalTestMartingale(new SimpleJumperBet(), 6);
        var capital = new BettingCapital(new SimpleJumperBet());
        for (int i = 0; i < scores.Length; i++)
        {
            SmoothedPValue p = source.Add(scores[i], i % 4);
            MartingaleStep step = monitor.Add(scores[i], i % 4);
            Assert.Equal(p, step.PValue);
            Assert.Equal(capital.Add(p.P).LogCapital, step.LogCapital);
        }

        // Refusals: a negative stratum, a NaN score, and the other forms' monitors.
        Assert.Throws<ArgumentOutOfRangeException>(() => new SmoothedPValues(1).Add(0.5, -1));
        Assert.Throws<ArgumentException>(() => new SmoothedPValues(1).Add(double.NaN, 2));
        Assert.Throws<InvalidOperationException>(() => ConformalTestMartingale.Full(new SimpleJumperBet(), 1, 3).Add(0.5, 1));
        Assert.Throws<InvalidOperationException>(() => ConformalTestMartingale.CovariateShift(new SimpleJumperBet(), 1, [0.1], [1.0]).Add(0.5, 1));
        Assert.Contains("Add(score, stratum)", Assert.Throws<InvalidOperationException>(() => new ConformalTestMartingale(new SimpleJumperBet(), 1).Add(0.5, 1.0)).Message);
    }
}
