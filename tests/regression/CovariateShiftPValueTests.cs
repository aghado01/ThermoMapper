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
/// E35's streams under covariate shift (<see cref="CovariateShiftPValues"/>): a reference sample from the source law, then a stream from
/// the target, each unit carrying its likelihood ratio.
/// <list type="bullet">
/// <item>The inclusion chances are their definition, enumerated over subsets, extreme and zero weights included.</item>
/// <item>The guarantee's finite core holds by enumeration: over every ordering of a bag, each weighed by its likelihood, the stream's
/// p-values are jointly uniform. The unweighted ranks are not.</item>
/// <item>On a stream under the stated shift the p-values are i.i.d. uniform and Ville holds, where the unweighted ranks pile low.</item>
/// <item>A change in y given x alarms.</item>
/// <item>Equal weights give the unweighted ranks, the monitor bets on the source's p-values, and bad inputs are refused.</item>
/// </list>
/// With the qualification switch, the cost per arrival.
/// </summary>
public sealed class CovariateShiftPValueTests(ITestOutputHelper output)
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

    // ---- The inclusion chances ----

    [Fact]
    public void TheInclusionChances_AreTheirDefinition_OverEverySubset()
    {
        // πⱼ = Σ over k-subsets holding j of Π w / Σ over k-subsets of Π w, by enumeration: 300 bags of 2 to 10 units with weights from
        // 10⁻⁸ to 10⁸, some zero, at every k the positive weights allow; both sides of the tables (k at most the reference, and above).
        var rng = new Xoshiro256PlusPlus(2026_1121);
        double worst = 0.0;
        for (int trial = 0; trial < 300; trial++)
        {
            int n = 2 + (int)(9 * rng.NextDouble());
            double spread = trial % 3 == 0 ? 16.0 : 2.0;
            double[] w = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble() < 0.1 ? 0.0 : Math.Pow(10.0, spread * (rng.NextDouble() - 0.5)))];
            int positive = w.Count(v => v > 0.0);
            for (int k = 1; k <= positive; k++)
            {
                double[] pi = CovariateShiftPValues.StreamShares(w, k);
                var held = new double[n];
                double total = 0.0;
                for (int mask = 0; mask < 1 << n; mask++)
                {
                    if (System.Numerics.BitOperations.PopCount((uint)mask) != k) continue;
                    double product = 1.0;
                    for (int j = 0; j < n; j++) if ((mask >> j & 1) == 1) product *= w[j];
                    total += product;
                    for (int j = 0; j < n; j++) if ((mask >> j & 1) == 1) held[j] += product;
                }
                for (int j = 0; j < n; j++)
                {
                    double expected = held[j] / total, error = Math.Abs(pi[j] - expected);
                    worst = Math.Max(worst, error / Math.Max(expected, 1e-300));
                    Assert.True(error <= 1e-10 * expected + 1e-290, $"trial {trial}, n {n}, k {k}, unit {j} (w {w[j]:E2}): {pi[j]:R} against {expected:R}");
                }
                Assert.Equal(k, pi.Sum(), 9);
            }
        }
        output.WriteLine($"the largest relative error over 300 bags: {worst:E2}");

        // Beyond enumeration: 150 units, 50 with weights 2x for x uniform and 100 for x of density 2x, against ln e_r by a log-space
        // recursion with each unit left out in turn, on both sides of the tables.
        double[] bag = [.. Enumerable.Range(0, 50).Select(_ => 2.0 * rng.NextDouble()), .. Enumerable.Range(0, 100).Select(_ => 2.0 * Math.Sqrt(rng.NextDouble()))];
        static double LogAdd(double a, double b) => a == double.NegativeInfinity ? b : Math.Max(a, b) + Math.Log(1.0 + Math.Exp(-Math.Abs(a - b)));
        double LogE(int k, int skip)
        {
            double[] e = [0.0, .. Enumerable.Repeat(double.NegativeInfinity, k)];
            for (int i = 0; i < bag.Length; i++)
            {
                if (i == skip) continue;
                for (int r = k; r >= 1; r--) e[r] = LogAdd(e[r], e[r - 1] + Math.Log(bag[i]));
            }
            return e[k];
        }
        double large = 0.0;
        foreach (int k in new[] { 1, 30, 50, 75, 100, 149 })
        {
            double[] pi = CovariateShiftPValues.StreamShares(bag, k);
            double whole = LogE(k, -1);
            for (int j = 0; j < bag.Length; j++)
            {
                double expected = Math.Exp(Math.Log(bag[j]) + LogE(k - 1, j) - whole);
                large = Math.Max(large, Math.Abs(pi[j] - expected) / expected);
            }
        }
        output.WriteLine($"150 units, k 1 to 149: the largest relative error {large:E2}");
        Assert.True(large < 1e-12);
    }

    // ---- The guarantee's finite core ----

    public static TheoryData<double[], double[], int> Bags => new()
    {
        { [0.3, 1.2, -0.4, 2.0, 0.9, -1.1], [0.5, 2.0, 1.0, 3.5, 0.2, 1.3], 3 },
        { [0.3, 1.2, -0.4, 2.0, 0.9, -1.1], [0.0, 2.0, 1.0, 3.5, 0.2, 1.3], 2 },
        { [1.0, 1.0, 2.0, 2.0, 0.5, 2.0], [0.7, 1.9, 0.4, 2.6, 1.1, 0.3], 4 },
        { [0.1, 0.2, 0.3, 0.4, 0.5, 0.6], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 3 },
    };

    [Theory]
    [MemberData(nameof(Bags))]
    public void OverEveryOrdering_WeighedByItsLikelihood_TheStreamsPValuesAreJointlyUniform(double[] scores, double[] weights, int reference)
    {
        // Six units, `reference` of them in the reference positions and the rest in the stream's. An ordering's chance is the product
        // of w over the stream positions, normalized (Tibshirani et al.'s weighted exchangeability, independent units). Given an
        // ordering, pₜ = Bₜ + θₜTₜ with θₜ uniform, so P(pₜ ≤ u) = clamp((u − Bₜ)/Tₜ, 0, 1); over the orderings, P(every pₜ ≤ uₜ) must
        // be Π uₜ on a grid of u. The unweighted ranks, under the same law, are not, unless the weights are equal.
        int n = scores.Length, k = n - reference;
        double[] grid = [0.1, 0.35, 0.5, 0.77, 0.95];
        var orderings = new List<(double Chance, double[] Below, double[] Tied, double[] PlainBelow, double[] PlainTied)>();
        double norm = 0.0;
        foreach (int[] order in Permutations(n))
        {
            double chance = 1.0;
            for (int t = reference; t < n; t++) chance *= weights[order[t]];
            if (chance == 0.0) continue;
            norm += chance;
            var stream = new CovariateShiftPValues([.. order.Take(reference).Select(i => scores[i])], [.. order.Take(reference).Select(i => weights[i])], seed: 1);
            var plain = new SmoothedPValues(seed: 1);
            foreach (int i in order.Take(reference)) plain.Add(scores[i]);
            var below = new double[k];
            var tied = new double[k];
            var plainBelow = new double[k];
            var plainTied = new double[k];
            for (int t = 0; t < k; t++)
            {
                int unit = order[reference + t];
                WeightedPValue p = stream.Add(scores[unit], weights[unit]);
                SmoothedPValue q = plain.Add(scores[unit]);
                (below[t], tied[t], plainBelow[t], plainTied[t]) = (p.Below, p.Tied, q.Below / (double)q.Count, q.Tied / (double)q.Count);
            }
            orderings.Add((chance, below, tied, plainBelow, plainTied));
        }
        double worst = 0.0, plainWorst = 0.0;
        foreach (int[] cell in Permutations(grid.Length).Select(perm => perm.Take(k).ToArray()).Distinct(new SequenceComparer()))
        {
            double[] u = [.. cell.Select(i => grid[i])];
            double target = u.Aggregate(1.0, (a, b) => a * b), weighted = 0.0, unweighted = 0.0;
            foreach ((double chance, double[] below, double[] tied, double[] plainBelow, double[] plainTied) in orderings)
            {
                double a = chance / norm, b = chance / norm;
                for (int t = 0; t < k; t++)
                {
                    a *= Math.Clamp((u[t] - below[t]) / tied[t], 0.0, 1.0);
                    b *= Math.Clamp((u[t] - plainBelow[t]) / plainTied[t], 0.0, 1.0);
                }
                weighted += a;
                unweighted += b;
            }
            worst = Math.Max(worst, Math.Abs(weighted - target));
            plainWorst = Math.Max(plainWorst, Math.Abs(unweighted - target));
        }
        output.WriteLine($"{orderings.Count} orderings of positive chance: the weighted p-values' largest departure from Π u {worst:E2}; the unweighted ranks' {plainWorst:E3}");
        Assert.True(worst < 1e-12);
        if (weights.Distinct().Count() > 1) Assert.True(plainWorst > 1e-3);
        else Assert.True(plainWorst < 1e-12);
    }

    private sealed class SequenceComparer : IEqualityComparer<int[]>
    {
        public bool Equals(int[]? a, int[]? b) => a!.SequenceEqual(b!);

        public int GetHashCode(int[] a) => a.Aggregate(17, (h, v) => h * 31 + v);
    }

    // ---- On a stream under the stated shift ----

    // Slice 14's fixture: x uniform for the source, density 2x for the target, so w(x) = 2x; y = sin(2πx) + (0.1 + 0.6x)ζ, the noise
    // growing with x; a one-scale fit's log density, so units with larger x conform less.
    private static (double Score, double Weight) Unit(Xoshiro256PlusPlus rng, bool target, double noise = 1.0)
    {
        double x = target ? Math.Sqrt(rng.NextDouble()) : rng.NextDouble();
        double y = Math.Sin(2.0 * Math.PI * x) + noise * (0.1 + 0.6 * x) * Normal.SampleStandard(rng);
        return (Normal.MeanSd(Math.Sin(2.0 * Math.PI * x), 0.35).LogPdf(y), 2.0 * x);
    }

    [Fact]
    public void UnderTheStatedShift_ThePValuesAreIidUniform_AndVilleHolds_WhereTheUnweightedRanksPileLow()
    {
        // 2,000 streams: 50 reference units from the source, then 100 from the target. The weighted p-values pooled, at three positions
        // and lag one apart, and Ville at c = 20 under the Sleepy Jumper; the unweighted ranks of the same units beside them. At 500
        // streams the pooled KS read 0.0075 against the 1% critical 0.0073, a draw of about p 0.008: a probe at 2,000 streams on another
        // seed read 0.0017 against 0.0036, every block of 25 arrivals flat, and the chances matched a log-space recursion at 150 units.
        const int streams = 2000, reference = 50, length = 100;
        const double c = 20.0;
        var rng = new Xoshiro256PlusPlus(2026_1122);
        var weighted = new List<double>();
        var unweighted = new List<double>();
        var at = new Dictionary<int, List<double>> { [1] = [], [50] = [], [100] = [] };
        var lagA = new List<double>();
        var lagB = new List<double>();
        int alarmed = 0, plainAlarmed = 0;
        for (int r = 0; r < streams; r++)
        {
            (double Score, double Weight)[] refs = [.. Enumerable.Range(0, reference).Select(_ => Unit(rng, target: false))];
            var monitor = ConformalTestMartingale.CovariateShift(new SleepyJumperBet(), r, [.. refs.Select(u => u.Score)], [.. refs.Select(u => u.Weight)],
                                                                 new MartingaleAlarm(AlarmProcedure.Ville, c));
            var plain = new ConformalTestMartingale(new SleepyJumperBet(), r, new MartingaleAlarm(AlarmProcedure.Ville, c));
            var plainP = new SmoothedPValues(seed: r + 7);
            foreach ((double score, _) in refs) { plain.Add(score); plainP.Add(score); }
            double previous = double.NaN;
            for (int k = 1; k <= length; k++)
            {
                (double score, double weight) = Unit(rng, target: true);
                double p = monitor.Add(score, weight).PValue.P;
                plain.Add(score);
                weighted.Add(p);
                unweighted.Add(plainP.Add(score).P);
                if (at.TryGetValue(k, out List<double>? list)) list.Add(p);
                if (k > 1) { lagA.Add(previous); lagB.Add(p); }
                previous = p;
            }
            if (monitor.Alarms[0].AlarmTimes.Count > 0) alarmed++;
            if (plain.Alarms[0].AlarmTimes.Any(t => t > reference)) plainAlarmed++;
        }
        double critical = 1.63 / Math.Sqrt(weighted.Count), positional = 1.63 / Math.Sqrt(streams);
        double ma = lagA.Average(), mb = lagB.Average();
        double lag = lagA.Zip(lagB, (a, b) => (a - ma) * (b - mb)).Sum() / Math.Sqrt(lagA.Sum(a => (a - ma) * (a - ma)) * lagB.Sum(b => (b - mb) * (b - mb)));
        double rate = alarmed / (double)streams, se = Math.Sqrt(1.0 / c * (1.0 - 1.0 / c) / streams);
        output.WriteLine($"weighted: KS pooled {KolmogorovSmirnov(weighted):F4} (1% critical {critical:F4}), at k 1, 50, 100 "
                         + string.Join(", ", at.Values.Select(v => KolmogorovSmirnov(v).ToString("F4"))) + $" (critical {positional:F4}), mean p {weighted.Average():F4}, "
                         + $"lag-1 correlation {lag:F4}, Ville's false alarms {rate:F4} against 1/c = {1.0 / c:F4}");
        output.WriteLine($"unweighted ranks: KS pooled {KolmogorovSmirnov(unweighted):F4}, mean p {unweighted.Average():F4}, Ville alarms in the stream {plainAlarmed / (double)streams:F4}");
        Assert.True(KolmogorovSmirnov(weighted) < critical);
        Assert.All(at.Values, v => Assert.True(KolmogorovSmirnov(v) < positional));
        Assert.InRange(lag, -4.0 / Math.Sqrt(lagA.Count), 4.0 / Math.Sqrt(lagA.Count));
        Assert.True(rate <= 1.0 / c + 3.0 * se);
        Assert.True(KolmogorovSmirnov(unweighted) > 3.0 * critical, $"the unweighted ranks' KS {KolmogorovSmirnov(unweighted):F4}");
    }

    [Fact]
    public void AChangeInYGivenX_Alarms_TheStatedShiftDoesNot()
    {
        // 100 streams: 50 reference units, 100 stream units under the stated shift, then 200 whose noise doubles. The weighted monitor
        // under the Sleepy Jumper with Ville at c = 20: false alarms before the change, and the share alarmed after it with the median
        // delay.
        const int streams = 100, reference = 50, before = 100, after = 200;
        const double c = 20.0;
        var rng = new Xoshiro256PlusPlus(2026_1123);
        int early = 0;
        var delays = new List<int>();
        for (int r = 0; r < streams; r++)
        {
            (double Score, double Weight)[] refs = [.. Enumerable.Range(0, reference).Select(_ => Unit(rng, target: false))];
            var monitor = ConformalTestMartingale.CovariateShift(new SleepyJumperBet(), r, [.. refs.Select(u => u.Score)], [.. refs.Select(u => u.Weight)],
                                                                 new MartingaleAlarm(AlarmProcedure.Ville, c));
            for (int k = 1; k <= before + after; k++)
            {
                (double score, double weight) = Unit(rng, target: true, noise: k <= before ? 1.0 : 2.0);
                monitor.Add(score, weight);
            }
            IReadOnlyList<int> times = monitor.Alarms[0].AlarmTimes;
            if (times.Any(t => t <= before)) early++;
            int first = times.FirstOrDefault(t => t > before);
            if (first > 0) delays.Add(first - before);
        }
        delays.Sort();
        output.WriteLine($"false alarms before the change {early / (double)streams:F3}; alarmed after it {delays.Count / (double)streams:F3}, "
                         + $"median delay {(delays.Count > 0 ? delays[delays.Count / 2].ToString() : "none")}");
        Assert.True(delays.Count >= 0.9 * streams);
    }

    // ---- Equal weights, the monitor, and the refusals ----

    [Fact]
    public void EqualWeights_GiveTheUnweightedRanks_TheMonitorBetsOnTheSource_AndBadInputsAreRefused()
    {
        var rng = new Xoshiro256PlusPlus(2026_1124);
        double[] refs = [.. Enumerable.Range(0, 20).Select(_ => Normal.SampleStandard(rng))];
        double[] stream = [.. Enumerable.Range(0, 30).Select(_ => Math.Round(2.0 * Normal.SampleStandard(rng)) / 2.0)];

        // Equal weights: each unit's chance is 1/n, so the masses are the unweighted counts over n; with no reference, p itself.
        var equal = new CovariateShiftPValues(refs, [.. refs.Select(_ => 3.0)], seed: 4);
        var counts = new SmoothedPValues(seed: 9);
        foreach (double s in refs) counts.Add(s);
        foreach (double s in stream)
        {
            WeightedPValue p = equal.Add(s, 3.0);
            SmoothedPValue q = counts.Add(s);
            Assert.Equal(q.Below / (double)q.Count, p.Below, 14);
            Assert.Equal(q.Tied / (double)q.Count, p.Tied, 14);
        }
        var alone = new CovariateShiftPValues([], [], seed: 5);
        var plain = new SmoothedPValues(seed: 5);
        foreach (double s in stream) Assert.Equal(plain.Add(s).P, alone.Add(s, 0.7).P, 14);

        // The monitor bets on the source's p-values.
        double[] weights = [.. refs.Select(_ => 0.2 + rng.NextDouble())];
        var source = new CovariateShiftPValues(refs, weights, seed: 6);
        var monitor = ConformalTestMartingale.CovariateShift(new SimpleJumperBet(), 6, refs, weights);
        var capital = new BettingCapital(new SimpleJumperBet());
        Assert.Equal(MonitorForm.CovariateShift, monitor.Form);
        foreach (double s in stream)
        {
            double w = 0.2 + rng.NextDouble();
            WeightedPValue p = source.Add(s, w);
            WeightedMartingaleStep step = monitor.Add(s, w);
            Assert.Equal(p, step.PValue);
            Assert.Equal(capital.Add(p.P).LogCapital, step.LogCapital);
        }
        Assert.Equal(stream.Length, monitor.Count);

        // Refusals: the reference's counts, NaN scores and bad weights; a stream unit outside the target law; each form's input.
        Assert.Throws<ArgumentException>(() => new CovariateShiftPValues([1.0, 2.0], [1.0], 1));
        Assert.Throws<ArgumentException>(() => new CovariateShiftPValues([double.NaN], [1.0], 1));
        Assert.Throws<ArgumentException>(() => new CovariateShiftPValues([1.0], [-1.0], 1));
        Assert.Throws<ArgumentException>(() => new CovariateShiftPValues([1.0], [double.PositiveInfinity], 1));
        Assert.True(new CovariateShiftPValues([1.0, 2.0], [0.0, 1.0], 1).Add(1.5, 1.0).P > 0.0);   // a reference unit outside the target law
        var refusing = new CovariateShiftPValues(refs, weights, 1);
        Assert.Throws<ArgumentException>(() => refusing.Add(0.5, 0.0));
        Assert.Throws<ArgumentException>(() => refusing.Add(0.5, double.NaN));
        Assert.Throws<ArgumentException>(() => refusing.Add(double.NaN, 1.0));
        Assert.Equal(0, refusing.Count);
        Assert.Throws<InvalidOperationException>(() => monitor.Add(0.5));
        Assert.Throws<InvalidOperationException>(() => monitor.Add([0.0, 0.0]));
        Assert.Throws<InvalidOperationException>(() => new ConformalTestMartingale(new SimpleJumperBet(), 1).Add(0.5, 1.0));
        Assert.Throws<InvalidOperationException>(() => ConformalTestMartingale.Full(new SimpleJumperBet(), 1, 2).Add(0.5, 1.0));
        Assert.Equal(MonitorForm.Inductive, new ConformalTestMartingale(new SimpleJumperBet(), 1).Form);
        Assert.Equal(MonitorForm.Full, ConformalTestMartingale.Full(new SimpleJumperBet(), 1, 2).Form);
    }

    // ---- The cost per arrival ----

    [Fact]
    public void CostPerArrival()
    {
        // One arrival's inclusion chances and rank, with 200 reference units, at stream lengths to 2,000. Wall or cycle clock, indicative.
        if (!BassQualification.Enabled) return;
        var rng = new Xoshiro256PlusPlus(2026_1125);
        (double Score, double Weight)[] refs = [.. Enumerable.Range(0, 200).Select(_ => Unit(rng, target: false))];
        var stream = new CovariateShiftPValues([.. refs.Select(u => u.Score)], [.. refs.Select(u => u.Weight)], seed: 1);
        (double Score, double Weight)[] units = [.. Enumerable.Range(0, 2000).Select(_ => Unit(rng, target: true))];
        var line = new List<string>();
        int done = 0;
        foreach (int upTo in new[] { 100, 200, 500, 1000, 2000 })
        {
            int from = done;
            long t0 = BassCost.Tick();
            for (; done < upTo; done++) stream.Add(units[done].Score, units[done].Weight);
            line.Add($"{(BassCost.Tick() - t0) / BassCost.TicksPerSecond * 1e6 / (upTo - from):F0} µs from k {from:N0} to {upTo:N0}");
        }
        output.WriteLine($"200 reference units ({(BassCost.Cycles ? "cycles" : "wall")} clock, indicative): {string.Join("; ", line)}");
    }
}
