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
/// E35's full form (<see cref="FullFormPValues"/>): every unit rescored at each arrival under the fit's draws reweighted by the stream so
/// far. M-1's finite core holds by enumeration with either score: over every ordering of a data set, each rank sequence occurs once. The
/// same scores, each kept from its own time, break it. The p-values are i.i.d. uniform on exchangeable streams, the draws' weights
/// degenerating as they grow, and Ville's alarm holds. The scores at each time are a function of the units so far, whatever their order,
/// computed exactly where the scaled sums underflow. The seam reads as the spans do, and what cannot be read is refused. Against the
/// inductive form on the same draws: the delay to an alarm after a change. With the qualification switch, the cost per arrival.
/// </summary>
public sealed class FullFormPValueTests(ITestOutputHelper output)
{
    private const double HalfLogTwoPi = 0.91893853320467274178;

    // The fixture: y ~ N(μ, 1), the fit's draws of μ.
    private static double LogDensity(double y, double mu) => -0.5 * (y - mu) * (y - mu) - HalfLogTwoPi;

    private static double[] Densities(double y, double[] mus) => [.. mus.Select(mu => LogDensity(y, mu))];

    private static double[] Cdfs(double y, double[] mus) => [.. mus.Select(mu => Normal.StandardLogCdf(y - mu))];

    // Draws of μ from the posterior given n0 training units drawn from N(0, 1), under a N(0, 10²) prior.
    private static double[] Draws(int n0, int draws, Xoshiro256PlusPlus rng)
    {
        double sum = 0.0;
        for (int i = 0; i < n0; i++) sum += Normal.SampleStandard(rng);
        double variance = 1.0 / (n0 + 0.01), mean = variance * sum;
        return [.. Enumerable.Range(0, draws).Select(_ => mean + Math.Sqrt(variance) * Normal.SampleStandard(rng))];
    }

    private static SmoothedPValue Add(FullFormPValues stream, double y, double[] mus)
        => stream.Add(Densities(y, mus), stream.Score == FullFormScore.Centrality ? Cdfs(y, mus) : default);

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

    // ---- M-1's finite core, by enumeration ----

    [Theory]
    [InlineData(FullFormScore.Density)]
    [InlineData(FullFormScore.Centrality)]
    public void M1_OverEveryOrdering_EachRankSequenceOnce_TheScoresKeptAtTheirOwnTimeDoNot(FullFormScore score)
    {
        // Six continuous units, so no ties: over the 720 orderings each rank sequence r₁, …, r₆ with rₙ in 1, …, n occurs exactly once
        // (seqRank_uniform with every T = 1). The scores at each time are recomputed under the stream's weights, a symmetric function of
        // the units so far. Each newest score kept from its own time and ranked among the earlier ones so kept is not: in-sample scores
        // drift, and some sequences repeat.
        var rng = new Xoshiro256PlusPlus(2026_1115);
        double[] mus = Draws(3, 40, rng);
        double[] data = [.. Enumerable.Range(0, 6).Select(_ => 0.4 + 1.3 * Normal.SampleStandard(rng))];
        var full = new Dictionary<string, int>();
        var kept = new Dictionary<string, int>();
        foreach (int[] order in Permutations(data.Length))
        {
            var stream = new FullFormPValues(mus.Length, seed: 1, score: score);
            SmoothedPValue[] steps = [.. order.Select(i => Add(stream, data[i], mus))];
            Assert.All(steps, s => Assert.Equal(1, s.Tied));
            string key = string.Join(",", steps.Select(s => s.Below + 1));
            full[key] = full.GetValueOrDefault(key) + 1;
            string drifted = string.Join(",", new SmoothedPValues(seed: 1).AddRange(steps.Select(s => s.Score)).Select(s => s.Below + 1));
            kept[drifted] = kept.GetValueOrDefault(drifted) + 1;
        }
        output.WriteLine($"{score}: {full.Count} rank sequences over 720 orderings recomputed; {kept.Count} kept from their own time, the most "
                         + $"frequent {kept.Values.Max()} times");
        Assert.Equal(720, full.Count);
        Assert.All(full.Values, v => Assert.Equal(1, v));
        Assert.True(kept.Count < 720);
    }

    // ---- i.i.d. uniform on exchangeable streams ----

    [Theory]
    [InlineData(FullFormScore.Density)]
    [InlineData(FullFormScore.Centrality)]
    public void ThePValues_AreIidUniform_OnExchangeableStreams_AndVilleHolds_TheWeightsDegenerating(FullFormScore score)
    {
        // 1,000 streams of 100 units from N(0.5, 1), each fit's 100 draws from five training units of N(0, 1): the fit is off, which
        // the guarantee does not mind. Pooled, at three positions (from n 2: p₁ is θ₁ by construction), and lag one apart; and Ville's
        // alarm at c = 20 under the Simple Jumper.
        const int streams = 1000, length = 100, draws = 100;
        const double c = 20.0;
        var rng = new Xoshiro256PlusPlus(2026_1116);
        var pooled = new List<double>();
        var at = new Dictionary<int, List<double>> { [2] = [], [50] = [], [100] = [] };
        var ess = new Dictionary<int, double> { [1] = 0, [10] = 0, [100] = 0 };
        var lagA = new List<double>();
        var lagB = new List<double>();
        int alarmed = 0;
        for (int r = 0; r < streams; r++)
        {
            double[] mus = Draws(5, draws, rng);
            var monitor = ConformalTestMartingale.Full(new SimpleJumperBet(), seed: r, draws, score: score,
                                                      alarms: [new MartingaleAlarm(AlarmProcedure.Ville, c)]);
            double previous = double.NaN;
            for (int n = 1; n <= length; n++)
            {
                double y = 0.5 + Normal.SampleStandard(rng);
                double p = monitor.Add(Densities(y, mus), score == FullFormScore.Centrality ? Cdfs(y, mus) : default).PValue.P;
                pooled.Add(p);
                if (at.TryGetValue(n, out List<double>? list)) list.Add(p);
                if (ess.ContainsKey(n)) ess[n] += monitor.EffectiveSampleSize / streams;
                if (n > 1) { lagA.Add(previous); lagB.Add(p); }
                previous = p;
            }
            if (monitor.Alarms[0].AlarmTimes.Count > 0) alarmed++;
        }
        double critical = 1.63 / Math.Sqrt(pooled.Count), positional = 1.63 / Math.Sqrt(streams);
        double ma = lagA.Average(), mb = lagB.Average();
        double lag = lagA.Zip(lagB, (a, b) => (a - ma) * (b - mb)).Sum() / Math.Sqrt(lagA.Sum(a => (a - ma) * (a - ma)) * lagB.Sum(b => (b - mb) * (b - mb)));
        double rate = alarmed / (double)streams, se = Math.Sqrt(1.0 / c * (1.0 - 1.0 / c) / streams);
        output.WriteLine($"{score}: KS pooled {KolmogorovSmirnov(pooled):F4} (1% critical {critical:F4}); at n 2, 50, 100: "
                         + string.Join(", ", at.Values.Select(v => KolmogorovSmirnov(v).ToString("F4"))) + $" (critical {positional:F4}); lag-1 correlation {lag:F4}; "
                         + $"Ville's false alarms {rate:F4} against 1/c = {1.0 / c:F4}; mean ESS of 100 draws at n 1, 10, 100: "
                         + string.Join(", ", ess.Values.Select(v => v.ToString("F1"))));
        Assert.True(KolmogorovSmirnov(pooled) < critical);
        Assert.All(at.Values, v => Assert.True(KolmogorovSmirnov(v) < positional));
        Assert.InRange(lag, -4.0 / Math.Sqrt(lagA.Count), 4.0 / Math.Sqrt(lagA.Count));
        Assert.True(rate <= 1.0 / c + 3.0 * se);
        Assert.True(ess[100] < ess[10] && ess[10] < ess[1]);
    }

    // ---- A function of the bag, computed exactly ----

    [Fact]
    public void TheScores_AreAFunctionOfTheUnitsSoFar_AndExactWhereTheScaledSumsUnderflow()
    {
        // 200 wide draws of μ from N(0, 10²); 300 units from N(0, 1) and one at 60. The 300 hold the stream's weights near μ = 0, where
        // the far unit's terms lie some e^−1,200 below their largest, at the widest draw: its scaled sum underflows, and its density is
        // read in logs. Every unit's density and centrality after the last arrival, in four orders of arrival, against the definitions.
        var rng = new Xoshiro256PlusPlus(2026_1117);
        double[] mus = [.. Enumerable.Range(0, 200).Select(_ => 10.0 * Normal.SampleStandard(rng))];
        double[] data = [.. Enumerable.Range(0, 300).Select(_ => Normal.SampleStandard(rng)), 60.0];
        int n = data.Length;
        double[] stream = new double[mus.Length];
        for (int s = 0; s < mus.Length; s++) stream[s] = data.Sum(y => LogDensity(y, mus[s]));
        static double LogSumExp(IEnumerable<double> v)
        {
            double top = v.Max();
            return top + Math.Log(v.Sum(x => Math.Exp(x - top)));
        }
        double logTotal = LogSumExp(stream), top = stream.Max();
        double[] w = [.. stream.Select(v => Math.Exp(v - top))];
        double[] density = [.. data.Select(y => LogSumExp(mus.Select((mu, s) => stream[s] + LogDensity(y, mu))) - logTotal)];
        double[] centrality = [.. data.Select(y => Math.Min(mus.Select((mu, s) => w[s] * Math.Exp(Normal.StandardLogCdf(y - mu))).Sum(),
                                                            mus.Select((mu, s) => w[s] * Math.Exp(Normal.StandardLogCdf(mu - y))).Sum()) / w.Sum())];
        double peak = mus.Max(mu => LogDensity(60.0, mu));
        double scaled = mus.Select((mu, s) => w[s] * Math.Exp(LogDensity(60.0, mu) - peak)).Sum();
        Assert.True(scaled < 1e-250, $"the far unit's scaled sum {scaled:E3} does not underflow");
        for (int trial = 0; trial < 4; trial++)
        {
            int[] order = [.. Enumerable.Range(0, n).OrderBy(_ => rng.NextDouble())];
            var byDensity = new FullFormPValues(mus.Length, seed: trial);
            var byCentrality = new FullFormPValues(mus.Length, seed: trial, score: FullFormScore.Centrality);
            foreach (int i in order)
            {
                Add(byDensity, data[i], mus);
                Add(byCentrality, data[i], mus);
            }
            for (int k = 0; k < n; k++)
            {
                double expected = density[order[k]], read = byDensity.Scores[k];
                Assert.True(Math.Abs(read - expected) <= 1e-8 * (1.0 + Math.Abs(expected)), $"order {trial}, unit {order[k]}'s density: {read} against {expected}");
                expected = centrality[order[k]];
                read = byCentrality.Scores[k];
                Assert.True(Math.Abs(read - expected) <= 1e-9 * expected + 1e-300, $"order {trial}, unit {order[k]}'s centrality: {read} against {expected}");
            }
        }
        output.WriteLine($"the far unit's scaled sum {scaled:E1}, its density {density[^1]:F2}; the others' densities from {density[..^1].Min():F3} to "
                         + $"{density[..^1].Max():F3}, centralities from {centrality[..^1].Min():E3} to {centrality[..^1].Max():F3}");
    }

    // ---- The seam, the monitor, and the refusals ----

    [Fact]
    public void TheSeam_ReadsAsTheSpans_TheMonitorBetsOnThem_AndWhatCannotBeRead_IsRefused()
    {
        var rng = new Xoshiro256PlusPlus(2026_1118);
        double[] mus = Draws(4, 30, rng);
        double[] ys = [.. Enumerable.Range(0, 8).Select(_ => Normal.SampleStandard(rng))];
        var values = new double[mus.Length, ys.Length];
        var logCdf = new double[mus.Length, ys.Length];
        for (int s = 0; s < mus.Length; s++)
            for (int t = 0; t < ys.Length; t++) { values[s, t] = LogDensity(ys[t], mus[s]); logCdf[s, t] = Normal.StandardLogCdf(ys[t] - mus[s]); }
        double[] weights = [.. mus.Select((_, s) => 1.0 + s % 3)];
        CandidateLogDensity Seam(ConformalForms forms = ConformalForms.Both, bool cdf = true, IReadOnlyList<double>? w = null)
            => new(values, [.. ys.Select(y => new[] { y })], []) { LogCdf = cdf ? logCdf : null, Weights = w ?? weights, Forms = forms, FormsReason = forms == ConformalForms.Both ? null : "a test's reason" };

        foreach (FullFormScore score in new[] { FullFormScore.Density, FullFormScore.Centrality })
        {
            SmoothedPValue[] bySeam = new FullFormPValues(mus.Length, 5, weights, score).AddRange(Seam());
            var spans = new FullFormPValues(mus.Length, 5, weights, score);
            SmoothedPValue[] bySpans = [.. ys.Select((y, t) => spans.Add(Densities(y, mus), Cdfs(y, mus)))];
            Assert.Equal(bySpans, bySeam);
            var monitor = ConformalTestMartingale.Full(new SimpleJumperBet(), 5, mus.Length, weights, score);
            MartingaleStep[] steps = monitor.AddRange(Seam());
            var capital = new BettingCapital(new SimpleJumperBet());
            for (int t = 0; t < steps.Length; t++)
            {
                Assert.Equal(bySeam[t], steps[t].PValue);
                Assert.Equal(capital.Add(bySeam[t].P).LogCapital, steps[t].LogCapital);
            }
            Assert.True(monitor.IsFullForm);
            Assert.Equal(spans.EffectiveSampleSize, monitor.EffectiveSampleSize);
            output.WriteLine($"{score}: p {string.Join(", ", bySeam.Select(p => p.P.ToString("F3")))}; ESS {spans.EffectiveSampleSize:F1} of {mus.Length}");
        }

        // The constructor's and each arrival's refusals; a refused arrival leaves the stream as it was.
        Assert.Throws<ArgumentOutOfRangeException>(() => new FullFormPValues(0, 1));
        Assert.Throws<ArgumentException>(() => new FullFormPValues(3, 1, [1.0, 1.0]));
        Assert.Throws<ArgumentException>(() => new FullFormPValues(3, 1, [1.0, -1.0, 1.0]));
        Assert.Throws<ArgumentException>(() => new FullFormPValues(3, 1, [0.0, 0.0, 0.0]));
        var stream = new FullFormPValues(3, 1);
        Assert.Throws<ArgumentException>(() => stream.Add([0.0, 0.0]));
        Assert.Throws<ArgumentException>(() => stream.Add([0.0, double.NaN, 0.0]));
        Assert.Throws<ArgumentException>(() => stream.Add([0.0, double.PositiveInfinity, 0.0]));
        Assert.Throws<ArgumentException>(() => stream.Add([double.NegativeInfinity, double.NegativeInfinity, double.NegativeInfinity]));
        Assert.Equal(0, stream.Count);
        stream.Add([-1.0, double.NegativeInfinity, double.NegativeInfinity]);
        Assert.Throws<ArgumentException>(() => stream.Add([double.NegativeInfinity, -1.0, -2.0]));
        Assert.Equal(1, stream.Count);
        Assert.Equal(2, stream.Add([-2.0, -1.0, -1.0]).Count);
        var central = new FullFormPValues(3, 1, score: FullFormScore.Centrality);
        Assert.Throws<ArgumentException>(() => central.Add([0.0, 0.0, 0.0]));
        Assert.Throws<ArgumentException>(() => central.Add([0.0, 0.0, 0.0], [-0.1, 0.1, -0.2]));
        Assert.Throws<ArgumentException>(() => central.Add([0.0, 0.0, 0.0], [-0.1, double.NaN, -0.2]));
        Assert.Equal(0, central.Count);

        // The seam's: a split-only seam (its reason quoted), other draws or weights, a grid of two values, no log CDF for the centrality.
        var fromSeam = new FullFormPValues(mus.Length, 1, weights, FullFormScore.Centrality);
        Assert.Contains("a test's reason", Assert.Throws<ArgumentException>(() => fromSeam.AddRange(Seam(ConformalForms.Split))).Message);
        Assert.Throws<ArgumentException>(() => new FullFormPValues(mus.Length, 1).AddRange(Seam()));
        Assert.Throws<ArgumentException>(() => fromSeam.AddRange(Seam(w: [.. weights.Select(v => v + 1.0)])));
        Assert.Throws<ArgumentException>(() => fromSeam.AddRange(Seam(cdf: false)));
        Assert.Throws<ArgumentException>(() => new FullFormPValues(mus.Length + 1, 1, [.. weights, 1.0]).AddRange(Seam()));
        var twoValues = new CandidateLogDensity(new double[mus.Length, 2], [[0.0, 1.0]], []);
        Assert.Throws<ArgumentException>(() => new FullFormPValues(mus.Length, 1).AddRange(twoValues));
        Assert.Equal(0, fromSeam.Count);

        // Each form refuses the other's input.
        var inductive = new ConformalTestMartingale(new SimpleJumperBet(), 1);
        Assert.False(inductive.IsFullForm);
        Assert.True(double.IsNaN(inductive.EffectiveSampleSize));
        Assert.Throws<InvalidOperationException>(() => inductive.Add([0.0, 0.0, 0.0]));
        Assert.Throws<InvalidOperationException>(() => inductive.AddRange(Seam()));
        Assert.Throws<InvalidOperationException>(() => ConformalTestMartingale.Full(new SimpleJumperBet(), 1, 3).Add(0.5));
    }

    // ---- The delay to an alarm, against the inductive form ----

    [Theory]
    [InlineData("mean", 3)]
    [InlineData("mean", 200)]
    [InlineData("noise", 3)]
    [InlineData("noise", 200)]
    public void AfterAChange_BothFormsAlarm_TheirDelaysSideBySide(string change, int trainingUnits)
    {
        // 200 streams: 100 units of N(0, 1), then 200 of N(1.5, 1) (mean) or N(0, 2.5²) (noise). One fit's 50 draws per stream, from
        // `trainingUnits` units of N(0, 1). The inductive form scores each unit by the fit's predictive density, the full form by the
        // in-sample density under the stream's weights; both under the Simple Jumper with Ville's alarm at c = 20.
        const int streams = 200, before = 100, after = 200, draws = 50;
        const double c = 20.0;
        var rng = new Xoshiro256PlusPlus(2026_1119 + trainingUnits + (change == "mean" ? 0 : 7));
        var delays = new Dictionary<string, List<int>> { ["inductive"] = [], ["full"] = [] };
        var early = new Dictionary<string, int> { ["inductive"] = 0, ["full"] = 0 };
        for (int r = 0; r < streams; r++)
        {
            double[] mus = Draws(trainingUnits, draws, rng);
            var inductive = new ConformalTestMartingale(new SimpleJumperBet(), r, new MartingaleAlarm(AlarmProcedure.Ville, c));
            var full = ConformalTestMartingale.Full(new SimpleJumperBet(), r, draws, alarms: [new MartingaleAlarm(AlarmProcedure.Ville, c)]);
            for (int n = 1; n <= before + after; n++)
            {
                double y = n <= before ? Normal.SampleStandard(rng) : change == "mean" ? 1.5 + Normal.SampleStandard(rng) : 2.5 * Normal.SampleStandard(rng);
                double[] terms = Densities(y, mus);
                double top = terms.Max();
                inductive.Add(top + Math.Log(terms.Average(t => Math.Exp(t - top))));
                full.Add(terms);
            }
            foreach ((string name, ConformalTestMartingale m) in new[] { ("inductive", inductive), ("full", full) })
            {
                IReadOnlyList<int> times = m.Alarms[0].AlarmTimes;
                if (times.Any(t => t <= before)) early[name]++;
                int first = times.FirstOrDefault(t => t > before);
                if (first > 0) delays[name].Add(first - before);
            }
        }
        foreach (string name in delays.Keys)
        {
            List<int> d = delays[name];
            output.WriteLine($"{change} change, {trainingUnits} training units, {name}: alarmed after the change {d.Count / (double)streams:F3}, "
                             + $"median delay {(d.Count > 0 ? d.Order().ElementAt(d.Count / 2) : double.NaN)}, false alarms before it {early[name] / (double)streams:F3}");
        }
        // Both forms alarm on a change in the noise, which a location model cannot absorb, and on any change under the tight fit. Under the
        // vague fit, the full form's posterior follows a shift in the mean, so the units after it score as conforming and it alarms less
        // often than the inductive form.
        double fullRate = delays["full"].Count / (double)streams, fixedFit = delays["inductive"].Count / (double)streams;
        if (change == "noise" || trainingUnits > 3) Assert.True(fullRate >= 0.9 && fixedFit >= 0.9, $"alarmed: full {fullRate}, inductive {fixedFit}");
        else Assert.True(fullRate < fixedFit && fullRate >= 0.5, $"alarmed: full {fullRate}, inductive {fixedFit}");
    }

    // ---- The cost per arrival ----

    [Fact]
    public void CostPerArrival()
    {
        // One arrival's reweighing, rescoring and rank with 1,000 draws, at stream lengths to 2,000, after a warm-up stream. Cycle clock,
        // indicative.
        if (!BassQualification.Enabled) return;
        var rng = new Xoshiro256PlusPlus(2026_1120);
        double[] mus = Draws(20, 1000, rng);
        var warm = new FullFormPValues(mus.Length, 1);
        for (int i = 0; i < 300; i++) Add(warm, Normal.SampleStandard(rng), mus);
        foreach (FullFormScore score in new[] { FullFormScore.Density, FullFormScore.Centrality })
        {
            var stream = new FullFormPValues(mus.Length, 2, score: score);
            double[] ys = [.. Enumerable.Range(0, 2000).Select(_ => Normal.SampleStandard(rng))];
            var line = new List<string>();
            int done = 0;
            foreach (int upTo in new[] { 250, 500, 1000, 2000 })
            {
                int from = done;
                long t0 = BassCost.Tick();
                for (; done < upTo; done++) Add(stream, ys[done], mus);
                double perArrival = (BassCost.Tick() - t0) / BassCost.TicksPerSecond * 1e6 / (upTo - from);
                line.Add($"{perArrival:F0} µs per arrival from n {from:N0} to {upTo:N0}");
            }
            output.WriteLine($"{score}, 1,000 draws ({(BassCost.Cycles ? "cycles" : "wall")} clock, indicative): {string.Join("; ", line)}; ESS at the end {stream.EffectiveSampleSize:F2}");
        }
    }
}
