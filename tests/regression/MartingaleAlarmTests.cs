using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The alarms of a conformal test martingale (E35's slice 6c).
/// <list type="bullet">
/// <item>Each statistic is its definition read from the capital's path, restarts included.</item>
/// <item>Ville's false alarms are at most 1/c on exchangeable streams, its guarantee proved.</item>
/// <item>CUSUM's and Shiryaev–Roberts's mean time to a false alarm is at least c, and their long-run frequency of false alarms at most
/// 1/c, measured and marked as posed.</item>
/// </list>
/// With the qualification switch: the delay to an alarm after a change in the mean or the noise, per bet and procedure; false alarms and
/// detection on a fixed BARS fit's and a fixed BASS fit's scores; and the cost per arrival.
/// </summary>
public sealed class MartingaleAlarmTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static IBettingFunction Bet(string name) => name switch
    {
        "jumper" => new SimpleJumperBet(),
        "power 0.5" => new PowerBet(0.5),
        "mixture" => new PowerMixtureBet(),
        "plug-in" => new PlugInBet(),
        "plug-in, 20 bins" => new PlugInBet(bins: 20),
        _ => throw new ArgumentOutOfRangeException(nameof(name)),
    };

    private static readonly AlarmProcedure[] Procedures = [AlarmProcedure.Ville, AlarmProcedure.Cusum, AlarmProcedure.ShiryaevRoberts];

    // ---- The statistics against their definitions ----

    [Theory]
    [InlineData(1e6)]
    [InlineData(3.0)]
    [InlineData(1.2)]
    public void EachStatistic_IsItsDefinition_ReadFromThePath(double threshold)
    {
        // Ville's ln(Sₙ/S_τ), CUSUM's ln max over τ ≤ i < n of Sₙ/Sᵢ and Shiryaev–Roberts's ln Σ Sₙ/Sᵢ over the same i, with τ the last
        // alarm (0 at the start), from the path of ln S. Low thresholds make them restart often.
        var rng = new Xoshiro256PlusPlus(seed: 6);
        foreach (string name in new[] { "jumper", "power 0.5", "mixture" })
        {
            var martingale = new ConformalTestMartingale(Bet(name), 1, [.. Procedures.Select(p => new MartingaleAlarm(p, threshold))]);
            var path = new List<double> { 0.0 };   // ln S₀
            int[] since = [0, 0, 0];
            for (int n = 1; n <= 600; n++)
            {
                MartingaleStep step = martingale.Add(Gaussian(rng) - (n > 300 ? 1.0 : 0.0));
                path.Add(step.LogCapital);
                for (int a = 0; a < 3; a++)
                {
                    MartingaleAlarm alarm = martingale.Alarms[a];
                    bool fired = step.Fired!.Contains(alarm.Procedure);
                    double definition = alarm.Procedure switch
                    {
                        AlarmProcedure.Ville => path[n] - path[since[a]],
                        AlarmProcedure.Cusum => path[n] - Enumerable.Range(since[a], n - since[a]).Min(i => path[i]),
                        _ => LogSumExp(Enumerable.Range(since[a], n - since[a]).Select(i => path[n] - path[i])),
                    };
                    Assert.Equal(definition >= Math.Log(threshold), fired);
                    if (fired) since[a] = n;
                    else Assert.True(Math.Abs(alarm.LogStatistic - definition) <= 1e-9 * Math.Max(1.0, Math.Abs(definition)),
                                     $"{name}, {alarm.Procedure} at {n}: {alarm.LogStatistic:R} against {definition:R}");
                }
            }
            output.WriteLine($"{name}, c {threshold}: alarms {string.Join(", ", martingale.Alarms.Select(a => $"{a.Procedure} {a.AlarmTimes.Count}"))}");
        }
    }

    private static double LogSumExp(IEnumerable<double> values)
    {
        double[] v = [.. values];
        double top = v.Max();
        return top + Math.Log(v.Sum(x => Math.Exp(x - top)));
    }

    // ---- Ville's guarantee, and CUSUM's and Shiryaev–Roberts's ----

    [Theory]
    [InlineData("jumper", 10.0)]
    [InlineData("jumper", 4.0)]
    [InlineData("power 0.5", 10.0)]
    [InlineData("mixture", 10.0)]
    [InlineData("plug-in, 20 bins", 10.0)]
    public void Ville_FalseAlarms_AtMostOneOverC(string name, double c)
    {
        // 10,000 streams of 300 i.i.d. uniform p-values, M-1's conclusion: the share of streams with any alarm.
        const int streams = 10_000, length = 300;
        var rng = new Xoshiro256PlusPlus(seed: 20261108);
        int alarmed = 0;
        for (int r = 0; r < streams; r++)
        {
            var capital = new BettingCapital(Bet(name));
            var alarm = new MartingaleAlarm(AlarmProcedure.Ville, c);
            for (int n = 0; n < length; n++)
                if (alarm.Step(capital.Add(1.0 - rng.NextDouble()).LogBet)) { alarmed++; break; }
        }
        double share = alarmed / (double)streams, se = Math.Sqrt(share * (1 - share) / streams);
        output.WriteLine($"{name}, c {c}: P(any alarm in {length}) {share:F4} (se {se:F4}) against 1/c = {1 / c:F4}");
        Assert.True(share <= 1.0 / c + 3.0 * se);
    }

    [Theory]
    [InlineData(AlarmProcedure.Cusum, "jumper", 20.0)]
    [InlineData(AlarmProcedure.ShiryaevRoberts, "jumper", 20.0)]
    [InlineData(AlarmProcedure.Cusum, "power 0.5", 20.0)]
    [InlineData(AlarmProcedure.ShiryaevRoberts, "power 0.5", 20.0)]
    [InlineData(AlarmProcedure.ShiryaevRoberts, "mixture", 50.0)]
    public void CusumAndShiryaevRoberts_AlarmFalselyNoSoonerThanC_Posed(AlarmProcedure procedure, string name, double c)
    {
        // Posed, not proved: the mean time to the first false alarm is at least c, and over one long stream the alarms, repeated, come at
        // most once in c arrivals. Runs capped at 100c count as the cap, which only lowers the mean.
        var rng = new Xoshiro256PlusPlus(seed: 77);
        const int runs = 1000;
        int cap = (int)(100 * c);
        double total = 0.0;
        for (int r = 0; r < runs; r++)
        {
            var capital = new BettingCapital(Bet(name));
            var alarm = new MartingaleAlarm(procedure, c);
            int n = 0;
            while (n < cap && !alarm.Step(capital.Add(1.0 - rng.NextDouble()).LogBet)) n++;
            total += Math.Min(n + 1, cap);
        }
        var longCapital = new BettingCapital(Bet(name));
        var longAlarm = new MartingaleAlarm(procedure, c);
        const int length = 200_000;
        for (int n = 0; n < length; n++) longAlarm.Step(longCapital.Add(1.0 - rng.NextDouble()).LogBet);
        double mean = total / runs, frequency = longAlarm.AlarmTimes.Count / (double)length;
        output.WriteLine($"{procedure}, {name}, c {c}: mean time to the first false alarm {mean:F1}; alarms per arrival over {length:N0}: {frequency:F5} against 1/c = {1 / c:F5}");
        Assert.True(mean >= c, $"mean {mean:F1}");
        Assert.True(frequency <= 1.0 / c, $"frequency {frequency:F5}");
        Assert.False(longAlarm.Guarantee.Proved);
    }

    [Fact]
    public void TheMartingale_ReadsOut_ItsAlarms_AndTheirGuarantees()
    {
        var martingale = new ConformalTestMartingale(new SimpleJumperBet(), 1,
                                                     new MartingaleAlarm(AlarmProcedure.Ville, 20.0), new MartingaleAlarm(AlarmProcedure.Cusum, 20.0));
        var rng = new Xoshiro256PlusPlus(seed: 2);
        var fired = new List<(int, AlarmProcedure)>();
        for (int n = 1; n <= 800; n++)
            foreach (AlarmProcedure p in martingale.Add(Gaussian(rng) - (n > 300 ? 2.0 : 0.0)).Fired!) fired.Add((n, p));
        Assert.Equal(martingale.Alarms[0].AlarmTimes.Count + martingale.Alarms[1].AlarmTimes.Count, fired.Count);
        Assert.Contains(fired, f => f.Item2 == AlarmProcedure.Ville && f.Item1 > 300);
        Assert.True(martingale.Alarms[0].Guarantee.Proved);
        Assert.Contains("1/c = 0.05", martingale.Alarms[0].Guarantee.Statement);
        Assert.Contains("posed", martingale.Alarms[1].Guarantee.Statement);
        output.WriteLine($"Ville at {string.Join(", ", martingale.Alarms[0].AlarmTimes)}; CUSUM at {string.Join(", ", martingale.Alarms[1].AlarmTimes)}");
        Assert.Throws<ArgumentOutOfRangeException>(() => new MartingaleAlarm(AlarmProcedure.Ville, 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new MartingaleAlarm((AlarmProcedure)7, 10.0));
        Assert.Empty(new ConformalTestMartingale(new PowerBet(0.5), 1).Add(0.0).Fired!);
    }

    [Fact]
    public void ALostCapital_StaysLostForVille_AndCusumAndShiryaevRobertsReadPastIt()
    {
        MartingaleAlarm[] alarms = [.. Procedures.Select(p => new MartingaleAlarm(p, 10.0))];
        foreach (MartingaleAlarm alarm in alarms) alarm.Step(double.NegativeInfinity);
        foreach (MartingaleAlarm alarm in alarms) for (int n = 0; n < 5; n++) alarm.Step(Math.Log(2.0));
        Assert.Equal(double.NegativeInfinity, alarms[0].LogStatistic);
        Assert.Single(alarms[1].AlarmTimes);   // 2⁴ = 16 ≥ 10 from the step after the loss
        Assert.Single(alarms[2].AlarmTimes);
        Assert.Throws<ArgumentException>(() => alarms[0].Step(double.NaN));
    }

    // ---- With the qualification switch: delays, a fixed fit's scores, and the cost ----

    [Theory]
    [InlineData("mean")]
    [InlineData("noise")]
    public void TheDelayToAnAlarm_AfterAChange(string change)
    {
        // 300 streams: 500 exchangeable arrivals, then 1,000 after a change. The scores are a fixed N(0, 1) fit's log density at y, with
        // y ~ N(0, 1) before; after, y ~ N(1, 1) (the mean) or N(0, 1.6²) (the noise). Each bet with each procedure at c 100: the share of
        // streams alarmed before the change (false), and after it the share detected and the median delay. A reading, indicative.
        if (!BassQualification.Enabled) return;
        const int streams = 300, before = 500, after = 1000;
        const double c = 100.0;
        foreach (string name in new[] { "jumper", "power 0.5", "mixture", "plug-in" })
        {
            var rows = new List<string>();
            foreach (AlarmProcedure procedure in Procedures)
            {
                var rng = new Xoshiro256PlusPlus(seed: 31);
                int falseAlarms = 0;
                var delays = new List<int>();
                for (int r = 0; r < streams; r++)
                {
                    var martingale = new ConformalTestMartingale(Bet(name), r, new MartingaleAlarm(procedure, c));
                    for (int n = 1; n <= before + after; n++)
                    {
                        double y = n <= before ? Gaussian(rng) : change == "mean" ? 1.0 + Gaussian(rng) : 1.6 * Gaussian(rng);
                        if (martingale.Add(-0.5 * y * y).Fired!.Count == 0) continue;
                        if (n <= before) falseAlarms++;
                        else delays.Add(n - before);
                        break;
                    }
                }
                delays.Sort();
                string median = delays.Count == 0 ? "none" : delays[delays.Count / 2].ToString();
                rows.Add($"{procedure}: false {falseAlarms / (double)streams:F3}, detected {delays.Count / (double)streams:F3}, median delay {median}");
            }
            output.WriteLine($"{change}, {name}: {string.Join("; ", rows)}");
        }
    }

    // f(x) = sin(2πx) + x/2 with noise 0.2ζ: the training set proper and the stream before the change; after it the noise is 0.5ζ.
    private static double Truth(double x, double sd, Xoshiro256PlusPlus rng) => Math.Sin(2.0 * Math.PI * x) + 0.5 * x + sd * Gaussian(rng);

    [Theory]
    [InlineData("bars")]
    [InlineData("bass")]
    public void OnAFixedFitsScores_VilleHolds_AndAChangeInTheNoiseAlarms(string engine)
    {
        // A fit on 200 points; 100 streams of 600 units scored by its predictive density, the noise rising from 0.2 to 0.5 at unit 300.
        // Simple Jumper with each procedure at c 20: the share of streams alarmed before the change, at most 1/20 for Ville, and after it
        // the share detected and the median delay.
        if (!BassQualification.Enabled) return;
        const int streams = 100, length = 600, change = 300;
        var rng = new Xoshiro256PlusPlus(seed: 20261108);
        double[] x = [.. Enumerable.Range(0, 200).Select(_ => rng.NextDouble())];
        Array.Sort(x);
        double[] y = [.. x.Select(v => Truth(v, 0.2, rng))];
        double[] xs = [.. Enumerable.Range(0, streams * length).Select(_ => rng.NextDouble())];
        double[] ys = [.. xs.Select((v, i) => Truth(v, i % length < change ? 0.2 : 0.5, rng))];
        double[][] grids = [.. ys.Select(v => new[] { v })];
        CandidateLogDensity units;
        if (engine == "bars")
        {
            var ensemble = new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0));
            BarsResult r = ensemble.Run(x, y, x, chains: 2, masterSeed: 1, burn: 500, samples: 100);
            units = ensemble.CandidateLogDensity(r, xs, grids);
        }
        else
        {
            var data = new BassData([new ContinuousPredictor("x", x)], new ScalarResponse(y));
            BassFit fit = Bass.Fit(BassStudentTTests.Config(null, burnIn: 1000, samples: 100), data);
            units = fit.CandidateLogDensity(data, PredictionCases.Of([new ContinuousPredictor("x", xs)]), grids);
        }
        double[] scores = StreamScores.Density(units);
        var rows = new List<string>();
        foreach (AlarmProcedure procedure in Procedures)
        {
            int falseAlarms = 0;
            var delays = new List<int>();
            for (int r = 0; r < streams; r++)
            {
                var martingale = new ConformalTestMartingale(new SimpleJumperBet(), r, new MartingaleAlarm(procedure, 20.0));
                for (int n = 0; n < length; n++)
                {
                    if (martingale.Add(scores[r * length + n]).Fired!.Count == 0) continue;
                    if (n < change) falseAlarms++;
                    else delays.Add(n - change + 1);
                    break;
                }
            }
            delays.Sort();
            rows.Add($"{procedure}: false {falseAlarms / (double)streams:F2}, detected {delays.Count / (double)streams:F2}, median delay {(delays.Count == 0 ? "none" : delays[delays.Count / 2].ToString())}");
            if (procedure == AlarmProcedure.Ville) Assert.True(falseAlarms / (double)streams <= 0.05 + 3.0 * Math.Sqrt(0.05 * 0.95 / streams));
        }
        output.WriteLine($"{engine}, S {units.Draws}: {string.Join("; ", rows)}");
    }

    [Fact]
    public void CostPerArrival()
    {
        // An arrival with its rank, Simple Jumper and the three alarms, over 100,000 Gaussian scores after a warm-up; and the three alarms'
        // recursions alone on the same bets, since the arrival's own cost swings with the machine's load. Cycle clock, indicative.
        if (!BassQualification.Enabled) return;
        var rng = new Xoshiro256PlusPlus(seed: 4);
        double[] scores = [.. Enumerable.Range(0, 100_000).Select(_ => Gaussian(rng))];
        for (int pass = 0; pass < 2; pass++)
        {
            var martingale = new ConformalTestMartingale(new SimpleJumperBet(), 1, [.. Procedures.Select(p => new MartingaleAlarm(p, 100.0))]);
            var bets = new double[scores.Length];
            long t0 = BassCost.Tick();
            for (int i = 0; i < scores.Length; i++) bets[i] = martingale.Add(scores[i]).LogBet;
            long t1 = BassCost.Tick();
            MartingaleAlarm[] alone = [.. Procedures.Select(p => new MartingaleAlarm(p, 100.0))];
            long t2 = BassCost.Tick();
            foreach (double bet in bets) foreach (MartingaleAlarm alarm in alone) alarm.Step(bet);
            long t3 = BassCost.Tick();
            if (pass == 0) continue;
            output.WriteLine($"({(BassCost.Cycles ? "cycles" : "wall")} clock, indicative): {(t1 - t0) / BassCost.TicksPerSecond * 1e9 / scores.Length:F0} ns per arrival "
                             + $"with the three alarms; the three alarms' recursions alone {(t3 - t2) / BassCost.TicksPerSecond * 1e9 / scores.Length:F0} ns");
        }
    }
}
