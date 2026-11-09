using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The bets and the capital of a conformal test martingale (E35's slice 6b, M-2).
/// <list type="bullet">
/// <item>Each bet integrates to 1, the plug-in to at most 1, by quadrature at states after random histories.</item>
/// <item>Each bet reads only the past: its curve before pₙ is a fresh bet's fed p₁, …, pₙ₋₁.</item>
/// <item>Simple Jumper is Vovk et al.'s Algorithm 1, and the power mixture is its integral over ε. The plug-in's incremental density is
/// its refit from the counts.</item>
/// <item>The capital has expectation 1 on exchangeable streams, end to end through the ranks.</item>
/// <item>A bet that peeks at the current p, or integrates above 1, breaks it.</item>
/// </list>
/// With the qualification switch, the cost per arrival.
/// </summary>
public sealed class BettingCapitalTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    public static TheoryData<string> Bets => ["power 0.2", "power 0.9", "mixture", "jumper", "jumper 0.3", "plug-in", "plug-in, floor 0.1"];

    private static IBettingFunction Bet(string name) => name switch
    {
        "power 0.2" => new PowerBet(0.2),
        "power 0.9" => new PowerBet(0.9),
        "mixture" => new PowerMixtureBet(),
        "jumper" => new SimpleJumperBet(),
        "jumper 0.3" => new SimpleJumperBet(0.3),
        "plug-in" => new PlugInBet(),
        "plug-in, floor 0.1" => new PlugInBet(floor: 0.1),
        _ => throw new ArgumentOutOfRangeException(nameof(name)),
    };

    // p-values from a history: uniform, or piled toward 0 (u²), which a drift gives.
    private static double HistoryP(Xoshiro256PlusPlus rng, bool skewed)
    {
        double u = 1.0 - rng.NextDouble();
        return skewed ? Math.Max(u * u, 1e-12) : u;
    }

    private static double Simpson(Func<double, double> f, double a, double b, int n)
    {
        double h = (b - a) / n, sum = f(a) + f(b);
        for (int k = 1; k < n; k++) sum += f(a + k * h) * (k % 2 == 1 ? 4.0 : 2.0);
        return sum * h / 3.0;
    }

    // ∫ from a = e^(−700) to 1 of a bet, in x = −ln p, where every bet here is smooth or piecewise smooth: fine steps over x ≤ 12, which
    // holds the plug-in's grid, coarse ones beyond.
    private static double MassAboveTail(IBettingFunction bet)
    {
        double G(double x) => Math.Exp(bet.LogBet(Math.Exp(-x)) - x);
        return Simpson(G, 0.0, 12.0, 240_000) + Simpson(G, 12.0, 700.0, 200_000);
    }

    [Theory]
    [MemberData(nameof(Bets))]
    public void EachBet_IntegratesToOne_OrAtMostOne_AfterAnyHistory(string name)
    {
        // The mass below a = e^(−700): a^ε for a power bet, (a − 1)/ln a for the mixture (∫₀¹ a^ε dε), at most a·max f for the bounded ones.
        var rng = new Xoshiro256PlusPlus(seed: 20261108 + name.Length);
        double a = Math.Exp(-700.0), worst = 0.0;
        foreach (int history in new[] { 0, 7, 10, 40, 333 })
            foreach (bool skewed in new[] { false, true })
            {
                IBettingFunction bet = Bet(name);
                for (int i = 0; i < history; i++) bet.Observe(HistoryP(rng, skewed));
                double tail = bet switch { PowerBet power => Math.Pow(a, power.Epsilon), PowerMixtureBet => (a - 1.0) / Math.Log(a), _ => 0.0 };
                double total = MassAboveTail(bet) + tail;
                worst = Math.Max(worst, Math.Abs(total - 1.0));
                if (bet.IntegratesToOne) Assert.True(Math.Abs(total - 1.0) < 1e-8, $"{name} after {history} ({(skewed ? "skewed" : "uniform")}): ∫ f = {total:R}");
                else Assert.True(total <= 1.0 + 1e-8 && total > 1.0 - 1e-6, $"{name} after {history} ({(skewed ? "skewed" : "uniform")}): ∫ f = {total:R}");
            }
        output.WriteLine($"{name}: the largest |∫ f − 1| {worst:E2}");
    }

    [Theory]
    [MemberData(nameof(Bets))]
    public void EachBet_ReadsOnlyThePast(string name)
    {
        // At each step the bet's curve, read before pₙ is observed, is a fresh bet's fed p₁, …, pₙ₋₁; reading it moves nothing.
        double[] grid = [1e-9, 1e-4, 0.01, 0.13, 0.5, 0.77, 0.999, 1.0];
        var rng = new Xoshiro256PlusPlus(seed: 7);
        double[] ps = [.. Enumerable.Range(0, 60).Select(i => HistoryP(rng, i > 30))];
        IBettingFunction running = Bet(name);
        for (int n = 0; n < ps.Length; n++)
        {
            double[] before = [.. grid.Select(running.LogBet)];
            double[] again = [.. grid.Reverse().Select(running.LogBet).Reverse()];
            Assert.Equal(before, again);
            if (n % 9 == 0)
            {
                IBettingFunction fresh = Bet(name);
                for (int i = 0; i < n; i++) fresh.Observe(ps[i]);
                Assert.Equal(before, grid.Select(fresh.LogBet).ToArray());
            }
            running.Observe(ps[n]);
        }
    }

    [Fact]
    public void SimpleJumper_IsVovksAlgorithmOne()
    {
        // Algorithm 1 as Vovk et al. write it: unnormalized capitals C_ε, a jump share J of their sum C, then each times its bet.
        foreach (double jump in new[] { 0.01, 0.3, 0.0, 1.0 })
        {
            var rng = new Xoshiro256PlusPlus(seed: 11);
            var capital = new BettingCapital(new SimpleJumperBet(jump));
            double[] c = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
            for (int n = 0; n < 400; n++)
            {
                double p = HistoryP(rng, n > 200), total = c.Sum();
                for (int e = 0; e < 3; e++) c[e] = (1.0 - jump) * c[e] + jump / 3.0 * total;
                for (int e = 0; e < 3; e++) c[e] *= 1.0 + (e - 1) * (p - 0.5);
                double logCapital = capital.Add(p).LogCapital, expected = Math.Log(c.Sum());
                Assert.True(Math.Abs(logCapital - expected) <= 1e-12 * Math.Max(1.0, Math.Abs(expected)), $"J {jump}, step {n}: {logCapital:R} against {expected:R}");
            }
        }
    }

    [Fact]
    public void ThePowerMixture_IsTheIntegralOfThePowerBets()
    {
        // ∫₀¹ εp^(ε−1) dε by Gauss–Legendre in ε, at p from 10⁻³⁰⁰ to 1, across the series' switch at p = e⁻¹.
        var rule = Maths.Distributions.GaussLegendre.Of(64);
        var mixture = new PowerMixtureBet();
        foreach (double p in new[] { 1e-300, 1e-30, 1e-6, 0.01, 0.3, Math.Exp(-1.0) * (1 - 1e-12), Math.Exp(-1.0), 0.5, 0.99, 1.0 - 1e-12, 1.0 })
        {
            // Panels geometric toward ε = 1, where εp^(ε−1) peaks at tiny p, each by the 64-point rule.
            double integral = 0.0;
            double[] edges = [0.0, 0.5, 0.75, 0.875, 0.9375, 0.96875, 0.984375, 0.9921875, 0.99609375, 0.998046875, 1.0];
            for (int k = 0; k + 1 < edges.Length; k++)
            {
                double lo = edges[k], hi = edges[k + 1], half = 0.5 * (hi - lo), mid = 0.5 * (hi + lo);
                for (int j = 0; j < rule.Count; j++)
                {
                    double eps = mid + half * rule.Nodes[j];
                    integral += half * rule.Weights[j] * eps * Math.Exp((eps - 1.0) * Math.Log(p));
                }
            }
            Assert.Equal(1.0, Math.Exp(mixture.LogBet(p)) / integral, 11);
        }
    }

    [Fact]
    public void ThePlugIn_IsItsKernelDensityRefittedFromTheCounts()
    {
        // The incremental density, after refits at 10, 20, 40, 80 and 160, against the binned kernel density rebuilt from scratch at its
        // last bandwidth, over its exact piecewise-linear integral and scaled down by 10⁻¹².
        var rng = new Xoshiro256PlusPlus(seed: 3);
        var bet = new PlugInBet(warmUp: 10, bins: 50);
        var counts = new int[50];
        for (int n = 0; n < 250; n++)
        {
            double p = HistoryP(rng, skewed: true);
            bet.Observe(p);
            counts[Math.Min((int)(p * 50), 49)]++;
        }
        double h = bet.Bandwidth;
        double Kernel(double x, double c) => Math.Exp(-0.5 * Math.Pow((x - c) / h, 2)) + Math.Exp(-0.5 * Math.Pow((x + c) / h, 2)) + Math.Exp(-0.5 * Math.Pow((x - 2 + c) / h, 2));
        double[] d = [.. Enumerable.Range(0, 50).Select(g => Enumerable.Range(0, 50).Sum(b => counts[b] * Kernel((g + 0.5) / 50, (b + 0.5) / 50)))];
        double integral = 0.5 / 50 * (d[0] + d[49]) + Enumerable.Range(0, 49).Sum(g => 0.5 / 50 * (d[g] + d[g + 1]));
        foreach (double p in new[] { 0.001, 0.01, 0.0137, 0.5, 0.987, 1.0 })
        {
            double u = p * 50 - 0.5;
            double dp = u <= 0 ? d[0] : u >= 49 ? d[49] : (1 - (u - Math.Floor(u))) * d[(int)u] + (u - Math.Floor(u)) * d[(int)u + 1];
            Assert.Equal(Math.Log((1.0 - 1e-12) * dp / integral), bet.LogBet(p), 11);
        }
        Assert.InRange(h, 1.0 / 50, 0.5);
    }

    [Theory]
    [InlineData("power 0.9")]
    [InlineData("jumper")]
    [InlineData("plug-in")]
    public void TheCapital_HasExpectationOne_OnExchangeableStreams(string name)
    {
        // End to end: Gaussian scores ranked into smoothed p-values, 20 arrivals, 100,000 streams. The mixture and the sharper power bets
        // have infinite variance, so their mean is checked by the integral above instead. The plug-in counts in 20 bins from a warm-up of 5.
        const int streams = 100_000, length = 20;
        var rng = new Xoshiro256PlusPlus(seed: 20261108);
        double sum = 0.0, square = 0.0;
        for (int r = 0; r < streams; r++)
        {
            IBettingFunction bet = name == "plug-in" ? new PlugInBet(warmUp: 5, bins: 20) : Bet(name);
            var martingale = new ConformalTestMartingale(bet, seed: r);
            for (int n = 0; n < length; n++) martingale.Add(Gaussian(rng));
            double s = martingale.EValue;
            sum += s;
            square += s * s;
        }
        double mean = sum / streams, se = Math.Sqrt((square / streams - mean * mean) / streams);
        output.WriteLine($"{name}: E S_{length} {mean:F4} (se {se:F4}) over {streams:N0} exchangeable streams");
        Assert.True(Math.Abs(mean - 1.0) < 4.0 * se, $"{name}: mean {mean:F4}, se {se:F4}");
    }

    // A bet that chooses its half after seeing p: f = 2 on the half holding pₙ. Each half integrates to 1, but the choice reads pₙ.
    private sealed class PeekingBet : IBettingFunction
    {
        public double LogBet(double p) => Math.Log(2.0);
        public void Observe(double p) { }
        public bool IntegratesToOne => true;   // as it claims; the choice is the breach
    }

    // A bet integrating to 1.1: f = 1.1 everywhere.
    private sealed class OverBet : IBettingFunction
    {
        public double LogBet(double p) => Math.Log(1.1);
        public void Observe(double p) { }
        public bool IntegratesToOne => true;
    }

    [Fact]
    public void ABetThatPeeks_OrIntegratesAboveOne_BreaksTheBound()
    {
        // Over 2,000 exchangeable streams of 40, an honest bet's capital reaches 10 in at most about a tenth of them (Ville); one that
        // reads the current p reaches it in every stream, and one integrating to 1.1 too.
        const int streams = 2000, length = 40;
        var rng = new Xoshiro256PlusPlus(seed: 5);
        int honest = 0, peeking = 0, over = 0;
        for (int r = 0; r < streams; r++)
        {
            var a = new ConformalTestMartingale(new SimpleJumperBet(), r);
            var b = new ConformalTestMartingale(new PeekingBet(), r);
            var c = new ConformalTestMartingale(new OverBet(), r);
            bool ra = false, rb = false, rc = false;
            for (int n = 0; n < length; n++)
            {
                double score = Gaussian(rng);
                ra |= a.Add(score).LogCapital >= Math.Log(10.0);
                rb |= b.Add(score).LogCapital >= Math.Log(10.0);
                rc |= c.Add(score).LogCapital >= Math.Log(10.0);
            }
            honest += ra ? 1 : 0;
            peeking += rb ? 1 : 0;
            over += rc ? 1 : 0;
        }
        output.WriteLine($"P(S reaches 10 within {length}): Simple Jumper {honest / (double)streams:F4}, peeking {peeking / (double)streams:F4}, ∫ f = 1.1 {over / (double)streams:F4}");
        Assert.True(honest / (double)streams <= 0.1);
        Assert.Equal(streams, peeking);
        Assert.Equal(streams, over);
    }

    [Fact]
    public void TheMartingale_ReadsOut_ItsPath_ItsBets_AndItsEValue()
    {
        var rng = new Xoshiro256PlusPlus(seed: 8);
        var martingale = new ConformalTestMartingale(new PowerBet(0.5), seed: 2);
        var capital = new BettingCapital(new PowerBet(0.5));
        double sum = 0.0;
        for (int n = 0; n < 50; n++)
        {
            MartingaleStep step = martingale.Add(Gaussian(rng));
            sum += step.LogBet;
            Assert.Equal(sum, step.LogCapital);
            Assert.Equal(capital.Add(step.PValue.P).LogCapital, step.LogCapital);
            Assert.Equal(Math.Log(0.5) - 0.5 * Math.Log(step.PValue.P), step.LogBet);
        }
        Assert.Equal(50, martingale.LogPath.Count);
        Assert.Equal(Math.Exp(sum), martingale.EValue);
        Assert.Throws<ArgumentOutOfRangeException>(() => capital.Add(0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => capital.Add(1.0 + 1e-16 * 3));
    }

    [Fact]
    public void EachBet_GrowsTheCapital_WhenTheStreamChanges()
    {
        // 200 exchangeable arrivals, then 200 whose scores fall by 1.5 sds: the late units conform less and their p-values pile toward 0.
        // The adaptive bets' capital grows by far more than 100. A fixed power bet has decayed under the null by then, by about
        // 1 − ε + ln ε per step, which is what CUSUM and Shiryaev–Roberts read past (6c). A reading, not a rate: 6c's delays measure it.
        var rise = new Dictionary<string, double>();
        foreach (string name in new[] { "power 0.2", "power 0.9", "mixture", "jumper", "plug-in" })
        {
            var rng = new Xoshiro256PlusPlus(seed: 12);
            var martingale = new ConformalTestMartingale(Bet(name), seed: 3);
            for (int n = 0; n < 400; n++) martingale.Add(Gaussian(rng) - (n >= 200 ? 1.5 : 0.0));
            double atChange = martingale.LogPath[199], peak = martingale.LogPath.Skip(200).Max();
            rise[name] = peak - atChange;
            output.WriteLine($"{name}: ln S {atChange:F2} at the change, its peak after {peak:F2}, {martingale.LogPath[^1]:F2} at 400");
        }
        Assert.True(rise["jumper"] > Math.Log(100.0) && rise["plug-in"] > Math.Log(100.0), $"Simple Jumper rose {rise["jumper"]:F2}, the plug-in {rise["plug-in"]:F2}");
    }

    [Fact]
    public void CostPerArrival()
    {
        // The rank and each bet, per arrival, over a stream of 100,000 Gaussian scores after a warm-up. Cycle clock, indicative.
        if (!BassQualification.Enabled) return;
        var rng = new Xoshiro256PlusPlus(seed: 4);
        double[] scores = [.. Enumerable.Range(0, 100_000).Select(_ => Gaussian(rng))];
        var warm = new ConformalTestMartingale(new PlugInBet(), seed: 1);
        for (int i = 0; i < 100_000; i++) warm.Add(scores[i]);
        foreach (string name in new[] { "power 0.2", "mixture", "jumper", "plug-in" })
        {
            var martingale = new ConformalTestMartingale(Bet(name), seed: 1);
            long t0 = BassCost.Tick();
            foreach (double score in scores) martingale.Add(score);
            output.WriteLine($"{name} ({(BassCost.Cycles ? "cycles" : "wall")} clock, indicative): {(BassCost.Tick() - t0) / BassCost.TicksPerSecond * 1e9 / scores.Length:F0} ns per arrival");
        }
    }
}
