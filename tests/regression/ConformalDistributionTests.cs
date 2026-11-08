using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The conformal predictive distribution (E16, slice 4). First the reading the slice turns on: whether the in-sample predictive
/// CDF, F̃(y) = Σ_s w̃_s(y)·F(y | x, θ_s) under the add-one-in weights w̃_s(y) ∝ W_s·p(y | x, θ_s), stays monotone in the candidate's
/// y, Vovk et al.'s condition for π(y) to be a CDF. It does not in general (ConformalPrediction.md, slice 4's reading).
/// </summary>
public sealed class ConformalDistributionTests(ITestOutputHelper output)
{
    /// <summary>A CDF score held by both tails: Σ w̃·F and Σ w̃·(1 − F), each to full relative precision where it is small.</summary>
    internal readonly record struct TailScore(double Lower, double Upper)
    {
        /// <summary>The sign of this score less <paramref name="other"/>, each side read on the tail that holds its digits.</summary>
        public int CompareTo(TailScore other)
        {
            bool low = Lower <= 0.5, otherLow = other.Lower <= 0.5;
            if (low && otherLow) return Lower.CompareTo(other.Lower);
            if (!low && !otherLow) return other.Upper.CompareTo(Upper);
            return low ? -1 : 1;
        }

        /// <summary>This score less <paramref name="other"/>, from the tail that holds both.</summary>
        public double Minus(TailScore other) => Lower <= 0.5 && other.Lower <= 0.5 ? Lower - other.Lower : other.Upper - Upper;
    }

    /// <summary>
    /// The CDF scores by direct sums: at each of <paramref name="count"/> candidate columns from <paramref name="first"/>, every unit's
    /// F̃ᵢ = Σ_s w̃_s F(y_i | θ_s) over <paramref name="unitLogCdf"/>'s columns <paramref name="members"/>, and the candidate's own
    /// F̃_c = Σ_s w̃_s F(y | θ_s) from <paramref name="logCdf"/>, the weights from the candidate's log density <paramref name="lc"/>.
    /// </summary>
    internal static (TailScore[][] Units, TailScore[] Own) Reference(double[,] unitLogCdf, int[] members, double[,] lc, double[,] logCdf, int first,
                                                                     int count, double[]? logWeights = null)
    {
        int s = lc.GetLength(0), n = members.Length;
        var lower = new double[s, n];
        var upper = new double[s, n];
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++)
            {
                double l = unitLogCdf[d, members[i]];
                lower[d, i] = Math.Exp(l);
                upper[d, i] = -SpecialFunctions.Expm1(l);
            }
        var units = new TailScore[count][];
        var own = new TailScore[count];
        var w = new double[s];
        var a = new double[n];
        var b = new double[n];
        for (int k = 0; k < count; k++)
        {
            int c = first + k;
            double top = double.NegativeInfinity;
            for (int d = 0; d < s; d++) top = Math.Max(top, (logWeights?[d] ?? 0.0) + lc[d, c]);
            double total = 0.0;
            for (int d = 0; d < s; d++) total += w[d] = Math.Exp((logWeights?[d] ?? 0.0) + lc[d, c] - top);
            Array.Clear(a);
            Array.Clear(b);
            double ownLower = 0.0, ownUpper = 0.0;
            for (int d = 0; d < s; d++)
            {
                double wd = w[d] / total, l = logCdf[d, c];
                ownLower += wd * Math.Exp(l);
                ownUpper += wd * -SpecialFunctions.Expm1(l);
                for (int i = 0; i < n; i++)
                {
                    a[i] += wd * lower[d, i];
                    b[i] += wd * upper[d, i];
                }
            }
            own[k] = new TailScore(ownLower, ownUpper);
            units[k] = [.. a.Zip(b, (x, y) => new TailScore(x, y))];
        }
        return (units, own);
    }

    /// <summary>
    /// Monotonicity along a grid: the candidate's score (Vovk's first condition), each unit's (his second: falling as the candidate
    /// rises), and the counts that make π(y), #{F̃ᵢ &lt; F̃_c}/m and (#{F̃ᵢ ≤ F̃_c} + 1)/m, which a CDF needs rising. Scores are
    /// compared on the tail that holds their digits; a fall of a score under 10⁻¹² of its tail is taken as rounding.
    /// </summary>
    internal sealed record Monotonicity(int OwnDrops, double OwnDrop, double OwnDropAt, int UnitsRising, double UnitRise, int CountDrops, int CountDrop,
                                        double CountDropAt, int Units)
    {
        public static Monotonicity Read(TailScore[][] units, TailScore[] own)
        {
            const double rounding = 1e-12;
            int g = own.Length, n = units[0].Length, ownDrops = 0, countDrops = 0, countDrop = 0;
            double ownDrop = 0.0, ownAt = double.NaN, unitRise = 0.0, countAt = double.NaN;
            var rising = new bool[n];
            int Below(int j) => units[j].Count(a => a.CompareTo(own[j]) < 0);
            int Upto(int j) => units[j].Count(a => a.CompareTo(own[j]) <= 0) + 1;
            int lastBelow = Below(0), lastUpto = Upto(0);
            for (int j = 1; j < g; j++)
            {
                double drop = own[j - 1].Minus(own[j]);
                if (drop > rounding * Math.Min(own[j].Lower, own[j].Upper) && drop > 0.0)
                {
                    ownDrops++;
                    if (drop > ownDrop) { ownDrop = drop; ownAt = own[j - 1].Lower; }
                }
                for (int i = 0; i < n; i++)
                {
                    double rise = units[j][i].Minus(units[j - 1][i]);
                    if (rise > rounding * Math.Min(units[j][i].Lower, units[j][i].Upper) && rise > 0.0) { rising[i] = true; unitRise = Math.Max(unitRise, rise); }
                }
                int below = Below(j), upto = Upto(j), fall = Math.Max(lastBelow - below, lastUpto - upto);
                if (fall > 0) { countDrops++; if (fall > countDrop) { countDrop = fall; countAt = own[j - 1].Lower; } }
                lastBelow = below;
                lastUpto = upto;
            }
            return new Monotonicity(ownDrops, ownDrop, ownAt, rising.Count(r => r), unitRise, countDrops, countDrop, countAt, n);
        }
    }

    [Fact]
    public void TailScore_ReadsEachSideOnTheTailThatHoldsItsDigits()
    {
        // Two scores within 10⁻²⁰ of 1 differ in their upper tails, which a probability would round to one value.
        var a = new TailScore(1.0, 1e-20);
        var b = new TailScore(1.0, 2e-20);
        Assert.True(b.CompareTo(a) < 0);
        Assert.Equal(1e-20, a.Minus(b), 1e-35);
        Assert.True(new TailScore(1e-30, 1.0).CompareTo(new TailScore(2e-30, 1.0)) < 0);
        Assert.True(new TailScore(0.4, 0.6).CompareTo(new TailScore(0.6, 0.4)) < 0);
        Assert.Equal(0, new TailScore(0.7, 0.3).CompareTo(new TailScore(0.7000000000000001, 0.3)));
    }

    [Fact]
    public void TwoDraws_TheScoreFallsBetweenThemOnceTheyStandMoreThan1_75σApart()
    {
        // Two equal draws N(±a, 1): at the midpoint the score's slope is φ(a) − a(Φ(a) − ½), the weights' shift outrunning the
        // density, negative past a ≈ 0.877, locations 1.75σ apart.
        static double Slope(double a)
        {
            double[,] lc = new double[2, 2], logCdf = new double[2, 2];
            double[] at = [-1e-4, 1e-4];
            for (int d = 0; d < 2; d++)
                for (int g = 0; g < 2; g++)
                {
                    double mu = d == 0 ? -a : a;
                    lc[d, g] = ConformalBayesTests.LogNormal(at[g], mu);
                    logCdf[d, g] = Normal.StandardLogCdf(at[g] - mu);
                }
            (_, TailScore[] own) = Reference(new double[2, 1], [0], lc, logCdf, 0, 2);
            return own[1].Minus(own[0]) / 2e-4;
        }
        static double Closed(double a) => Math.Exp(-0.5 * a * a) / Math.Sqrt(2 * Math.PI) - a * (Math.Exp(Normal.StandardLogCdf(a)) - 0.5);
        foreach (double a in new[] { 0.5, 0.8, 0.86, 0.89, 1.0, 2.0 })
        {
            Assert.Equal(Closed(a), Slope(a), 6);
            Assert.Equal(a < 0.877, Slope(a) > 0);
        }
    }

    // ---- the reading ----------------------------------------------------------------------------------------------------------------

    [Fact]
    public void Monotonicity_OnTheConjugateModel_UnimodalDrawsHold_TwoChainsApartDoNot()
    {
        // The reading's quick form, 20 cases a kind: under the exact posterior and a chain stopped early, which sit in one mode, neither
        // the score nor π falls; under two chains stuck apart both do.
        List<Monotonicity>[] readings = Reading(20, ["exact, S 50", "exact, S 400", "a chain stopped early", "two chains stuck apart"]);
        for (int k = 0; k < 3; k++)
        {
            Assert.All(readings[k], r => Assert.Equal(0, r.OwnDrops));
            Assert.All(readings[k], r => Assert.Equal(0, r.CountDrops));
        }
        Assert.Contains(readings[3], r => r.OwnDrops > 0 && r.CountDrops > 0);
        // Vovk's second condition fails all the same: the units on the far side of x = 0 from the candidate rise with it.
        Assert.Contains(readings[1], r => r.UnitsRising > 0);
    }

    [Fact]
    public void Monotonicity_OnTheConjugateModel_TheFullReading()
    {
        if (!BassQualification.Enabled) return;
        Reading(100, ["exact, S 50", "exact, S 400", "exact, S 4000", "a chain stopped early", "two chains stuck apart"]);
    }

    // y | x, β ~ N(βx, 1), β ~ N(0, 10), n = 19, data from a truth the model does not share, a test case at x_t; the grid's step 0.05
    // over [−25, 25]. Draws: the exact posterior, a chain stopped early, or two chains stuck apart (the exact draws shifted by ∓1 in turn).
    private List<Monotonicity>[] Reading(int replications, string[] kinds)
    {
        const int n = 19;
        double[] grid = [.. Enumerable.Range(0, 1001).Select(g => -25.0 + g / 20.0)];
        int[] members = [.. Enumerable.Range(0, n)];
        int[] seeds = SeedTree.Derive(2026_1109, replications);
        var all = new List<Monotonicity>[kinds.Length];
        for (int k = 0; k < kinds.Length; k++)
        {
            string kind = kinds[k];
            var readings = all[k] = [];
            for (int r = 0; r < replications; r++)
            {
                var rng = new Xoshiro256PlusPlus(seeds[r]);
                var x = new double[n];
                var y = new double[n];
                for (int i = 0; i < n; i++) { x[i] = 4 * rng.NextDouble() - 2; y[i] = ConformalBayesTests.Truth(x[i], 0.5, rng); }
                double xt = 4 * rng.NextDouble() - 2;
                double[] beta = kind switch
                {
                    "exact, S 50" => ConformalBayesTests.ExactDraws(x, y, 50, rng),
                    "exact, S 400" => ConformalBayesTests.ExactDraws(x, y, 400, rng),
                    "exact, S 4000" => ConformalBayesTests.ExactDraws(x, y, 4000, rng),
                    "a chain stopped early" => SplitConformalTests.StoppedEarly(x, y, rng),
                    _ => [.. ConformalBayesTests.ExactDraws(x, y, 400, rng).Select((b, d) => b + (d % 2 == 0 ? -1.0 : 1.0))],
                };
                int s = beta.Length;
                var unitLogCdf = new double[s, n];
                var lc = new double[s, grid.Length];
                var logCdf = new double[s, grid.Length];
                for (int d = 0; d < s; d++)
                {
                    for (int i = 0; i < n; i++) unitLogCdf[d, i] = Normal.StandardLogCdf(y[i] - beta[d] * x[i]);
                    for (int g = 0; g < grid.Length; g++)
                    {
                        lc[d, g] = ConformalBayesTests.LogNormal(grid[g], beta[d] * xt);
                        logCdf[d, g] = Normal.StandardLogCdf(grid[g] - beta[d] * xt);
                    }
                }
                (TailScore[][] units, TailScore[] own) = Reference(unitLogCdf, members, lc, logCdf, 0, grid.Length);
                readings.Add(Monotonicity.Read(units, own));
            }
            Report(output, kind, readings, n + 1);
        }
        return all;
    }

    internal static void Report(ITestOutputHelper output, string kind, List<Monotonicity> readings, int m)
    {
        int cases = readings.Count;
        Monotonicity worstOwn = readings.MaxBy(r => r.OwnDrop)!, worstCount = readings.MaxBy(r => r.CountDrop)!;
        output.WriteLine($"{kind}: {cases} cases, m {m}");
        output.WriteLine($"  own score falls on {readings.Count(r => r.OwnDrops > 0)} cases, {readings.Sum(r => r.OwnDrops)} steps; largest {worstOwn.OwnDrop:G3} from F {worstOwn.OwnDropAt:F3}");
        output.WriteLine($"  counts fall on {readings.Count(r => r.CountDrops > 0)} cases, {readings.Sum(r => r.CountDrops)} steps; largest {worstCount.CountDrop} of m {m} from F {worstCount.CountDropAt:F3}");
        output.WriteLine($"  units rising with the candidate {readings.Sum(r => r.UnitsRising)} of {readings.Sum(r => r.Units)}; largest step {readings.Max(r => r.UnitRise):G3}");
    }
}
