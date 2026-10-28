using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// Wallenius' variable-set probability against brute force over the K! draw orders, at R BASS's default weights and at
/// weights spanning five decades, K = 1 … 5; the probabilities over every K-subset summing to one; F-O5's values; the
/// two evaluations agreeing; R's alternating sum corrected; the boundary cases; and the sampler drawing the law.
/// </summary>
public sealed class WalleniusTests
{
    // F-O5's weights: R's w₂ = 5 plus appearance counts, and weights spanning five decades.
    private static readonly double[] RLike = { 5, 8, 12, 5, 6, 17, 5, 9 };
    private static readonly double[] FiveDecades = { 1e-3, 1e-2, 1e-1, 1, 10, 100 };

    public static TheoryData<string, int> WeightsAndOrders()
    {
        var data = new TheoryData<string, int>();
        foreach (string name in new[] { nameof(RLike), nameof(FiveDecades) })
            for (int k = 1; k <= 5; k++) data.Add(name, k);
        return data;
    }

    private static double[] Weights(string name) => name == nameof(RLike) ? RLike : FiveDecades;

    // F-O5's reference: Σ over the orders of S of Π z / the weight still in the urn, each denominator a sum of
    // positives, so every term keeps its relative precision.
    private static double BruteForce(double[] z, int[] set)
    {
        double total = 0.0;
        foreach (int[] order in Permutations(set))
        {
            double product = 1.0;
            var drawn = new HashSet<int>();
            foreach (int u in order)
            {
                double remaining = 0.0;
                for (int k = 0; k < z.Length; k++)
                    if (!drawn.Contains(k)) remaining += z[k];
                product *= z[u] / remaining;
                drawn.Add(u);
            }
            total += product;
        }
        return total;
    }

    private static IEnumerable<int[]> Permutations(int[] items)
    {
        if (items.Length <= 1) { yield return items; yield break; }
        for (int i = 0; i < items.Length; i++)
            foreach (int[] rest in Permutations(items.Where((_, j) => j != i).ToArray()))
                yield return new[] { items[i] }.Concat(rest).ToArray();
    }

    private static IEnumerable<int[]> Subsets(int p, int k)
    {
        if (k == 0) { yield return Array.Empty<int>(); yield break; }
        for (int last = k - 1; last < p; last++)
            foreach (int[] head in Subsets(last, k - 1))
                yield return head.Append(last).ToArray();
    }

    private static void Relative(double expected, double actual, double tolerance)
        => Assert.True(Math.Abs(actual / expected - 1.0) <= tolerance, $"expected {expected:R}, got {actual:R}");

    [Theory]
    [MemberData(nameof(WeightsAndOrders))]
    public void Probability_MatchesBruteForce_OverEveryDrawOrder(string weights, int k)
    {
        double[] z = Weights(weights);
        foreach (int[] set in Subsets(z.Length, k))
        {
            double brute = BruteForce(z, set);
            Relative(brute, Wallenius.Probability(z, set), 1e-14);
            Relative(brute, Wallenius.Probability(z, set, WalleniusMethod.Quadrature), 1e-13);
        }
    }

    [Theory]
    [MemberData(nameof(WeightsAndOrders))]
    public void Probabilities_SumToOne_OverEveryKSubset(string weights, int k)
    {
        double[] z = Weights(weights);
        foreach (WalleniusMethod method in new[] { WalleniusMethod.DrawOrders, WalleniusMethod.Quadrature })
            Assert.Equal(1.0, Subsets(z.Length, k).Sum(set => Wallenius.Probability(z, set, method)), 13);
    }

    [Fact]
    public void EqualWeights_GiveOneOverTheSubsetCount()
    {
        // p = 50, where R's alternating sum has lost three digits by K = 3, and p = 14 at K = 11, past Auto's
        // recursion (which then takes the quadrature) but within the explicit one's range.
        double[] z = Enumerable.Repeat(5.0, 50).ToArray();
        foreach ((int[] set, double count) in new[] { (new[] { 0, 7 }, 1225.0), (new[] { 3, 11, 40 }, 19600.0), (new[] { 1, 2, 30, 49 }, 230300.0) })
            foreach (WalleniusMethod method in Enum.GetValues<WalleniusMethod>())
                Relative(1.0 / count, Wallenius.Probability(z, set, method), 1e-13);

        double[] fourteen = Enumerable.Repeat(1.0, 14).ToArray();
        int[] eleven = Enumerable.Range(2, 11).ToArray();
        Assert.Equal(Wallenius.LogProbability(fourteen, eleven, WalleniusMethod.Quadrature), Wallenius.LogProbability(fourteen, eleven));
        foreach (WalleniusMethod method in Enum.GetValues<WalleniusMethod>())
            Relative(1.0 / 364.0, Wallenius.Probability(fourteen, eleven, method), 1e-13);
    }

    [Fact]
    public void TheTwoEvaluations_Agree_OnWeightsSpanningTwelveDecades()
    {
        double[] z = Enumerable.Range(0, 12).Select(i => Math.Pow(10.0, i - 6)).ToArray();
        for (int k = 2; k <= 6; k++)
            foreach (int[] set in Subsets(z.Length, k).Where((_, i) => i % 7 == 0))
                Assert.Equal(Wallenius.LogProbability(z, set, WalleniusMethod.Quadrature),
                    Wallenius.LogProbability(z, set, WalleniusMethod.DrawOrders), 12);
    }

    [Fact]
    public void FO5Values_AreReproduced_WithRsAlternatingSumCorrected()
    {
        // bass_wallenius_oracle.R's brute force (subsets 0-based here), at each order's subset where R's dmwnchBass
        // lost the most digits.
        Relative(1.6200834368726581e-09, Wallenius.Probability(FiveDecades, new[] { 0, 1 }), 1e-14);
        Relative(4.3783879363938001e-12, Wallenius.Probability(FiveDecades, new[] { 0, 1, 2 }), 1e-14);
        Relative(1.8268244113956862e-12, Wallenius.Probability(FiveDecades, new[] { 0, 1, 2, 4 }), 1e-14);
        Relative(8.8143129500747435e-14, Wallenius.Probability(FiveDecades, new[] { 0, 1, 2, 3, 4 }), 1e-14);
        Relative(0.014562302185529254, Wallenius.Probability(RLike, new[] { 0, 4 }), 1e-14);
        Relative(0.0061293753481082083, Wallenius.Probability(RLike, new[] { 0, 3, 7 }), 1e-14);
        Relative(0.0026213344330012957, Wallenius.Probability(RLike, new[] { 0, 1, 3, 4 }), 1e-14);
        Relative(0.0014654355017197087, Wallenius.Probability(RLike, new[] { 0, 1, 3, 4, 6 }), 1e-14);

        // dmwnchBass there: 1.6200836228108528e-09 (1.1·10⁻⁷ off), 4.3782755199117673e-12 (2.6·10⁻⁵),
        // 1.8274270985330077e-12 (3.3·10⁻⁴), 8.8817841970012523e-14 (7.7·10⁻³). In R's range at R's weights it is
        // right, and the recursion agrees with it to its own error: 0.014562302185529008.
        Relative(0.014562302185529008, Wallenius.Probability(RLike, new[] { 0, 4 }), 1e-13);

        // At K = 1 dmwnchBass returns twice the probability (0.35820895522388074), a path R never takes; z/Σz here.
        Relative(0.17910447761194029, Wallenius.Probability(RLike, new[] { 2 }), 1e-15);
        foreach (WalleniusMethod method in Enum.GetValues<WalleniusMethod>())
            Assert.Equal(Wallenius.LogProbability(RLike, new[] { 2 }), Wallenius.LogProbability(RLike, new[] { 2 }, method));
    }

    [Fact]
    public void Boundaries_TheWholeSet_ProbabilitiesNearOne_AndProbabilitiesBelowUnderflow()
    {
        Assert.Equal(0.0, Wallenius.LogProbability(FiveDecades, new[] { 5, 0, 3, 1, 4, 2 }));
        Assert.Equal(0.0, Wallenius.LogProbability(new[] { 2.0 }, new[] { 0 }));

        // K = 1 near one keeps the log's digits: ln(1e10 / (1e10 + 3)) = −ln(1 + 3·10⁻¹⁰), by its series.
        Relative(-(3e-10 - 4.5e-20 + 9e-30), Wallenius.LogProbability(new[] { 1.0, 1e10, 2.0 }, new[] { 1 }), 1e-14);

        // Two dominant weights drawn first with probability near one.
        double[] heavy = { 1e8, 2e8, 1.0, 3.0 };
        foreach (WalleniusMethod method in Enum.GetValues<WalleniusMethod>())
            Relative(BruteForce(heavy, new[] { 0, 1 }), Wallenius.Probability(heavy, new[] { 0, 1 }, method), 1e-14);

        // P(S) ≈ 2·10⁻²⁷¹ against brute force, and P(S) ≈ 10⁻³⁵⁰, below double range, against K!·Πα, whose next term
        // is O(α) = 10⁻⁷¹: both evaluations keep the log's digits.
        double[] tiny = { 1e-90, 1e-90, 1e-90, 1.0, 2.0 };
        double brute = BruteForce(tiny, new[] { 0, 1, 2 });
        Assert.InRange(brute, 1e-272, 1e-270);
        double[] tinier = { 1e-70, 1e-70, 1e-70, 1e-70, 1e-70, 1.0, 2.0 };
        double leading = Math.Log(120.0) + 5.0 * Math.Log(1e-70 / 3.0);
        foreach (WalleniusMethod method in Enum.GetValues<WalleniusMethod>())
        {
            Relative(Math.Log(brute), Wallenius.LogProbability(tiny, new[] { 0, 1, 2 }, method), 1e-14);
            Relative(leading, Wallenius.LogProbability(tinier, new[] { 0, 1, 2, 3, 4 }, method), 1e-14);
        }
    }

    [Fact]
    public void Draws_FollowTheLaw()
    {
        var rng = new Xoshiro256PlusPlus(seed: 11);
        foreach ((double[] z, int k) in new[] { (RLike, 2), (RLike, 3), (FiveDecades, 2) })
        {
            const int draws = 200_000;
            var counts = new Dictionary<string, int>();
            var set = new int[k];
            for (int d = 0; d < draws; d++)
            {
                Wallenius.Sample(rng, z, set);
                Assert.Equal(k, set.Distinct().Count());
                string key = string.Join(",", set.OrderBy(v => v));
                counts[key] = counts.GetValueOrDefault(key) + 1;
            }
            foreach (int[] s in Subsets(z.Length, k))
            {
                double p = Wallenius.Probability(z, s);
                double frequency = counts.GetValueOrDefault(string.Join(",", s)) / (double)draws;
                Assert.True(Math.Abs(frequency - p) <= 4.5 * Math.Sqrt(p * (1.0 - p) / draws) + 1.0 / draws,
                    $"{{{string.Join(",", s)}}}: frequency {frequency} against {p}");
            }
        }
    }

    [Fact]
    public void InvalidInputs_AreRefused()
    {
        Assert.Throws<ArgumentException>(() => Wallenius.LogProbability(RLike, new[] { 1, 1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => Wallenius.LogProbability(RLike, new[] { 1, 8 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => Wallenius.LogProbability(RLike, Array.Empty<int>()));
        Assert.Throws<ArgumentOutOfRangeException>(() => Wallenius.LogProbability(new[] { 1.0, 0.0, 2.0 }, new[] { 0, 2 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => Wallenius.LogProbability(new[] { 1.0, double.PositiveInfinity }, new[] { 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => Wallenius.Sample(new Xoshiro256PlusPlus(seed: 1), RLike, new int[9]));
    }
}
