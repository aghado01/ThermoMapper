using System;
using System.Linq;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Samplers.Tests;

/// <summary>
/// A self-sizing ladder's policy (<see cref="TemperingLadder"/>): rungs placed on a known barrier profile take equal
/// increments of it at any count, and at the current count the placement is <see cref="TemperingLadder.Equalize"/>; the
/// sizing rule's grow, deadband and clamp cases; each placed rung's source, against a brute-force reading of the
/// barrier; and the split rule (F-M2's lever 4b) on synthetic rejection profiles with a known steep band.
/// </summary>
public sealed class TemperingLadderTests
{
    private readonly ITestOutputHelper _out;
    public TemperingLadderTests(ITestOutputHelper output) => _out = output;

    // A barrier profile with a ramp in ln b and a step near b = 0.05, as a transition gives: Λ(b) from the hot end 1e-3.
    private static double Profile(double b)
        => 0.8 * Math.Log(b / 1e-3) + 1.0 * (Erf((Math.Log(b) - Math.Log(0.05)) / 0.4) - Erf((Math.Log(1e-3) - Math.Log(0.05)) / 0.4));

    private static double Erf(double x) => 2.0 * Distributions.Normal.Standard.Cdf(x * Math.Sqrt(2.0)) - 1.0;

    [Theory]
    [InlineData(5)]
    [InlineData(9)]
    [InlineData(17)]
    [InlineData(33)]
    public void Place_OnAKnownBarrier_GivesEqualIncrements(int levels)
    {
        // The current ladder: 64 geometric rungs, each pair's rejection the profile's increment across it.
        double[] inverseTemperatures = ParallelTempering<int>.GeometricLadder(64, 1e-3);
        double[] rejection = [.. Enumerable.Range(0, 63).Select(i => Profile(inverseTemperatures[i]) - Profile(inverseTemperatures[i + 1]))];
        Assert.All(rejection, r => Assert.InRange(r, 0.0, 1.0));
        double total = Profile(1.0);

        double[] placed = TemperingLadder.Place(inverseTemperatures, rejection, levels);
        Assert.Equal(levels, placed.Length);
        Assert.Equal((1.0, 1e-3), (placed[0], placed[^1]));
        double step = total / (levels - 1), worst = 0.0;
        for (int i = 0; i + 1 < levels; i++)
        {
            Assert.True(placed[i + 1] < placed[i]);
            worst = Math.Max(worst, Math.Abs(Profile(placed[i]) - Profile(placed[i + 1]) - step) / step);
        }
        _out.WriteLine($"{levels} rungs: Λ = {total:F3}, equal steps of {step:F3} to within {100 * worst:F2}%");
        Assert.True(worst < 0.02, $"{100 * worst:F2}%");
    }

    [Fact]
    public void Place_AtTheCurrentCount_IsEqualize()
    {
        double[] inverseTemperatures = [1.0, 0.6, 0.3, 0.1, 0.03, 0.0];
        double[][] rates = [[0.1, 0.2, 0.9, 0.3, 0.6], [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.4, 0.0, 0.0, 0.2], [1.0, 1.0, 1.0, 1.0, 1.0]];
        foreach (double[] r in rates)
        {
            double[] equalized = TemperingLadder.Equalize(inverseTemperatures, r);
            Assert.Equal(equalized, TemperingLadder.Place(inverseTemperatures, r, inverseTemperatures.Length));
            Assert.Equal((1.0, 0.0), (equalized[0], equalized[^1]));   // the ends stay, b = 0 among them
        }
        // No rejection: the current count comes back unchanged.
        Assert.Equal(inverseTemperatures, TemperingLadder.Place(inverseTemperatures, rates[1], inverseTemperatures.Length));
        // Two rungs, at their count, are left alone; at another, they are placed between their ends.
        Assert.Equal([1.0, 0.1], TemperingLadder.Place([1.0, 0.1], [0.5], 2));
        double[] three = TemperingLadder.Place([1.0, 0.1], [0.5], 3);
        Assert.Equal((1.0, 0.1), (three[0], three[2]));
        Assert.Equal(0.55, three[1], 12);   // one pair: the cubic is the straight line
    }

    [Fact]
    public void Place_WithoutAReading_PlacesAnotherCountByTheRungsIndex()
    {
        double[] inverseTemperatures = [1.0, 0.5, 0.2, 0.05];
        double[] grown = TemperingLadder.Place(inverseTemperatures, [0.0, 0.0, 0.0], 7);
        Assert.Equal(7, grown.Length);
        for (int i = 0; i < inverseTemperatures.Length; i++) Assert.Equal(inverseTemperatures[i], grown[2 * i], 12);   // every old rung, a new one between
        for (int i = 1; i < 7; i++) Assert.True(grown[i] < grown[i - 1]);
        double[] shrunk = TemperingLadder.Place(inverseTemperatures, [0.0, 0.0, 0.0], 2);
        Assert.Equal([1.0, 0.05], shrunk);
    }

    [Fact]
    public void Place_RefusesWhatItCannotPlace()
    {
        Assert.Throws<ArgumentNullException>(() => TemperingLadder.Place(null!, [], 3));
        Assert.Equal("rejection", Assert.Throws<ArgumentException>(() => TemperingLadder.Place([1.0, 0.5], [0.1, 0.2], 3)).ParamName);
        Assert.Equal("rejection", Assert.Throws<ArgumentOutOfRangeException>(() => TemperingLadder.Place([1.0, 0.5], [1.5], 3)).ParamName);
        Assert.Equal("rejection", Assert.Throws<ArgumentOutOfRangeException>(() => TemperingLadder.Place([1.0, 0.5], [double.NaN], 3)).ParamName);
        Assert.Equal("levels", Assert.Throws<ArgumentOutOfRangeException>(() => TemperingLadder.Place([1.0, 0.5], [0.1], 1)).ParamName);
        Assert.Equal("inverseTemperatures", Assert.Throws<ArgumentException>(() => TemperingLadder.Place([1.0], [], 3)).ParamName);
        Assert.Equal([1.0], TemperingLadder.Equalize([1.0], []));   // a lone rung, at its count, as Equalize always took it
    }

    // (Λ̂, current, minimum, maximum) → the count: N ≈ 2Λ asks for ⌈2Λ̂⌉ + 1 rungs.
    [Theory]
    [InlineData(7.0, 8, 3, 64, 15)]      // a saturated ladder grows at once
    [InlineData(7.0, 15, 3, 64, 15)]     // and stops when the reading stops asking
    [InlineData(7.2, 15, 3, 64, 16)]     // one more asked, one more taken
    [InlineData(3.0, 8, 3, 64, 8)]       // 7 asked: within the deadband
    [InlineData(2.5, 8, 3, 64, 6)]       // 6 asked: past it, shrinks to the reading
    [InlineData(2.0, 8, 3, 64, 5)]       // ⌈4⌉ + 1
    [InlineData(2.0000001, 8, 3, 64, 6)] // ⌈4.0000002⌉ + 1
    [InlineData(40.0, 8, 3, 64, 64)]     // the cap binds
    [InlineData(40.0, 80, 3, 64, 64)]    // a ladder above the cap comes down to it
    [InlineData(0.0, 8, 3, 64, 3)]       // no barrier: down to the minimum
    [InlineData(0.0, 3, 3, 64, 3)]
    [InlineData(0.0, 3, 2, 64, 3)]       // 2 asked, one fewer: within the deadband
    [InlineData(0.0, 4, 2, 64, 2)]
    [InlineData(1.0, 2, 2, 64, 3)]
    [InlineData(10.0, 8, 4, 4, 4)]       // a fixed range
    public void Size_GrowsAtOnce_ShrinksPastTheDeadband_AndClamps(double barrier, int current, int minimum, int maximum, int expected)
        => Assert.Equal(expected, TemperingLadder.Size(barrier, current, minimum, maximum));

    [Fact]
    public void Size_RefusesWhatIsNoLadder()
    {
        Assert.Equal("barrier", Assert.Throws<ArgumentOutOfRangeException>(() => TemperingLadder.Size(double.NaN, 8, 3, 64)).ParamName);
        Assert.Equal("barrier", Assert.Throws<ArgumentOutOfRangeException>(() => TemperingLadder.Size(-0.1, 8, 3, 64)).ParamName);
        Assert.Equal("barrier", Assert.Throws<ArgumentOutOfRangeException>(() => TemperingLadder.Size(double.PositiveInfinity, 8, 3, 64)).ParamName);
        Assert.Equal("current", Assert.Throws<ArgumentOutOfRangeException>(() => TemperingLadder.Size(1.0, 1, 2, 64)).ParamName);
        Assert.Equal("minimum", Assert.Throws<ArgumentOutOfRangeException>(() => TemperingLadder.Size(1.0, 8, 1, 64)).ParamName);
        Assert.Equal("maximum", Assert.Throws<ArgumentOutOfRangeException>(() => TemperingLadder.Size(1.0, 8, 5, 4)).ParamName);
    }

    [Fact]
    public void Sources_AreTheNearestRungsOnTheBarrier()
    {
        double[] inverseTemperatures = [1.0, 0.5, 0.2, 0.05];
        double[] rejection = [0.15, 0.4, 0.5];   // Λ from the hot end: 0, 0.5, 0.9, 1.05
        double[] cumulative = [1.05, 0.9, 0.5, 0.0];   // each current rung's, cold first
        foreach (int levels in new[] { 2, 3, 6, 11 })
        {
            double[] placed = TemperingLadder.Place(inverseTemperatures, rejection, levels);
            int[] sources = TemperingLadder.Sources(inverseTemperatures, rejection, placed);
            Assert.Equal(0, sources[0]);
            Assert.Equal(inverseTemperatures.Length - 1, sources[^1]);
            for (int k = 0; k < levels; k++)
            {
                // The placed rung k sits at Λ = Λ̂ (K − 1 − k)/(K − 1) by construction.
                double at = 1.05 * (levels - 1 - k) / (levels - 1);
                int nearest = Enumerable.Range(0, inverseTemperatures.Length).MinBy(i => Math.Abs(at - cumulative[i]));
                Assert.Equal(nearest, sources[k]);
            }
        }
        Assert.Equal([0, 1, 2, 3], TemperingLadder.Sources(inverseTemperatures, rejection, inverseTemperatures));

        // Without a reading the nearest in b decides, and a tie there goes to the colder rung.
        Assert.Equal([0, 2, 3, 3], TemperingLadder.Sources(inverseTemperatures, [0.0, 0.0, 0.0], [1.0, 0.3, 0.1, 0.05]));
        Assert.Equal([0, 1, 2], TemperingLadder.Sources([1.0, 0.5, 0.25], [0.0, 0.0], [1.0, 0.375, 0.25]));
        // Beyond the current ends the barrier holds its end's value: a rung below the hot end starts from the hottest.
        Assert.Equal([0, 2], TemperingLadder.Sources([1.0, 0.5, 0.2], [0.3, 0.4], [1.0, 0.01]));
        // A flat pair ties on the barrier: the nearer in b takes it.
        Assert.Equal([0, 0, 1, 2], TemperingLadder.Sources([1.0, 0.5, 0.25], [0.0, 0.5], [1.0, 0.875, 0.625, 0.25]));
        Assert.Throws<ArgumentException>(() => TemperingLadder.Sources(inverseTemperatures, [0.1], inverseTemperatures));
        Assert.Throws<ArgumentException>(() => TemperingLadder.Sources(inverseTemperatures, rejection, [1.0, double.NaN]));
        Assert.Throws<ArgumentNullException>(() => TemperingLadder.Sources(inverseTemperatures, rejection, null!));
    }

    [Fact]
    public void Split_PlacesEachSteepPairsMidpoint_AndKeepsEveryRung()
    {
        // Equal spacing and equal rejections make the cumulative barrier linear in b, so each midpoint is the pair's own.
        double[] ladder = [1.0, 0.8, 0.6, 0.4, 0.2];
        double[] split = TemperingLadder.Split(ladder, [0.95, 0.95, 0.95, 0.95], 0.9, 64);
        Assert.Equal(9, split.Length);
        double[] expected = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2];
        for (int i = 0; i < expected.Length; i++) Assert.Equal(expected[i], split[i], 1e-12);
        Assert.All(ladder, b => Assert.Contains(b, split));

        // A steep pair between gentle ones splits inside itself, on the monotone cubic, every current rung kept.
        double[] uneven = [1.0, 0.7, 0.5, 0.3, 0.1];
        double[] one = TemperingLadder.Split(uneven, [0.2, 0.97, 0.3, 0.1], 0.9, 64);
        Assert.Equal(6, one.Length);
        Assert.All(uneven, b => Assert.Contains(b, one));
        double added = one.Except(uneven).Single();
        Assert.InRange(added, 0.5 + 1e-9, 0.7 - 1e-9);
        for (int i = 1; i < one.Length; i++) Assert.True(one[i] < one[i - 1]);

        // With a reference rung the pair above b = 0 splits as any other, and b = 0 stays the last rung.
        double[] reference = TemperingLadder.Split([1.0, 0.5, 0.1, 0.0], [0.3, 0.4, 0.99], 0.9, 64);
        Assert.Equal(5, reference.Length);
        Assert.Equal(0.0, reference[^1]);
        Assert.InRange(reference[3], 1e-12, 0.1 - 1e-12);
    }

    [Fact]
    public void Split_TheSteepestFirst_NeverPastTheMaximum()
    {
        double[] ladder = [1.0, 0.7, 0.5, 0.3, 0.1];
        double[] rejection = [0.2, 0.95, 0.99, 0.92];
        double[] six = TemperingLadder.Split(ladder, rejection, 0.9, 6);   // room for one: the 0.99 pair
        Assert.Equal(6, six.Length);
        Assert.InRange(six.Except(ladder).Single(), 0.3 + 1e-9, 0.5 - 1e-9);
        double[] seven = TemperingLadder.Split(ladder, rejection, 0.9, 7);   // then the 0.95 pair
        Assert.Equal(7, seven.Length);
        Assert.Single(seven.Except(ladder), b => b > 0.5 && b < 0.7);
        Assert.Equal(8, TemperingLadder.Split(ladder, rejection, 0.9, 64).Length);
        // Equal rejections tie: the colder pair first.
        double[] tie = TemperingLadder.Split([1.0, 0.8, 0.6, 0.4], [0.95, 0.95, 0.95], 0.9, 5);
        Assert.InRange(tie.Except([1.0, 0.8, 0.6, 0.4]).Single(), 0.8 + 1e-9, 1.0 - 1e-9);
        // At the maximum, or above it, nothing splits.
        Assert.Equal(ladder, TemperingLadder.Split(ladder, rejection, 0.9, 5));
        Assert.Equal(ladder, TemperingLadder.Split(ladder, rejection, 0.9, 3));
    }

    [Fact]
    public void Split_LeavesALadderWithNoSteepPairUnchanged()
    {
        double[] ladder = [1.0, 0.6, 0.3, 0.1, 0.03];
        foreach (double[] r in new[] { new[] { 0.1, 0.2, 0.9, 0.3 }, [0.0, 0.0, 0.0, 0.0], [0.89, 0.9, 0.5, 0.2] })
        {
            double[] split = TemperingLadder.Split(ladder, r, 0.9, 64);
            Assert.NotSame(ladder, split);
            Assert.Equal(ladder.Select(BitConverter.DoubleToInt64Bits), split.Select(BitConverter.DoubleToInt64Bits));
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => TemperingLadder.Split(ladder, [0.1, 0.2, 0.9, 0.3], 1.0, 64));
        Assert.Throws<ArgumentOutOfRangeException>(() => TemperingLadder.Split(ladder, [0.1, 0.2, 0.9, 0.3], 0.0, 64));
        Assert.Throws<ArgumentOutOfRangeException>(() => TemperingLadder.Split(ladder, [0.1, 0.2, 0.9, 0.3], 0.9, 1));
        Assert.Throws<ArgumentException>(() => TemperingLadder.Split(ladder, [0.1, 0.2], 0.9, 64));
    }

    // A barrier with a steep band: a gentle ramp in ln b from the hot end 0.01, and a step of 6 near b = 0.2, which a pair's
    // censored rejection r = 1 − e^(−ΔΛ) reads as at most 1.
    private static double SteepBand(double b)
        => 0.4 * Math.Log(b / 0.01) + 3.0 * (Erf((Math.Log(b) - Math.Log(0.2)) / 0.05) - Erf((Math.Log(0.01) - Math.Log(0.2)) / 0.05));

    private static double[] Censored(double[] ladder)
        => [.. Enumerable.Range(0, ladder.Length - 1).Select(i => 1.0 - Math.Exp(-(SteepBand(ladder[i]) - SteepBand(ladder[i + 1]))))];

    [Theory]
    [InlineData(64)]
    [InlineData(12)]
    public void Split_OverRounds_ResolvesAKnownSteepBand_WithinTheMaximum(int maximum)
    {
        double[] ladder = ParallelTempering<int>.GeometricLadder(8, 0.01);
        double[] r = Censored(ladder);
        double before = TemperingLadder.Barrier(r), truth = SteepBand(1.0);
        Assert.True(r.Max() > 0.99);
        int rounds = 0;
        while (r.Max() > 0.9 && rounds < 30)
        {
            double[] next = TemperingLadder.Split(ladder, r, 0.9, maximum);
            if (next.Length == ladder.Length) break;
            Assert.All(ladder, b => Assert.Contains(b, next));
            Assert.True(next.Length <= maximum);
            ladder = next;
            r = Censored(ladder);
            rounds++;
        }
        double after = TemperingLadder.Barrier(r);
        int inBand = ladder.Count(b => b > 0.18 && b < 0.22);
        _out.WriteLine($"max {maximum}: {rounds} rounds to {ladder.Length} rungs, {inBand} in the band; Σr {before:F2} → {after:F2} against Λ = {truth:F2}; " +
                       $"largest r {r.Max():F3}");
        Assert.True(ladder.Length <= maximum);
        Assert.True(after > before, "the finer pairs read more of the barrier");
        if (maximum == 64)
        {
            Assert.True(r.Max() <= 0.9, $"the band still rejects {r.Max():F3}");
            Assert.True(inBand >= 4, $"{inBand} rungs in the band");
        }
        else Assert.Equal(maximum, ladder.Length);
        // A gentle profile, every pair below the threshold, is left as it was.
        double[] gentle = ParallelTempering<int>.GeometricLadder(8, 0.5);
        double[] g = Censored(gentle);
        Assert.True(g.Max() <= 0.9);
        Assert.Equal(gentle, TemperingLadder.Split(gentle, g, 0.9, 64));
    }
}
