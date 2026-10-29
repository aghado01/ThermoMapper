using System;
using System.Linq;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Samplers.Tests;

/// <summary>
/// A self-sizing ladder's policy (<see cref="TemperingLadder"/>): rungs placed on a known barrier profile take equal
/// increments of it at any count, and at the current count the placement is <see cref="TemperingLadder.Equalize"/>; the
/// sizing rule's grow, deadband and clamp cases; and each placed rung's source, against a brute-force reading of the
/// barrier.
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
        double[] betas = ParallelTempering<int>.GeometricLadder(64, 1e-3);
        double[] rejection = [.. Enumerable.Range(0, 63).Select(i => Profile(betas[i]) - Profile(betas[i + 1]))];
        Assert.All(rejection, r => Assert.InRange(r, 0.0, 1.0));
        double total = Profile(1.0);

        double[] placed = TemperingLadder.Place(betas, rejection, levels);
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
        double[] betas = [1.0, 0.6, 0.3, 0.1, 0.03, 0.0];
        double[][] rates = [[0.1, 0.2, 0.9, 0.3, 0.6], [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.4, 0.0, 0.0, 0.2], [1.0, 1.0, 1.0, 1.0, 1.0]];
        foreach (double[] r in rates)
        {
            double[] equalized = TemperingLadder.Equalize(betas, r);
            Assert.Equal(equalized, TemperingLadder.Place(betas, r, betas.Length));
            Assert.Equal((1.0, 0.0), (equalized[0], equalized[^1]));   // the ends stay, b = 0 among them
        }
        // No rejection: the current count comes back unchanged.
        Assert.Equal(betas, TemperingLadder.Place(betas, rates[1], betas.Length));
        // Two rungs, at their count, are left alone; at another, they are placed between their ends.
        Assert.Equal([1.0, 0.1], TemperingLadder.Place([1.0, 0.1], [0.5], 2));
        double[] three = TemperingLadder.Place([1.0, 0.1], [0.5], 3);
        Assert.Equal((1.0, 0.1), (three[0], three[2]));
        Assert.Equal(0.55, three[1], 12);   // one pair: the cubic is the straight line
    }

    [Fact]
    public void Place_WithoutAReading_PlacesAnotherCountByTheRungsIndex()
    {
        double[] betas = [1.0, 0.5, 0.2, 0.05];
        double[] grown = TemperingLadder.Place(betas, [0.0, 0.0, 0.0], 7);
        Assert.Equal(7, grown.Length);
        for (int i = 0; i < betas.Length; i++) Assert.Equal(betas[i], grown[2 * i], 12);   // every old rung, a new one between
        for (int i = 1; i < 7; i++) Assert.True(grown[i] < grown[i - 1]);
        double[] shrunk = TemperingLadder.Place(betas, [0.0, 0.0, 0.0], 2);
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
        Assert.Equal("betas", Assert.Throws<ArgumentException>(() => TemperingLadder.Place([1.0], [], 3)).ParamName);
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
        double[] betas = [1.0, 0.5, 0.2, 0.05];
        double[] rejection = [0.15, 0.4, 0.5];   // Λ from the hot end: 0, 0.5, 0.9, 1.05
        double[] cumulative = [1.05, 0.9, 0.5, 0.0];   // each current rung's, cold first
        foreach (int levels in new[] { 2, 3, 6, 11 })
        {
            double[] placed = TemperingLadder.Place(betas, rejection, levels);
            int[] sources = TemperingLadder.Sources(betas, rejection, placed);
            Assert.Equal(0, sources[0]);
            Assert.Equal(betas.Length - 1, sources[^1]);
            for (int k = 0; k < levels; k++)
            {
                // The placed rung k sits at Λ = Λ̂ (K − 1 − k)/(K − 1) by construction.
                double at = 1.05 * (levels - 1 - k) / (levels - 1);
                int nearest = Enumerable.Range(0, betas.Length).MinBy(i => Math.Abs(at - cumulative[i]));
                Assert.Equal(nearest, sources[k]);
            }
        }
        Assert.Equal([0, 1, 2, 3], TemperingLadder.Sources(betas, rejection, betas));

        // Without a reading the nearest in b decides, and a tie there goes to the colder rung.
        Assert.Equal([0, 2, 3, 3], TemperingLadder.Sources(betas, [0.0, 0.0, 0.0], [1.0, 0.3, 0.1, 0.05]));
        Assert.Equal([0, 1, 2], TemperingLadder.Sources([1.0, 0.5, 0.25], [0.0, 0.0], [1.0, 0.375, 0.25]));
        // A flat pair ties on the barrier: the nearer in b takes it.
        Assert.Equal([0, 0, 1, 2], TemperingLadder.Sources([1.0, 0.5, 0.25], [0.0, 0.5], [1.0, 0.875, 0.625, 0.25]));
        Assert.Throws<ArgumentException>(() => TemperingLadder.Sources(betas, [0.1], betas));
        Assert.Throws<ArgumentException>(() => TemperingLadder.Sources(betas, rejection, [1.0, double.NaN]));
        Assert.Throws<ArgumentNullException>(() => TemperingLadder.Sources(betas, rejection, null!));
    }
}
