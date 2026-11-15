using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The residual scores' grid-free sets (E16, slice 17's follow-up; RoBAS's Algorithm 1): each interval holds the grid form's accepted values,
/// its ends within the tolerance of the decision's change, read on both sides of each end; distance to origin's ends are the k-th largest
/// |r|; the intervals cover under a mean shift; k = 0 is the whole line, and RoBAS-EB below n = 4 is refused.
/// </summary>
public sealed class ResidualIntervalTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    [Fact]
    public void TheIntervalHoldsTheGridSet_AndEachEndIsTheDecisionsChange()
    {
        const double tolerance = 1e-9, centre = 0.7;
        var rng = new Xoshiro256PlusPlus(2026_1146);
        double[] grid = [.. Enumerable.Range(0, 16001).Select(g => -8.0 + 0.001 * g)];
        foreach (int n in new[] { 19, 49, 99 })
            foreach (double shift in new[] { 0.0, 1.5 })
            {
                double[] residuals = [.. Enumerable.Range(0, n).Select(_ => shift + Gaussian(rng))];
                foreach (ResidualScore score in Enum.GetValues<ResidualScore>())
                {
                    ResidualInterval interval = ResidualShrinkage.Intervals(residuals, [centre], score, 0.1, tolerance).Intervals[0];
                    ConformalSet set = ResidualShrinkage.Prediction(residuals, [centre], [grid], score, 0.1).Sets[0];
                    int first = -1, last = -1;
                    for (int g = 0; g < grid.Length; g++)
                        if (set.Accepts(g)) { if (first < 0) first = g; last = g; }
                    Assert.True(first > 0 && last < grid.Length - 1, "the grid should reach past both ends");
                    for (int g = first; g <= last; g++) Assert.True(set.Accepts(g), "the grid's accepted values should be one interval");
                    Assert.InRange(grid[first], interval.Lower, interval.Upper);
                    Assert.InRange(grid[last], interval.Lower, interval.Upper);
                    Assert.True(interval.Upper <= grid[last + 1] + tolerance && interval.Lower >= grid[first - 1] - tolerance);

                    // Each end read on both sides: accepted a tolerance inside, refused at the end.
                    ConformalSet ends = ResidualShrinkage.Prediction(residuals, [centre],
                        [[interval.Lower, interval.Lower + 2 * tolerance, interval.Upper - 2 * tolerance, interval.Upper]], score, 0.1).Sets[0];
                    Assert.Equal([false, true, true, false], [ends.Accepts(0), ends.Accepts(1), ends.Accepts(2), ends.Accepts(3)]);
                }
            }
    }

    [Fact]
    public void DistanceToOrigin_EndsAtTheKthLargestAbsoluteResidual()
    {
        var rng = new Xoshiro256PlusPlus(2026_1147);
        double[] residuals = [.. Enumerable.Range(0, 59).Select(_ => 0.3 + Gaussian(rng))];
        ResidualInterval interval = ResidualShrinkage.Intervals(residuals, [2.0], ResidualScore.DistanceToOrigin, 0.1, 1e-10).Intervals[0];
        double q = residuals.Select(Math.Abs).OrderDescending().ElementAt(5);   // k = ⌊0.1 · 60⌋ = 6: the 6th largest
        Assert.Equal((60, 6), (interval.Units, interval.Level));
        Assert.InRange(interval.Upper - 2.0 - q, 0.0, 1e-10);
        Assert.InRange(2.0 - q - interval.Lower, 0.0, 1e-10);
    }

    [Fact]
    public void TheIntervals_CoverUnderAMeanShift()
    {
        // Residuals N(1.5, 1), the predictor biased: 49 calibration residuals and a test residual, 2,000 replications, α 0.1.
        const int replications = 2000;
        var rng = new Xoshiro256PlusPlus(2026_1148);
        ResidualScore[] scores = Enum.GetValues<ResidualScore>();
        var covered = new int[scores.Length];
        var width = new double[scores.Length];
        for (int r = 0; r < replications; r++)
        {
            double[] residuals = [.. Enumerable.Range(0, 49).Select(_ => 1.5 + Gaussian(rng))];
            double test = 1.5 + Gaussian(rng);
            for (int f = 0; f < scores.Length; f++)
            {
                ResidualInterval interval = ResidualShrinkage.Intervals(residuals, [0.0], scores[f]).Intervals[0];
                if (interval.Contains(test)) covered[f]++;
                width[f] += interval.Width / replications;
            }
        }
        double se = Math.Sqrt(0.09 / replications);
        output.WriteLine(string.Join("; ", scores.Select((s, f) => $"{s}: coverage {covered[f] / (double)replications:F4}, mean width {width[f]:F3}")) + $"; se {se:F4}");
        for (int f = 0; f < scores.Length; f++) Assert.True(covered[f] / (double)replications >= 0.9 - 3.5 * se, $"{scores[f]}");
        Assert.True(width[2] < width[0], "under the shift RoBAS-EB should be narrower than distance to origin");
    }

    [Fact]
    public void TheWholeLine_AndWhatTheIntervalsDoNotServe()
    {
        ResidualInterval whole = ResidualShrinkage.Intervals([0.1, -0.2, 0.3], [1.0]).Intervals[0];
        Assert.True(whole.WholeSpace);
        Assert.Equal((double.NegativeInfinity, double.PositiveInfinity), (whole.Lower, whole.Upper));
        Assert.Contains("n ≥ 4", Assert.Throws<ArgumentException>(() => ResidualShrinkage.Intervals([0.1, -0.2, 0.3], [1.0], ResidualScore.Shrinkage, 0.5)).Message);
        Assert.False(ResidualShrinkage.Intervals([0.1, -0.2, 0.3], [1.0], ResidualScore.DistanceToAverage, 0.5).Intervals[0].WholeSpace);
        Assert.Throws<ArgumentOutOfRangeException>(() => ResidualShrinkage.Intervals([0.1, 0.2], [0.0], tolerance: 0.0));
        Assert.Throws<ArgumentException>(() => ResidualShrinkage.Intervals([0.1, double.NaN], [0.0]));
        Assert.Throws<ArgumentException>(() => ResidualShrinkage.Intervals([.. Enumerable.Range(0, 20).Select(i => 0.1 * i)], [double.NaN]));
        Assert.Throws<ArgumentException>(() => ResidualShrinkage.Intervals([0.1, 0.2], [double.NaN]));
        Assert.Throws<ArgumentOutOfRangeException>(() => ResidualShrinkage.Intervals([0.1, 0.2], [0.0], (ResidualScore)7));
    }
}
