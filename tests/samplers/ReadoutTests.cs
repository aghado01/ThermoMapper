using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Samplers.Tests;

/// <summary>
/// Readout statistics (step 3.4): additive and mergeable, so pooling does not depend on order, with an optional weight
/// per fold. Unit weights must reproduce the unweighted statistics exactly; merging in any order must equal folding
/// everything into one readout.
/// </summary>
public sealed class ReadoutTests
{
    [Fact]
    public void Moments_UnitWeights_MatchTheUnweightedUpdateExactly()
    {
        double[] x = [0.3, -1.7, 2.25, 0.1, 5.5, -0.02];
        RunningMoments unweighted = RunningMoments.Of(x), weighted = new();
        foreach (double v in x) weighted.Add(v, 1.0);
        Assert.Equal(unweighted.Mean, weighted.Mean);
        Assert.Equal(unweighted.SumSquaredDeviations, weighted.SumSquaredDeviations);
        Assert.Equal(x.Length, weighted.SumWeights);
    }

    [Fact]
    public void Moments_Weighted_MatchTheDirectFormula_AndMergeInAnyOrder()
    {
        double[] x = [1.0, 4.0, -2.0, 3.5, 0.25, 7.0];
        double[] w = [0.5, 2.0, 1.0, 0.25, 3.0, 1.5];
        double sw = 0.0, swx = 0.0;
        for (int i = 0; i < x.Length; i++) { sw += w[i]; swx += w[i] * x[i]; }
        double mean = swx / sw, ss = 0.0;
        for (int i = 0; i < x.Length; i++) ss += w[i] * (x[i] - mean) * (x[i] - mean);

        RunningMoments Fold(int from, int to)
        {
            var m = new RunningMoments();
            for (int i = from; i < to; i++) m.Add(x[i], w[i]);
            return m;
        }
        RunningMoments whole = Fold(0, 6), a = Fold(0, 2), b = Fold(2, 4), c = Fold(4, 6);
        Assert.Equal(mean, whole.Mean, 12);
        Assert.Equal(ss / sw, whole.PopulationVariance, 12);
        foreach (RunningMoments m in new[] { RunningMoments.Merge(RunningMoments.Merge(a, b), c),
                                             RunningMoments.Merge(a, RunningMoments.Merge(b, c)),
                                             RunningMoments.Merge(RunningMoments.Merge(c, a), b) })
        {
            Assert.Equal(whole.Mean, m.Mean, 12);
            Assert.Equal(whole.SumSquaredDeviations, m.SumSquaredDeviations, 11);
            Assert.Equal(6, m.Count);
        }
    }

    [Fact]
    public void WeightedSums_FoldDenseAndSparse_AndMergeInAnyOrder()
    {
        WeightedSums Make(int seed)
        {
            var s = new WeightedSums(4);
            s.Add([1.0, 2.0, 0.0, 0.5], weight: seed);
            s.AddSparse([3, 1, 3], [1.0, 1.0, 1.0], weight: 0.5);   // repeats accumulate
            s.AddSparse([], []);                                     // an empty draw still counts
            return s;
        }
        WeightedSums a = Make(1), b = Make(2), c = Make(3);
        var ab = Make(1); ab.MergeFrom(b); ab.MergeFrom(c);
        var cb = Make(3); cb.MergeFrom(b); cb.MergeFrom(a);
        Assert.Equal(ab.Sums.ToArray(), cb.Sums.ToArray());
        Assert.Equal(ab.TotalWeight, cb.TotalWeight);
        Assert.Equal(9, ab.Folds);

        Assert.Equal([1.0, 2.5, 0.0, 1.5], a.Sums.ToArray());   // 1·[1, 2, 0, 0.5] + 0.5·(1 at 1, 2 at 3)
        Assert.Equal(2.5, a.TotalWeight);
        Assert.Equal(1.0, a.Mean(1));
        Assert.True(double.IsNaN(new WeightedSums(2).Mean(0)));
        Assert.Throws<ArgumentException>(() => a.MergeFrom(new WeightedSums(3)));
    }

    [Fact]
    public void DrawSample_UnitWeights_GiveTheType7Quantile()
    {
        var s = new DrawSample();
        foreach (double v in new[] { 5.0, 1.0, 4.0, 2.0, 3.0 }) s.Add(v);
        Assert.Equal(2.0, s.Quantile(0.25));
        Assert.Equal(3.0, s.Quantile(0.5));
        Assert.Equal(4.6, s.Quantile(0.9), 12);
        Assert.True(double.IsNaN(new DrawSample().Quantile(0.5)));
    }

    [Fact]
    public void DrawSample_Weighted_InterpolatesTheWeightedEcdf_AndMergesWithUnitWeights()
    {
        var s = new DrawSample();
        s.Add(0.0, 1.0);
        s.Add(1.0, 3.0);   // midpoints of the weights at 0.5 and 2.5 out of 4
        Assert.Equal(0.75, s.Quantile(0.5), 12);
        Assert.Equal(0.0, s.Quantile(0.1));
        Assert.Equal(1.0, s.Quantile(0.9));

        var unit = new DrawSample();
        unit.Add(2.0);
        s.MergeFrom(unit);
        Assert.Equal(3, s.Count);
        Assert.Equal(5.0, s.TotalWeight);
    }

    [Fact]
    public void DrawSample_Weighted_PoolsTiesAndDropsZeroMass()
    {
        // The external review's cases. Ties: (0, 1), (0, 3), (1, 1) in either order is (0, 4), (1, 1), whose weight
        // midpoints are 2 and 4.5 of 5, so the 70th percentile is 0.6.
        foreach (double[] order in new[] { new[] { 1.0, 3.0 }, new[] { 3.0, 1.0 } })
        {
            var tied = new DrawSample();
            tied.Add(0.0, order[0]);
            tied.Add(0.0, order[1]);
            tied.Add(1.0, 1.0);
            Assert.Equal(0.6, tied.Quantile(0.7), 12);
        }
        // A zero weight carries no mass: every quantile lies within the positive mass's range.
        var zero = new DrawSample();
        zero.Add(10.0, 0.5);
        zero.Add(20.0, 0.5);
        zero.Add(-100.0, 0.0);
        Assert.Equal(10.0, zero.Quantile(0.05));
        Assert.Equal(15.0, zero.Quantile(0.5), 12);
        // No mass at all: no quantile.
        var none = new DrawSample();
        none.Add(123.0, 0.0);
        Assert.True(double.IsNaN(none.Quantile(0.5)));
    }

    [Fact]
    public void DrawSample_Weighted_IsInvariantToOrderMergeAndScale()
    {
        // A weighted multiset with ties and zero weights, as SMC's duplicated particles give. Its quantiles must not
        // depend on insertion order or on how it was split and merged (exactly), nor on a common rescaling of the weights.
        var rng = new Maths.Rng.Xoshiro256PlusPlus(seed: 17);
        var draws = new List<(double Value, double Weight)>();
        for (int i = 0; i < 60; i++)
        {
            double value = Math.Round(rng.NextDouble() * 8.0) / 2.0;   // many ties
            double weight = i % 7 == 0 ? 0.0 : 0.1 + rng.NextDouble();
            draws.Add((value, weight));
        }
        double[] ps = [0.0, 0.03, 0.25, 0.5, 0.61, 0.9, 1.0];
        DrawSample From(IEnumerable<(double Value, double Weight)> items, double scale = 1.0)
        {
            var s = new DrawSample();
            foreach ((double v, double w) in items) s.Add(v, w * scale);
            return s;
        }
        double[] reference = ps.Select(From(draws).Quantile).ToArray();

        for (int trial = 0; trial < 20; trial++)
        {
            var shuffled = draws.OrderBy(_ => rng.NextDouble()).ToList();
            Assert.Equal(reference, ps.Select(From(shuffled).Quantile).ToArray());

            int cut = 1 + (int)(rng.NextDouble() * (draws.Count - 2));
            DrawSample left = From(shuffled.Take(cut)), right = From(shuffled.Skip(cut));
            DrawSample leftThenRight = left.Copy(), rightThenLeft = right.Copy();
            leftThenRight.MergeFrom(right);
            rightThenLeft.MergeFrom(left);
            Assert.Equal(reference, ps.Select(leftThenRight.Quantile).ToArray());
            Assert.Equal(reference, ps.Select(rightThenLeft.Quantile).ToArray());
        }
        double[] scaled = ps.Select(From(draws, scale: 3.7).Quantile).ToArray();
        for (int i = 0; i < ps.Length; i++) Assert.Equal(reference[i], scaled[i], 12);
    }

    private enum Outcome { Kept, LostA, LostB }

    [Fact]
    public void Tally_CountsByCategory_AndMerges()
    {
        var a = new Tally<Outcome>();
        a.Add(Outcome.Kept); a.Add(Outcome.LostB); a.Add(Outcome.Kept);
        var b = new Tally<Outcome>();
        b.Add(Outcome.LostA);
        a.MergeFrom(b);
        Assert.Equal(2, a.Count(Outcome.Kept));
        Assert.Equal(1, a.Count(Outcome.LostA));
        Assert.Equal(4, a.Total);
        Assert.Equal([(Outcome.Kept, 2L), (Outcome.LostA, 1L), (Outcome.LostB, 1L)], a.Entries);
    }
}
