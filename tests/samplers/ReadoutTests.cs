using System;
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
