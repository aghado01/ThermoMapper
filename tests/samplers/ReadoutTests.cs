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
    public void DrawSample_Weighted_IsAkinshinsType7_AndMergesWithUnitWeights()
    {
        // (0, 1) and (1, 3): n* = 4²/10 = 1.6, and the window's CDF at the first draw's cumulative weight ¼ is
        // F = 1.6·¼ − 0.6p, the share of the window the value 0 takes.
        var s = new DrawSample();
        s.Add(0.0, 1.0);
        s.Add(1.0, 3.0);
        Assert.Equal(0.9, s.Quantile(0.5), 12);    // F = 0.1
        Assert.Equal(0.66, s.Quantile(0.1), 12);   // F = 0.34
        Assert.Equal(1.0, s.Quantile(0.9), 12);    // F = 0
        Assert.Equal(0.6, s.Quantile(0.0), 12);    // not anchored: the least draw holds less than 1/n* of the weight

        var unit = new DrawSample();
        unit.Add(2.0);
        s.MergeFrom(unit);
        Assert.Equal(3, s.Count);
        Assert.Equal(5.0, s.TotalWeight);
    }

    [Fact]
    public void DrawSample_Weighted_PoolsTiesAndDropsZeroMass()
    {
        // Ties: (0, 1), (0, 3), (1, 1) in either order is the quantile function of (0, 4), (1, 1), with n* = 5²/11 over
        // the draws. At p = 0.7, F at the cumulative weight 4/5 is (25/11)(4/5) − (14/11)(0.7) = 10.2/11, so the quantile is
        // the value 1's share, 0.8/11.
        foreach (double[] order in new[] { new[] { 1.0, 3.0 }, new[] { 3.0, 1.0 } })
        {
            var tied = new DrawSample();
            tied.Add(0.0, order[0]);
            tied.Add(0.0, order[1]);
            tied.Add(1.0, 1.0);
            Assert.Equal(0.8 / 11.0, tied.Quantile(0.7), 12);
        }
        // A zero weight carries no mass: every quantile lies within the positive mass's range, and the equal weights
        // left give the type-7 quantile of their draws.
        var zero = new DrawSample();
        zero.Add(10.0, 0.5);
        zero.Add(20.0, 0.5);
        zero.Add(-100.0, 0.0);
        Assert.Equal(10.5, zero.Quantile(0.05), 12);
        Assert.Equal(15.0, zero.Quantile(0.5), 12);
        // No mass at all: no quantile.
        var none = new DrawSample();
        none.Add(123.0, 0.0);
        Assert.True(double.IsNaN(none.Quantile(0.5)));
    }

    [Fact]
    public void DrawSample_EqualWeights_AtAnyScale_GiveTheType7Quantile()
    {
        // Equal weights are unit weights up to a common scale: normalized MCMC weights (1/n each) must not move the
        // quantile off type 7, with or without zero-weight draws beside them.
        double[] values = [5.0, 1.0, 4.0, 2.0, 3.0, 4.0];
        double[] ps = [0.0, 0.1, 0.25, 0.5, 0.9, 1.0];
        var unit = new DrawSample();
        foreach (double v in values) unit.Add(v);
        foreach (double scale in new[] { 1.0 / values.Length, 0.2, 3.7 })
        {
            var scaled = new DrawSample();
            foreach (double v in values) scaled.Add(v, scale);
            scaled.Add(-50.0, 0.0);
            foreach (double p in ps) Assert.Equal(unit.Quantile(p), scaled.Quantile(p));
        }
        Assert.Equal(4.0, unit.Quantile(0.6), 12);   // sorted 1, 2, 3, 4, 4, 5 at position 0.6 · 5 = 3
    }

    [Fact]
    public void DrawSample_Weighted_IsContinuousAsOneWeightLeavesEqual()
    {
        // One draw's weight moved from the common 1/4 by a relative ε, for the least, an inner and the largest draw: every
        // quantile stays within O(ε) of type 7 (n·ε·range bounds the window's shift), so the rule has no jump between equal
        // and near-equal weights.
        double[] values = [5.0, 1.0, 4.0, 2.0, 3.0, 4.0, 7.5];
        double[] ps = [0.0, 0.1, 0.25, 0.5, 0.6, 0.9, 1.0];
        var equal = new DrawSample();
        foreach (double v in values) equal.Add(v, 0.25);
        foreach (int moved in new[] { 1, 0, 6 })
            foreach (double epsilon in new[] { 1e-2, -1e-2, 1e-5, 1e-8, -1e-11 })
            {
                var s = new DrawSample();
                for (int i = 0; i < values.Length; i++) s.Add(values[i], i == moved ? 0.25 * (1.0 + epsilon) : 0.25);
                foreach (double p in ps)
                    Assert.True(Math.Abs(s.Quantile(p) - equal.Quantile(p)) <= values.Length * Math.Abs(epsilon) * 6.5 + 1e-12,
                                $"draw {moved}, ε {epsilon}, p {p}: {s.Quantile(p)} against {equal.Quantile(p)}");
            }
    }

    [Fact]
    public void DrawSample_Weighted_StaysInRange_AndANegligibleDrawMovesNoQuantile()
    {
        // The midpoint rule's failures (stage 8's brief): with (0, ½) and (1, ½), a 100 at weight 10⁻⁹ read Q(0.99) = 97.03,
        // and with (0, ½) and (10, ½), a 1 at 10⁻⁹ read Q(0.4) = 0.80. Without the negligible draw both are type 7: 0.99 and 4.
        DrawSample With(double extra, double weight, double low, double high)
        {
            var s = new DrawSample();
            s.Add(low, 0.5);
            s.Add(high, 0.5);
            if (weight > 0.0) s.Add(extra, weight);
            return s;
        }
        Assert.Equal(0.99, With(100.0, 1e-9, 0.0, 1.0).Quantile(0.99), 6);
        Assert.Equal(4.0, With(1.0, 1e-9, 0.0, 10.0).Quantile(0.4), 6);

        // As the draw's weight goes to 0, every quantile approaches the one without it, within O(w).
        double[] ps = [0.0, 0.01, 0.4, 0.5, 0.99, 1.0];
        foreach ((double extra, double low, double high) in new[] { (100.0, 0.0, 1.0), (1.0, 0.0, 10.0), (-30.0, 2.0, 3.0) })
        {
            DrawSample without = With(extra, 0.0, low, high);
            foreach (double w in new[] { 1e-3, 1e-6, 1e-9, 1e-12, 1e-300 })
                foreach (double p in ps)
                {
                    double moved = Math.Abs(With(extra, w, low, high).Quantile(p) - without.Quantile(p));
                    Assert.True(moved <= 10.0 * w * (Math.Abs(extra) + Math.Abs(high)) + 1e-12, $"w {w}, p {p}: moved {moved}");
                }
        }

        // Random unequal weights, some zero: every quantile lies within the weighted draws' range and rises with p.
        var rng = new Maths.Rng.Xoshiro256PlusPlus(seed: 23);
        for (int trial = 0; trial < 50; trial++)
        {
            var s = new DrawSample();
            double min = double.PositiveInfinity, max = double.NegativeInfinity;
            int n = 1 + (int)(rng.NextDouble() * 12);
            for (int i = 0; i < n; i++)
            {
                double value = rng.NextDouble() * 20.0 - 10.0, weight = rng.NextDouble() < 0.2 ? 0.0 : Math.Pow(rng.NextDouble(), 4);
                s.Add(value, weight);
                if (weight > 0.0) (min, max) = (Math.Min(min, value), Math.Max(max, value));
            }
            double previous = double.NegativeInfinity;
            for (int k = 0; k <= 40; k++)
            {
                double q = s.Quantile(k / 40.0);
                if (double.IsInfinity(min)) Assert.True(double.IsNaN(q));
                else
                {
                    Assert.InRange(q, min, max);
                    Assert.True(q >= previous - 1e-12, $"Q falls at p = {k / 40.0}");
                    previous = q;
                }
            }
        }
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
