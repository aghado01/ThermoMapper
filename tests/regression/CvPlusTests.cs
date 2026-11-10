using System;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;
using Unit = Maths.Regression.Tests.CrossConformalTests.Unit;

namespace Maths.Regression.Tests;

public sealed class CvPlusTests(ITestOutputHelper output)
{
    internal static FoldReadings Readings(Unit[] train, int[] foldOf, int folds, double[] tests,
                                          Func<Unit[], Func<double, double>> fit)
    {
        var held = new double[train.Length];
        var predictions = new double[folds, tests.Length];
        for (int fold = 0; fold < folds; fold++)
        {
            Func<double, double> f = fit([.. train.Where((_, i) => foldOf[i] != fold)]);
            for (int i = 0; i < train.Length; i++)
                if (foldOf[i] == fold) held[i] = f(train[i].X);
            for (int t = 0; t < tests.Length; t++) predictions[fold, t] = f(tests[t]);
        }
        return new(train.Select(u => u.Y).ToArray(), held, foldOf, predictions);
    }

    [Fact]
    public void Bounds_AreSelectedOnlyUnderTheirOwnAssumptions()
    {
        const int n = 1000, folds = 5;
        var r = new FoldReadings(new double[n], new double[n],
                                 Enumerable.Range(0, n).Select(i => i % folds).ToArray(), new double[folds, 1]);
        foreach (CvPlusAssumptions assumptions in Enumerable.Range(0, 16).Select(i => (CvPlusAssumptions)i))
        {
            var result = CrossConformal.CvPlus(r, assumptions);
            var b = result.Bounds;
            Assert.Equal(100, b.Level);
            Assert.Equal(200, b.FoldSize);
            Assert.Equal(1 - 208.0 / 1005, b.Merged.Value);
            Assert.Equal(1 - 100.0 / 201, b.OneFold.Value);
            Assert.Equal(1 - 398.0 / 1200, b.Tournament.Value);
            Assert.Equal((assumptions & CvPlusAssumptions.MergedRanks) == CvPlusAssumptions.MergedRanks, b.Merged.Applies);
            Assert.Equal(b.Merged.Applies, b.OneFold.Applies);
            Assert.Equal((assumptions & CvPlusAssumptions.FoldTournament) == CvPlusAssumptions.FoldTournament, b.Tournament.Applies);
            Assert.Equal("cvPlus_merged", b.Merged.Theorem);
            Assert.Equal("cvPlus_oneFold", b.OneFold.Theorem);
            Assert.Equal("cvPlus_tournament", b.Tournament.Theorem);
            double? expected = null;
            if (b.Merged.Applies) expected = b.Merged.Value;
            if (b.OneFold.Applies) expected = Math.Max(expected ?? 0, b.OneFold.Value);
            if (b.Tournament.Applies) expected = Math.Max(expected ?? 0, b.Tournament.Value);
            Assert.Equal(expected, b.Coverage);
            Assert.Equal(expected, result.Intervals[0].Coverage);
        }
        var whole = CrossConformal.CvPlus(r, CvPlusAssumptions.None, double.Epsilon);
        Assert.Equal(1, whole.Bounds.Coverage);
        Assert.True(whole.Intervals[0].WholeSpace);
        Assert.False(whole.Bounds.Merged.Applies || whole.Bounds.OneFold.Applies || whole.Bounds.Tournament.Applies);
        Assert.Throws<ArgumentOutOfRangeException>(() => CrossConformal.CvPlus(r, (CvPlusAssumptions)16));
        var smallerFolds = new FoldReadings(new double[40], new double[40],
                                            Enumerable.Range(0, 40).Select(i => i % 10).ToArray(), new double[10, 1]);
        var mergedOnly = CrossConformal.CvPlus(smallerFolds, CvPlusAssumptions.MergedRanks).Bounds;
        Assert.True(mergedOnly.Tournament.Value > mergedOnly.Merged.Value);
        Assert.Equal(mergedOnly.Merged.Value, mergedOnly.Coverage);
        Assert.Equal(0.48, mergedOnly.Coverage!.Value, 15);
        var twoFolds = CrossConformal.CvPlus(r with
        {
            FoldOf = Enumerable.Range(0, n).Select(i => i % 2).ToArray(), Predictions = new double[2, 1],
        }, CvPlusAssumptions.Both).Bounds;
        Assert.True(twoFolds.OneFold.Value > Math.Max(twoFolds.Merged.Value, twoFolds.Tournament.Value));
        Assert.Equal(twoFolds.OneFold.Value, twoFolds.Coverage);
    }

    [Fact]
    public void CVRefuses_EqualFoldsClosedEndsStrictCountsAndFarScales()
    {
        var rng = new Xoshiro256PlusPlus(20261117);
        int probes = 0, empty = 0;
        for (int trial = 0; trial < 1000; trial++)
        {
            int folds = 1 + rng.NextInt(10), size = 1 + rng.NextInt(10), n = folds * size;
            double scale = Math.ScaleB(1, rng.NextInt(1001) - 500);
            double[] Values() => [.. Enumerable.Range(0, n).Select(_ => scale * (rng.NextInt(13) - 6))];
            double[] y = Values(), held = Values();
            int[] labels = [.. Enumerable.Range(0, n).Select(i => i % folds)];
            var predictions = new double[folds, 2];
            for (int f = 0; f < folds; f++) { predictions[f, 0] = scale * (rng.NextInt(13) - 6); predictions[f, 1] = -predictions[f, 0]; }
            double alpha = 0.01 + 0.98 * rng.NextDouble();
            var result = CrossConformal.CvPlus(new(y, held, labels, predictions), CvPlusAssumptions.Both, alpha);
            int level = (int)Math.Floor(alpha * (n + 1));
            foreach (var d in result.Intervals)
            {
                double[] lower = [.. Enumerable.Range(0, n).Select(i => predictions[labels[i], d.Case] - Math.Abs(y[i] - held[i]))];
                double[] upper = [.. Enumerable.Range(0, n).Select(i => predictions[labels[i], d.Case] + Math.Abs(y[i] - held[i]))];
                foreach (double value in lower.Concat(upper).Append(-100 * scale).Append(0).Append(100 * scale))
                {
                    bool refused = upper.Count(v => v < value) + level >= n + 1 || lower.Count(v => v > value) + level >= n + 1;
                    Assert.Equal(!refused, d.Contains(value));
                    probes++;
                }
                if (d.Empty) { empty++; Assert.Equal(0, d.Width); }
                if (level == 0) { Assert.Equal(double.NegativeInfinity, d.Lower); Assert.Equal(double.PositiveInfinity, d.Upper); }
            }
        }
        Assert.True(empty > 0);
        output.WriteLine($"CVRefuses: {probes} comparisons including ties, endpoints, empty sets and dyadic scales near 1e-150 to 1e150");
    }

    [Fact]
    public void SingletonFolds_ReduceToJackknifePlusUnderTheTournamentAssumptions()
    {
        var rng = new Xoshiro256PlusPlus(20261118);
        for (int n = 1; n <= 50; n++)
        {
            double[] y = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            double[] held = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
            int[] labels = [.. Enumerable.Range(0, n).Reverse()];
            var predictions = new double[n, 2];
            var rows = new double[n, 2];
            for (int i = 0; i < n; i++)
                for (int t = 0; t < 2; t++) rows[i, t] = predictions[labels[i], t] = rng.NextDouble();
            foreach (double alpha in new[] { 0.01, 0.1, 0.3, 0.8 })
            {
                var cv = CrossConformal.CvPlus(new(y, held, labels, predictions), CvPlusAssumptions.FoldTournament, alpha);
                var jack = CrossConformal.JackknifePlus(new(y, held, rows), alpha);
                for (int t = 0; t < 2; t++)
                {
                    Assert.Equal(jack.Intervals[t].Lower, cv.Intervals[t].Lower);
                    Assert.Equal(jack.Intervals[t].Upper, cv.Intervals[t].Upper);
                    Assert.Equal(jack.Intervals[t].Coverage, cv.Intervals[t].Coverage);
                }
                Assert.False(cv.Bounds.Merged.Applies);
            }
        }
        output.WriteLine("K=n, including permuted singleton fold labels: CV+ endpoints and the applicable tournament bound match jackknife+");
    }

    private static Func<double, double> Mean(Unit[] units)
    {
        double mean = units.Average(u => u.Y);
        return _ => mean;
    }

    [Theory]
    [InlineData(2, 2)]
    [InlineData(3, 2)]
    public void CvPlusMerged_AllOrderingsIncludingOrderReadingFits(int folds, int size)
    {
        // cvPlus_merged and sum_rankLe_lt_of_cvRefuses. The fits read only units outside their fold.
        int n = folds * size;
        double[] values = [0, 2, 1, 3, 0, 1, 7];
        Unit[] data = [.. Enumerable.Range(0, n + 1).Select(i => new Unit(i, values[i]))];
        var fits = new Func<Unit[], Func<double, double>>[]
        {
            Mean,
            u => { double first = u[0].Y; return _ => first; },
            u => { double weighted = u.Select((v, i) => (i + 1) * v.Y).Sum() / 2; return _ => weighted; },
        };
        int[] labels = [.. Enumerable.Range(0, n).Select(i => i / size)];
        int orderings = 0, refusedCases = 0;
        foreach (var fit in fits)
        {
            var covered = new int[n + 1];
            int total = 0;
            foreach (int[] p in CrossConformalTests.Permutations(data.Length))
            {
                Unit[] z = [.. p.Select(i => data[i])];
                Unit test = z[^1];
                FoldReadings readings = Readings(z[..n], labels, folds, [test.X], fit);
                int rankSum = 0;
                int largestRank = 0;
                for (int f = 0; f < folds; f++)
                {
                    double candidate = -Math.Abs(test.Y - readings.Predictions[f, 0]);
                    int rank = 1 + Enumerable.Range(0, n).Count(i => labels[i] == f &&
                        -Math.Abs(z[i].Y - readings.HeldOut[i]) <= candidate);
                    rankSum += rank;
                    largestRank = Math.Max(largestRank, rank);
                }
                for (int level = 0; level <= n; level++)
                {
                    var result = CrossConformal.CvPlus(readings, CvPlusAssumptions.MergedRanks, (level + 0.25) / (n + 1));
                    Assert.Equal(level, result.Bounds.Level);
                    Assert.False(result.Bounds.Tournament.Applies);
                    if (result.Intervals[0].Contains(test.Y)) covered[level]++;
                    else { Assert.True(rankSum < level + folds); Assert.True(largestRank <= level); refusedCases++; }
                }
                total++;
            }
            for (int level = 0; level <= n; level++)
            {
                double bound = Math.Max(Math.Max(0, 1 - 2.0 * (level + folds - 1) / (n + folds)),
                                        1 - level / (size + 1.0));
                Assert.True(covered[level] / (double)total >= bound - 1e-15);
            }
            orderings += total;
        }
        output.WriteLine($"cvPlus_merged and cvPlus_oneFold: {orderings} full orderings with symmetric, first-unit and position-weighted fits; {refusedCases} refusals obey both rank implications");
    }

    [Fact]
    public void OneFold_EmptyTrainingFitHasTheSplitRankBound()
    {
        double[] bag = [0, 1, 1, 3, 8];
        var covered = new int[bag.Length];
        int total = 0;
        foreach (int[] p in CrossConformalTests.Permutations(bag.Length))
        {
            double[] y = [.. p[..^1].Select(i => bag[i])];
            for (int level = 0; level < bag.Length; level++)
            {
                double alpha = (level + 0.25) / bag.Length;
                var result = CrossConformal.CvPlus(new(y, new double[y.Length], new int[y.Length], new double[1, 1]),
                                                    CvPlusAssumptions.MergedRanks, alpha);
                var split = SplitConformal.Residuals(new(new double[y.Length]), y, new(new double[1]), alpha);
                Assert.Equal(split.Intervals[0].Lower, result.Intervals[0].Lower);
                Assert.Equal(split.Intervals[0].Upper, result.Intervals[0].Upper);
                Assert.Equal(split.Intervals[0].Coverage, result.Bounds.Coverage);
                if (result.Intervals[0].Contains(bag[p[^1]])) covered[level]++;
            }
            total++;
        }
        for (int level = 0; level < bag.Length; level++)
            Assert.True(covered[level] / (double)total >= 1 - level / (double)bag.Length);
        output.WriteLine($"cvPlus_oneFold: {total} orderings with ties, empty-training predictions, and split endpoint/bound parity");
    }

    [Fact]
    public void CvPlusTournament_ExtendedFoldOrderingsAndSymmetricFits()
    {
        // cvPlus_tournament over n+m units. Only the n training units and the selected test unit enter the public readings.
        Unit[] data = [new(0, 0, 3), new(1, 2, 0), new(2, 1, 5), new(3, 3, 1), new(4, 8, 4), new(5, 1, 2)];
        const int n = 4, size = 2, folds = 2;
        int[] labels = [0, 0, 1, 1];
        var fits = new Func<Unit[], Func<double, double>>[]
        {
            Mean,
            u => { double median = u.Select(v => v.Y).Order().ElementAt(u.Length / 2); return _ => median; },
            u => x => { double distance = u.Min(v => Math.Abs(v.X - x)); return u.Where(v => Math.Abs(v.X - x) == distance).Average(v => v.Y); },
            u => { double firstBySharedPriority = u.MinBy(v => v.Priority)!.Y; return _ => firstBySharedPriority; },
        };
        int orderings = 0;
        foreach (var fit in fits)
        {
            var covered = new int[n + 1];
            int total = 0;
            foreach (int[] p in CrossConformalTests.Permutations(data.Length))
            {
                Unit[] z = [.. p.Select(i => data[i])];
                FoldReadings readings = Readings(z[..n], labels, folds, [z[n].X], fit);
                for (int level = 0; level <= n; level++)
                {
                    var result = CrossConformal.CvPlus(readings, CvPlusAssumptions.FoldTournament, (level + 0.25) / (n + 1));
                    Assert.False(result.Bounds.Merged.Applies);
                    if (result.Intervals[0].Contains(z[n].Y)) covered[level]++;
                }
                total++;
            }
            for (int level = 0; level <= n; level++)
                Assert.True(covered[level] / (double)total >= Math.Max(0, 1 - Math.Max(0, 2.0 * level + size - 2) / (n + size)) - 1e-15);
            orderings += total;
        }
        output.WriteLine($"cvPlus_tournament: {orderings} orderings of n+m units, including the test-fold filler omitted from all inputs");
    }

    [Fact]
    public void FoldValidation_RejectsEmptyUnequalMissingAndMismatchedRows()
    {
        FoldReadings R(int[] labels, int folds) => new(new double[labels.Length], new double[labels.Length], labels, new double[folds, 1]);
        foreach (FoldReadings bad in new[] { R([], 0), R([0, 1, 1], 2), R([0, 0, 0, 1], 2), R([0, 0], 2), R([-1, 1], 2), R([0, 2], 2) })
            Assert.Throws<ArgumentException>(() => CrossConformal.CvPlus(bad, CvPlusAssumptions.Both));
        Assert.Throws<ArgumentException>(() => CrossConformal.CvPlus(R([0, 1], 2) with { HeldOut = [0] }, CvPlusAssumptions.Both));
        Assert.Throws<ArgumentException>(() => CrossConformal.CvPlus(R([0, 1], 2) with { Response = [double.NaN, 0] }, CvPlusAssumptions.Both));
        Assert.Throws<ArgumentException>(() => CrossConformal.CvPlus(R([0, 1], 2) with { Predictions = new double[,] { { double.PositiveInfinity }, { 0 } } }, CvPlusAssumptions.Both));
        // A single fold is mathematically allowed when its empty-training fit supplies the readings.
        Assert.Single(CrossConformal.CvPlus(R([0, 0], 1), CvPlusAssumptions.MergedRanks).Intervals);
    }
}
