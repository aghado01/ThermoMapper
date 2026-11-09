using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class CrossConformalTests(ITestOutputHelper output)
{
    internal sealed record Unit(double X, double Y, int Priority = 0);

    internal static IEnumerable<int[]> Permutations(int n)
    {
        int[] order = [.. Enumerable.Range(0, n)];
        IEnumerable<int[]> Walk(int at)
        {
            if (at == n) { yield return (int[])order.Clone(); yield break; }
            for (int j = at; j < n; j++)
            {
                (order[at], order[j]) = (order[j], order[at]);
                foreach (int[] p in Walk(at + 1)) yield return p;
                (order[at], order[j]) = (order[j], order[at]);
            }
        }
        return Walk(0);
    }

    internal static LeaveOneOutReadings ResidualReadings(
        Unit[] train, double[] tests, Func<Unit[], Func<double, double>> fit)
    {
        var held = new double[train.Length];
        var predictions = new double[train.Length, tests.Length];
        for (int i = 0; i < train.Length; i++)
        {
            Func<double, double> f = fit([.. train.Where((_, j) => j != i)]);
            held[i] = f(train[i].X);
            for (int t = 0; t < tests.Length; t++) predictions[i, t] = f(tests[t]);
        }
        return new([.. train.Select(u => u.Y)], held, predictions);
    }

    private static ConformalSet PairAt(double[,] matrix, int t, double alpha)
    {
        int[] others = [.. Enumerable.Range(0, matrix.GetLength(0)).Where(i => i != t)];
        double[] own = [.. others.Select(i => matrix[i, t])];
        var candidate = new double[others.Length, 1];
        for (int i = 0; i < others.Length; i++) candidate[i, 0] = matrix[t, others[i]];
        return CrossConformal.Pairwise(own, candidate, [[0.0]], alpha).Sets[0];
    }

    [Fact]
    public void CardRefusedPairwiseLe_ExhaustiveTiesAndSeededFarScores()
    {
        // card_refused_pairwise_le and acceptsPairwise_zero. Each edge has three outcomes, including a tie.
        int outcomes = 0;
        for (int m = 1; m <= 4; m++)
        {
            var edges = (from i in Enumerable.Range(0, m) from j in Enumerable.Range(i + 1, m - i - 1) select (i, j)).ToArray();
            int count = (int)Math.Pow(3, edges.Length);
            for (int code = 0; code < count; code++)
            {
                var score = new double[m, m];
                int c = code;
                foreach (var (i, j) in edges)
                {
                    int outcome = c % 3;
                    c /= 3;
                    score[i, j] = outcome == 0 ? double.NegativeInfinity : outcome == 1 ? double.PositiveInfinity : -1e300;
                    score[j, i] = outcome == 0 ? double.PositiveInfinity : outcome == 1 ? double.NegativeInfinity : -1e300;
                }
                for (int k = 0; k < m; k++)
                {
                    double alpha = (k + 0.25) / m;
                    int refused = 0;
                    for (int t = 0; t < m; t++)
                    {
                        ConformalSet set = PairAt(score, t, alpha);
                        int strange = Enumerable.Range(0, m).Count(j => score[t, j] < score[j, t]);
                        Assert.Equal(strange + k < m, set.Accepts(0));
                        Assert.Equal(ConformalForm.CrossFitted, set.Form);
                        Assert.Equal(k == 0 ? 1 : Math.Max(0, 1 - (2.0 * k - 1) / m), set.Coverage);
                        Assert.Equal(set.Coverage, set.CoverageUnderError);
                        if (!set.Accepts(0)) refused++;
                    }
                    Assert.True(refused <= Math.Max(0, 2 * k - 1));
                }
                outcomes++;
            }
        }
        var rng = new Xoshiro256PlusPlus(20261110);
        for (int trial = 0; trial < 400; trial++)
        {
            int m = 2 + rng.NextInt(29), k = rng.NextInt(m);
            var score = new double[m, m];
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++) score[i, j] = rng.NextInt(5) - 2;
            int refused = Enumerable.Range(0, m).Count(t => !PairAt(score, t, (k + 0.25) / m).Accepts(0));
            Assert.True(refused <= Math.Max(0, 2 * k - 1));
        }
        output.WriteLine($"{outcomes} complete comparison patterns, ties and infinite scores included; 400 seeded matrices to m=30");
    }

    [Fact]
    public void JackRefuses_ClosedEndsAgreeWithStrictCounts()
    {
        var rng = new Xoshiro256PlusPlus(20261111);
        int checkedPoints = 0;
        for (int trial = 0; trial < 1000; trial++)
        {
            int n = rng.NextInt(35);
            double scale = Math.ScaleB(1.0, rng.NextInt(1001) - 500);
            double[] y = [.. Enumerable.Range(0, n).Select(_ => scale * (rng.NextInt(7) - 3))];
            double[] held = [.. Enumerable.Range(0, n).Select(_ => scale * (rng.NextInt(7) - 3))];
            var pred = new double[n, 2];
            for (int i = 0; i < n; i++) { pred[i, 0] = scale * (rng.NextInt(7) - 3); pred[i, 1] = -pred[i, 0]; }
            double alpha = 0.01 + 0.95 * rng.NextDouble();
            var result = CrossConformal.JackknifePlus(new(y, held, pred), alpha);
            for (int t = 0; t < 2; t++)
            {
                JackknifePlusInterval d = result.Intervals[t];
                double[] lower = [.. Enumerable.Range(0, n).Select(i => pred[i, t] - Math.Abs(y[i] - held[i]))];
                double[] upper = [.. Enumerable.Range(0, n).Select(i => pred[i, t] + Math.Abs(y[i] - held[i]))];
                double[] probes = [.. lower, .. upper, -50 * scale, 0, 50 * scale];
                foreach (double v in probes)
                {
                    bool refused = upper.Count(x => x < v) + d.Level >= n + 1 || lower.Count(x => x > v) + d.Level >= n + 1;
                    Assert.Equal(!refused, d.Contains(v));
                    checkedPoints++;
                }
                if (!d.Empty && !d.WholeSpace) { Assert.True(d.Contains(d.Lower)); Assert.True(d.Contains(d.Upper)); }
            }
        }
        output.WriteLine($"JackRefuses: {checkedPoints} comparisons, including exact endpoints and scales 1e-150 to 1e150");
    }

    [Fact]
    public void JackknifePlusOrderings_AndPairScoreEquivariance_WithTravelingRandomness()
    {
        // jackknifePlus_orderings, pairScore_equivariant, not_acceptsPairwise_of_jackRefuses.
        // Priority is the run's common random ordering attached to each unit and relabeled with it.
        Unit[] data = [new(0, 0, 3), new(1, 2, 0), new(2, 1, 4), new(3, -1, 1), new(4, 3, 2)];
        var fits = new Func<Unit[], Func<double, double>>[]
        {
            u => { double mean = u.OrderBy(v => v.X).Average(v => v.Y); return _ => mean; },
            u => { double median = u.Select(v => v.Y).Order().ElementAt(u.Length / 2); return _ => median; },
            u => { double y = u.MinBy(v => v.Priority)!.Y; return x => y + 0.5 * x; },
            u => { double mean = u.OrderBy(v => v.X).Average(v => v.Y) + (u.Length >= 4 ? 100 : 0); return _ => mean; },
        };
        int orderings = 0;
        foreach (var fit in fits)
            foreach (int k in new[] { 1, 2 })
            {
                int accepted = 0, pairAccepted = 0, total = 0;
                double[,] Pair(Unit[] z)
                {
                    var matrix = new double[z.Length, z.Length];
                    for (int i = 0; i < z.Length; i++)
                        for (int j = i + 1; j < z.Length; j++)
                        {
                            var f = fit([.. z.Where((_, a) => a != i && a != j)]);
                            matrix[i, j] = -Math.Abs(z[i].Y - f(z[i].X));
                            matrix[j, i] = -Math.Abs(z[j].Y - f(z[j].X));
                        }
                    return matrix;
                }
                double[,] reference = Pair(data);
                foreach (int[] p in Permutations(data.Length))
                {
                    Unit[] z = [.. p.Select(i => data[i])];
                    double[,] pair = Pair(z);
                    for (int i = 0; i < z.Length; i++)
                        for (int j = 0; j < z.Length; j++) Assert.Equal(reference[p[i], p[j]], pair[i, j]);
                    Unit test = z[^1];
                    var r = ResidualReadings(z[..^1], [test.X], fit);
                    var interval = CrossConformal.JackknifePlus(r, (k + 0.25) / z.Length).Intervals[0];
                    bool pairHit = PairAt(pair, z.Length - 1, (k + 0.25) / z.Length).Accepts(0);
                    Assert.True(!pairHit || interval.Contains(test.Y));
                    if (pairHit) pairAccepted++;
                    if (interval.Contains(test.Y)) accepted++;
                    total++;
                    Assert.Equal(Math.Max(0, 1 - (2.0 * k - 1) / z.Length), interval.Coverage);
                }
                Assert.True(accepted * data.Length >= total * (data.Length - (2 * k - 1)));
                Assert.True(pairAccepted * data.Length >= total * (data.Length - (2 * k - 1)));
                orderings += total;
            }
        // A fixed index order can violate the premise. The same unit's seed/priority has to travel with it.
        double IndexScore(Unit[] z, int i, int j) => -Math.Abs(z[i].Y - z.Where((_, a) => a != i && a != j).First().Y);
        Assert.NotEqual(IndexScore(data, 4, 3), IndexScore(data.Reverse().ToArray(), 0, 1));
        output.WriteLine($"{orderings} public jackknife+ intervals over every ordering, including an unstable fit and a shared random ordering");
    }

    [Theory]
    [InlineData(3, 1)]
    [InlineData(5, 2)]
    [InlineData(7, 3)]
    [InlineData(20, 2)]
    public void JackknifePlusCoverage_TheFiniteBoundIsAttained(int m, int k)
    {
        // Symmetric lookup fit identified by its two missing input labels. Its residual matrix is a regular tournament.
        int size = 2 * k - 1;
        var residual = new double[m, m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++) residual[i, j] = 1;
        for (int i = 0; i < size; i++)
            for (int j = 0; j < m; j++)
                if (j >= size || ((j - i + size) % size >= 1 && (j - i + size) % size <= (size - 1) / 2))
                { residual[i, j] = 2; residual[j, i] = 1; }
        int covered = 0;
        for (int t = 0; t < m; t++)
        {
            int[] others = [.. Enumerable.Range(0, m).Where(i => i != t)];
            double[] y = new double[m - 1], held = [.. others.Select(i => -residual[i, t])];
            var pred = new double[m - 1, 1];
            for (int i = 0; i < others.Length; i++) pred[i, 0] = -residual[t, others[i]];
            double alpha = m == 20 ? 0.1 : (k + 0.25) / m;
            var interval = CrossConformal.JackknifePlus(new(y, held, pred), alpha).Intervals[0];
            Assert.Equal(k, interval.Level);
            if (interval.Contains(0)) covered++;
        }
        Assert.Equal(m - size, covered);
        if (m == 5)
        {
            Unit[] data = [.. Enumerable.Range(0, m).Select(i => new Unit(i, 0))];
            Func<double, double> Lookup(Unit[] training)
            {
                int[] missing = [.. Enumerable.Range(0, m).Where(i => !training.Any(u => u.X == i))];
                int a = missing[0], b = missing[1];
                return x => x == a ? -residual[a, b] : x == b ? -residual[b, a] : 0;
            }
            int hits = 0, total = 0;
            foreach (int[] p in Permutations(m))
            {
                Unit[] z = [.. p.Select(i => data[i])];
                var readings = ResidualReadings(z[..^1], [z[^1].X], Lookup);
                if (CrossConformal.JackknifePlus(readings, (k + 0.25) / m).Intervals[0].Contains(0)) hits++;
                total++;
            }
            Assert.Equal(total * (m - size), hits * m);
        }
        output.WriteLine($"m={m}, k={k}: coverage {covered}/{m} attains 1-(2k-1)/m = {1 - (double)size / m}");
    }

    [Fact]
    public void Density_UsesEachFitsOwnDrawsAndStrictComparisons()
    {
        CandidateLogDensity Seam(double[,] values, double[][] grids) => new(values, grids, []);
        CandidateLogDensity[] held =
        [
            Seam(new double[,] { { -1000 }, { -1002 } }, [[1.0]]) with { Weights = [1, 0] },
            Seam(new double[,] { { -3 } }, [[2.0]]),
            Seam(new double[,] { { -1 }, { -2 }, { -3 } }, [[3.0]]) with { ChainLengths = [1, 2] },
        ];
        CandidateLogDensity[] tests =
        [
            Seam(new double[,] { { -1001, -1000, -999 }, { -1002, -1002, -1002 } }, [[-1.0, 0, 1]]) with { Weights = [1, 0] },
            Seam(new double[,] { { -4, -3, -2 } }, [[-1.0, 0, 1]]),
            Seam(new double[,] { { -2, -1, 0 }, { -3, -2, -1 }, { -4, -3, -2 } }, [[-1.0, 0, 1]]) with { ChainLengths = [1, 2] },
        ];
        var d = CrossConformal.Density(held, tests, 0.5).Sets[0];
        Assert.Equal([false, true, true], Enumerable.Range(0, 3).Select(d.Accepts));
        Assert.Equal(3, d.Tied![1]);
        Assert.Equal(0.25, d.Coverage);
        Assert.Null(d.CandidateScores);
        Assert.Throws<ArgumentException>(() => CrossConformal.Density(held, [tests[0] with { Forms = ConformalForms.Full }, tests[1], tests[2]]));
        Assert.Throws<ArgumentException>(() => CrossConformal.Density(held, [tests[0] with { Weights = [1, 1] }, tests[1], tests[2]]));
        Assert.Throws<ArgumentException>(() => CrossConformal.Pairwise([1], new double[,] { { double.NaN } }, [[0.0]]));
        Assert.Throws<ArgumentException>(() => CrossConformal.Pairwise([1], new double[1, 2], [[0.0, 0.0]]));
        Assert.Throws<ArgumentException>(() => CrossConformal.JackknifePlus(new([0], [double.NaN], new double[1, 1])));
    }
}
