using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Sensitivity;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The sampled-sensitivity brief's slice 0 gates on the pinned effect estimator (<see cref="PinnedEffects"/>), against the closed-form
/// ANOVA components of <see cref="SobolReferences"/>: every main effect and two-way effects on five-point grids, replicated over seeds,
/// each cell's mean within its spread of the truth and its stated error within 15% of its spread; the error falling as N^(−1/2); a
/// zero-variance function; shift invariance; and a seeded sweep over far ranges.
/// </summary>
public sealed class PinnedEffectsTests(ITestOutputHelper output)
{
    public static TheoryData<string> Names()
    {
        var names = new TheoryData<string>();
        foreach (SobolReference r in SobolReferences.All()) names.Add(r.Name);
        return names;
    }

    private static SobolReference Named(string name) => SobolReferences.All().Single(r => r.Name == name);

    private static double[] Grid(SobolReference f, int v, int points) =>
        [.. Enumerable.Range(0, points).Select(i => f.Law.Lower[v] + (f.Law.Upper[v] - f.Law.Lower[v]) * i / (points - 1))];

    // Every varying input's main effect, and the two-way effects of the first two pairs of varying inputs.
    private static List<PinnedEffectDesign> Designs(SobolReference f, int points)
    {
        int[] varying = [.. Enumerable.Range(0, f.Dimension).Where(f.Varies)];
        var designs = varying.Select(v => new PinnedEffectDesign(f.Dimension, [v], [Grid(f, v, points)])).ToList();
        foreach (int[] pair in varying.SelectMany((i, a) => varying.Skip(a + 1).Select(j => new[] { i, j })).Take(2))
            designs.Add(new PinnedEffectDesign(f.Dimension, pair, [Grid(f, pair[0], points), Grid(f, pair[1], points)]));
        return designs;
    }

    // Each cell's truth: the reference's ANOVA component at the cell, the last input fastest.
    private static double[] Truth(SobolReference f, PinnedEffectDesign design)
    {
        var truth = new double[design.Cells];
        int q = design.Subset.Length;
        for (int c = 0; c < design.Cells; c++)
        {
            var at = new double[q];
            int rest = c;
            for (int t = q - 1; t >= 0; t--)
            {
                at[t] = design.Grids[t][rest % design.Grids[t].Length];
                rest /= design.Grids[t].Length;
            }
            truth[c] = f.Effect(design.Subset, at);
        }
        return truth;
    }

    private static double Scale(SobolReference f) => Math.Abs(f.Mean) + Math.Sqrt(f.Variance);

    [Theory]
    [MemberData(nameof(Names))]
    public void Cells_AreUnbiasedWithinTheirStatedError(string name)
    {
        const int Replicates = 200, N = 1024;
        SobolReference f = Named(name);
        List<PinnedEffectDesign> designs = Designs(f, 5);
        var runs = new List<List<PickFreezeTests.Piece>>();
        foreach (int seed in SeedTree.Derive(20261115, Replicates))
        {
            (double[,] x, _) = PickFreezeTests.Bases(f, N, seed);
            var pieces = new List<PickFreezeTests.Piece>();
            foreach (PinnedEffectDesign design in designs)
            {
                PinnedEffectEstimate e = PinnedEffects.Estimate(design, f.Value, x);
                double[] truth = Truth(f, design);
                string label = string.Join(",", design.Subset.Select(i => i + 1));
                for (int c = 0; c < design.Cells; c++) pieces.Add(new($"f_{label}[{c}]", truth[c], e.Values[c], e.Errors[c]));
            }
            runs.Add(pieces);
        }
        PickFreezeTests.Summary s = PickFreezeTests.Summarize(runs, Scale(f));
        output.WriteLine($"{name}: {s.Tested} cells tested, {s.Exact} exact zeros; worst |bias| / its error {s.WorstBias:F2} ({s.WorstBiasPiece}); "
                         + $"stated error / spread in [{s.LowestCalibration:F3}, {s.HighestCalibration:F3}] ({s.LowestPiece}, {s.HighestPiece})");
        Assert.True(s.WorstBias <= 4.0, $"{name}: {s.WorstBiasPiece} is biased by {s.WorstBias:F2} of its error");
        Assert.True(s.LowestCalibration >= 0.85 && s.HighestCalibration <= 1.15,
                    $"{name}: stated error / spread runs from {s.LowestCalibration:F3} ({s.LowestPiece}) to {s.HighestCalibration:F3} ({s.HighestPiece})");
    }

    [Fact]
    public void Error_FallsAsTheInverseRootOfN()
    {
        const int Replicates = 100;
        int[] sizes = [256, 1024, 4096, 16384];
        SobolReference f = SobolReferences.Ishigami();
        PinnedEffectDesign[] designs =
        [
            new(3, [0], [[-2.0, 0.5, 3.0]]),
            new(3, [0, 2], [[-2.0, 1.0], [-3.0, 0.5, 2.5]]),
        ];
        var spread = new List<double[]>();
        var stated = new List<double[]>();
        foreach (int n in sizes)
        {
            var estimates = new List<double[]>();
            var errors = new List<double[]>();
            foreach (int seed in SeedTree.Derive(n, Replicates))
            {
                (double[,] x, _) = PickFreezeTests.Bases(f, n, seed);
                PinnedEffectEstimate[] e = [.. designs.Select(d => PinnedEffects.Estimate(d, f.Value, x))];
                estimates.Add([.. e.SelectMany(t => t.Values)]);
                errors.Add([.. e.SelectMany(t => t.Errors)]);
            }
            int cells = estimates[0].Length;
            spread.Add([.. Enumerable.Range(0, cells).Select(c => PickFreezeTests.Deviation([.. estimates.Select(t => t[c])]))]);
            stated.Add([.. Enumerable.Range(0, cells).Select(c => errors.Average(t => t[c]))]);
        }
        double[] logN = [.. sizes.Select(n => Math.Log(n))];
        for (int c = 0; c < spread[0].Length; c++)
        {
            double slope = Slope(logN, [.. spread.Select(s => Math.Log(s[c]))]), statedSlope = Slope(logN, [.. stated.Select(s => Math.Log(s[c]))]);
            output.WriteLine($"cell {c}: spread slope {slope:F3}, stated slope {statedSlope:F3}");
            Assert.True(Math.Abs(slope + 0.5) <= 0.08, $"cell {c}: the spread falls as N^{slope:F3}");
            Assert.True(Math.Abs(statedSlope + 0.5) <= 0.05, $"cell {c}: the stated error falls as N^{statedSlope:F3}");
        }
    }

    [Fact]
    public void ZeroVariance_GivesZeroEffects()
    {
        foreach (double c in new[] { 0.0, -2.5, 1e200 })
        {
            var f = new ConstantReference(c, 3);
            (double[,] x, _) = PickFreezeTests.Bases(f, 64, 5);
            foreach (PinnedEffectDesign design in new[] { new PinnedEffectDesign(3, [1], [[0.0, 0.5, 1.0]]), new PinnedEffectDesign(3, [0, 2], [[0.2, 0.9], [0.1, 0.4]]) })
            {
                PinnedEffectEstimate e = PinnedEffects.Estimate(design, f.Value, x);
                Assert.All(e.Values.Concat(e.Errors), v => Assert.Equal(0.0, v));
            }
        }
    }

    [Fact]
    public void Estimates_AreShiftInvariant()
    {
        SobolReference f = SobolReferences.G();
        var design = new PinnedEffectDesign(6, [0, 3], [[0.1, 0.6], [0.0, 0.3, 0.9]]);
        (double[,] x, _) = PickFreezeTests.Bases(f, 2048, 6);
        PinnedEffectEstimate e0 = PinnedEffects.Estimate(design, f.Value, x);
        foreach (double shift in new[] { 10.0, 1e4, 1e8 })
        {
            PinnedEffectEstimate e1 = PinnedEffects.Estimate(design, z => f.Value(z) + shift, x);
            double worst = Enumerable.Range(0, design.Cells).Max(c => Math.Max(Math.Abs(e1.Values[c] - e0.Values[c]), Math.Abs(e1.Errors[c] - e0.Errors[c])) / e0.Errors[c]);
            output.WriteLine($"shift {shift:E0}: worst change / error {worst:E2}");
            Assert.True(worst <= 1e-4, $"shift {shift:E0}: {worst:E2}");
        }
    }

    [Fact]
    public void SeededSweep_FarRanges()
    {
        // 200 cases from the estimator sweep's families, scales and offsets: one effect of one or two inputs at three random points each,
        // N = 2048 base rows; each cell's z = (estimate − truth)/stated error.
        var rng = new Xoshiro256PlusPlus(20261116);
        var z = new List<(double Z, string Where)>();
        int exact = 0;
        for (int c = 0; c < 200; c++)
        {
            SobolReference f = PickFreezeTests.SweepCase(rng);
            int p = f.Dimension, first = rng.NextInt(p);
            int[] subset = p >= 2 && rng.NextDouble() < 0.5 ? [.. new[] { first, (first + 1 + rng.NextInt(p - 1)) % p }.Order()] : [first];
            double[][] grids = [.. subset.Select(v => Enumerable.Range(0, 3).Select(_ => f.Law.Lower[v] + (f.Law.Upper[v] - f.Law.Lower[v]) * rng.NextDouble()).ToArray())];
            var design = new PinnedEffectDesign(p, subset, grids);
            (double[,] x, _) = PickFreezeTests.Bases(f, 2048, (int)rng.NextUInt64());
            PinnedEffectEstimate e = PinnedEffects.Estimate(design, f.Value, x);
            double[] truth = Truth(f, design);
            double scale = Math.Sqrt(f.Variance) + Math.Abs(f.Mean);
            for (int cell = 0; cell < design.Cells; cell++)
            {
                if (Math.Abs(truth[cell]) <= 1e-12 * scale && Math.Abs(e.Values[cell]) <= 1e-10 * scale)
                {
                    exact++;
                    continue;
                }
                Assert.True(e.Errors[cell] > 0.0, $"case {c} ({f.Name}) cell {cell}: estimate {e.Values[cell]:R}, truth {truth[cell]:R}, no error");
                z.Add(((e.Values[cell] - truth[cell]) / e.Errors[cell], $"case {c} ({f.Name}) [{string.Join(",", subset)}] cell {cell}"));
            }
        }
        double mean = z.Average(t => t.Z), sd = Math.Sqrt(z.Sum(t => (t.Z - mean) * (t.Z - mean)) / (z.Count - 1));
        double beyond3 = z.Count(t => Math.Abs(t.Z) > 3.0) / (double)z.Count;
        var worst = z.MaxBy(t => Math.Abs(t.Z));
        output.WriteLine($"{z.Count} z-scores ({exact} exact zeros): mean {mean:F3}, sd {sd:F3}, beyond 3: {beyond3:P2}, worst {worst.Z:F2} at {worst.Where}");
        Assert.True(Math.Abs(mean) <= 0.12, $"mean z {mean:F3}");
        Assert.True(sd >= 0.88 && sd <= 1.12, $"sd of z {sd:F3}");
        Assert.True(beyond3 <= 0.012, $"{beyond3:P2} beyond 3");
        Assert.True(Math.Abs(worst.Z) <= 6.0, $"worst z {worst.Z:F2} at {worst.Where}");
    }

    [Fact]
    public void Design_BlocksCellsAndPoints()
    {
        var design = new PinnedEffectDesign(3, [0, 2], [[1.0, 2.0], [7.0, 8.0, 9.0]]);
        Assert.Equal(6, design.Cells);
        Assert.Equal(1 + 2 + 3 + 6, design.Blocks.Count);
        Assert.Empty(design.Blocks[0].Inputs);
        Assert.Equal([0], design.Blocks[1].Inputs);
        Assert.Equal([2], design.Blocks[3].Inputs);
        Assert.Equal([0, 2], design.Blocks[6].Inputs);
        Assert.Equal([1.0, 8.0], design.Blocks[7].Values);   // the pair's cells run with the last input fastest
        double[,] x = { { 0.1, 0.2, 0.3 }, { 0.4, 0.5, 0.6 } };
        Assert.Equal(new double[,] { { 2.0, 0.2, 9.0 }, { 2.0, 0.5, 9.0 } }, design.Points(x, 11));

        // f = x₀·x₂ + x₁: per row the alternating sum is (g₀ − x₀)(g₂ − x₂), x₁ cancelling; the cell is their mean, its error their sd over √N.
        PinnedEffectEstimate e = PinnedEffects.Estimate(design, z => z[0] * z[2] + z[1], x);
        for (int c = 0; c < 6; c++)
        {
            double g0 = design.Grids[0][c / 3], g2 = design.Grids[1][c % 3];
            double[] rows = [.. Enumerable.Range(0, 2).Select(r => g0 * g2 - g0 * x[r, 2] - x[r, 0] * g2 + x[r, 0] * x[r, 2])];
            Assert.Equal(rows.Average(), e.Values[c], 12);
            Assert.Equal(Math.Abs(rows[0] - rows[1]) / 2.0, e.Errors[c], 12);   // sd over √2 of two values
        }

        Assert.Throws<ArgumentException>(() => new PinnedEffectDesign(3, [], []));
        Assert.Throws<ArgumentException>(() => new PinnedEffectDesign(3, [2, 0], [[1.0], [1.0]]));
        Assert.Throws<ArgumentException>(() => new PinnedEffectDesign(3, [0], [[]]));
        Assert.Throws<ArgumentException>(() => new PinnedEffectDesign(3, [0], [[double.NaN]]));
        Assert.Throws<ArgumentException>(() => new PinnedEffectDesign(3, [0, 1], [[1.0, 2.0], [1.0, 2.0]], budget: 5));
        Assert.Throws<ArgumentException>(() => PinnedEffects.Estimate(new PinnedEffectDesign(1, [0], [[0.5]]), [[1.0], [2.0]]));
    }

    [Fact]
    public void Estimate_FromValuesMeetsTheFunctionsForm()
    {
        SobolReference f = SobolReferences.Ishigami();
        var design = new PinnedEffectDesign(3, [0, 2], [[-1.0, 2.0], [0.5, 3.0]]);
        (double[,] x, _) = PickFreezeTests.Bases(f, 256, 2);
        double[] Values(double[,] points) => [.. Enumerable.Range(0, points.GetLength(0)).Select(j => f.Value(PickFreezeTests.Row(points, j)))];
        PinnedEffectEstimate byValues = PinnedEffects.Estimate(design, [.. Enumerable.Range(0, design.Blocks.Count).Select(k => Values(design.Points(x, k)))]);
        PinnedEffectEstimate byFunction = PinnedEffects.Estimate(design, f.Value, x);
        Assert.Equal(byFunction.Values, byValues.Values);
        Assert.Equal(byFunction.Errors, byValues.Errors);
    }

    private static double Slope(double[] x, double[] y)
    {
        double mx = x.Average(), my = y.Average();
        return x.Zip(y).Sum(t => (t.First - mx) * (t.Second - my)) / x.Sum(t => (t - mx) * (t - mx));
    }
}
