using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The cut SMC's ordering of the emulator draws (stage 8, C5; D09): the distances against the spread of the
/// log-likelihood differences computed directly, unmoved by a draw's constant; the tour visiting every draw once from
/// the first leg's; 2-opt never lengthening the path, against brute force on small instances; the tour and the random
/// order deterministic in their inputs and seed; and, with weighted draws, the tour ordering the resampled set, each
/// duplicate beside its twin.
/// </summary>
public sealed class BassCalibrationOrderingTests
{
    private readonly ITestOutputHelper _out;
    public BassCalibrationOrderingTests(ITestOutputHelper output) => _out = output;

    // Euclidean distances between `count` random points in the plane.
    private static double[,] Plane(int count, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[][] points = [.. Enumerable.Range(0, count).Select(_ => new[] { rng.NextDouble(), rng.NextDouble() })];
        var d = new double[count, count];
        for (int a = 0; a < count; a++)
            for (int b = 0; b < count; b++) d[a, b] = Math.Sqrt(Math.Pow(points[a][0] - points[b][0], 2) + Math.Pow(points[a][1] - points[b][1], 2));
        return d;
    }

    [Fact]
    public void Tour_VisitsEveryDrawOnce_FromTheFirstLegs()
    {
        foreach (int count in new[] { 1, 2, 3, 5, 30, 120 })
        {
            int start = count / 2;
            int[] tour = CutSmcOrdering.Tour(Plane(count, count), start);
            Assert.Equal(start, tour[0]);
            Assert.Equal(Enumerable.Range(0, count), tour.Order());
            int[] random = CutSmcOrdering.Random(count, start, new Xoshiro256PlusPlus(count));
            Assert.Equal(start, random[0]);
            Assert.Equal(Enumerable.Range(0, count), random.Order());
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => CutSmcOrdering.Tour(Plane(4, 1), 4));
        Assert.Throws<ArgumentException>(() => CutSmcOrdering.Tour(new double[0, 0], 0));
    }

    [Fact]
    public void TwoOpt_NeverLengthensThePath_AndNearsTheOptimum()
    {
        var ratios = new List<double>();
        for (int seed = 1; seed <= 200; seed++)
        {
            int count = 3 + seed % 40;
            double[,] d = Plane(count, 1000 + seed);
            int[] nearest = CutSmcOrdering.NearestNeighbour(d, 0);
            double before = CutSmcOrdering.Length(d, nearest), after = CutSmcOrdering.Length(d, CutSmcOrdering.TwoOpt(d, nearest));
            Assert.True(after <= before + 1e-12, $"2-opt lengthened {before} to {after}");
            int[] random = CutSmcOrdering.Random(count, 0, new Xoshiro256PlusPlus(seed));
            Assert.True(CutSmcOrdering.Length(d, CutSmcOrdering.TwoOpt(d, random)) <= CutSmcOrdering.Length(d, random) + 1e-12);
        }
        // Against brute force over every open path from node 0, at seven nodes.
        for (int seed = 1; seed <= 30; seed++)
        {
            double[,] d = Plane(7, 5000 + seed);
            double best = Permutations([1, 2, 3, 4, 5, 6]).Min(rest => CutSmcOrdering.Length(d, [0, .. rest]));
            double tour = CutSmcOrdering.Length(d, CutSmcOrdering.Tour(d, 0));
            Assert.True(tour >= best - 1e-12);
            ratios.Add(tour / best);
        }
        _out.WriteLine($"tour against the optimum over 30 seven-node instances: mean ratio {ratios.Average():F4}, worst {ratios.Max():F4}, optimal in {ratios.Count(r => r < 1 + 1e-12)}");
        // On a line from an end, the tour is the sorted order.
        var line = new double[6, 6];
        for (int a = 0; a < 6; a++)
            for (int b = 0; b < 6; b++) line[a, b] = Math.Abs(Math.Pow(a, 1.5) - Math.Pow(b, 1.5));
        Assert.Equal([0, 1, 2, 3, 4, 5], CutSmcOrdering.Tour(line, 0));
    }

    private static IEnumerable<int[]> Permutations(int[] items)
    {
        if (items.Length <= 1) { yield return items; yield break; }
        for (int i = 0; i < items.Length; i++)
            foreach (int[] rest in Permutations([.. items.Where((_, k) => k != i)]))
                yield return [items[i], .. rest];
    }

    [Fact]
    public void Distances_AreTheWeightedSpreadOfTheLogLikelihoodDifferences()
    {
        var rng = new Xoshiro256PlusPlus(9);
        const int s = 5, n = 300;
        var l = new double[s, n];
        double[] raw = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        double[] w = [.. raw.Select(x => x / raw.Sum())];
        w[7] = 0;   // a particle of no weight does not count
        double total = w.Sum();
        w = [.. w.Select(x => x / total)];
        for (int a = 0; a < s; a++)
            for (int i = 0; i < n; i++) l[a, i] = -1e6 * (a + 1) + 3 * Math.Sin(i * 0.1 + a) + 0.5 * a * rng.NextDouble();
        double[,] d = CutSmcOrdering.Distances(l, w);
        for (int a = 0; a < s; a++)
        {
            Assert.Equal(0.0, d[a, a]);
            for (int b = 0; b < s; b++)
            {
                Assert.Equal(d[a, b], d[b, a]);
                // Directly, on the differences with their large constants removed exactly (the rows' −10⁶·(a + 1)).
                double[] diff = [.. Enumerable.Range(0, n).Select(i => (l[a, i] + 1e6 * (a + 1)) - (l[b, i] + 1e6 * (b + 1)))];
                double mean = diff.Zip(w, (x, v) => x * v).Sum();
                double expected = Math.Sqrt(diff.Zip(w, (x, v) => v * (x - mean) * (x - mean)).Sum());
                Assert.Equal(expected, d[a, b], 1e-9 * (1 + expected));
            }
        }
        // A draw whose likelihood moves by a constant, a normalizer, stays where it was.
        var shifted = (double[,])l.Clone();
        for (int i = 0; i < n; i++) shifted[2, i] += 1234.5;
        double[,] ds = CutSmcOrdering.Distances(shifted, w);
        for (int a = 0; a < s; a++) Assert.Equal(d[2, a], ds[2, a], 1e-9);
        Assert.Throws<ArgumentException>(() => CutSmcOrdering.Distances(l, w.AsSpan(1)));
    }

    [Fact]
    public void Tour_AndRandomOrder_AreDeterministicInTheirInputsAndSeed()
    {
        double[,] d = Plane(40, 77);
        Assert.Equal(CutSmcOrdering.Tour(d, 3), CutSmcOrdering.Tour((double[,])d.Clone(), 3));
        Assert.Equal(CutSmcOrdering.Random(40, 3, new Xoshiro256PlusPlus(5)), CutSmcOrdering.Random(40, 3, new Xoshiro256PlusPlus(5)));
        Assert.NotEqual(CutSmcOrdering.Random(40, 3, new Xoshiro256PlusPlus(5)), CutSmcOrdering.Random(40, 3, new Xoshiro256PlusPlus(6)));
        // The random order is uniform over the other nodes: node 1's position among 1 … 4.
        var rng = new Xoshiro256PlusPlus(11);
        var counts = new int[5];
        for (int i = 0; i < 20000; i++) counts[Array.IndexOf(CutSmcOrdering.Random(5, 0, rng), 1)]++;
        Assert.Equal(0, counts[0]);
        for (int k = 1; k < 5; k++) Assert.InRange(counts[k] / 20000.0, 0.25 - 0.015, 0.25 + 0.015);
    }

    /// <summary>An emulator whose draws carry weights, over a <see cref="TestEmulator"/>.</summary>
    private sealed class WeightedEmulator(TestEmulator inner, double[] weights) : ICalibrationEmulator
    {
        public IReadOnlyList<CalibrationParameter> Parameters => inner.Parameters;
        public int DrawCount => inner.DrawCount;
        public IReadOnlyList<double>? DrawWeights => weights;
        public int Outputs => inner.Outputs;
        public EmulatorNoise Noise(int draw) => inner.Noise(draw);
        public double[,]? TruncationResiduals() => inner.TruncationResiduals();
        public void Predict(int draw, ReadOnlySpan<double> continuous, ReadOnlySpan<int> levels, Span<double> output) => inner.Predict(draw, continuous, levels, output);
    }

    [Fact]
    public void WithWeightedDraws_TheTourOrdersTheResampledSet()
    {
        TestEmulator inner = CalibrationFixture.Emulator();
        var emulator = new WeightedEmulator(inner, [0.2, 0.5, 0.3]);
        var model = new CalibrationModel(CalibrationFixture.Config(inner), CalibrationFixture.Observations(inner), emulator);
        int[] outer = NestedCalibrationSampler.OuterDraws(emulator, 9, new Xoshiro256PlusPlus(4));
        Assert.True(outer.Distinct().Count() < outer.Length);   // a resample, with duplicates
        var rng = new Xoshiro256PlusPlus(8);
        CalibrationState[] particles = [.. Enumerable.Range(0, 200).Select(_ => model.DrawPrior(rng))];
        double[] weights = [.. Enumerable.Repeat(1.0 / 200, 200)];
        CalibrationTarget[] targets = [.. Enumerable.Range(0, 3).Select(model.Target)];
        double[,] l = CutSmcOrdering.LogLikelihoods(targets, particles, maxDegreeOfParallelism: 4);
        // The matrix is each target's log-likelihood, and leaves the particles' caches as they were.
        Assert.Equal(targets[1].EvaluateState(new CalibrationState(particles[5].Unit.ToArray(), [], particles[5].Variances.ToArray(), particles[5].Coefficients.ToArray())).LogLikelihoodValue, l[1, 5]);
        Assert.Null(particles[5].Prediction);
        double[,] byDraw = CutSmcOrdering.Distances(l, weights);
        var d = new double[outer.Length, outer.Length];
        for (int a = 0; a < outer.Length; a++)
            for (int b = 0; b < outer.Length; b++) d[a, b] = byDraw[outer[a], outer[b]];
        int[] tour = CutSmcOrdering.Tour(d, 0);
        Assert.Equal(Enumerable.Range(0, outer.Length), tour.Order());
        int[] visited = [.. tour.Select(k => outer[k])];
        _out.WriteLine($"resampled set {string.Join(",", outer)}; the tour visits draws {string.Join(",", visited)}");
        // Each draw's copies are visited in one run, at distance 0 from each other.
        for (int draw = 0; draw < 3; draw++)
        {
            int[] at = [.. Enumerable.Range(0, visited.Length).Where(k => visited[k] == draw)];
            if (at.Length > 0) Assert.Equal(at.Length - 1, at[^1] - at[0]);
        }
    }
}
