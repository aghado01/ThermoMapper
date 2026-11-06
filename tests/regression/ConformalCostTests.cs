using System;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The conformal layer's cost row (E16, slices 0 and 2), run with the qualification switch (<c>BASS_QUALIFICATION</c> = 1) on the
/// clock <c>BASS_COST_CLOCK</c> names: per grid of G = 200 candidates and per candidate, the weights with the errors' sums and
/// each candidate's relative efficiency, the three products (the scores and the errors' two) and the ranks with the margin's
/// band and k̂, with the training side's once-per-stratum likelihoods beside them, at the cost harness's base (n 1,000, S 1,000)
/// and at n 1,000 with S 4,000, for independent draws and for two chains; the product beside the scalar loop it replaces.
/// Indicative: the best of five on a shared machine.
/// </summary>
public sealed class ConformalCostTests(ITestOutputHelper output)
{
    private const int Grid = 200, Repeats = 5;

    private static double Best(Action call)
    {
        call();
        double best = double.PositiveInfinity;
        for (int r = 0; r < Repeats; r++)
        {
            long t0 = BassCost.Tick();
            call();
            best = Math.Min(best, (BassCost.Tick() - t0) / BassCost.TicksPerSecond * 1e3);
        }
        return best;   // ms
    }

    [Theory]
    [InlineData(1000, 1000, 1)]
    [InlineData(1000, 1000, 2)]
    [InlineData(1000, 4000, 1)]
    [InlineData(1000, 4000, 2)]
    public void CostRow(int n, int s, int chains)
    {
        if (!BassQualification.Enabled) return;
        var rng = new Xoshiro256PlusPlus(n + s);
        var ll = new double[s, n];
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++) ll[d, i] = -1 - 4 * rng.NextDouble();
        var lc = new double[s, Grid];
        for (int d = 0; d < s; d++)
            for (int g = 0; g < Grid; g++) lc[d, g] = -1 - 6 * rng.NextDouble() * (1 + g / 50.0);
        LatentReport[] latents = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];
        int[]? lengths = chains == 1 ? null : [s / 2, s - s / 2];
        var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), latents, lengths);
        var candidates = new CandidateLogDensity(lc, [[.. Enumerable.Range(0, Grid).Select(g => (double)g)]], latents, chainLengths: lengths);
        int[] rows = [.. Enumerable.Range(0, s)], members = [.. Enumerable.Range(0, n)];

        ConformalBayes.StratumScorer scorer = null!;
        int level = (int)Math.Floor(0.1 * (n + 1));
        double training1 = Best(() => scorer = new ConformalBayes.StratumScorer(ll, rows, null, members, lengths, level, 3.5));
        var reading = new ConformalBayes.Reading(Grid);
        double weigh = Best(() => scorer.Weigh(lc, 0, Grid, reading));
        double product = Best(() => scorer.Multiply(Grid));
        // Count overwrites the score product with the scores, so each timed call reads a fresh product.
        double count = Best(() => { scorer.Multiply(Grid); scorer.Count(lc, 0, Grid, reading); }) - product;
        double step = Best(() => ConformalBayes.FullStep(training, candidates, new ConformalOptions { Alpha = 0.1 }));

        // The scalar loop the product replaces, on the same operands: one fused multiply-add at a time, rows by columns.
        var left = new double[n * s];
        var right = new double[s * Grid];
        for (int i = 0; i < left.Length; i++) left[i] = rng.NextDouble();
        for (int i = 0; i < right.Length; i++) right[i] = rng.NextDouble();
        var c = new double[n * Grid];
        PackedColumns packed = DenseProduct.Pack(right, s, Grid);
        double kernel = Best(() => DenseProduct.Multiply(left, n, packed, c));
        double scalar = Best(() =>
        {
            for (int i = 0; i < n; i++)
                for (int j = 0; j < Grid; j++)
                {
                    double sum = 0.0;
                    for (int k = 0; k < s; k++) sum = Math.FusedMultiplyAdd(left[i * s + k], right[k * Grid + j], sum);
                    c[i * Grid + j] = sum;
                }
        });
        double gflops = 2.0 * n * s * Grid / (kernel * 1e-3) / 1e9;
        string clock = BassCost.Cycles ? "cycles" : "wall";
        output.WriteLine($"n {n}, S {s}, {(chains == 1 ? "independent draws" : "two chains")}, G {Grid} ({clock} clock, indicative, best of {Repeats}):");
        output.WriteLine($"  training side (once per stratum): {training1:F2} ms");
        output.WriteLine($"  per grid: weights {weigh:F2} ms, the three products {product:F2} ms, ranks {count:F2} ms; the full step {step:F2} ms");
        output.WriteLine($"  per candidate: weights {1e3 * weigh / Grid:F1} µs, products {1e3 * product / Grid:F1} µs, ranks {1e3 * count / Grid:F1} µs");
        output.WriteLine($"  the product alone {kernel:F2} ms ({gflops:F1} GFLOP/s), the scalar loop {scalar:F2} ms: {scalar / kernel:F1}× faster");
    }
}
