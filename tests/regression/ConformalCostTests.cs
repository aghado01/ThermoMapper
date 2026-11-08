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
/// and at n 1,000 with S 4,000, for independent draws and for two chains; the product beside the scalar loop it replaces. The
/// split form's row (slice 3): the calibration's predictive once per stratum, then per grid the predictive and the ranks, and CQR's
/// and the residual's calibration and intervals. Indicative: the best of five on a shared machine.
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
        // The split-half diagnostic: two plain-mode scorers on the halves' rows, their likelihoods built once, then a grid each.
        int[] halfA, halfB;
        if (lengths is null) { halfA = [.. rows.Where(d => d % 2 == 0)]; halfB = [.. rows.Where(d => d % 2 == 1)]; }
        else { halfA = [.. Enumerable.Range(0, s / 4).Concat(Enumerable.Range(s / 2, s / 4))]; halfB = [.. Enumerable.Range(s / 4, s / 4).Concat(Enumerable.Range(3 * s / 4, s / 4))]; }
        ConformalBayes.StratumScorer a = null!, b = null!;
        double halvesTraining = Best(() =>
        {
            a = new ConformalBayes.StratumScorer(ll, halfA, null, members, null, level, 1.0, plain: true);
            b = new ConformalBayes.StratumScorer(ll, halfB, null, members, null, level, 1.0, plain: true);
        });
        var ra = new ConformalBayes.Reading(Grid);
        var rb = new ConformalBayes.Reading(Grid);
        double halves = Best(() => { a.Rank(lc, 0, Grid, ra); b.Rank(lc, 0, Grid, rb); });

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
        output.WriteLine($"  training side (once per stratum): {training1:F2} ms, the halves' {halvesTraining:F2} ms");
        output.WriteLine($"  the split-half diagnostic per grid: {halves:F2} ms; the plain step's own parts, the weights, the score product and the ranks, are slice 0's row");
        output.WriteLine($"  per grid: weights {weigh:F2} ms, the three products {product:F2} ms, ranks {count:F2} ms; the full step {step:F2} ms");
        output.WriteLine($"  per candidate: weights {1e3 * weigh / Grid:F1} µs, products {1e3 * product / Grid:F1} µs, ranks {1e3 * count / Grid:F1} µs");
        output.WriteLine($"  the product alone {kernel:F2} ms ({gflops:F1} GFLOP/s), the scalar loop {scalar:F2} ms: {scalar / kernel:F1}× faster");
    }

    [Theory]
    [InlineData(1000, 1000)]
    [InlineData(1000, 4000)]
    public void SplitCostRow(int n, int s)
    {
        // The split form at n calibration units and S draws: the density score's calibration side (n predictives, once per stratum),
        // then per grid of G its predictive and ranks, the whole step on one grid; CQR's and the residual's calibration (n scores and
        // a sort) and their intervals per case.
        if (!BassQualification.Enabled) return;
        var rng = new Xoshiro256PlusPlus(n + s + 1);
        var lcal = new double[s, n];
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++) lcal[d, i] = -1 - 4 * rng.NextDouble();
        var lc = new double[s, Grid];
        for (int d = 0; d < s; d++)
            for (int g = 0; g < Grid; g++) lc[d, g] = -1 - 6 * rng.NextDouble() * (1 + g / 50.0);
        LatentReport[] latents = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];
        var calibration = new CandidateLogDensity(lcal, [.. Enumerable.Range(0, n).Select(i => new[] { (double)i })], latents);
        var candidates = new CandidateLogDensity(lc, [[.. Enumerable.Range(0, Grid).Select(g => (double)g)]], latents);
        int[] rows = [.. Enumerable.Range(0, s)];
        var scores = new double[n];
        double calibrationSide = Best(() =>
        {
            SplitConformal.LogPredictive(lcal, rows, null, 0, n, scores);
            SplitConformal.Calibrate(scores);
        });
        SplitCalibration fit = SplitConformal.Calibrate(scores);
        var own = new double[Grid];
        double perGrid = Best(() =>
        {
            SplitConformal.LogPredictive(lc, rows, null, 0, Grid, own);
            for (int g = 0; g < Grid; g++) fit.Rank(0, own[g]);
        });
        double step = Best(() => SplitConformal.Density(calibration, candidates));
        double[] lo = [.. Enumerable.Range(0, n).Select(_ => -1 - rng.NextDouble())], hi = [.. lo.Select(v => v + 2)];
        double[] y = [.. Enumerable.Range(0, n).Select(_ => 2 * rng.NextDouble() - 1)];
        var test = new QuantileReadings(lo[..Grid], hi[..Grid]);
        var calibrationReadings = new QuantileReadings(lo, hi);
        double cqr = Best(() => SplitConformal.Quantiles(calibrationReadings, y, test));
        var centres = new ResidualReadings(y, [.. lo.Select(v => -v)]);
        double residual = Best(() => SplitConformal.Residuals(centres, y, centres with { Center = y[..Grid], Scale = [.. lo[..Grid].Select(v => -v)] }));
        string clock = BassCost.Cycles ? "cycles" : "wall";
        output.WriteLine($"the split form, n {n} calibration units, S {s}, G {Grid} ({clock} clock, indicative, best of {Repeats}):");
        output.WriteLine($"  the density score: calibration side (once per stratum) {calibrationSide:F2} ms; per grid {perGrid:F2} ms, "
                         + $"{1e3 * perGrid / Grid:F1} µs per candidate; the whole step on one grid {step:F2} ms");
        output.WriteLine($"  CQR, calibration and {Grid} intervals: {cqr:F3} ms; the residual's: {residual:F3} ms");
    }
}
