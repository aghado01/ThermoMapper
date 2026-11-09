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
/// and the residual's calibration and intervals. The predictive distribution's row (slice 4): the split form's calibration PITs once per
/// stratum and its whole call on one grid, beside the full step's opt-in, two products a batch. Indicative: the best of five on a shared
/// machine.
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
        // The default step reads no margin; asked for, the margin adds each unit's error under every unit's weights, three n × S × n
        // products once per stratum, a fourth product per grid and the bands.
        var reading = new ConformalBayes.Reading(Grid);
        double training1 = Best(() => scorer = new ConformalBayes.StratumScorer(ll, rows, null, members, lengths, level, 3.5));
        double weigh = Best(() => scorer.Weigh(lc, 0, Grid, reading));
        double product = Best(() => scorer.Multiply(Grid));
        // Count overwrites the score product with the scores, so each timed call reads a fresh product.
        double count = Best(() => { scorer.Multiply(Grid); scorer.Count(lc, 0, Grid, reading); }) - product;
        double step = Best(() => ConformalBayes.FullStep(training, candidates, new ConformalOptions { Alpha = 0.1 }));
        double orderings = Best(() => scorer = new ConformalBayes.StratumScorer(ll, rows, null, members, lengths, level, 3.5, MarginBounds.Carried));
        scorer.Weigh(lc, 0, Grid, reading);
        double marginProduct = Best(() => scorer.Multiply(Grid));
        double marginCount = Best(() => { scorer.Multiply(Grid); scorer.Count(lc, 0, Grid, reading); }) - marginProduct;
        double marginStep = Best(() => ConformalBayes.FullStep(training, candidates, new ConformalOptions { Alpha = 0.1, MarginBounds = MarginBounds.Carried }));
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
        output.WriteLine($"  training side (once per stratum): {training1:F2} ms, with the margin's orderings {orderings:F2} ms, the halves' {halvesTraining:F2} ms");
        output.WriteLine($"  the split-half diagnostic per grid: {halves:F2} ms; the plain step's own parts, the weights, the score product and the ranks, are slice 0's row");
        output.WriteLine($"  per grid: weights {weigh:F2} ms, the three products {product:F2} ms, ranks {count:F2} ms; the full step {step:F2} ms");
        output.WriteLine($"  with the margin: the four products {marginProduct:F2} ms, ranks and bands {marginCount:F2} ms; the full step {marginStep:F2} ms");
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

    [Theory]
    [InlineData(1000, 1000)]
    [InlineData(1000, 4000)]
    public void DistributionCostRow(int n, int s)
    {
        // The predictive distribution at n units and S draws, one grid of G: the split form's calibration PITs (n by both tails and a
        // sort, once per stratum), then per grid the candidate's CDF by both tails and its ranks, and its whole call; the full step's
        // whole call, two (n × S)·(S × G) products and the counts.
        if (!BassQualification.Enabled) return;
        var rng = new Xoshiro256PlusPlus(n + s + 2);
        var ll = new double[s, n];
        var lcdf = new double[s, n];
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++) { ll[d, i] = -1 - 4 * rng.NextDouble(); lcdf[d, i] = Math.Log(rng.NextDouble() + 1e-300); }
        var lc = new double[s, Grid];
        var cc = new double[s, Grid];
        for (int d = 0; d < s; d++)
            for (int g = 0; g < Grid; g++) { lc[d, g] = -1 - 6 * rng.NextDouble() * (1 + g / 50.0); cc[d, g] = Math.Log((g + rng.NextDouble()) / (Grid + 1)); }
        LatentReport[] latents = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];
        var calibration = new CandidateLogDensity(ll, [.. Enumerable.Range(0, n).Select(i => new[] { (double)i })], latents) { LogCdf = lcdf };
        var candidates = new CandidateLogDensity(lc, [[.. Enumerable.Range(0, Grid).Select(g => (double)g)]], latents) { LogCdf = cc };
        var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), latents) { LogCdf = lcdf };
        int[] rows = [.. Enumerable.Range(0, s)];
        double[] weights = CdfScores.Weights(rows, null);
        var scores = new CdfScore[n];
        double calibrationSide = Best(() =>
        {
            CdfScores.Fixed(lcdf, rows, weights, 0, n, scores);
            Array.Sort(scores);
        });
        CdfScore[] sorted = [.. scores];
        var own = new CdfScore[Grid];
        double perGrid = Best(() =>
        {
            CdfScores.Fixed(cc, rows, weights, 0, Grid, own);
            for (int g = 0; g < Grid; g++) CdfScores.Rank(sorted, own[g]);
        });
        double split = Best(() => SplitConformal.Distribution(calibration, candidates));
        double full = Best(() => ConformalBayes.Distribution(training, candidates));
        string clock = BassCost.Cycles ? "cycles" : "wall";
        output.WriteLine($"the predictive distribution, n {n} units, S {s}, G {Grid} ({clock} clock, indicative, best of {Repeats}):");
        output.WriteLine($"  the split form: calibration PITs (once per stratum) {calibrationSide:F2} ms; per grid {perGrid:F2} ms, "
                         + $"{1e3 * perGrid / Grid:F1} µs per candidate; the whole call on one grid {split:F2} ms");
        output.WriteLine($"  the full step's opt-in: the whole call on one grid {full:F2} ms, {1e3 * full / Grid:F1} µs per candidate");
    }

    [Theory]
    [InlineData(1000, 1000)]
    [InlineData(1000, 4000)]
    public void ExactEndsCostRow(int n, int s)
    {
        if (!BassQualification.Enabled) return;
        var rng = new Xoshiro256PlusPlus(20261109);
        double[] beta = [.. Enumerable.Range(0, s).Select(_ => rng.NextDouble())];
        double[] x = [.. Enumerable.Range(0, n).Select(_ => 2 * rng.NextDouble() - 1)];
        double[][] y = [.. x.Select(v => new[] { ConformalBayesTests.Truth(v, 0.5, rng) })];
        var cal = ConformalPredictiveDistributionTests.Seam(beta, 1, null, x, y);
        // At x = 0 the common N(0,1) predictive has a native inverse. Time the actual candidate seam, rank and readouts on both paths.
        CandidateLogDensity Test(double[] g) => ConformalPredictiveDistributionTests.Seam(beta, 1, null, [0], [g]);
        SplitPitCalibration fit = null!;
        double sort = Best(() => fit = SplitConformal.CalibrateDistribution(cal));
        double[] grid = [.. Enumerable.Range(0, Grid).Select(i => -6 + 12.0 * i / (Grid - 1))];
        double gridded = Best(() =>
        {
            var d = fit.OnGrid(Test(grid)).Distributions[0];
            _ = d.Quantile(0.5);
            _ = d.Intervals(0.1);
        });
        int calls = 0;
        double Measure(Func<double, double>? inverse) => Best(() =>
        {
            calls = 0;
            var d = fit.AtCase(v => { calls++; return Test([v]); }, predictiveQuantile: inverse);
            _ = d.Interval(0.1);
            _ = d.Quantile(0.5);
        });
        double native = Measure(Maths.Distributions.Normal.StandardQuantile);
        int nativeProbes = calls;
        double root = Measure(null);
        output.WriteLine($"n {n}, S {s}, G {Grid}, {(BassCost.Cycles ? "cycles" : "wall")}, best of five: calibration {sort:F2} ms; " +
                         $"median and 90% ends: grid {gridded:F2} ms ({1000 * gridded / Grid:F2} µs/candidate), native {native:F2} ms ({nativeProbes} probes), CDF root {root:F2} ms ({calls} probes)");
    }

    [Fact]
    public void RefinedBoundariesCostRow()
    {
        if (!BassQualification.Enabled) return;
        const int n = 39, s = 40;
        var rng = new Xoshiro256PlusPlus(20261109);
        double[] xs = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        double[][] ys = [.. xs.Select(x => new[] { ConformalBayesTests.Truth(x, 0.5, rng) })];
        double[] beta = [.. Enumerable.Range(0, s).Select(_ => rng.NextDouble())];
        var cal = ConformalPredictiveDistributionTests.Seam(beta, 1, null, xs, ys);
        ConformalSet Set(double[] grid) => SplitConformal.Density(cal, ConformalPredictiveDistributionTests.Seam(beta, 1, null, [0.2], [grid])).Sets[0];
        double[] coarse = [.. Enumerable.Range(0, 17).Select(i => -4.0 + i * 0.5)];
        double[] fine = [.. Enumerable.Range(0, 8001).Select(i => -4.0 + i * 0.001)];
        int calls = 0;
        double gridCost = Best(() => _ = Set(fine).Intervals);
        double cost = Best(() =>
        {
            calls = 0;
            var bounds = Set(coarse).RefineBoundaries(y => { calls++; return Set([y]); }, 0.001);
            Assert.Equal(2, bounds.Count);
        });
        output.WriteLine($"split density, n {n}, S {s}, resolution 0.001, {(BassCost.Cycles ? "cycles" : "wall")}, best of five: " +
                         $"grid G={fine.Length} {gridCost:F2} ms; grid G={coarse.Length} plus {calls} singleton probes {cost:F2} ms (calibration reread at each probe)");
    }
}
