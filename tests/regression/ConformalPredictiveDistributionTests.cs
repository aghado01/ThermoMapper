using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The conformal predictive distribution (E16, slice 4; the user's decision of 2026-11-08): the split form's, on the fit's predictive
/// CDF, the default (<see cref="SplitConformal.Distribution"/>), and the full step's on the in-sample predictive CDF, an opt-in
/// that flags the cases whose score falls (<see cref="ConformalBayes.Distribution"/>). On the conjugate model: the counts against
/// references in seeded sweeps, the readouts on counts given by hand, the PIT uniform on held-out data for any fit
/// (<c>coverage_tieBreak</c>, <c>rank_uniform</c>, <c>split_coverage</c>), Vovk's conditions on every grid, the full step's PIT and
/// flag, agreement with slice 3's sets on a fit that is right, and CRPS against a quadrature reference. With the qualification switch,
/// the measured default: both forms side by side on the same fits.
/// </summary>
public sealed class ConformalPredictiveDistributionTests(ITestOutputHelper output)
{
    private static readonly LatentReport[] Theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];
    private const double Z95 = 1.6448536269514722;

    private static double LogNormal(double y, double mean, double sd)
        => -0.5 * Math.Log(2 * Math.PI) - Math.Log(sd) - 0.5 * (y - mean) * (y - mean) / (sd * sd);

    /// <summary>A conjugate fit's seam at cases <paramref name="xs"/> on <paramref name="grids"/>: ln N(y | βx, sd²) and ln Φ((y − βx)/sd)
    /// per draw.</summary>
    internal static CandidateLogDensity Seam(double[] beta, double sd, int[]? chains, double[] xs, double[][] grids)
    {
        int s = beta.Length, columns = grids.Sum(g => g.Length);
        var values = new double[s, columns];
        var logCdf = new double[s, columns];
        for (int d = 0; d < s; d++)
            for (int t = 0, c = 0; t < xs.Length; t++)
                foreach (double v in grids[t])
                {
                    double mean = beta[d] * xs[t];
                    values[d, c] = LogNormal(v, mean, sd);
                    logCdf[d, c++] = Normal.StandardLogCdf((v - mean) / sd);
                }
        return new CandidateLogDensity(values, grids, Theta, chainLengths: chains) { LogCdf = logCdf };
    }

    /// <summary>The training density with its log CDF, each unit one observation.</summary>
    internal static PointwiseLogDensity Training(double[] beta, double sd, int[]? chains, double[] xs, double[] ys)
    {
        CandidateLogDensity at = Seam(beta, sd, chains, xs, [.. ys.Select(v => new[] { v })]);
        return new PointwiseLogDensity(at.Values, UnitPartition.EachObservation(xs.Length), Theta, chains) { LogCdf = at.LogCdf };
    }

    /// <summary>The fit's own predictive CDF, Σ_s F(y | θ_s)/S, at each grid value.</summary>
    private static double[] OwnCdf(double[] beta, double sd, double x, double[] grid)
        => [.. grid.Select(v => beta.Average(b => Math.Exp(Normal.StandardLogCdf((v - b * x) / sd))))];

    /// <summary>A CDF's CRPS on a grid, read as the layer reads Q̄: linear between grid values, tails cut at the grid's ends.</summary>
    internal static double GridCrps(double[] grid, double[] cdf, double y)
    {
        static double Piece(double a, double b, double ua, double ub) => (b - a) * (ua * ua + ua * ub + ub * ub) / 3.0;
        int last = grid.Length - 1;
        double total = Math.Max(0.0, grid[0] - y) + Math.Max(0.0, y - grid[last]);
        for (int g = 0; g < last; g++)
        {
            double a = grid[g], b = grid[g + 1], qa = cdf[g], qb = cdf[g + 1];
            if (y <= a) total += Piece(a, b, qa - 1, qb - 1);
            else if (y >= b) total += Piece(a, b, qa, qb);
            else
            {
                double qy = qa + (y - a) / (b - a) * (qb - qa);
                total += Piece(a, y, qa, qy) + Piece(y, b, qy - 1, qb - 1);
            }
        }
        return total;
    }

    /// <summary>
    /// Ten bins of the PITs against a tenth each. Given the calibration (or the training set), a bin's chance is a sum of the
    /// calibration's spacings, Beta(wm, (1 − w)m) for width w, so its variance w(1 − w)/(m + 1) over r replications, plus the cases'
    /// binomial remainder over the PITs. Returns each bin's share, that standard error and the Kolmogorov–Smirnov distance.
    /// </summary>
    internal static (double[] Share, double Se, double Ks) Uniformity(IReadOnlyList<double> pits, int m, int r)
    {
        const double w = 0.1;
        var share = new double[10];
        foreach (double p in pits) share[Math.Min(9, (int)(p * 10))] += 1.0 / pits.Count;
        double between = w * (1 - w) / (m + 1);
        double se = Math.Sqrt(between / r + (w * (1 - w) - between) / pits.Count);
        double[] sorted = [.. pits.Order()];
        double ks = 0.0;
        for (int i = 0; i < sorted.Length; i++) ks = Math.Max(ks, Math.Max((i + 1.0) / sorted.Length - sorted[i], sorted[i] - (double)i / sorted.Length));
        return (share, se, ks);
    }

    private void ReportUniformity(string label, IReadOnlyList<double> pits, int m, int r, bool assert)
    {
        (double[] share, double se, double ks) = Uniformity(pits, m, r);
        output.WriteLine($"  {label}: {pits.Count} PITs, bins {string.Join(" ", share.Select(v => v.ToString("F3")))} (se {se:F4}), KS {ks:F4}");
        if (assert) Assert.All(share, v => Assert.True(Math.Abs(v - 0.1) <= 4 * se, $"{label}: a bin holds {v:F4}"));
    }

    // ---- the readouts ---------------------------------------------------------------------------------------------------------------

    [Fact]
    public void Readouts_OnCountsGivenByHand()
    {
        // m = 5, one unit below per grid step: Q⁻ = B/5, Q⁺ = (B + 1)/5, Q̄ = (B + ½)/5.
        var d = new ConformalDistribution(0, 0, 5, [0, 1, 2, 3, 4], [0, 1, 2, 3, 4], [0, 0, 0, 0, 0], ConformalForm.Split, 0, 0.0);
        Assert.Equal([0.0, 0.2, 0.4, 0.6, 0.8], Enumerable.Range(0, 5).Select(d.Lower));
        Assert.Equal([0.2, 0.4, 0.6, 0.8, 1.0], Enumerable.Range(0, 5).Select(d.Upper));
        Assert.Equal(0.3, d.At(1, 0.5), 15);
        Assert.Throws<ArgumentOutOfRangeException>(() => d.At(0, 1.5));
        // α 0.4: k = ⌊0.4·5/2⌋ = 1 per tail, so B + T + 1 > 1 and 5 − B > 1: grid values 1 to 3, coverage 1 − 2/5.
        Assert.Equal(1, d.Level(0.4));
        Assert.Equal(0.6, d.Coverage(0.4), 15);
        Assert.Equal([false, true, true, true, false], Enumerable.Range(0, 5).Select(g => d.Accepts(g, 0.4)));
        Assert.Equal([new GridInterval(1, 3)], d.Intervals(0.4));
        Assert.Equal(3.0, d.Width(0.4), 15);
        Assert.True(d.WholeSpace(0.3));   // k = ⌊0.75⌋ = 0
        Assert.Equal(2.0, d.Quantile(0.5));
        Assert.Equal(double.NegativeInfinity, d.Quantile(0.05));
        Assert.Equal(double.PositiveInfinity, d.Quantile(0.95));
        // CRPS against a midpoint sum of the same integrand, 10⁶ pieces, at a y inside, on, and outside the grid.
        foreach (double y in new[] { 2.3, 3.0, -1.5, 6.0 })
        {
            double reference = 0.0, h = 10.0 / 1e6;
            for (int i = 0; i < 1_000_000; i++)
            {
                double z = -3.0 + (i + 0.5) * h;
                double q, spread;
                if (z < 0) { q = 0; spread = 0; }
                else if (z > 4) { q = 1; spread = 0; }
                else
                {
                    int g = Math.Min(3, (int)z);
                    double t = z - g;
                    q = (1 - t) * d.Mid(g) + t * d.Mid(g + 1);
                    spread = 0.2;
                }
                double u = q - (z >= y ? 1 : 0);
                reference += h * (u * u + spread * spread / 12);
            }
            Assert.Equal(reference, d.Crps(y), 6);
        }
        // A case whose counts fall has no distribution: its quantiles and CRPS are NaN, its set stands.
        int[] below = [0, 2, 1, 3, 4], tied = [0, 0, 0, 0, 0];
        Assert.Equal(1, CdfScores.Fall(below, tied));
        Assert.Equal(2, CdfScores.Fall([0, 3, 1], [0, 0, 1]));
        var falling = new ConformalDistribution(0, 0, 5, [0, 1, 2, 3, 4], below, tied, ConformalForm.Full, 1, 0.1);
        Assert.False(falling.Monotone);
        Assert.True(double.IsNaN(falling.Quantile(0.5)));
        Assert.True(double.IsNaN(falling.Crps(2.0)));
        Assert.NotEmpty(falling.Intervals(0.4));
    }

    // ---- the counts against references ----------------------------------------------------------------------------------------------

    // A log CDF for the sweeps: the body, near 1 (ln F a few 10⁻²⁰ below 0), the far lower tail, or −∞.
    private static double SweepLogCdf(Xoshiro256PlusPlus rng) => rng.NextDouble() switch
    {
        < 0.5 => Math.Log(rng.NextDouble() + 1e-300),
        < 0.75 => -1e-20 * (1 + 9 * rng.NextDouble()),
        < 0.95 => -700 * rng.NextDouble(),
        _ => double.NegativeInfinity,
    };

    [Fact]
    public void Split_CountsAreTheReference_InASeededSweep()
    {
        // 400 trials: weighted draws with zeros or chains, strata with no calibration unit, both tails far out, copied columns. The
        // reference sums each tail in row order under the normalized weights and counts by the tail that holds the digits.
        int[] seeds = SeedTree.Derive(2026_1110, 400);
        int compared = 0;
        foreach (int seed in seeds)
        {
            var rng = new Xoshiro256PlusPlus(seed);
            int s = 1 + (int)(rng.NextDouble() * 40), n = (int)(rng.NextDouble() * 12), cases = 1 + (int)(rng.NextDouble() * 3);
            double[]? weights = rng.NextDouble() < 0.4 ? [.. Enumerable.Range(0, s).Select(_ => rng.NextDouble() < 0.2 ? 0.0 : rng.NextDouble())] : null;
            if (weights is not null && weights.Sum() == 0.0) weights[0] = 1.0;
            int[]? chains = weights is null && s > 1 && rng.NextDouble() < 0.5 ? [s / 2, s - s / 2] : null;
            var lcal = new double[s, Math.Max(n, 1)];
            var cal = new double[s, Math.Max(n, 1)];
            for (int d = 0; d < s; d++)
                for (int i = 0; i < n; i++) { lcal[d, i] = -rng.NextDouble(); cal[d, i] = SweepLogCdf(rng); }
            if (n == 0) { for (int d = 0; d < s; d++) cal[d, 0] = -1.0; }
            int units = Math.Max(n, 1);
            int[] calStrata = [.. Enumerable.Range(0, units).Select(_ => (int)(rng.NextDouble() * 2))];
            double[][] grids = [.. Enumerable.Range(0, cases).Select(_ => Enumerable.Range(0, 1 + (int)(rng.NextDouble() * 4)).Select(g => (double)g).ToArray())];
            int columns = grids.Sum(g => g.Length);
            var lc = new double[s, columns];
            var cc = new double[s, columns];
            for (int d = 0; d < s; d++)
                for (int c = 0; c < columns; c++)
                {
                    lc[d, c] = -rng.NextDouble();
                    // A third of the columns copy a calibration unit's, so they tie with it.
                    cc[d, c] = n > 0 && c % 3 == 0 ? cal[d, c % n] : SweepLogCdf(rng);
                }
            int[] testStrata = [.. Enumerable.Range(0, cases).Select(_ => (int)(rng.NextDouble() * 3))];
            var calibration = new CandidateLogDensity(lcal, [.. Enumerable.Range(0, units).Select(i => new[] { (double)i })], Theta, calStrata) { LogCdf = cal };
            var candidates = new CandidateLogDensity(lc, grids, Theta, testStrata);
            calibration = weights is null ? calibration with { ChainLengths = chains } : calibration with { Weights = weights };
            candidates = weights is null ? candidates with { ChainLengths = chains, LogCdf = cc } : candidates with { Weights = weights, LogCdf = cc };
            ConformalDistributions result = SplitConformal.Distribution(calibration, candidates);

            double total = 0.0;   // in order, as PredictiveChecks.Rows sums it
            if (weights is not null) foreach (double w in weights) total += w;
            ConformalDistributionTests.TailScore Score(double[,] logCdf, int column)
            {
                double lower = 0.0, upper = 0.0;
                for (int d = 0; d < s; d++)
                {
                    // The layer's weights: the normalized logs of the draws of positive weight, exponentiated.
                    if (weights is not null && !(weights[d] > 0.0)) continue;
                    double w = weights is null ? 1.0 / s : Math.Exp(Math.Log(weights[d] / total));
                    // 1 − F directly at or below the median, −expm1 above it, as the layer reads both tails.
                    double l = logCdf[d, column];
                    lower += w * Math.Exp(l);
                    upper += w * (l <= -0.69314718055994530942 ? 1.0 - Math.Exp(l) : -SpecialFunctions.Expm1(l));
                }
                return new(lower, upper);
            }
            for (int t = 0, c = 0; t < cases; t++)
            {
                ConformalDistribution dist = result.Distributions[t];
                int[] members = [.. Enumerable.Range(0, units).Where(i => calStrata[i] == testStrata[t])];
                Assert.Equal(members.Length + 1, dist.Units);
                for (int g = 0; g < grids[t].Length; g++, c++)
                {
                    ConformalDistributionTests.TailScore own = Score(cc, c);
                    int below = members.Count(i => Score(cal, i).CompareTo(own) < 0), tied = members.Count(i => Score(cal, i).CompareTo(own) == 0);
                    Assert.Equal((below, tied), (dist.Below[g], dist.Tied[g]));
                    compared++;
                }
            }
        }
        output.WriteLine($"{compared} grid values against the reference, every count equal");
    }

    [Fact]
    public void FullStep_CountsAreTheReference_InASeededSweep()
    {
        // 300 trials of the full step's two products against the direct sums: the counts agree, a quarter of the columns copied from a
        // unit so they tie with it exactly, except where a unit's score and the candidate's differ within 10⁻¹² on their tail, where
        // the products' order of summation decides.
        int[] seeds = SeedTree.Derive(2026_1111, 300);
        int compared = 0, ambiguous = 0;
        foreach (int seed in seeds)
        {
            var rng = new Xoshiro256PlusPlus(seed);
            int s = 1 + (int)(rng.NextDouble() * 50), n = 1 + (int)(rng.NextDouble() * 14), g = 1 + (int)(rng.NextDouble() * 20);
            var ll = new double[s, n];
            var lcdf = new double[s, n];
            for (int d = 0; d < s; d++)
                for (int i = 0; i < n; i++) { ll[d, i] = -3 * rng.NextDouble(); lcdf[d, i] = SweepLogCdf(rng); }
            int[] strata = [.. Enumerable.Range(0, n).Select(_ => (int)(rng.NextDouble() * 2))];
            double[]? weights = rng.NextDouble() < 0.4 ? [.. Enumerable.Range(0, s).Select(_ => rng.NextDouble() < 0.2 ? 0.0 : rng.NextDouble())] : null;
            if (weights is not null && weights.Sum() == 0.0) weights[0] = 1.0;
            var training = new PointwiseLogDensity(ll, new UnitPartition([.. Enumerable.Range(0, n)], strata), Theta) { LogCdf = lcdf };
            var lc = new double[s, g];
            var cc = new double[s, g];
            for (int d = 0; d < s; d++)
                for (int k = 0; k < g; k++)
                {
                    lc[d, k] = -3 * rng.NextDouble();
                    cc[d, k] = k % 4 == 0 ? lcdf[d, k % n] : SweepLogCdf(rng);
                }
            int stratum = (int)(rng.NextDouble() * 3);
            var candidates = new CandidateLogDensity(lc, [[.. Enumerable.Range(0, g).Select(k => (double)k)]], Theta, [stratum]) { LogCdf = cc };
            if (weights is not null)
            {
                training = training with { Weights = weights };
                candidates = candidates with { Weights = weights };
            }
            ConformalDistribution dist = ConformalBayes.Distribution(training, candidates).Distributions[0];
            int[] members = [.. Enumerable.Range(0, n).Where(i => strata[i] == stratum)];
            double[]? logWeights = weights is null ? null : [.. weights.Select(w => Math.Log(w))];
            (ConformalDistributionTests.TailScore[][] units, ConformalDistributionTests.TailScore[] own) =
                ConformalDistributionTests.Reference(lcdf, members, lc, cc, 0, g, logWeights);
            for (int k = 0; k < g; k++)
            {
                bool close = units[k].Any(u => u.CompareTo(own[k]) != 0
                                               && Math.Abs(u.Minus(own[k])) <= 1e-12 * Math.Max(1e-300, Math.Min(Math.Min(u.Lower, u.Upper), Math.Min(own[k].Lower, own[k].Upper))));
                int below = units[k].Count(u => u.CompareTo(own[k]) < 0), tied = units[k].Count(u => u.CompareTo(own[k]) == 0);
                if (close) { ambiguous++; continue; }
                Assert.Equal((below, tied), (dist.Below[k], dist.Tied[k]));
                compared++;
            }
        }
        output.WriteLine($"{compared} grid values against the reference, every count equal; {ambiguous} with a unit within rounding of the candidate and not equal to it, not compared");
    }

    [Fact]
    public void Refuses_WhatIsNotAPredictiveDistributionsInput()
    {
        double[] beta = [0.5, 0.6, 0.7];
        CandidateLogDensity calibration = Seam(beta, 1.0, null, [0.1, 0.2], [[0.0], [1.0]]);
        CandidateLogDensity candidates = Seam(beta, 1.0, null, [0.3], [[0.0, 1.0]]);
        SplitConformal.Distribution(calibration, candidates);
        Assert.Throws<ArgumentException>(() => SplitConformal.Distribution(calibration with { LogCdf = null }, candidates));
        Assert.Throws<ArgumentException>(() => SplitConformal.Distribution(calibration, candidates with { LogCdf = null }));
        Assert.Throws<ArgumentException>(() => SplitConformal.Distribution(Seam(beta, 1.0, null, [0.1], [[0.0, 1.0]]), candidates));
        Assert.Throws<ArgumentException>(() => SplitConformal.Distribution(calibration, Seam([0.5, 0.6], 1.0, null, [0.3], [[0.0]])));
        PointwiseLogDensity training = Training(beta, 1.0, null, [0.1, 0.2], [0.0, 1.0]);
        ConformalBayes.Distribution(training, candidates);
        Assert.Throws<ArgumentException>(() => ConformalBayes.Distribution(training with { LogCdf = null }, candidates));
        Assert.Throws<ArgumentException>(() => ConformalBayes.Distribution(training, candidates with { LogCdf = null }));
        Assert.Throws<ArgumentException>(() => ConformalBayes.Distribution(
            training with { Latents = [new("θ", LatentScope.Global, LatentTreatment.IntegratedGivenOthers)] }, candidates));
    }

    // ---- the split form: PIT and Vovk's conditions ----------------------------------------------------------------------------------

    private sealed record Fit(string Name, double[] Beta, double NoiseSd, int[]? Chains);

    // Slice 3's three fits on a training split: the exact posterior, a wrong fit (the slope's sign flipped under an overconfident noise
    // sd 0.25) and a chain stopped early (forty Metropolis steps from β = 5); and two chains stuck apart (the exact draws shifted ∓1).
    private static Fit[] Fits(double[] x, double[] y, int s, Xoshiro256PlusPlus rng, bool apart)
    {
        double[] exact = ConformalBayesTests.ExactDraws(x, y, s, rng);
        var fits = new List<Fit>
        {
            new("the exact posterior", exact, 1.0, null),
            new("a wrong fit", [.. exact.Select(b => -b)], 0.25, null),
            new("a chain stopped early", SplitConformalTests.StoppedEarly(x, y, rng), 1.0, [40]),
        };
        if (apart) fits.Add(new("two chains stuck apart", [.. exact.Select((b, d) => b + (d % 2 == 0 ? -1.0 : 1.0))], 1.0, null));
        return [.. fits];
    }

    [Fact]
    public void Split_PitIsUniform_ForAnyFit_AndVovksConditionsHold()
    {
        // split_coverage with random ties (coverage_tieBreak): 30 training units, 19 calibration units (m 20), 2,000 replications of 20
        // test cases under slice 0's truth, the PIT at each test case's response with τ uniform. On 50 replications' first two test
        // cases, a grid over [−50, 50]: the candidate's score never falls, nor Q, from Q⁻ = 0 below to Q⁺ = 1 above (Vovk's conditions:
        // the calibration's scores do not move). Beside it the fits' own PITs, F̂(y), and the 90% set's width and CRPS against the fit's.
        const int train = 30, calibrationUnits = 19, tests = 20, replications = 2000, gridReps = 50;
        double[] grid = [.. Enumerable.Range(0, 801).Select(g => -50.0 + 0.125 * g)];
        var pits = new List<double>[3];
        var own = new List<double>[3];
        var width = new double[3];
        var crps = new double[3];
        var ownCrps = new double[3];
        for (int f = 0; f < 3; f++) { pits[f] = []; own[f] = []; }
        int[] seeds = SeedTree.Derive(2026_1112, replications);
        for (int r = 0; r < replications; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            (double[] x, double[] y) Units(int count)
            {
                var xs = new double[count];
                var ys = new double[count];
                for (int i = 0; i < count; i++) { xs[i] = 4 * rng.NextDouble() - 2; ys[i] = ConformalBayesTests.Truth(xs[i], 0.5, rng); }
                return (xs, ys);
            }
            (double[] xTrain, double[] yTrain) = Units(train);
            (double[] xCal, double[] yCal) = Units(calibrationUnits);
            (double[] xTest, double[] yTest) = Units(tests);
            Fit[] fits = Fits(xTrain, yTrain, 200, rng, apart: false);
            var tau = new Xoshiro256PlusPlus(seeds[r] ^ 0x5bd1e995);
            for (int f = 0; f < 3; f++)
            {
                Fit fit = fits[f];
                CandidateLogDensity calibration = Seam(fit.Beta, fit.NoiseSd, fit.Chains, xCal, [.. yCal.Select(v => new[] { v })]);
                ConformalDistributions atTruth = SplitConformal.Distribution(calibration, Seam(fit.Beta, fit.NoiseSd, fit.Chains, xTest, [.. yTest.Select(v => new[] { v })]));
                for (int t = 0; t < tests; t++)
                {
                    pits[f].Add(atTruth.Distributions[t].At(0, tau.NextDouble()));
                    own[f].Add(OwnCdf(fit.Beta, fit.NoiseSd, xTest[t], [yTest[t]])[0]);
                }
                if (r >= gridReps) continue;
                ConformalDistributions onGrid = SplitConformal.Distribution(calibration, Seam(fit.Beta, fit.NoiseSd, fit.Chains, xTest[..2], [grid, grid]));
                Assert.Equal(0, onGrid.Flagged);
                for (int t = 0; t < 2; t++)
                {
                    ConformalDistribution d = onGrid.Distributions[t];
                    Assert.Equal(0, d.Fall);
                    Assert.Equal(0.0, d.ScoreFall);
                    Assert.Equal(0.0, d.Lower(0));
                    Assert.Equal(1.0, d.Upper(grid.Length - 1));
                    width[f] += d.Width(0.1) / (2.0 * gridReps);
                    crps[f] += d.Crps(yTest[t]) / (2.0 * gridReps);
                    ownCrps[f] += GridCrps(grid, OwnCdf(fit.Beta, fit.NoiseSd, xTest[t], grid), yTest[t]) / (2.0 * gridReps);
                }
            }
        }
        string[] names = ["the exact posterior", "a wrong fit", "a chain stopped early"];
        for (int f = 0; f < 3; f++)
        {
            output.WriteLine($"{names[f]}: the 90% set's mean width {width[f]:F2}; CRPS {crps[f]:F3} against the fit's own {ownCrps[f]:F3}");
            ReportUniformity("the split form's PIT", pits[f], calibrationUnits + 1, replications, assert: true);
            ReportUniformity("the fit's own PIT, F̂(y)", own[f], calibrationUnits + 1, replications, assert: false);
        }
    }

    // ---- the full step's opt-in: PIT and the flag -----------------------------------------------------------------------------------

    [Fact]
    public void FullStep_PitIsUniform_OnTheExactPosterior_FlaggedOrNot_AndTwoChainsApartAreFlagged()
    {
        // The full step on all 49 units (m 50), S 200 exact draws, 400 replications of 5 test cases: the PIT at each response, and each
        // case's flag read on a grid of step 0.5 over [−25, 25]. A conformal p-value needs no monotone score, so the PIT is uniform
        // over every case and over the unflagged ones. Two chains stuck apart give falling scores: some cases are flagged, with no
        // quantiles or CRPS, and keep their sets.
        const int n = 49, tests = 5, replications = 400;
        double[] grid = [.. Enumerable.Range(0, 101).Select(g => -25.0 + 0.5 * g)];
        var all = new List<double>();
        var unflagged = new List<double>();
        int flagged = 0, apartFlagged = 0, apartCases = 0, scoreFalls = 0;
        int[] seeds = SeedTree.Derive(2026_1113, replications);
        for (int r = 0; r < replications; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            var x = new double[n + tests];
            var y = new double[n + tests];
            for (int i = 0; i < x.Length; i++) { x[i] = 4 * rng.NextDouble() - 2; y[i] = ConformalBayesTests.Truth(x[i], 0.5, rng); }
            double[] beta = ConformalBayesTests.ExactDraws(x[..n], y[..n], 200, rng);
            PointwiseLogDensity training = Training(beta, 1.0, null, x[..n], y[..n]);
            ConformalDistributions atTruth = ConformalBayes.Distribution(training, Seam(beta, 1.0, null, x[n..], [.. y[n..].Select(v => new[] { v })]));
            ConformalDistributions onGrid = ConformalBayes.Distribution(training, Seam(beta, 1.0, null, x[n..], [.. Enumerable.Repeat(grid, tests)]));
            var tau = new Xoshiro256PlusPlus(seeds[r] ^ 0x5bd1e995);
            for (int t = 0; t < tests; t++)
            {
                double pit = atTruth.Distributions[t].At(0, tau.NextDouble());
                all.Add(pit);
                if (onGrid.Distributions[t].Monotone) unflagged.Add(pit);
                else flagged++;
            }
            if (r >= 100) continue;
            double[] apart = [.. beta.Select((b, d) => b + (d % 2 == 0 ? -1.0 : 1.0))];
            ConformalDistributions apartGrid = ConformalBayes.Distribution(Training(apart, 1.0, null, x[..n], y[..n]),
                                                                           Seam(apart, 1.0, null, x[n..], [.. Enumerable.Repeat(grid, tests)]));
            foreach (ConformalDistribution d in apartGrid.Distributions)
            {
                apartCases++;
                if (d.Monotone) continue;
                apartFlagged++;
                // π can fall while the candidate's own score rises, where a unit's score rises faster and overtakes it.
                Assert.True(d.Fall > 0);
                if (d.ScoreFall > 0.0) scoreFalls++;
                Assert.True(double.IsNaN(d.Quantile(0.5)) && double.IsNaN(d.Crps(0.0)));
                Assert.NotEmpty(d.Intervals(0.1));
            }
        }
        output.WriteLine($"the exact posterior: {flagged} of {all.Count} cases flagged; two chains stuck apart: {apartFlagged} of {apartCases}, "
                         + $"{scoreFalls} of them with the candidate's own score falling");
        ReportUniformity("the full step's PIT, every case", all, n + 1, replications, assert: true);
        ReportUniformity("the full step's PIT, unflagged cases", unflagged, n + 1, replications, assert: true);
        Assert.True(apartFlagged > 0);
    }

    // ---- agreement with slice 3's sets, and CRPS against quadrature ------------------------------------------------------------------

    [Fact]
    public void Agreement_TheEquiTailedSetIsTheSplitRulePerTail_AndMeetsSlice3sSetsOnAFitThatIsRight()
    {
        // Data from the model itself, y = 1.3x + N(0, 1), the exact posterior on 30 training units, 1,999 calibration units (m 2,000):
        // the calibration PITs are uniform, so the distribution's 90% set, the density's (slice 3) and CQR's on the fit's 5% and 95%
        // quantiles all tend to the fit's own central 90% interval. Their ends' mean differences over 20 replications of 5 test cases
        // lie within four of their standard errors and a grid step. And exactly: the equi-tailed set is slice 3's rule on the CDF score
        // in each tail at α/2 (SplitCalibration on F̂ and on −F̂), grid value for grid value.
        const int train = 30, calibrationUnits = 1999, tests = 5, replications = 20, s = 200;
        const double step = 0.02;
        double[] grid = [.. Enumerable.Range(0, 1201).Select(g => -12.0 + step * g)];
        var gaps = new Dictionary<string, List<double>> { ["distribution − fit's own, lower"] = [], ["distribution − fit's own, upper"] = [],
                                                          ["density − distribution, lower"] = [], ["density − distribution, upper"] = [],
                                                          ["CQR − distribution, lower"] = [], ["CQR − distribution, upper"] = [] };
        int checkedValues = 0, ambiguous = 0;
        int[] seeds = SeedTree.Derive(2026_1114, replications);
        for (int r = 0; r < replications; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            (double[] x, double[] y) Units(int count)
            {
                var xs = new double[count];
                var ys = new double[count];
                for (int i = 0; i < count; i++) { xs[i] = 4 * rng.NextDouble() - 2; ys[i] = 1.3 * xs[i] + Normal.SampleStandard(rng); }
                return (xs, ys);
            }
            (double[] xTrain, double[] yTrain) = Units(train);
            (double[] xCal, double[] yCal) = Units(calibrationUnits);
            (double[] xTest, _) = Units(tests);
            double[] beta = ConformalBayesTests.ExactDraws(xTrain, yTrain, s, rng);
            CandidateLogDensity calibration = Seam(beta, 1.0, null, xCal, [.. yCal.Select(v => new[] { v })]);
            CandidateLogDensity candidates = Seam(beta, 1.0, null, xTest, [.. Enumerable.Repeat(grid, tests)]);
            ConformalDistributions distributions = SplitConformal.Distribution(calibration, candidates);
            ConformalPrediction density = SplitConformal.Density(calibration, candidates, 0.1);
            double mean = beta.Average(), variance = beta.Sum(b => (b - mean) * (b - mean)) / s;
            QuantileReadings Quantiles(double[] xs) => new([.. xs.Select(v => mean * v - Z95 * Math.Sqrt(1 + v * v * variance))],
                                                           [.. xs.Select(v => mean * v + Z95 * Math.Sqrt(1 + v * v * variance))]);
            ConformalIntervals cqr = SplitConformal.Quantiles(Quantiles(xCal), yCal, Quantiles(xTest), 0.1);
            // The split rule per tail on the CDF score, from sums in the layer's order.
            double Score(double v, double xv) => beta.Sum(b => 1.0 / s * Math.Exp(Normal.StandardLogCdf(v - b * xv)));
            double[] calScores = [.. yCal.Select((v, i) => Score(v, xCal[i]))];
            SplitCalibration lowerTail = SplitConformal.Calibrate(calScores, alpha: 0.05), upperTail = SplitConformal.Calibrate([.. calScores.Select(v => -v)], alpha: 0.05);
            double[] sortedCal = [.. calScores.Order()];
            for (int t = 0; t < tests; t++)
            {
                ConformalDistribution d = distributions.Distributions[t];
                Assert.Equal(2000, d.Units);
                Assert.Equal(100, d.Level(0.1));
                GridInterval set = Assert.Single(d.Intervals(0.1));
                GridInterval byDensity = Assert.Single(density.Sets[t].Intervals);
                double sd = Math.Sqrt(1 + xTest[t] * xTest[t] * variance), centre = mean * xTest[t];
                gaps["distribution − fit's own, lower"].Add(set.Lower - (centre - Z95 * sd));
                gaps["distribution − fit's own, upper"].Add(set.Upper - (centre + Z95 * sd));
                gaps["density − distribution, lower"].Add(byDensity.Lower - set.Lower);
                gaps["density − distribution, upper"].Add(byDensity.Upper - set.Upper);
                gaps["CQR − distribution, lower"].Add(cqr.Intervals[t].Lower - set.Lower);
                gaps["CQR − distribution, upper"].Add(cqr.Intervals[t].Upper - set.Upper);
                for (int g = 0; g < grid.Length; g++)
                {
                    if (Math.Abs(grid[g] - centre) > 6 * sd) continue;   // the bulk, where F̂ as a probability keeps its order
                    double score = Score(grid[g], xTest[t]);
                    int at = Array.BinarySearch(sortedCal, score);
                    bool near = at >= 0 || (~at < sortedCal.Length && sortedCal[~at] - score < 1e-14) || (~at > 0 && score - sortedCal[~at - 1] < 1e-14);
                    if (near) { ambiguous++; continue; }
                    Assert.Equal(lowerTail.Accepts(0, score) && upperTail.Accepts(0, -score), d.Accepts(g, 0.1));
                    checkedValues++;
                }
            }
        }
        output.WriteLine($"the split rule per tail: {checkedValues} grid values in the bulk agree; {ambiguous} within rounding of a calibration score");
        foreach ((string label, List<double> values) in gaps)
        {
            double mean = values.Average(), se = Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / (values.Count - 1) / values.Count);
            output.WriteLine($"  {label}: mean {mean:+0.0000;-0.0000} (se {se:F4}), spread {Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / (values.Count - 1)):F4}");
            Assert.True(Math.Abs(mean) <= 4 * se + step, $"{label}: {mean:F4}");
        }
    }

    [Fact]
    public void Crps_AgainstAQuadratureReference()
    {
        // The split form's CRPS on a grid of step h against the exact distribution's, the step function Q̄ read from the exact F̂ at
        // every point of a grid 16 times finer by the midpoint rule, tails cut at the same ends. Reading Q̄ linearly between grid
        // values moves each of its m − 1 steps of 1/m within a cell, which costs at most h(1 + 1/(3m)) over the whole integral; the
        // typical gap is far smaller, its signs mixed.
        const int train = 30, tests = 3, replications = 5, s = 50;
        const double h = 0.01;
        double[] grid = [.. Enumerable.Range(0, 1601).Select(g => -8.0 + h * g)];
        foreach (int m in new[] { 20, 200 })
        {
            double worst = 0.0, sum = 0.0;
            int count = 0;
            int[] seeds = SeedTree.Derive(2026_1115 + m, replications);
            for (int r = 0; r < replications; r++)
            {
                var rng = new Xoshiro256PlusPlus(seeds[r]);
                (double[] x, double[] y) Units(int k)
                {
                    var xs = new double[k];
                    var ys = new double[k];
                    for (int i = 0; i < k; i++) { xs[i] = 4 * rng.NextDouble() - 2; ys[i] = ConformalBayesTests.Truth(xs[i], 0.5, rng); }
                    return (xs, ys);
                }
                (double[] xTrain, double[] yTrain) = Units(train);
                (double[] xCal, double[] yCal) = Units(m - 1);
                (double[] xTest, double[] yTest) = Units(tests);
                double[] beta = ConformalBayesTests.ExactDraws(xTrain, yTrain, s, rng);
                ConformalDistributions d = SplitConformal.Distribution(Seam(beta, 1.0, null, xCal, [.. yCal.Select(v => new[] { v })]),
                                                                       Seam(beta, 1.0, null, xTest, [.. Enumerable.Repeat(grid, tests)]));
                double Cdf(double v, double xv) => beta.Sum(b => Math.Exp(Normal.StandardLogCdf(v - b * xv))) / s;
                double[] cal = [.. yCal.Select((v, i) => Cdf(v, xCal[i])).Order()];
                for (int t = 0; t < tests; t++)
                {
                    double yt = yTest[t], fine = h / 16, reference = Math.Max(0.0, grid[0] - yt) + Math.Max(0.0, yt - grid[^1]);
                    for (int i = 0; i < 16 * (grid.Length - 1); i++)
                    {
                        double z = grid[0] + (i + 0.5) * fine, f = Cdf(z, xTest[t]);
                        int below = 0, tied = 0;
                        foreach (double c in cal) { if (c < f) below++; else if (c == f) tied++; }
                        double q = (below + 0.5 * (tied + 1)) / m, spread = (tied + 1.0) / m, u = q - (z >= yt ? 1 : 0);
                        reference += fine * (u * u + spread * spread / 12);
                    }
                    double gap = Math.Abs(d.Distributions[t].Crps(yt) - reference);
                    worst = Math.Max(worst, gap);
                    sum += gap;
                    count++;
                    Assert.True(gap <= h * (1 + 1.0 / (3 * m)), $"m {m}: CRPS {d.Distributions[t].Crps(yt):F6} against {reference:F6}");
                }
            }
            output.WriteLine($"m {m}: {count} cases, the grid's CRPS within {worst:E2} of the reference (mean {sum / count:E2}); the bound h(1 + 1/(3m)) {h * (1 + 1.0 / (3 * m)):E2}");
        }
    }

    // ---- the measured default (qualification switch) -----------------------------------------------------------------------------

    [Fact]
    public void MeasuredDefault_OnTheConjugateModel()
    {
        // Both forms side by side on the same data: the split form's fit on 30 training units with 19 calibration units (m 20), the full
        // step's on all 49 (m 50). Four fits: the exact posterior, a wrong fit, a chain stopped early, two chains stuck apart. Per fit:
        // each form's PIT, the full step's PIT over every case and its unflagged ones, its flag rate, the 90% width and CRPS where both
        // give a distribution, and each form's cost per case (the layer alone, the seams built before the clock).
        if (!BassQualification.Enabled) return;
        const int train = 30, calibrationUnits = 19, tests = 10, replications = 500, gridCases = 2, s = 200;
        double[] grid = [.. Enumerable.Range(0, 401).Select(g => -25.0 + 0.125 * g)];
        string[] names = ["the exact posterior", "a wrong fit", "a chain stopped early", "two chains stuck apart"];
        var splitPit = names.Select(_ => new List<double>()).ToArray();
        var fullPit = names.Select(_ => new List<double>()).ToArray();
        var fullUnflagged = names.Select(_ => new List<double>()).ToArray();
        var flags = new int[4];
        var gridCount = new int[4];
        var both = new int[4];
        var splitWidth = new double[4];
        var fullWidth = new double[4];
        var splitCrps = new double[4];
        var fullCrps = new double[4];
        var splitTicks = new long[4];
        var fullTicks = new long[4];
        int[] seeds = SeedTree.Derive(2026_1116, replications);
        for (int r = 0; r < replications; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            (double[] x, double[] y) Units(int count)
            {
                var xs = new double[count];
                var ys = new double[count];
                for (int i = 0; i < count; i++) { xs[i] = 4 * rng.NextDouble() - 2; ys[i] = ConformalBayesTests.Truth(xs[i], 0.5, rng); }
                return (xs, ys);
            }
            (double[] xTrain, double[] yTrain) = Units(train);
            (double[] xCal, double[] yCal) = Units(calibrationUnits);
            (double[] xTest, double[] yTest) = Units(tests);
            double[] xAll = [.. xTrain, .. xCal], yAll = [.. yTrain, .. yCal];
            Fit[] splitFits = Fits(xTrain, yTrain, s, rng, apart: true);
            Fit[] fullFits = Fits(xAll, yAll, s, rng, apart: true);
            var tau = new Xoshiro256PlusPlus(seeds[r] ^ 0x5bd1e995);
            for (int f = 0; f < 4; f++)
            {
                Fit sf = splitFits[f], ff = fullFits[f];
                CandidateLogDensity calibration = Seam(sf.Beta, sf.NoiseSd, sf.Chains, xCal, [.. yCal.Select(v => new[] { v })]);
                PointwiseLogDensity training = Training(ff.Beta, ff.NoiseSd, ff.Chains, xAll, yAll);
                double[][] atTruth = [.. yTest.Select(v => new[] { v })];
                ConformalDistributions splitTruth = SplitConformal.Distribution(calibration, Seam(sf.Beta, sf.NoiseSd, sf.Chains, xTest, atTruth));
                ConformalDistributions fullTruth = ConformalBayes.Distribution(training, Seam(ff.Beta, ff.NoiseSd, ff.Chains, xTest, atTruth));
                CandidateLogDensity splitGrid = Seam(sf.Beta, sf.NoiseSd, sf.Chains, xTest, [.. Enumerable.Repeat(grid, tests)]);
                CandidateLogDensity fullGrid = Seam(ff.Beta, ff.NoiseSd, ff.Chains, xTest, [.. Enumerable.Repeat(grid, tests)]);
                long t0 = Stopwatch.GetTimestamp();
                ConformalDistributions splitOn = r < 100 ? SplitConformal.Distribution(calibration, splitGrid) : splitTruth;
                long t1 = Stopwatch.GetTimestamp();
                ConformalDistributions fullOn = ConformalBayes.Distribution(training, fullGrid);
                long t2 = Stopwatch.GetTimestamp();
                if (r < 100) { splitTicks[f] += t1 - t0; fullTicks[f] += t2 - t1; }
                for (int t = 0; t < tests; t++)
                {
                    splitPit[f].Add(splitTruth.Distributions[t].At(0, tau.NextDouble()));
                    double pit = fullTruth.Distributions[t].At(0, tau.NextDouble());
                    fullPit[f].Add(pit);
                    ConformalDistribution full = fullOn.Distributions[t];
                    gridCount[f]++;
                    if (!full.Monotone) { flags[f]++; continue; }
                    fullUnflagged[f].Add(pit);
                    if (r >= 100 || t >= gridCases) continue;
                    ConformalDistribution split = splitOn.Distributions[t];
                    both[f]++;
                    splitWidth[f] += split.Width(0.1);
                    fullWidth[f] += full.Width(0.1);
                    splitCrps[f] += split.Crps(yTest[t]);
                    fullCrps[f] += full.Crps(yTest[t]);
                }
            }
        }
        double perCase = 1e3 / Stopwatch.Frequency / (100 * tests);
        for (int f = 0; f < 4; f++)
        {
            output.WriteLine($"{names[f]}: the full step flags {flags[f]} of {gridCount[f]} cases ({100.0 * flags[f] / gridCount[f]:F1}%); "
                             + $"where both give a distribution ({both[f]} cases), 90% width {splitWidth[f] / both[f]:F2} split against {fullWidth[f] / both[f]:F2} full, "
                             + $"CRPS {splitCrps[f] / both[f]:F3} against {fullCrps[f] / both[f]:F3}; per case on a grid of {grid.Length}, "
                             + $"{splitTicks[f] * perCase:F3} ms split, {fullTicks[f] * perCase:F3} ms full; the split form's data cost {calibrationUnits} of 49 units");
            ReportUniformity("split PIT", splitPit[f], calibrationUnits + 1, replications, assert: f < 3);
            ReportUniformity("full PIT, every case", fullPit[f], 50, replications, assert: false);
            ReportUniformity("full PIT, unflagged", fullUnflagged[f], 50, replications, assert: false);
        }
    }
}
