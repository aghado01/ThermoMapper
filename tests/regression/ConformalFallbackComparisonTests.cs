using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Two readings of the fallback signal side by side (E16, slice 2), with the qualification switch (<c>BASS_QUALIFICATION</c> = 1):
/// the layer's, the gap between the inflated and plain sets (<see cref="ConformalSet.Gap"/>, with k̂'s trust), and a split-half
/// comparator, the main session's construction and not from a paper: the plain set from each half of the draws, and the mass,
/// under the full draws' predictive at x, of the region where the two half-sets disagree. Each case reads both signals on a grid
/// at x and the plain set's decision at the true value at the same x, so a rule's firing is read against the coverage the plain
/// set achieved: a useful rule fires where the plain set falls short and stays quiet where it holds. Fixtures: exact draws at S 40
/// to 2,560; slowly mixing chains, where the plain set may fall short; slice 1's Gaussian and t BASS fits, with the degenerate
/// Gaussian weights read apart.
/// </summary>
public sealed class ConformalFallbackComparisonTests(ITestOutputHelper output)
{
    private static readonly LatentReport[] Theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];

    private sealed record Case(bool Covered, bool Fallback, double Gap, double Mass, double Ess);

    // The two halves of the draws: each chain's first and second halves where there are chains, alternate draws otherwise.
    private static (int[] A, int[] B, int[]? ChainsA, int[]? ChainsB) Halves(int draws, IReadOnlyList<int>? chains)
    {
        if (chains is null)
            return ([.. Enumerable.Range(0, draws).Where(d => d % 2 == 0)], [.. Enumerable.Range(0, draws).Where(d => d % 2 == 1)], null, null);
        var a = new List<int>();
        var b = new List<int>();
        var la = new int[chains.Count];
        var lb = new int[chains.Count];
        int start = 0;
        for (int c = 0; c < chains.Count; c++)
        {
            int half = chains[c] / 2;
            a.AddRange(Enumerable.Range(start, half));
            b.AddRange(Enumerable.Range(start + half, half));
            la[c] = lb[c] = half;
            start += chains[c];
        }
        return ([.. a], [.. b], la, lb);
    }

    private static double[,] Rows(double[,] values, int[] rows)
    {
        var sub = new double[rows.Length, values.GetLength(1)];
        for (int d = 0; d < rows.Length; d++)
            for (int j = 0; j < values.GetLength(1); j++) sub[d, j] = values[rows[d], j];
        return sub;
    }

    /// <summary>The split-half reading of each case: the predictive mass where the half-draws' plain sets disagree.</summary>
    internal static double[] SplitHalfMass(PointwiseLogDensity training, CandidateLogDensity candidates, ConformalOptions? options = null)
    {
        (int[] a, int[] b, int[]? ca, int[]? cb) = Halves(training.Draws, training.ChainLengths);
        double[]? w = training.Weights is null ? null : [.. training.Weights];
        PointwiseLogDensity Half(int[] rows, int[]? chains)
            => new PointwiseLogDensity(Rows(training.Values, rows), training.Units, training.Latents, chains)
            { Weights = w is null ? null : [.. rows.Select(d => w[d])] };
        CandidateLogDensity HalfC(int[] rows, int[]? chains)
            => new CandidateLogDensity(Rows(candidates.Values, rows), candidates.Grids, candidates.Latents, candidates.StratumOf, chains)
            { Weights = w is null ? null : [.. rows.Select(d => w[d])] };
        ConformalPrediction pa = ConformalBayes.FullStep(Half(a, ca), HalfC(a, ca), options);
        ConformalPrediction pb = ConformalBayes.FullStep(Half(b, cb), HalfC(b, cb), options);
        var mass = new double[candidates.Cases];
        int s = candidates.Draws;
        double logTotal = w is null ? Math.Log(s) : Math.Log(w.Sum());
        for (int t = 0; t < candidates.Cases; t++)
        {
            ConformalSet sa = pa.Sets[t], sb = pb.Sets[t];
            if (sa.WholeSpace) continue;
            double[] grid = candidates.Grids[t];
            int off = candidates.Offset(t), n = grid.Length;
            for (int g = 0; g < n; g++)
            {
                if (sa.AcceptsPlain(g) == sb.AcceptsPlain(g)) continue;
                // The full draws' predictive density at y_g, times the value's cell.
                double max = double.NegativeInfinity;
                for (int d = 0; d < s; d++) max = Math.Max(max, (w is null ? 0.0 : Math.Log(w[d])) + candidates.Values[d, off + g]);
                double sum = 0.0;
                for (int d = 0; d < s; d++) sum += Math.Exp((w is null ? 0.0 : Math.Log(w[d])) + candidates.Values[d, off + g] - max);
                double density = Math.Exp(max + Math.Log(sum) - logTotal);
                double cell = n == 1 ? 0.0 : 0.5 * ((g > 0 ? grid[g] - grid[g - 1] : grid[1] - grid[0]) + (g < n - 1 ? grid[g + 1] - grid[g] : grid[n - 1] - grid[n - 2]));
                mass[t] += density * cell;
            }
        }
        return mass;
    }

    // Each rule's firing rate, and the plain set's coverage where it fires and where it does not.
    private string Summary(string label, IReadOnlyList<Case> cases, double alpha = 0.1)
    {
        string Rule(string name, Func<Case, bool> fires)
        {
            Case[] on = [.. cases.Where(fires)], off = [.. cases.Where(c => !fires(c))];
            string Cov(Case[] set) => set.Length == 0 ? "–" : $"{set.Count(c => c.Covered) / (double)set.Length:F3} of {set.Length}";
            return $"{name} fires {on.Length / (double)cases.Count:F3} (covered where it fires {Cov(on)}, where quiet {Cov(off)})";
        }
        double[] mass = [.. cases.Select(c => c.Mass).Order()];
        return $"{label}: plain coverage {cases.Count(c => c.Covered) / (double)cases.Count:F4} over {cases.Count}; "
               + Rule("the layer's signal", c => c.Fallback) + "; " + Rule("the gap alone above 0.5", c => c.Gap > 0.5) + "; "
               + Rule("split-half mass above 0.1α", c => c.Mass > 0.1 * alpha) + "; " + Rule("above 0.25α", c => c.Mass > 0.25 * alpha) + "; "
               + $"the mass median {mass[mass.Length / 2]:F4}, 90th percentile {mass[(int)(0.9 * (mass.Length - 1))]:F4}";
    }

    // Slice 0's conjugate fixture: n 19 under noise the model does not share, a grid and the true value at one x, the draws of β
    // given by `draws` (exact, or a chain).
    private static Case ConjugateCase(Xoshiro256PlusPlus rng, Func<double, double, double[]> draws, int[]? chains)
    {
        const int n = 19;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++) { x[i] = 4 * rng.NextDouble() - 2; y[i] = ConformalBayesTests.Truth(x[i], 0.5, rng); }
        double xt = 4 * rng.NextDouble() - 2, yt = ConformalBayesTests.Truth(xt, 0.5, rng);
        double precision = 0.1 + x.Sum(v => v * v), center = x.Zip(y, (a, b) => a * b).Sum() / precision;
        double[] beta = draws(center, 1.0 / Math.Sqrt(precision));
        int s = beta.Length;
        var ll = new double[s, n];
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++) ll[d, i] = ConformalBayesTests.LogNormal(y[i], beta[d] * x[i]);
        double[] grid = [.. Enumerable.Range(0, 201).Select(g => -10.0 + g / 10.0)];
        var lc = new double[s, grid.Length + 1];
        for (int d = 0; d < s; d++)
        {
            for (int g = 0; g < grid.Length; g++) lc[d, g] = ConformalBayesTests.LogNormal(grid[g], beta[d] * xt);
            lc[d, grid.Length] = ConformalBayesTests.LogNormal(yt, beta[d] * xt);
        }
        var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), Theta, chains);
        var candidates = new CandidateLogDensity(lc, [grid, [yt]], Theta, chainLengths: chains);
        ConformalPrediction p = ConformalBayes.FullStep(training, candidates);
        double[] mass = SplitHalfMass(training, candidates);
        ConformalSet set = p.Sets[0], truth = p.Sets[1];
        return new Case(truth.AcceptsPlain(0), set.FallBack, set.Gap, mass[0], truth.EffectiveSampleSizes![0]);
    }

    [Fact]
    public void ExactDraws_AndSlowChains_BothRulesAgainstThePlainSetsCoverage()
    {
        if (!BassQualification.Enabled) return;
        const int sets = 1000;
        foreach (int s in new[] { 40, 160, 640, 2560 })
        {
            int[] seeds = SeedTree.Derive(2026_1110 + s, sets);
            var cases = new List<Case>();
            for (int r = 0; r < sets; r++)
            {
                var rng = new Xoshiro256PlusPlus(seeds[r]);
                cases.Add(ConjugateCase(rng, (c, sd) => [.. Enumerable.Range(0, s).Select(_ => c + sd * Normal.SampleStandard(rng))], null));
            }
            output.WriteLine(Summary($"exact draws, S {s}", cases));
        }
        // Two chains of 500 moving as an AR(1) about the exact posterior: they mix more slowly as ρ grows, so the draws represent the
        // posterior less well, which is where the plain set can fall short.
        foreach (double rho in new[] { 0.9, 0.99, 0.995, 0.999 })
        {
            int[] seeds = SeedTree.Derive(2026_1111 + (int)(1000 * rho), sets);
            var cases = new List<Case>();
            for (int r = 0; r < sets; r++)
            {
                var rng = new Xoshiro256PlusPlus(seeds[r]);
                cases.Add(ConjugateCase(rng, (c, sd) =>
                {
                    var chain = new double[1000];
                    for (int k = 0; k < 2; k++)
                    {
                        double z = Normal.SampleStandard(rng);
                        for (int d = 0; d < 500; d++)
                        {
                            z = rho * z + Math.Sqrt(1 - rho * rho) * Normal.SampleStandard(rng);
                            chain[k * 500 + d] = c + sd * z;
                        }
                    }
                    return chain;
                }, [500, 500]));
            }
            output.WriteLine(Summary($"two chains, ρ {rho}", cases));
        }
    }

    [Theory]
    [InlineData("gaussian")]
    [InlineData("t")]
    public void BassFits_BothRulesAgainstThePlainSetsCoverage(string law)
    {
        // Slice 1's fixture: 20 fits of n 99 (k 10 of m 100, S 1,000 in two chains), every held-out case read on a grid and at its
        // true value; the cases whose weights at the true value have an ESS below 50 are summarized apart.
        if (!BassQualification.Enabled) return;
        const int n = 99, tests = 50, fits = 20;
        double[] grid = [.. Enumerable.Range(0, 101).Select(g => -3.0 + 0.08 * g)];
        int[] seeds = SeedTree.Derive(2026_1107, fits);
        var cases = new List<Case>();
        for (int r = 0; r < fits; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            (double[] x0, double[] x1, double[] y) = BassCandidateTests.Draw(n, rng);
            (double[] t0, double[] t1, double[] yt) = BassCandidateTests.Draw(tests, rng);
            var data = new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
            BassFit fit = Bass.Fit(BassStudentTTests.Config(BassCandidateTests.Noise(law), burnIn: 1000, samples: 500), data);
            PredictionCases both = PredictionCases.Of([new ContinuousPredictor("x0", [.. t0, .. t0]), new ContinuousPredictor("x1", [.. t1, .. t1])]);
            double[][] grids = [.. Enumerable.Range(0, tests).Select(_ => grid), .. yt.Select(v => new[] { v })];
            PointwiseLogDensity training = fit.PointwiseLogDensity(data);
            CandidateLogDensity candidates = fit.CandidateLogDensity(data, both, grids);
            ConformalPrediction p = ConformalBayes.FullStep(training, candidates);
            double[] mass = SplitHalfMass(training, candidates);
            for (int t = 0; t < tests; t++)
            {
                ConformalSet set = p.Sets[t], truth = p.Sets[tests + t];
                cases.Add(new Case(truth.AcceptsPlain(0), set.FallBack, set.Gap, mass[t], truth.EffectiveSampleSizes![0]));
            }
        }
        output.WriteLine(Summary($"{law}, all cases", cases));
        Case[] degenerate = [.. cases.Where(c => c.Ess < 50)];
        output.WriteLine($"{law}, ESS below 50 at the true value ({degenerate.Length}): covered {degenerate.Count(c => c.Covered)}, "
                         + $"the layer's signal {degenerate.Count(c => c.Fallback)}, the gap above 0.5 {degenerate.Count(c => c.Gap > 0.5)}, "
                         + $"split-half mass above 0.1α {degenerate.Count(c => c.Mass > 0.01)}; masses {string.Join(", ", degenerate.Select(c => c.Mass.ToString("F4")))}");
        output.WriteLine(Summary($"{law}, ESS 50 or more", [.. cases.Where(c => c.Ess >= 50)]));
    }
}
