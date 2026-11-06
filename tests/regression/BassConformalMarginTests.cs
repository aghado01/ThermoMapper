using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// What CB-3's reading says where it matters (E16, slice 2), with the qualification switch (<c>BASS_QUALIFICATION</c> = 1): the
/// Gaussian BASS fits of slice 1's held-out check, whose add-one-in weights fall to an ESS of 1.5 at outlying candidates, beside
/// Student t fits of the same data; and how the margin tightens as the draws grow, on exact draws.
/// </summary>
public sealed class BassConformalMarginTests(ITestOutputHelper output)
{
    private sealed record Candidate(double Ess, double Error, double Uniform, int Fragile, int Carried, double ParetoK, double PValue, bool Plain, bool Inflated,
                                    bool FallBack, bool CarriedFallBack);

    [Theory]
    [InlineData("gaussian")]
    [InlineData("t")]
    public void TheDegenerateWeights_WhatTheLayerSays(string law)
    {
        // Slice 1's fixture, fit for fit: n 99, k 10 of m 100, twenty fits, fifty held-out cases each, and the width read on twenty.
        // At each true value: its weights' ESS, its bound, the uniform bound, the band, k̂ where the inflated rule accepts, both
        // rules' decisions and the case's fallback signal.
        if (!BassQualification.Enabled) return;
        const int n = 99, tests = 50, widths = 20, fits = 20;
        double[] grid = [.. Enumerable.Range(0, 101).Select(g => -3.0 + 0.08 * g)];
        int[] seeds = SeedTree.Derive(2026_1107, fits);
        var at = new List<Candidate>();
        double plainWidth = 0.0, inflatedWidth = 0.0, gridMargin = 0.0, gridUnder = 0.0, carriedMargin = 0.0, carriedUnder = 0.0;
        int gridFallBacks = 0, untrustedGrids = 0, carriedFallBacks = 0;
        for (int r = 0; r < fits; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            (double[] x0, double[] x1, double[] y) = BassCandidateTests.Draw(n, rng);
            (double[] t0, double[] t1, double[] yt) = BassCandidateTests.Draw(tests, rng);
            var data = new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
            BassFit fit = Bass.Fit(BassStudentTTests.Config(BassCandidateTests.Noise(law), burnIn: 1000, samples: 500), data);
            PredictionCases cases = PredictionCases.Of([new ContinuousPredictor("x0", [.. t0, .. t0[..widths]]), new ContinuousPredictor("x1", [.. t1, .. t1[..widths]])]);
            double[][] grids = [.. yt.Select(v => new[] { v }), .. Enumerable.Range(0, widths).Select(_ => grid)];
            PointwiseLogDensity pointwise = fit.PointwiseLogDensity(data);
            CandidateLogDensity candidates = fit.CandidateLogDensity(data, cases, grids);
            ConformalPrediction p = ConformalBayes.FullStep(pointwise, candidates);
            ConformalPrediction pc = ConformalBayes.FullStep(pointwise, candidates, new ConformalOptions { MarginBounds = MarginBounds.Carried });
            for (int t = 0; t < tests; t++)
            {
                ConformalSet set = p.Sets[t];
                ConformalMargin m = set.Margin!;
                at.Add(new Candidate(set.EffectiveSampleSizes![0], m.CandidateErrors[0], m.UniformErrors[0], m.Fragile[0], pc.Sets[t].Margin!.Fragile[0],
                                     m.ParetoK[0], set.PValues![0], set.Accepts(0), (set with { Inflated = true }).Accepts(0), set.FallBack, pc.Sets[t].FallBack));
            }
            for (int t = 0; t < widths; t++)
            {
                ConformalSet wide = p.Sets[tests + t];
                plainWidth += wide.Width;
                inflatedWidth += (wide with { Inflated = true }).Width;
                gridMargin += wide.Margin!.Count;
                gridUnder += wide.CoverageUnderError;
                if (wide.FallBack) gridFallBacks++;
                if (wide.Margin.Untrusted(wide.Level) > 0) untrustedGrids++;
                ConformalSet carried = pc.Sets[tests + t];
                carriedMargin += carried.Margin!.Count;
                carriedUnder += carried.CoverageUnderError;
                if (carried.FallBack) carriedFallBacks++;
            }
        }
        int count = at.Count, gridCount = fits * widths;
        void Line(string label, IReadOnlyList<Candidate> c)
        {
            if (c.Count == 0) { output.WriteLine($"{law}, {label}: none"); return; }
            double[] khat = [.. c.Where(x => !double.IsNaN(x.ParetoK)).Select(x => x.ParetoK)];
            output.WriteLine($"{law}, {label} ({c.Count}): ESS median {Median(c.Select(x => x.Ess))}, its bound median {Median(c.Select(x => x.Error))}, "
                             + $"the uniform bound median {Median(c.Select(x => x.Uniform))}, the band m̂ median {Median(c.Select(x => (double)x.Fragile))} (largest {c.Max(x => x.Fragile)}), "
                             + $"carried {Median(c.Select(x => (double)x.Carried))} (largest {c.Max(x => x.Carried)}); "
                             + $"plain accepts {c.Count(x => x.Plain)}, inflated {c.Count(x => x.Inflated)}; k̂ read at {khat.Length}, above the threshold at {khat.Count(v => v > PredictiveChecks.ParetoKThreshold(1000))}; "
                             + $"fallback at {c.Count(x => x.FallBack)}, carried {c.Count(x => x.CarriedFallBack)}");
        }
        output.WriteLine($"{law}: coverage at the true values, plain {at.Count(x => x.Plain) / (double)count:F4}, inflated {at.Count(x => x.Inflated) / (double)count:F4} (1 − k/m 0.900)");
        Line("ESS below 50", [.. at.Where(x => x.Ess < 50)]);
        Line("ESS 50 or more", [.. at.Where(x => x.Ess >= 50)]);
        foreach (Candidate c in at.OrderBy(x => x.Ess).Take(8))
            output.WriteLine($"  ESS {c.Ess:F1}: π {c.PValue:F2}, bound {c.Error:F3}, uniform {c.Uniform:F3}, m̂ {c.Fragile}, carried {c.Carried}, k̂ {c.ParetoK:F2}, plain {c.Plain}, inflated {c.Inflated}, fallback {c.FallBack}, carried {c.CarriedFallBack}");
        output.WriteLine($"{law}, the width grids: plain width {plainWidth / gridCount:F3}, inflated {inflatedWidth / gridCount:F3}; m̂ over the grid averaged {gridMargin / gridCount:F1} of 100, "
                         + $"coverage under error averaged {gridUnder / gridCount:F3}; fallback on {gridFallBacks} of {gridCount}, untrusted candidates in the inflated set on {untrustedGrids}; "
                         + $"carried: m̂ {carriedMargin / gridCount:F1}, coverage under error {carriedUnder / gridCount:F3}, fallback on {carriedFallBacks}");
        Assert.All(at.Where(x => x.Ess >= 1.0), x => Assert.True(double.IsFinite(x.Error) && x.Fragile >= 0));
    }

    private static string Median(IEnumerable<double> values)
    {
        double[] v = [.. values.Order()];
        return v.Length == 0 ? "–" : v[v.Length / 2].ToString("G3");
    }

    [Fact]
    public void TheMargin_TightensAsTheDrawsGrow()
    {
        // Slice 0's conjugate fixture on exact draws (n 19, k 2 of m 20): the band over a grid, the fallback signal and both rules'
        // widths at S 40, 160, 640 and 2,560, on the same data sets.
        if (!BassQualification.Enabled) return;
        const int n = 19, sets = 100;
        LatentReport[] theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];
        double[] grid = [.. Enumerable.Range(0, 401).Select(g => -10.0 + g / 20.0)];
        int[] seeds = SeedTree.Derive(2026_1109, sets);
        foreach (int s in new[] { 40, 160, 640, 2560 })
        {
            double margin = 0.0, under = 0.0, plain = 0.0, inflated = 0.0, uniform = 0.0, carried = 0.0, carriedUnder = 0.0;
            int fallBacks = 0, carriedFallBacks = 0;
            for (int r = 0; r < sets; r++)
            {
                var rng = new Xoshiro256PlusPlus(seeds[r]);
                var x = new double[n];
                var y = new double[n];
                for (int i = 0; i < n; i++) { x[i] = 4 * rng.NextDouble() - 2; y[i] = ConformalBayesTests.Truth(x[i], 0.5, rng); }
                double xt = 4 * rng.NextDouble() - 2;
                double[] beta = ConformalBayesTests.ExactDraws(x, y, s, rng);
                var ll = new double[s, n];
                for (int d = 0; d < s; d++)
                    for (int i = 0; i < n; i++) ll[d, i] = ConformalBayesTests.LogNormal(y[i], beta[d] * x[i]);
                var lc = new double[s, grid.Length];
                for (int d = 0; d < s; d++)
                    for (int g = 0; g < grid.Length; g++) lc[d, g] = ConformalBayesTests.LogNormal(grid[g], beta[d] * xt);
                var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), theta);
                var candidates = new CandidateLogDensity(lc, [grid], theta);
                ConformalSet set = ConformalBayes.FullStep(training, candidates).Sets[0];
                ConformalSet own = ConformalBayes.FullStep(training, candidates, new ConformalOptions { MarginBounds = MarginBounds.Carried }).Sets[0];
                carried += own.Margin!.Count;
                carriedUnder += own.CoverageUnderError;
                if (own.FallBack) carriedFallBacks++;
                margin += set.Margin!.Count;
                under += set.CoverageUnderError;
                plain += set.Width;
                inflated += (set with { Inflated = true }).Width;
                double accepted = 0.0;
                int within = 0;
                for (int g = 0; g < grid.Length; g++)
                    if (set.Accepts(g)) { accepted += set.Margin.UniformErrors[g]; within++; }
                uniform += within > 0 ? accepted / within : 0.0;
                if (set.FallBack) fallBacks++;
            }
            output.WriteLine($"S {s}: m̂ {margin / sets:F2} of 20, coverage under error {under / sets:F3}, fallback {fallBacks} of {sets}; "
                             + $"the uniform bound in the set {uniform / sets:F3}; widths plain {plain / sets:F3}, inflated {inflated / sets:F3}; "
                             + $"carried: m̂ {carried / sets:F2}, coverage under error {carriedUnder / sets:F3}, fallback {carriedFallBacks}");
        }
    }
}
