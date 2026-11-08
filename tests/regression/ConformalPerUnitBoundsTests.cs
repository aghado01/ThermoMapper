using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// CB-3's per-unit bounds in the conformal layer (E16, the slice between 4 and 5), against <c>lean/Enthymemata/ConformalValidity.lean</c>:
/// the margin on bounds carried by the units (<c>coverage_close</c>, <c>card_fragile_le_estimated</c>) gains where a unit stands
/// apart, a unit in a tail, and falls back to the uniform bound where none does, a Beta–Bernoulli posterior symmetric about ½
/// (<c>fragileEstimated_mono</c>); and, with the qualification switch, how far slice 2's carried bound, read from two orderings,
/// fell short of each unit's error over every ordering the draws represent.
/// </summary>
public sealed class ConformalPerUnitBoundsTests(ITestOutputHelper output)
{
    private static readonly LatentReport[] Theta = [new("θ", LatentScope.Global, LatentTreatment.Conditioned)];

    private static readonly ConformalOptions UniformBounds = new() { MarginBounds = MarginBounds.Uniform, InflationBounds = InflationBounds.Uniform };

    private static readonly ConformalOptions Carried = new() { MarginBounds = MarginBounds.Carried };

    // The conjugate fixture of slice 0: n units of y = βx + noise the model does not share, exact posterior draws of β, and, with a
    // tail, the last unit at x 3.5 with a residual of `tail` beyond the slope 0.7, so its density varies far more over the draws
    // than any other unit's.
    private static (PointwiseLogDensity Training, double[] Beta) Fixture(int n, int s, double tail, Xoshiro256PlusPlus rng)
    {
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++) { x[i] = 4 * rng.NextDouble() - 2; y[i] = ConformalBayesTests.Truth(x[i], 0.5, rng); }
        if (tail > 0) { x[n - 1] = 3.5; y[n - 1] = 0.7 * 3.5 + tail; }
        double[] beta = ConformalBayesTests.ExactDraws(x, y, s, rng);
        var ll = new double[s, n];
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++) ll[d, i] = ConformalBayesTests.LogNormal(y[i], beta[d] * x[i]);
        return (new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), Theta), beta);
    }

    private static CandidateLogDensity Candidates(double[] beta, double xt, double[] grid)
    {
        var lc = new double[beta.Length, grid.Length];
        for (int d = 0; d < beta.Length; d++)
            for (int g = 0; g < grid.Length; g++) lc[d, g] = ConformalBayesTests.LogNormal(grid[g], beta[d] * xt);
        return new CandidateLogDensity(lc, [grid], Theta);
    }

    [Theory]
    [InlineData(40)]
    [InlineData(640)]
    public void TheMargin_GainsWhereAUnitStandsInATail(int s)
    {
        // n 20 (19 units and one in a tail), α 0.1, k 2. Over 40 grids: the carried band never above the uniform one at any grid
        // value, and below it on average; printed, how far the units' bounds over the orderings spread, and where the tail unit's
        // falls among them. Over 2,000 true values: the plain set covers at least the carried margin's reading.
        const int n = 20, sets = 40, truths = 2000;
        double[] grid = [.. Enumerable.Range(0, 161).Select(g => -8.0 + g / 10.0)];
        int[] seeds = SeedTree.Derive(2026_1109 + s, sets + truths);
        double carried = 0.0, uniform = 0.0, tailShare = 0.0, spread = 0.0, underCarried = 0.0, underUniform = 0.0;
        int covered = 0, tailWidest = 0;
        for (int r = 0; r < sets + truths; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            (PointwiseLogDensity training, double[] beta) = Fixture(n, s, 6.0, rng);
            double xt = 4 * rng.NextDouble() - 2, yt = ConformalBayesTests.Truth(xt, 0.5, rng);
            CandidateLogDensity candidates = Candidates(beta, xt, r < sets ? grid : [yt]);
            ConformalSet set = ConformalBayes.FullStep(training, candidates, Carried).Sets[0];
            ConformalSet same = ConformalBayes.FullStep(training, candidates, UniformBounds).Sets[0];
            ConformalMargin m = set.Margin!;
            for (int g = 0; g < m.Fragile!.Length; g++)
                if (m.Fragile[g] >= 0) Assert.True(m.Fragile[g] <= same.Margin!.Fragile![g], "fragileEstimated_mono");
            if (r < sets)
            {
                carried += m.Count!.Value;
                uniform += same.Margin!.Count!.Value;
                double[] units = m.UnitCarriedErrors!;
                double median = units.Order().ElementAt(n / 2);
                tailShare += units[n - 1] / median;
                spread += units.Max() / median;
                if (units[n - 1] == units.Max()) tailWidest++;
                continue;
            }
            if (set.AcceptsPlain(0)) covered++;
            underCarried += set.CoverageUnderError;
            underUniform += same.CoverageUnderError;
        }
        double se = Math.Sqrt(0.09 / truths), rate = covered / (double)truths;
        output.WriteLine($"S {s}, a unit in a tail: m̂ carried {carried / sets:F2}, uniform {uniform / sets:F2} of {n + 1}; the units' bounds over the orderings: the widest "
                         + $"{spread / sets:F1} times the median, the tail unit's {tailShare / sets:F1} times, the widest in {tailWidest} of {sets}; "
                         + $"plain coverage {rate:F4} (se {se:F4}) against the margin's reading {underCarried / truths:F4} carried, {underUniform / truths:F4} uniform");
        Assert.True(carried < uniform, "the carried margin gains");
        Assert.True(rate >= underCarried / truths - 3.5 * se, "coverage_close, carried");
    }

    [Fact]
    public void TheMargin_FallsBackToTheUniformBound_WhereNoUnitStandsApart()
    {
        // The Lean checks' conformal-Bayes case: a Beta–Bernoulli posterior, ten ones and ten zeros under a flat prior, so the
        // posterior Beta(11, 11) is symmetric about ½, drawn in pairs θ and 1 − θ. Each unit's density is θ or 1 − θ, so every
        // unit's error over the orderings is the same, the candidate's too: the carried bounds all equal the uniform one, to
        // rounding, and the bands are equal at both candidates, 0 and 1.
        const int n = 20, s = 4000;
        var rng = new Xoshiro256PlusPlus(2026_1110);
        double[] theta = new double[s];
        for (int d = 0; d < s / 2; d++)
        {
            double a = Gamma.SampleStandard(rng, 11.0), b = Gamma.SampleStandard(rng, 11.0);
            theta[d] = a / (a + b);
            theta[s / 2 + d] = 1.0 - theta[d];
        }
        double Log(double y, double t) => y == 1.0 ? Math.Log(t) : Math.Log(1.0 - t);
        var ll = new double[s, n];
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++) ll[d, i] = Log(i < n / 2 ? 1.0 : 0.0, theta[d]);
        var lc = new double[s, 2];
        for (int d = 0; d < s; d++) { lc[d, 0] = Log(0.0, theta[d]); lc[d, 1] = Log(1.0, theta[d]); }
        var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), Theta);
        var candidates = new CandidateLogDensity(lc, [[0.0, 1.0]], Theta);
        foreach (double alpha in new[] { 0.1, 0.3, 0.5 })
        {
            ConformalSet set = ConformalBayes.FullStep(training, candidates, Carried with { Alpha = alpha }).Sets[0];
            ConformalSet same = ConformalBayes.FullStep(training, candidates, UniformBounds with { Alpha = alpha }).Sets[0];
            ConformalMargin m = set.Margin!;
            for (int g = 0; g < 2; g++)
            {
                double top = m.CarriedUniformErrors![g];
                foreach (double bound in m.UnitCarriedErrors!.Append(m.CandidateCarriedErrors![g]))
                    Assert.True(Math.Abs(bound - top) <= 1e-9 * top, $"a carried bound {bound:R} against the uniform {top:R}");
                Assert.Equal(same.Margin!.Fragile![g], m.Fragile![g]);
            }
            output.WriteLine($"α {alpha}: k {set.Level}, the uniform bound {m.CarriedUniformErrors![0]:F4} and {m.CarriedUniformErrors[1]:F4}, "
                             + $"bands {m.Fragile![0]} and {m.Fragile[1]}, carried and uniform alike");
        }
    }

    [Theory]
    [InlineData(40, 0.0)]
    [InlineData(640, 0.0)]
    [InlineData(40, 6.0)]
    [InlineData(640, 6.0)]
    public void SliceTwosCarriedBound_AgainstEveryOrderingTheDrawsRepresent(int s, double tail)
    {
        // A reading, with the qualification switch: slice 2 carried each unit's bound from two orderings, the held one and its own
        // turn as the candidate. Here every unit's error is read under every unit's weights from public steps, the grid's
        // candidates as extra units, and the layer's carried band is checked against the band those give. Printed: how often and
        // by how much some ordering's error exceeded slice 2's bound, and the bands it gives beside the layer's and the uniform.
        if (!BassQualification.Enabled) return;
        const int n = 20, sets = 20;
        double[] grid = [.. Enumerable.Range(0, 161).Select(g => -8.0 + g / 10.0)];
        int[] seeds = SeedTree.Derive(2026_1108, sets);
        int pairs = 0, shortfalls = 0;
        double sliceTwo = 0.0, layer = 0.0, uniform = 0.0, largest = 1.0;
        var ratios = new List<double>();
        for (int r = 0; r < sets; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            (PointwiseLogDensity training, double[] beta) = Fixture(n, s, tail, rng);
            double xt = 4 * rng.NextDouble() - 2;
            CandidateLogDensity candidates = Candidates(beta, xt, grid);
            ConformalSet set = ConformalBayes.FullStep(training, candidates, Carried).Sets[0];
            ConformalSet same = ConformalBayes.FullStep(training, candidates, UniformBounds).Sets[0];
            double z = set.Margin!.Multiple;
            var aug = new double[s, n + grid.Length];
            for (int d = 0; d < s; d++)
            {
                for (int i = 0; i < n; i++) aug[d, i] = training.Values[d, i];
                for (int g = 0; g < grid.Length; g++) aug[d, n + g] = candidates.Values[d, g];
            }
            var augmented = new PointwiseLogDensity(aug, UnitPartition.EachObservation(n + grid.Length), Theta);
            var worst = new double[n];
            var self = new double[n];
            var candidateWorst = new double[grid.Length];
            for (int i = 0; i < n; i++)
            {
                var column = new double[s, 1];
                for (int d = 0; d < s; d++) column[d, 0] = training.Values[d, i];
                ConformalScores turn = ConformalBayes.Scores(augmented, new CandidateLogDensity(column, [[0.0]], Theta), 0);
                for (int j = 0; j < n; j++) worst[j] = Math.Max(worst[j], turn.Errors[j]);
                for (int g = 0; g < grid.Length; g++) candidateWorst[g] = Math.Max(candidateWorst[g], turn.Errors[n + g]);
                self[i] = turn.CandidateError;
            }
            int mTwo = 0, mAll = 0;
            for (int g = 0; g < grid.Length; g++)
            {
                ConformalScores sc = ConformalBayes.Scores(training, candidates, g);
                double[] all = [.. sc.Scores, sc.Candidate];
                double[] two = [.. Enumerable.Range(0, n).Select(j => z * Math.Max(sc.Errors[j], self[j])), z * sc.CandidateError];
                double[] every = [.. Enumerable.Range(0, n).Select(j => z * Math.Max(sc.Errors[j], worst[j])), z * Math.Max(sc.CandidateError, candidateWorst[g])];
                for (int j = 0; j <= n; j++)
                {
                    pairs++;
                    double ratio = every[j] / two[j];
                    if (ratio > 1.0 + 1e-9) { shortfalls++; ratios.Add(ratio); }
                    largest = Math.Max(largest, ratio);
                }
                mTwo = Math.Max(mTwo, ConformalBayes.EstimatedBand(all, two, set.Level));
                mAll = Math.Max(mAll, ConformalBayes.EstimatedBand(all, every, set.Level));
            }
            Assert.Equal(mAll, set.Margin.Count);
            sliceTwo += mTwo;
            layer += set.Margin.Count!.Value;
            uniform += same.Margin!.Count!.Value;
        }
        double[] q = [.. ratios.Order()];
        output.WriteLine($"S {s}, tail {tail}: some ordering's error exceeds slice 2's carried bound at {shortfalls} of {pairs} (unit, grid value) pairs, "
                         + $"by a median {(q.Length > 0 ? q[q.Length / 2] : 1.0):F2} times, at most {largest:F1}; m̂ slice 2's {sliceTwo / sets:F2}, "
                         + $"carried over every ordering {layer / sets:F2}, uniform {uniform / sets:F2} of {n + 1}");
    }
}
