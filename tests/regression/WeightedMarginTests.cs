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
/// CB-3's margin and inflated rule under shift weights (E16, slice 14's follow-up 3): the weighted band is the Lean's
/// <c>WFragileEstimated</c> share on a seeded sweep, and the count's band at equal weights; the full step's weighted sets carry the band's
/// share and the weighted inflated rule's p-values, reading 1 − α − M̂ − ν and 1 − α − ν (<c>coverage_close_weighted</c>,
/// <c>coverage_inflated_weighted</c>). With the qualification switch (<c>BASS_QUALIFICATION</c> = 1): coverage under covariate shift
/// from few draws.
/// </summary>
public sealed class WeightedMarginTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    [Fact]
    public void TheWeightedBand_IsTheLeansDefinition_OnASweep()
    {
        // 3,000 cases: 2 to 9 units (the last the candidate), scores on a lattice so they tie, bounds of 0 to 0.3 with zeros, weights
        // with zeros, α in (0.05, 0.5). The band's share is Σ q_u over the units the weighted inflated rule accepts, more than α of the
        // weight with s_j − ε_j ≤ s_u + ε_u, with at most α at j = u or s_j + 2ε_j ≤ s_u − 2ε_u. At equal weights and α between k/m and
        // (k + 1)/m it is the count's band over m.
        var rng = new Xoshiro256PlusPlus(1601);
        int nonzero = 0;
        for (int r = 0; r < 3000; r++)
        {
            int m = 2 + rng.NextInt(8);
            double alpha = 0.05 + 0.45 * rng.NextDouble();
            double[] s = [.. Enumerable.Range(0, m).Select(_ => 0.1 * rng.NextInt(12))];
            double[] e = [.. Enumerable.Range(0, m).Select(_ => rng.NextDouble() < 0.3 ? 0.0 : 0.05 * rng.NextInt(7))];
            double[] w = [.. Enumerable.Range(0, m).Select(_ => rng.NextDouble() < 0.15 ? 0.0 : rng.NextDouble())];
            double total = w.Sum(), expected = 0.0;
            if (total > 0.0)
                for (int u = 0; u < m; u++)
                {
                    double accepted = 0.0, shift = 0.0;
                    for (int j = 0; j < m; j++)
                    {
                        if (s[j] - e[j] <= s[u] + e[u]) accepted += w[j] / total;
                        if (j == u || s[j] + 2 * e[j] <= s[u] - 2 * e[u]) shift += w[j] / total;
                    }
                    if (accepted > alpha && shift <= alpha) expected += w[u] / total;
                }
            double band = ConformalBayes.WeightedBand(s, e, w[..^1], w[^1], alpha);
            if (band > 0.0) nonzero++;
            Assert.True(Math.Abs(band - expected) < 1e-12, $"case {r}: the band {band}, the definition's {expected}");

            int k = rng.NextInt(m);
            double share = ConformalBayes.WeightedBand(s, e, [.. Enumerable.Repeat(0.4, m - 1)], 0.4, (k + 0.5) / m);
            Assert.Equal(ConformalBayes.EstimatedBand(s, e, k), (int)Math.Round(share * m));
        }
        output.WriteLine($"3,000 cases, {nonzero} with a band");
    }

    // The conjugate model's terms, y = βx + N(0, 1).
    private static (PointwiseLogDensity Training, CandidateLogDensity Candidates) Conjugate(double[] x, double[] y, double[] beta, double[] cases, double[] grid)
    {
        int s = beta.Length, n = x.Length;
        var tv = new double[s, n];
        var cv = new double[s, cases.Length * grid.Length];
        for (int d = 0; d < s; d++)
        {
            for (int i = 0; i < n; i++) tv[d, i] = ConformalBayesTests.LogNormal(y[i], beta[d] * x[i]);
            for (int t = 0; t < cases.Length; t++)
                for (int g = 0; g < grid.Length; g++) cv[d, t * grid.Length + g] = ConformalBayesTests.LogNormal(grid[g], beta[d] * cases[t]);
        }
        return (new PointwiseLogDensity(tv, UnitPartition.EachObservation(n), []),
                new CandidateLogDensity(cv, [.. Enumerable.Repeat(grid, cases.Length)], []));
    }

    [Fact]
    public void TheFullStepsWeightedSets_CarryTheirMargin()
    {
        var rng = new Xoshiro256PlusPlus(1602);
        double[] x = [.. Enumerable.Range(0, 29).Select(_ => Gaussian(rng))];
        double[] y = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
        double[] grid = [.. Enumerable.Range(0, 41).Select(g => -4.0 + 0.2 * g)];
        (PointwiseLogDensity training, CandidateLogDensity candidates) = Conjugate(x, y, ConformalBayesTests.ExactDraws(x, y, 40, rng), [0.4, 1.2], grid);
        var weights = ShiftWeights.ByCase([.. x.Select(v => Math.Exp(v - 0.5))], [Math.Exp(-0.1), Math.Exp(0.7)], candidates.Grids);
        foreach (MarginBounds bounds in new[] { MarginBounds.Uniform, MarginBounds.Carried })
        {
            ConformalPrediction plain = ConformalBayes.FullStep(training, candidates, new ConformalOptions { Weights = weights, MarginBounds = bounds });
            ConformalPrediction inflated = ConformalBayes.FullStep(training, candidates, new ConformalOptions { Weights = weights, MarginBounds = bounds, Inflate = true });
            for (int t = 0; t < 2; t++)
            {
                ConformalSet set = plain.Sets[t];
                ConformalMargin margin = set.Margin!;
                Assert.Equal(0.1, set.WeightedAlpha);
                Assert.NotNull(margin.WeightedFragile);
                double band = margin.WeightedCount!.Value;
                Assert.Equal(0.9 - band - set.ErrorLevel, set.CoverageUnderError, 12);
                Assert.Equal(0.9 - inflated.Sets[t].ErrorLevel, inflated.Sets[t].CoverageUnderError, 12);
                for (int g = 0; g < grid.Length; g++)
                {
                    Assert.True(margin.WeightedInflated![g] >= set.PValues![g], "the inflated rule counts every unit the plain one does");
                    if (set.AcceptsPlain(g)) Assert.True(set.AcceptsInflated(g));
                    Assert.Equal(set.AcceptsInflated(g), inflated.Sets[t].Accepts(g));
                    Assert.InRange(margin.WeightedFragile![g], 0.0, 1.0);
                }
                output.WriteLine($"{bounds}, case {t}: M̂ {band:F4}, the plain set's coverage under error {set.CoverageUnderError:F4}, its gap {set.Gap:F3}");
            }
        }

        // Equal weights: CB-1's set and margin, counts and all.
        var flat = ShiftWeights.ByCase([.. x.Select(_ => 0.3)], [0.3, 0.3], candidates.Grids);
        ConformalPrediction counts = ConformalBayes.FullStep(training, candidates, new ConformalOptions { MarginBounds = MarginBounds.Carried });
        ConformalPrediction equal = ConformalBayes.FullStep(training, candidates, new ConformalOptions { Weights = flat, MarginBounds = MarginBounds.Carried });
        for (int t = 0; t < 2; t++)
        {
            Assert.Null(equal.Sets[t].WeightedAlpha);
            Assert.Equal(counts.Sets[t].Margin!.Fragile, equal.Sets[t].Margin!.Fragile);
            Assert.Equal(counts.Sets[t].CoverageUnderError, equal.Sets[t].CoverageUnderError);
        }
    }

    // ---- Behind the qualification switch ----

    [Fact]
    public void Measurement_CoverageUnderShift_FromFewDraws()
    {
        // n 49, 40 exact draws of the conjugate model under the conformal tests' heteroscedastic truth, training x ~ N(0, 1) and test x ~
        // N(1, 1), w(x) = e^(x − 1/2); 2,000 replications at the truth, the carried margin: the plain weighted set covers at least
        // 1 − α − E[M̂] − ν and the inflated one at least 1 − α − ν, up to the sampling error.
        if (!BassQualification.Enabled) return;
        const int n = 49, replications = 2000;
        var rng = new Xoshiro256PlusPlus(2026_1125);
        int plain = 0, inflated = 0;
        double bands = 0.0, nu = 0.0;
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))];
            double[] y = [.. x.Select(v => ConformalBayesTests.Truth(v, 1.0, rng))];
            double xt = 1.0 + Gaussian(rng), yt = ConformalBayesTests.Truth(xt, 1.0, rng);
            (PointwiseLogDensity training, CandidateLogDensity candidates) = Conjugate(x, y, ConformalBayesTests.ExactDraws(x, y, 40, rng), [xt], [yt]);
            var weights = new ShiftWeights([.. x.Select(v => Math.Exp(v - 0.5))], [Math.Exp(xt - 0.5)]);
            ConformalSet set = ConformalBayes.FullStep(training, candidates, new ConformalOptions { Weights = weights, MarginBounds = MarginBounds.Carried }).Sets[0];
            if (set.AcceptsPlain(0)) plain++;
            if (set.AcceptsInflated(0)) inflated++;
            bands += set.Margin!.WeightedCount ?? 0.0;
            nu = set.ErrorLevel;
        }
        double se = Math.Sqrt(0.09 / replications), p = plain / (double)replications, q = inflated / (double)replications, band = bands / replications;
        output.WriteLine($"from 40 draws under covariate shift: the plain weighted set covers {p:F4} (mean M̂ {band:F4}, ν {nu:F3}), the inflated {q:F4}; se {se:F4}");
        Assert.True(p >= 0.9 - band - nu - 3.5 * se);
        Assert.True(q >= 0.9 - nu - 3.5 * se);
    }
}
