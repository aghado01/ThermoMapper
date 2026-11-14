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
/// The full step's predictive distribution under shift weights with its Monte Carlo error readouts (E16, slice 14's follow-up 3b): the
/// weighted band at α and the weighted inflated shares are the counts' at equal weights; a weighted case reads 1 − α − M̂ − ν for the plain
/// set and 1 − α − ν for the inflated one, keeps every value the plain set does, and reads no PIT counts.
/// </summary>
public sealed class WeightedCdfMarginTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    [Fact]
    public void TheWeightedBandAndShares_AreTheCounts_AtEqualWeights()
    {
        // 2,000 orderings of 2 to 9 CDF scores on a lattice, so they tie, under bounds of 0 to 0.06 with zeros. At equal weights the
        // weighted band at α/2 = (k + ½)/m, times m, is the counts' band at level k, and the inflated shares times m the inflated ranks.
        var rng = new Xoshiro256PlusPlus(1701);
        for (int r = 0; r < 2000; r++)
        {
            int m = 2 + rng.NextInt(8);
            CdfScore[] scores = [.. Enumerable.Range(0, m).Select(_ => { double f = 0.05 * (1 + rng.NextInt(19)); return new CdfScore(f, 1.0 - f); })];
            double[] bounds = [.. Enumerable.Range(0, m).Select(_ => rng.NextDouble() < 0.3 ? 0.0 : 0.01 * rng.NextInt(7))];
            double[] equal = [.. Enumerable.Repeat(1.0, m)];
            CdfMarginReading counts = CdfMarginCounts.Read(scores, bounds, bounds, false, null, MarginBounds.Carried, InflationBounds.PerPosition);
            CdfMarginReading weighed = CdfMarginCounts.Read(scores, bounds, bounds, false, null, MarginBounds.Carried, InflationBounds.PerPosition, equal);
            Assert.Equal(counts.InflatedLowerRank, (int)Math.Round(weighed.WeightedLower * m));
            Assert.Equal(counts.InflatedUpperRank, (int)Math.Round(weighed.WeightedUpper * m));
            Assert.Null(weighed.Fragile);
            for (int k = 0; k < m; k++)
                Assert.Equal(counts.Fragile![k], (int)Math.Round(CdfMarginCounts.WeightedBand(scores, bounds, equal, 2.0 * (k + 0.5) / m) * m));
        }
    }

    [Fact]
    public void TheFullStepsDistribution_CarriesItsWeightedMargin()
    {
        var rng = new Xoshiro256PlusPlus(1702);
        double[] beta = [.. Enumerable.Range(0, 40).Select(_ => 0.5 + 0.15 * Gaussian(rng))];
        double[] x = [.. Enumerable.Range(0, 29).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
        double[] grid = [.. Enumerable.Range(0, 41).Select(g => -4.0 + 0.2 * g)];
        PointwiseLogDensity training = ConformalPredictiveDistributionTests.Training(beta, 1.0, null, x, y);
        CandidateLogDensity candidates = ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, [0.4, 1.2], [grid, grid]);
        var weights = ShiftWeights.ByCase([.. x.Select(v => Math.Exp(v - 0.5))], [Math.Exp(-0.1), Math.Exp(0.7)], candidates.Grids);
        foreach (MarginBounds bounds in new[] { MarginBounds.Uniform, MarginBounds.Carried })
        {
            ConformalDistributions plain = ConformalBayes.Distribution(training, candidates, new ConformalOptions { Weights = weights, MarginBounds = bounds });
            ConformalDistributions inflated = ConformalBayes.Distribution(training, candidates, new ConformalOptions { Weights = weights, MarginBounds = bounds, Inflate = true });
            for (int t = 0; t < 2; t++)
            {
                ConformalDistribution d = plain.Distributions[t];
                Assert.NotNull(d.Weighted);
                Assert.True(d.Margin!.Weighted);
                double band = d.Margin.WeightedBand(0.1)!.Value;
                Assert.Equal(0.9 - band - d.Margin.ErrorLevel, d.CoverageUnderError(0.1), 12);
                Assert.Equal(0.9 - d.Margin.ErrorLevel, inflated.Distributions[t].CoverageUnderError(0.1), 12);
                Assert.True(double.IsNaN(d.Margin.Pit(1).Lower), "under weights the PIT counts are not read");
                for (int g = 0; g < grid.Length; g++)
                {
                    if (d.AcceptsPlain(g, 0.1)) Assert.True(d.AcceptsInflated(g, 0.1));
                    Assert.Equal(d.AcceptsInflated(g, 0.1), inflated.Distributions[t].Accepts(g, 0.1));
                }
                output.WriteLine($"{bounds}, case {t}: M̂ at α 0.1 {band:F4}, the plain set's coverage under error {d.CoverageUnderError(0.1):F4}");
            }
        }

        // Equal weights: the counts' margin.
        var flat = ShiftWeights.ByCase([.. x.Select(_ => 0.3)], [0.3, 0.3], candidates.Grids);
        ConformalDistributions counts = ConformalBayes.Distribution(training, candidates, new ConformalOptions { MarginBounds = MarginBounds.Carried });
        ConformalDistributions equal = ConformalBayes.Distribution(training, candidates, new ConformalOptions { Weights = flat, MarginBounds = MarginBounds.Carried });
        for (int t = 0; t < 2; t++)
        {
            Assert.Null(equal.Distributions[t].Weighted);
            Assert.False(equal.Distributions[t].Margin!.Weighted);
            Assert.Equal(counts.Distributions[t].Margin!.Counts, equal.Distributions[t].Margin!.Counts);
            Assert.Equal(counts.Distributions[t].CoverageUnderError(0.1), equal.Distributions[t].CoverageUnderError(0.1));
        }
    }
}
