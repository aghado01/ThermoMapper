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
/// The split predictive distribution's smoothed cell law (E16, slice 18's follow-up): G̃(F̂(y | x)) lies inside the conformal band
/// [Q⁻, Q⁺] and rises, every cell carries mass, the log score is finite where the step law's is +∞ and still ranks a right fit ahead of a
/// wrong one, its PIT is uniform within 1/(n + 1), the grid-free readout agrees with the grid, and the full step and weighed cases have none.
/// </summary>
public sealed class SmoothedCellLawTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    // A fit's draws of the slope, the calibration units at their responses, and the test cases on a grid.
    private static (CandidateLogDensity Calibration, CandidateLogDensity Test) Seams(double slope, double[] x, double[] y, double[] xt, double[] grid)
    {
        double[] beta = [.. Enumerable.Repeat(slope, 20)];
        return (ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, x, [.. y.Select(v => new[] { v })]),
                ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, xt, [.. xt.Select(_ => grid)]));
    }

    [Fact]
    public void TheSmoothedLaw_LiesInTheBand_Rises_AndGivesEveryCellMass()
    {
        var rng = new Xoshiro256PlusPlus(2101);
        double[] x = [.. Enumerable.Range(0, 49).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
        double[] xt = [.. Enumerable.Range(0, 20).Select(_ => Gaussian(rng))];
        double[] grid = [.. Enumerable.Range(0, 161).Select(g => -8.0 + 0.1 * g)];
        (CandidateLogDensity calibration, CandidateLogDensity test) = Seams(0.5, x, y, xt, grid);
        ConformalDistributions d = SplitConformal.Distribution(calibration, test);
        foreach (ConformalDistribution c in d.Distributions)
        {
            SmoothedCdf s = c.Smoothed!;
            for (int g = 0; g < grid.Length; g++)
            {
                Assert.InRange(s.Lower[g], c.Lower(g) - 1e-12, c.Upper(g) + 1e-12);
                Assert.Equal(1.0, s.Lower[g] + s.Upper[g], 12);
                if (g > 0) Assert.True(s.Lower[g] >= s.Lower[g - 1] - 1e-15);
            }
            double[] masses = c.CellMasses(CellLaw.Smoothed)!;
            Assert.Equal(1.0, masses.Sum(), 12);
            Assert.All(masses, p => Assert.True(p > 0.0));
            int empty = c.CellMasses()!.Count(p => p <= 0.0);
            Assert.True(empty > 0, "the fixture's step law should leave cells empty");
        }
    }

    [Fact]
    public void TheSmoothedLogScore_IsFinite_AndRanksTheFits()
    {
        // Slice 18's fixture: 199 calibration units, 400 test cases, cells of width 1; the right fit (slope ½) against one ignoring x.
        var rng = new Xoshiro256PlusPlus(1902);
        double[] grid = [.. Enumerable.Range(0, 9).Select(g => -4.0 + 1.0 * g)];
        double[] x = [.. Enumerable.Range(0, 199).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
        double[] xt = [.. Enumerable.Range(0, 400).Select(_ => Gaussian(rng))], yt = [.. xt.Select(v => 0.5 * v + Gaussian(rng))];
        (double Smoothed, int StepInfinite) Score(double slope)
        {
            (CandidateLogDensity calibration, CandidateLogDensity test) = Seams(slope, x, y, xt, grid);
            ConformalDistributions d = SplitConformal.Distribution(calibration, test);
            double total = 0.0;
            int infinite = 0;
            for (int t = 0; t < xt.Length; t++)
            {
                double smoothed = d.Distributions[t].LogScore(yt[t], CellLaw.Smoothed);
                Assert.True(double.IsFinite(smoothed));
                total += smoothed;
                if (double.IsPositiveInfinity(d.Distributions[t].LogScore(yt[t]))) infinite++;
            }
            return (total / xt.Length, infinite);
        }
        (double right, int rightInfinite) = Score(0.5);
        (double wrong, int wrongInfinite) = Score(0.0);
        output.WriteLine($"the smoothed log score, the right fit against one ignoring x: {right:F4} and {wrong:F4}; the step law's +∞ in {rightInfinite} and {wrongInfinite} of {xt.Length}");
        Assert.True(right < wrong);
    }

    [Fact]
    public void TheSmoothedPit_IsUniformWithinOneOverNPlusOne()
    {
        // 4,000 replications of 99 calibration units and a test case: G̃(F̂(y)) at the true y against uniform at the deciles.
        const int n = 99, replications = 4000;
        var rng = new Xoshiro256PlusPlus(2026_1135);
        var pits = new List<double>();
        for (int r = 0; r < replications; r++)
        {
            double[] x = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.5 * v + 1.3 * Gaussian(rng))];
            double xt = Gaussian(rng), yt = 0.5 * xt + 1.3 * Gaussian(rng);
            double[] beta = [.. Enumerable.Repeat(0.5, 10)];
            CandidateLogDensity calibration = ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, x, [.. y.Select(v => new[] { v })]);
            CandidateLogDensity test = ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, [xt], [[yt]]);
            pits.Add(SplitConformal.Distribution(calibration, test).Distributions[0].Smoothed!.Lower[0]);
        }
        double worst = 0.0;
        for (int k = 1; k <= 9; k++)
        {
            double level = k / 10.0;
            worst = Math.Max(worst, Math.Abs(pits.Count(v => v <= level) / (double)replications - level));
        }
        double se = Math.Sqrt(0.25 / replications);
        output.WriteLine($"the smoothed PIT's largest departure from uniform at the deciles: {worst:F4} (1/(n + 1) = {1.0 / (n + 1):F4}, se at most {se:F4})");
        Assert.True(worst <= 1.0 / (n + 1) + 4.0 * se);
    }

    [Fact]
    public void TheGridFreeReadout_IsTheGrids_AndTheFullStepHasNone()
    {
        var rng = new Xoshiro256PlusPlus(2102);
        double[] x = [.. Enumerable.Range(0, 30).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
        double[] grid = [.. Enumerable.Range(0, 41).Select(g => -4.0 + 0.2 * g)];
        double[] beta = [.. Enumerable.Range(0, 20).Select(_ => 0.5 + 0.1 * Gaussian(rng))];
        CandidateLogDensity calibration = ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, x, [.. y.Select(v => new[] { v })]);
        CandidateLogDensity test = ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, [0.4], [grid]);
        SplitPitCalibration fit = SplitConformal.CalibrateDistribution(calibration);
        ConformalDistribution d = fit.OnGrid(test).Distributions[0];
        SplitPredictiveDistribution free = fit.AtCase(v => ConformalPredictiveDistributionTests.Seam(beta, 1.0, null, [0.4], [[v]]));
        for (int g = 0; g < grid.Length; g++) Assert.Equal(d.Smoothed!.Lower[g], free.SmoothedAt(grid[g]));

        // The full step and a weighed case: no smoothed law, the scores on it NaN.
        PointwiseLogDensity training = ConformalPredictiveDistributionTests.Training(beta, 1.0, null, x, y);
        ConformalDistribution full = ConformalBayes.Distribution(training, test).Distributions[0];
        Assert.Null(full.Smoothed);
        Assert.True(double.IsNaN(full.LogScore(0.3, CellLaw.Smoothed)));
        ConformalDistribution weighed = fit.OnGrid(test, ShiftWeights.ByCase([.. x.Select(v => 1.0 + v * v)], [2.0], [grid])).Distributions[0];
        Assert.Null(weighed.Smoothed);
        Assert.Null(weighed.CellMasses(CellLaw.Smoothed));
        Assert.Throws<ArgumentOutOfRangeException>(() => d.CellMasses((CellLaw)5));
    }
}
