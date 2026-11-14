using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The conformal layer's cell rule (E16, GBASS's follow-up for E21's slice 5): where both seams give cells, a case in a known group ranks only
/// among units of its own cell, in the full step and the split forms; a new group's case is exempt, the group-level form reads no cells,
/// a seam without cells is not checked, and the cells change no result. BASS's builders give each seam its cells.
/// </summary>
public sealed class ConformalCellRuleTests
{
    private static readonly double[] Beta = [0.4, 0.5, 0.6, 0.5, 0.45, 0.55];
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);
    private static int[] Identity(int n) => [.. Enumerable.Range(0, n)];

    private static CandidateLogDensity Cases(CandidateLogDensity at, int[] strata, int[]? cells, bool[]? fresh = null)
        => new CandidateLogDensity(at.Values, at.Grids, at.Latents, strata) { Cells = cells, NewGroup = fresh };

    private static void SamePValues(ConformalPrediction a, ConformalPrediction b)
    {
        Assert.Equal(a.Sets.Count, b.Sets.Count);
        for (int t = 0; t < a.Sets.Count; t++) Assert.Equal(a.Sets[t].PValues, b.Sets[t].PValues);
    }

    [Fact]
    public void TheFullStep_RanksAKnownGroupsCase_OnlyAmongItsCell()
    {
        var rng = new Xoshiro256PlusPlus(2301);
        const int n = 24;
        double[] x = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
        int[] cellOf = [.. Enumerable.Range(0, n).Select(i => i % 2)];
        double[] grid = [.. Enumerable.Range(0, 21).Select(g => -3.0 + 0.3 * g)];
        PointwiseLogDensity bare = ConformalPredictiveDistributionTests.Training(Beta, 1.0, null, x, y);
        var byCell = new PointwiseLogDensity(bare.Values, new UnitPartition(Identity(n), cellOf), bare.Latents) { Cells = cellOf };
        var pooled = new PointwiseLogDensity(bare.Values, UnitPartition.EachObservation(n), bare.Latents) { Cells = cellOf };
        CandidateLogDensity at = ConformalPredictiveDistributionTests.Seam(Beta, 1.0, null, [0.3, -0.2], [grid, grid]);

        // Each case in its cell's stratum: served, and the cells change nothing.
        SamePValues(ConformalBayes.FullStep(byCell with { Cells = null }, Cases(at, [0, 1], null)), ConformalBayes.FullStep(byCell, Cases(at, [0, 1], [0, 1])));

        // A case of cell 0 among cell 1's units, and a stratum whose units span both cells: refused.
        Assert.Contains("of cell 0, ranks among stratum 1", Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(byCell, Cases(at, [1, 1], [0, 1]))).Message);
        Assert.Contains("span cells 0 and 1", Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(pooled, Cases(at, [0, 0], [0, 0]))).Message);
        // The per-case readouts go through the same check.
        Assert.Throws<ArgumentException>(() => ConformalBayes.Scores(pooled, Cases(at, [0, 0], [0, 0]), 0));

        // A seam without cells: nothing to check, and a stratum no training unit holds reads no cell.
        ConformalBayes.FullStep(pooled, Cases(at, [0, 0], null));
        ConformalBayes.FullStep(pooled with { Cells = null }, Cases(at, [0, 0], [0, 1]));
        ConformalBayes.FullStep(byCell, Cases(at, [2, 3], [0, 1]));

        Assert.Throws<ArgumentException>(() => byCell with { Cells = [0, 1] });
        Assert.Throws<ArgumentException>(() => at with { Cells = [0, 1, 2] });
    }

    [Fact]
    public void TheSplitForms_RankAKnownGroupsCase_OnlyAmongItsCell_AndExemptNewGroups()
    {
        var rng = new Xoshiro256PlusPlus(2302);
        const int n = 30;
        double[] x = [.. Enumerable.Range(0, n).Select(_ => Gaussian(rng))], y = [.. x.Select(v => 0.5 * v + Gaussian(rng))];
        int[] cellOf = [.. Enumerable.Range(0, n).Select(i => i % 3)];
        double[] grid = [.. Enumerable.Range(0, 21).Select(g => -3.0 + 0.3 * g)];
        CandidateLogDensity heldOut = ConformalPredictiveDistributionTests.Seam(Beta, 1.0, null, x, [.. y.Select(v => new[] { v })]);
        CandidateLogDensity byCell = Cases(heldOut, cellOf, cellOf), pooled = Cases(heldOut, new int[n], cellOf);
        CandidateLogDensity at = ConformalPredictiveDistributionTests.Seam(Beta, 1.0, null, [0.3, -0.2, 1.0], [grid, grid, grid]);

        // Each case in its cell's stratum: served, and the cells change nothing.
        SamePValues(SplitConformal.Density(byCell with { Cells = null }, Cases(at, [0, 1, 2], null)), SplitConformal.Density(byCell, Cases(at, [0, 1, 2], [0, 1, 2])));

        // A case among another cell's units, by the density form and by the calibrated distribution; a stratum spanning cells.
        Assert.Contains("of cell 1, ranks among stratum 2", Assert.Throws<ArgumentException>(() => SplitConformal.Density(byCell, Cases(at, [0, 1, 2], [0, 1, 1]))).Message);
        Assert.Throws<ArgumentException>(() => SplitConformal.CalibrateDistribution(byCell).OnGrid(Cases(at, [0, 1, 2], [2, 1, 0])));
        Assert.Contains("span cells 0 and 1", Assert.Throws<ArgumentException>(() => SplitConformal.Density(pooled, Cases(at, [0, 0, 0], [0, 0, 0]))).Message);

        // A new group's cases are exempt, and the group-level form reads groups, which span cells.
        SamePValues(SplitConformal.Density(pooled with { Cells = null }, Cases(at, [0, 0, 0], null)),
                    SplitConformal.Density(pooled, Cases(at, [0, 0, 0], [7, 8, 9], [true, true, true])));
        SamePValues(GroupConformal.Density(pooled with { Cells = null }, cellOf, Cases(at, [0, 0, 0], null)),
                    GroupConformal.Density(pooled, cellOf, Cases(at, [0, 0, 0], [7, 8, 9])));
    }

    [Fact]
    public void BassSeams_CarryTheirCells_AndACrossedStratumIsRefused()
    {
        var (config, data) = RandomEffectsMultipleReference.Fixture("crossed", drawn: false);
        config = config with { Noise = BassCandidateTests.Noise("gaussian"), Intercept = GPriorCentre.Mean,
                               Sampling = config.Sampling with { BurnIn = 40, Samples = 16 } };
        BassFit fit = Bass.Fit(config, data);
        var cases = PredictionCases.Of(data);
        BassGroupCells cells = fit.GroupCells(data, cases);
        PointwiseLogDensity training = fit.PointwiseLogDensity(data, new UnitPartition(Identity(cells.Training.Length), cells.Training));
        double[][] grids = [.. data.Response.Values.ToArray().Select(v => new[] { v })];
        CandidateLogDensity candidates = fit.CandidateLogDensity(data, cases, grids, strata: cells.Candidates);
        Assert.Equal(cells.Training, training.Cells);
        Assert.Equal(cells.Candidates, candidates.Cells);
        Assert.Equal(cases.Cases, ConformalBayes.FullStep(training, candidates).Sets.Count);

        // Strata relabelled among the candidates lie inside cells for the builder, but put each case among another cell's training units.
        Assert.True(cells.TrainingCells > 1);
        int[] shifted = [.. cells.Candidates!.Select(c => (c + 1) % cells.TrainingCells)];
        CandidateLogDensity crossed = fit.CandidateLogDensity(data, cases, grids, strata: shifted);
        Assert.Contains("ranks among stratum", Assert.Throws<ArgumentException>(() => ConformalBayes.FullStep(training, crossed)).Message);
    }
}
