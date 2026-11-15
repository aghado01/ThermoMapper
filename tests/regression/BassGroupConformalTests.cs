using System;
using System.Linq;
using Maths.Regression.Conformal;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's opt-in group-units call (<see cref="BassFit.GroupConformal"/>; the user's choice (b) of 2026-11-14), end to end on the crossed
/// school × reader fixture. Two held-out readers, each seen at two schools, span four cells. Named by their grouping they are two groups,
/// reproducing the explicit labels' p-values, while BASS's case-level call still refuses them pooled in one stratum. A fit with one
/// grouping reads its groups from the cells.
/// </summary>
public sealed class BassGroupConformalTests
{
    private const UnseenLevel New = UnseenLevel.OutsideEverySubset;

    private static (BassFit Fit, BassData Data) Crossed(bool readerOnly)
    {
        var (config, data) = RandomEffectsMultipleReference.Fixture("crossed", drawn: false);
        config = config with { Noise = BassCandidateTests.Noise("gaussian"), Intercept = GPriorCentre.Mean,
                               Sampling = config.Sampling with { BurnIn = 40, Samples = 16 } };
        if (readerOnly)
        {
            config = config with { RandomIntercepts = [config.RandomIntercepts!.Single(s => s.Grouping == "reader")] };
            data = new BassData(data.Predictors, data.Response, [data.Groupings.Single(g => g.Name == "reader")]);
        }
        return (Bass.Fit(config, data), data);
    }

    // Two held-out readers, N1 and N2, each with repeated cases at the training schools L0 and L1; a third new reader's two test cases.
    // The one-grouping fit's cases carry no school column, its data no school grouping.
    private static PredictionCases HeldOut(bool schools = true, bool readers = true) => new([new ContinuousColumn("x", [0.2, 0.5, 0.7, 0.3, 0.6, 0.4, 0.8]),
        .. schools ? new PredictionColumn[] { new CategoricalColumn("school", ["L0", "L0", "L1", "L1", "L0", "L0", "L1"]) } : [],
        .. readers ? new PredictionColumn[] { new CategoricalColumn("reader", ["N1", "N1", "N1", "N1", "N2", "N2", "N2"]) } : []]);
    private static readonly double[] Responses = [0.4, 0.9, 0.7, 1.1, 0.2, 0.5, 0.8];
    private static readonly PredictionCases Test = new([new ContinuousColumn("x", [0.45, 0.55]),
        new CategoricalColumn("school", ["L0", "L1"]), new CategoricalColumn("reader", ["N3", "N3"])]);
    private static readonly PredictionCases ReaderTest = new([new ContinuousColumn("x", [0.45, 0.55]), new CategoricalColumn("reader", ["N3", "N3"])]);
    private static readonly double[] Grid = [.. Enumerable.Range(0, 26).Select(g => -2.0 + 0.2 * g)];
    private static readonly int[] Labels = [0, 0, 0, 0, 1, 1, 1];

    [Fact]
    public void HeldOutSubjectsAcrossSites_AreTwoGroups_AndTheCaseLevelCallStillRefuses()
    {
        var (fit, data) = Crossed(readerOnly: false);
        Assert.Equal(4, fit.GroupCells(data, HeldOut(), New).Candidates!.Distinct().Count());
        foreach (GroupScore score in Enum.GetValues<GroupScore>())
        {
            ConformalPrediction named = fit.GroupConformal(data, HeldOut(), Responses, Test, [Grid, Grid], grouping: "reader", alpha: 0.5,
                                                           score: score, seed: 3, unseenLevels: New);
            ConformalPrediction labelled = fit.GroupConformal(data, HeldOut(), Responses, Test, [Grid, Grid], groupOf: Labels, alpha: 0.5,
                                                              score: score, seed: 3, unseenLevels: New);
            for (int t = 0; t < 2; t++)
            {
                Assert.Equal(3, named.Sets[t].Units);
                Assert.NotNull(named.Sets[t].PValues);
                Assert.Equal(labelled.Sets[t].PValues, named.Sets[t].PValues);
            }
        }

        // The cells would count four groups, so a crossed fit names its groups.
        Assert.Contains("name the grouping", Assert.Throws<ArgumentException>(
            () => fit.GroupConformal(data, HeldOut(), Responses, Test, [Grid, Grid], alpha: 0.5, unseenLevels: New)).Message);
        // BASS's case-level call still refuses the pooled stratum, left to its default or given as one.
        double[][] at = [.. Responses.Select(y => new[] { y })];
        Assert.Contains("GroupCells", Assert.Throws<ArgumentException>(() => fit.CandidateLogDensity(data, HeldOut(), at, unseenLevels: New)).Message);
        Assert.Contains("spans cells", Assert.Throws<ArgumentException>(() => fit.CandidateLogDensity(data, HeldOut(), at, new int[7], New)).Message);

        // Both readings at once; a grouping no effect reads; a held-out reader the fit knows, at a new school; a case of known groups only.
        Assert.Throws<ArgumentException>(() => fit.GroupConformal(data, HeldOut(), Responses, Test, [Grid, Grid], Labels, "reader", unseenLevels: New));
        Assert.Contains("No effect", Assert.Throws<ArgumentException>(
            () => fit.GroupConformal(data, HeldOut(), Responses, Test, [Grid, Grid], grouping: "class", unseenLevels: New)).Message);
        var knownReader = new PredictionCases([new ContinuousColumn("x", [0.3]), new CategoricalColumn("school", ["S9"]),
                                               new CategoricalColumn("reader", ["L1"])]);
        Assert.Contains("one the fit knows", Assert.Throws<ArgumentException>(
            () => fit.GroupConformal(data, knownReader, [0.5], Test, [Grid, Grid], grouping: "reader", unseenLevels: New)).Message);
        var known = new PredictionCases([new ContinuousColumn("x", [0.3]), new CategoricalColumn("school", ["L0"]),
                                         new CategoricalColumn("reader", ["L1"])]);
        Assert.Contains("all known", Assert.Throws<ArgumentException>(
            () => fit.GroupConformal(data, known, [0.5], Test, [Grid, Grid], groupOf: [0], unseenLevels: New)).Message);
    }

    [Fact]
    public void OneGrouping_ReadsItsGroupsFromTheCells()
    {
        var (fit, data) = Crossed(readerOnly: true);
        ConformalPrediction byCells = fit.GroupConformal(data, HeldOut(schools: false), Responses, ReaderTest, [Grid, Grid], alpha: 0.5, unseenLevels: New);
        ConformalPrediction named = fit.GroupConformal(data, HeldOut(schools: false), Responses, ReaderTest, [Grid, Grid], grouping: "reader",
                                                       alpha: 0.5, unseenLevels: New);
        for (int t = 0; t < 2; t++)
        {
            Assert.Equal(3, byCells.Sets[t].Units);
            Assert.Equal(named.Sets[t].PValues, byCells.Sets[t].PValues);
        }
        // With the reader's column left out each held-out case is a cell of its own, which names no group.
        Assert.Contains("cell of its own", Assert.Throws<ArgumentException>(
            () => fit.GroupConformal(data, HeldOut(schools: false, readers: false), Responses, ReaderTest, [Grid, Grid], alpha: 0.5, unseenLevels: New)).Message);
    }

    [Fact]
    public void TheConfidence_PassesThroughToTheGroupForm()
    {
        // Twelve held-out readers, each at schools L0 and L1: m 13, α 0.5, so CB-1's k is 6 and the (α, 0.9) level is smaller.
        var (fit, data) = Crossed(readerOnly: false);
        const int groups = 12;
        var heldOut = new PredictionCases([new ContinuousColumn("x", [.. Enumerable.Range(0, 2 * groups).Select(i => (i % 9 + 0.5) / 9.0)]),
            new CategoricalColumn("school", [.. Enumerable.Range(0, 2 * groups).Select(i => i % 2 == 0 ? "L0" : "L1")]),
            new CategoricalColumn("reader", [.. Enumerable.Range(0, 2 * groups).Select(i => $"R{i / 2}")])]);
        double[] responses = [.. Enumerable.Range(0, 2 * groups).Select(i => 0.6 + 0.4 * Math.Sin(1.3 * i))];
        int[] labels = [.. Enumerable.Range(0, 2 * groups).Select(i => i / 2)];

        // The same seams through the case-level call, each case in its cell's stratum, then pooled: the strata do not touch the values.
        static CandidateLogDensity Pooled(CandidateLogDensity c)
            => new(c.Values, c.Grids, c.Latents, null, c.ChainLengths) { NewGroup = c.NewGroup, Forms = c.Forms, FormsReason = c.FormsReason };
        CandidateLogDensity calibration = Pooled(fit.CandidateLogDensity(data, heldOut, [.. responses.Select(y => new[] { y })],
                                                                         fit.GroupCells(data, heldOut, New).Candidates, New));
        CandidateLogDensity test = Pooled(fit.CandidateLogDensity(data, Test, [Grid, Grid], fit.GroupCells(data, Test, New).Candidates, New));

        ConformalPrediction omitted = fit.GroupConformal(data, heldOut, responses, Test, [Grid, Grid], labels, alpha: 0.5, unseenLevels: New);
        ConformalPrediction none = fit.GroupConformal(data, heldOut, responses, Test, [Grid, Grid], labels, alpha: 0.5, unseenLevels: New, confidence: null);
        ConformalPrediction direct = GroupConformal.Density(calibration, labels, test, 0.5);
        ConformalPrediction built = fit.GroupConformal(data, heldOut, responses, Test, [Grid, Grid], labels, alpha: 0.5, unseenLevels: New, confidence: 0.9);
        ConformalPrediction builtDirect = GroupConformal.Density(calibration, labels, test, 0.5, confidence: 0.9);
        Assert.Null(omitted.Confidence);
        Assert.Equal(0.9, built.Confidence);
        for (int t = 0; t < 2; t++)
        {
            Assert.Equal((13, 6), (omitted.Sets[t].Units, omitted.Sets[t].Level));
            Assert.Equal(omitted.Sets[t].PValues, none.Sets[t].PValues);
            Assert.Equal(direct.Sets[t].PValues, omitted.Sets[t].PValues);
            Assert.Equal(TrainingConditional.Level(13, 0.5, 0.9), built.Sets[t].Level);
            Assert.True(built.Sets[t].Level < omitted.Sets[t].Level);
            Assert.Equal(builtDirect.Sets[t].Level, built.Sets[t].Level);
            Assert.Equal(builtDirect.Sets[t].PValues, built.Sets[t].PValues);
            Assert.True(built.Sets[t].TrainingConditionalMiscoverage(0.9) <= 0.5);
        }
    }
}
