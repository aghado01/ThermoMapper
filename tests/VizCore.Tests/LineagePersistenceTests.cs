using System;
using System.Linq;
using Clustering.Dendrograms;
using Clustering.Primitives;
using Graphs.Primitives;
using UserRepl.Commands;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// Track 2 — the lineage-persistence resolver (<see cref="LineagePersistence"/>,
/// wave_clus's insight lifted): lineage tracking across the T-stack, split
/// detection, the score in ln T, the two readings and the representatives.
/// Oracles hand-computed (validation independence). Temperatures are powers of
/// e, so a level's cell is a sum of half-integers in ln T.
/// </summary>
public sealed class LineagePersistenceTests
{
    private static double E(double x) => Math.Exp(x);

    // Path graph 0—1—2—3.
    private static CsrGraph PathGraph() => CsrGraph.FromEdges(
        new[] { new Edge(0, 1, 1.0), new Edge(1, 2, 1.0), new Edge(2, 3, 1.0) }, 4);

    private static double[] Column(CsrGraph graph, double v01, double v12, double v23)
    {
        var col = new double[graph.Targets.Length];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
        {
            double v = (edge.Source, edge.Target) switch
            {
                (0, 1) or (1, 0) => v01,
                (1, 2) or (2, 1) => v12,
                _                => v23,   // (2,3)
            };
            col[edge.Slot] = v;
        }
        return col;
    }

    /// <summary>
    /// A stack that splits once and the children persist, given the giant's and the halves' ln T:
    ///   giant levels: all bonds on → {0,1,2,3}
    ///   half levels:  the middle bond off → {0,1} | {2,3} (a clean split)
    /// </summary>
    private static (CsrGraph Graph, PartitionHierarchy Stack, double[][] Cols) SplitStack(double[] giantLnT, double[] halfLnT)
    {
        var graph = PathGraph();
        double[] temps = giantLnT.Concat(halfLnT).Select(E).ToArray();
        double[][] cols = giantLnT.Select(_ => Column(graph, 0.9, 0.9, 0.9))
            .Concat(halfLnT.Select(_ => Column(graph, 0.9, 0.1, 0.9))).ToArray();
        return (graph, DenseTStack.Build(graph, temps, cols, theta: 0.5), cols);
    }

    /// <summary>The giant at ln T = 0, 1 and the halves at ln T = 2, 3, 4.</summary>
    private static (CsrGraph Graph, PartitionHierarchy Stack, double[][] Cols) SplitStack()
        => SplitStack(new[] { 0.0, 1.0 }, new[] { 2.0, 3.0, 4.0 });

    private static LineagePersistenceConfig Min2 => new() { MinClusterSize = 2 };

    [Fact]
    public void SelectsThePersistentHalves_NotTheColdGiant()
    {
        var (graph, stack, cols) = SplitStack();

        LineagePersistenceResult result = LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2);

        // The giant's cells are 0.5 + 1 = 1.5 e-folds, and each half's 1 + 1 + 0.5 = 2.5. δ̄ is 0.9 inside
        // both, so the scores are 1.35 and 2.25. The halves outrank the giant, which overlaps them and drops.
        var giant = result.AllLineages.Single(l => l.Members.Length == 4);
        Assert.Equal(1.5, giant.Lifetime, precision: 12);
        Assert.Equal(1.35, giant.Score, precision: 12);
        Assert.Equal(2, result.Selected.Count);
        Assert.All(result.Selected, l =>
        {
            Assert.Equal(2.5, l.Lifetime, precision: 12);
            Assert.Equal(2.25, l.Score, precision: 12);
            Assert.Equal(0.9, l.Decidedness, precision: 12);
        });

        int[] labels = result.Assignment.Labels;
        Assert.Equal(2, result.Assignment.Count);
        Assert.Equal(labels[0], labels[1]);           // {0,1} together
        Assert.Equal(labels[2], labels[3]);           // {2,3} together
        Assert.NotEqual(labels[0], labels[2]);        // distinct halves
        Assert.DoesNotContain(Assignment.Unassigned, labels);   // all four covered
    }

    [Fact]
    public void TheGiantLineageEndsAtTheSplit()
    {
        var (graph, stack, cols) = SplitStack();
        LineagePersistenceResult result = LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2);

        // The size-4 lineage spans only the two cold levels: it ends at the split
        // rather than chaining into one of the halves.
        var giant = result.AllLineages.Single(l => l.Members.Length == 4);
        Assert.Equal(1.0, giant.TBirth, precision: 12);
        Assert.Equal(E(1.0), giant.TDeath, precision: 12);
        Assert.Equal(2, giant.LevelCount);
        Assert.Equal(3, result.AllLineages.Count);
    }

    [Fact]
    public void TemperatureWindow_BoundsTheAnalyzedLevels()
    {
        var (graph, stack, cols) = SplitStack();

        // Window excludes the cold giant levels: only the split levels remain, so the two
        // halves are tracked and the giant never appears. Their cells shrink with the window:
        // 0.5 + 1 + 0.5 = 2 e-folds.
        LineagePersistenceResult result = LineagePersistence.Resolve(graph, stack, cols, q: 20,
            Min2 with { TemperatureWindow = (E(1.5), E(4.5)) });

        Assert.DoesNotContain(result.AllLineages, l => l.Members.Length == 4);
        Assert.Equal(2, result.Selected.Count);
        Assert.All(result.Selected, l => Assert.Equal(2.0, l.Lifetime, precision: 12));
    }

    [Fact]
    public void Score_IsInvariantToTheGridsDensity_WhereALevelCountIsNot()
    {
        // One more level inside the giant's range and one inside the halves': the counts go from 2 to 3
        // and from 3 to 4, and the cells still sum to 1.5 and 2.5.
        var (graph, stack, cols) = SplitStack(new[] { 0.0, 0.5, 1.0 }, new[] { 2.0, 2.5, 3.0, 4.0 });
        LineagePersistenceResult result = LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2);

        var giant = result.AllLineages.Single(l => l.Members.Length == 4);
        Assert.Equal(3, giant.LevelCount);
        Assert.Equal(1.5, giant.Lifetime, precision: 12);
        Assert.Equal(1.35, giant.Score, precision: 12);
        Assert.Equal(2, result.Selected.Count);
        Assert.All(result.Selected, l =>
        {
            Assert.Equal(4, l.LevelCount);
            Assert.Equal(2.5, l.Lifetime, precision: 12);
        });
    }

    // Two components: 0—1 and 2—3. Each list says which pair's bond is on, cold to hot.
    private static (CsrGraph Graph, PartitionHierarchy Stack, double[][] Cols) TwoPairs(bool[] first, bool[] second)
    {
        var graph = CsrGraph.FromEdges(new[] { new Edge(0, 1, 1.0), new Edge(2, 3, 1.0) }, 4);
        var cols = new double[first.Length][];
        for (int i = 0; i < first.Length; i++)
        {
            cols[i] = new double[graph.Targets.Length];
            foreach (UndirectedEdge edge in graph.UndirectedEdges())
                cols[i][edge.Slot] = (edge.Source is 0 or 1 ? first[i] : second[i]) ? 0.9 : 0.05;
        }
        double[] temps = Enumerable.Range(0, first.Length).Select(i => E(i)).ToArray();
        return (graph, DenseTStack.Build(graph, temps, cols, theta: 0.5), cols);
    }

    [Fact]
    public void ADisjointLineage_IsKept_HoweverMuchLongerAnotherLives()
    {
        // {0,1} at every level; {2,3} melts after the first. Lifetimes 4 and 0.5 e-folds: a ratio of 8,
        // which a cut at the largest ratio in the sorted scores would drop.
        var (graph, stack, cols) = TwoPairs(new[] { true, true, true, true, true }, new[] { true, false, false, false, false });
        LineagePersistenceResult result = LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2);

        Assert.Equal(2, result.Selected.Count);
        Assert.Equal(4.0, result.Selected[0].Lifetime, precision: 12);
        Assert.Equal(0.5, result.Selected[1].Lifetime, precision: 12);
        Assert.DoesNotContain(Assignment.Unassigned, result.Assignment.Labels);
    }

    [Fact]
    public void SelectTopK_KeepsTheHighestScores()
    {
        var (graph, stack, cols) = TwoPairs(new[] { true, true, true }, new[] { true, false, false });
        LineagePersistenceResult result = LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2 with { SelectTopK = 1 });

        var only = Assert.Single(result.Selected);
        Assert.Equal(new[] { 0, 1 }, only.Members);
        Assert.Equal(Assignment.Unassigned, result.Assignment.Labels[2]);
    }

    // A path of 30 nodes, with a chip of `chip` nodes cut from one end above ln T = 0.
    private static (CsrGraph Graph, PartitionHierarchy Stack, double[][] Cols) ChipStack(int chip)
    {
        const int n = 30;
        var graph = CsrGraph.FromEdges(Enumerable.Range(0, n - 1).Select(i => new Edge(i, i + 1, 1.0)).ToArray(), n);
        double[] whole = new double[graph.Targets.Length], chipped = new double[graph.Targets.Length];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
        {
            whole[edge.Slot] = 0.9;
            chipped[edge.Slot] = Math.Min(edge.Source, edge.Target) == chip - 1 ? 0.1 : 0.9;
        }
        double[][] cols = { whole, chipped, chipped };
        return (graph, DenseTStack.Build(graph, new[] { E(0.0), E(1.0), E(2.0) }, cols, theta: 0.5), cols);
    }

    private static LineagePersistenceResult ChipRun(int chip, double? share = null)
    {
        var (graph, stack, cols) = ChipStack(chip);
        var config = new LineagePersistenceConfig { MinClusterSize = 2 };
        if (share is double s) config = config with { SplitShare = s };
        return LineagePersistence.Resolve(graph, stack, cols, q: 20, config);
    }

    [Fact]
    public void ShareBelowTheDefault_IsErosion_AndTheLineageContinues()
    {
        // A chip of 2 on 30 is a share of 0.067: the main lineage runs through all three levels.
        LineagePersistenceResult result = ChipRun(chip: 2);
        var main = result.AllLineages.Single(l => l.Members.Length == 30);
        Assert.Equal(3, main.LevelCount);
        Assert.Equal(2.0, main.Lifetime, precision: 12);
    }

    [Fact]
    public void ShareAtTheDefault_IsASplit_AndTheLineageEnds()
    {
        // A chip of 3 on 30 is a share of exactly 0.10.
        LineagePersistenceResult result = ChipRun(chip: 3);
        Assert.Equal(LineagePersistenceConfig.DefaultSplitShare, result.Config.SplitShare);
        var cold = result.AllLineages.Single(l => l.Members.Length == 30);
        Assert.Equal(1, cold.LevelCount);
        Assert.Equal(0.5, cold.Lifetime, precision: 12);
    }

    [Fact]
    public void Share_IsExposed_AndATinyOneIsTheSizeRule()
    {
        Assert.Equal(1, ChipRun(chip: 2, share: 1e-9).AllLineages.Single(l => l.Members.Length == 30).LevelCount);
        Assert.Equal(3, ChipRun(chip: 3, share: 0.2).AllLineages.Single(l => l.Members.Length == 30).LevelCount);
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(-0.1)]
    [InlineData(1.5)]
    [InlineData(double.NaN)]
    public void Share_OutsideItsRange_Throws(double share)
    {
        var (graph, stack, cols) = SplitStack();
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2 with { SplitShare = share }));
    }

    /// <summary>
    /// One lineage eroding across four levels, ln T = 0, 1.5, 2, 3 (a path of 6 nodes):
    ///   6 points, then 5, then 4, then 4, with δ̄ falling and then recovering. The dropped nodes are
    /// singletons, below the size floor, so the lineage is one chain. Its midpoint, ln T = 1.5, is a level.
    /// </summary>
    private static (CsrGraph Graph, PartitionHierarchy Stack, double[][] Cols) ErodingStack()
    {
        const int n = 6;
        var graph = CsrGraph.FromEdges(Enumerable.Range(0, n - 1).Select(i => new Edge(i, i + 1, 1.0)).ToArray(), n);
        double[] Col(int keptEdges, double inside)
        {
            var col = new double[graph.Targets.Length];
            foreach (UndirectedEdge edge in graph.UndirectedEdges())
                col[edge.Slot] = Math.Min(edge.Source, edge.Target) < keptEdges ? inside : 0.1;
            return col;
        }
        double[][] cols = { Col(5, 0.95), Col(4, 0.85), Col(3, 0.60), Col(3, 0.75) };
        double[] temps = { E(0.0), E(1.5), E(2.0), E(3.0) };
        return (graph, DenseTStack.Build(graph, temps, cols, theta: 0.5), cols);
    }

    [Theory]
    [InlineData(LineageRepresentative.Largest, 6, 0.0)]        // the coldest level holds all 6
    [InlineData(LineageRepresentative.Midpoint, 5, 1.5)]
    [InlineData(LineageRepresentative.LeastDecided, 4, 2.0)]   // δ̄ = 0.60, at ln T = 2 and not at the hot end
    public void Representative_PicksItsLevel(LineageRepresentative kind, int size, double lnT)
    {
        var (graph, stack, cols) = ErodingStack();
        LineagePersistenceResult result = LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2 with { Representative = kind });

        var lineage = Assert.Single(result.Selected);
        Assert.Equal(size, lineage.Members.Length);
        Assert.Equal(size, result.Assignment.Labels.Count(l => l != Assignment.Unassigned));
        Assert.Equal(E(lnT), lineage.Representative(kind).Temperature, precision: 12);
    }

    [Fact]
    public void Lineage_KeepsEveryLevel_SoOtherReductionsAreViews()
    {
        var (graph, stack, cols) = ErodingStack();
        var lineage = Assert.Single(LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2).Selected);

        Assert.Equal(new[] { 6, 5, 4, 4 }, lineage.Levels.Select(l => l.Members.Length));
        double[] cells = { 0.75, 1.0, 0.75, 0.5 }, decided = { 0.95, 0.85, 0.60, 0.75 };
        for (int i = 0; i < cells.Length; i++)
        {
            Assert.Equal(cells[i], lineage.Levels[i].CellWidth, precision: 12);
            Assert.Equal(decided[i], lineage.Levels[i].Decidedness, precision: 12);
        }
        // ∫ δ̄ d ln T = 0.75·0.95 + 0.85 + 0.75·0.60 + 0.5·0.75 = 2.3875, over 3 e-folds.
        Assert.Equal(2.3875, lineage.Score, precision: 12);
        Assert.Equal(3.0, lineage.Lifetime, precision: 12);
        Assert.Equal(2.3875 / 3.0, lineage.Decidedness, precision: 12);
        Assert.Equal(4, lineage.LevelCount);
    }

    [Fact]
    public void ALevelHoldingTwoOfALineagesClusters_CountsOnce()
    {
        // δ̄ rises with T on the middle bond, so this stack doesn't nest: {0,1} and {2,3} at ln T = 0
        // both continue into {0,1,2,3} at ln T = 1, and the one lineage holds two clusters at ln T = 0.
        var graph = PathGraph();
        double[][] cols = { Column(graph, 0.8, 0.1, 0.6), Column(graph, 0.9, 0.9, 0.9) };
        var stack = DenseTStack.Build(graph, new[] { E(0.0), E(1.0) }, cols, theta: 0.5);

        var lineage = Assert.Single(LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2).AllLineages);

        Assert.Equal(3, lineage.Levels.Count);
        Assert.Equal(2, lineage.LevelCount);
        Assert.Equal(1.0, lineage.Lifetime, precision: 12);    // 0.5 + 0.5, not 0.5 + 0.5 + 0.5
        // ∫ δ̄ d ln T = 0.5·mean(0.8, 0.6) + 0.5·0.9: the giant's internal δ̄ is (0.9 + 0.9 + 0.9) / 3.
        Assert.Equal(0.8, lineage.Score, precision: 12);
    }

    [Fact]
    public void Slice_ReadsTheLevelNearestInLnT_WithEachLineageWhole()
    {
        var (graph, stack, cols) = SplitStack();

        // T* = 1.7: ln T* = 0.53, so level ln T = 1 is nearest in ln T, where in T it would be level 0.
        LineagePersistenceResult cold = LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2 with { SliceTemperature = 1.7 });
        Assert.Equal(E(1.0), cold.SliceLevel!.Value, precision: 12);
        Assert.Equal(4, Assert.Single(cold.Selected).Members.Length);
        Assert.DoesNotContain(Assignment.Unassigned, cold.Assignment.Labels);

        LineagePersistenceResult hot = LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2 with { SliceTemperature = E(3.2) });
        Assert.Equal(E(3.0), hot.SliceLevel!.Value, precision: 12);
        Assert.Equal(2, hot.Selected.Count);
        Assert.Equal(hot.Assignment.Labels[0], hot.Assignment.Labels[1]);
        Assert.NotEqual(hot.Assignment.Labels[0], hot.Assignment.Labels[2]);
    }

    [Fact]
    public void Slice_TakesSelectTopK_ByScore()
    {
        // {0,1} at all three levels and {2,3} at the first: both alive at ln T = 0, and the longer-lived kept.
        var (graph, stack, cols) = TwoPairs(new[] { true, true, true }, new[] { true, false, false });
        LineagePersistenceResult result = LineagePersistence.Resolve(graph, stack, cols, q: 20,
            Min2 with { SliceTemperature = E(0.0), SelectTopK = 1 });

        Assert.Equal(new[] { 0, 1 }, Assert.Single(result.Selected).Members);
    }

    [Fact]
    public void ALabel_IsItsLineagesIndexInSelected()
    {
        var (graph, stack, cols) = TwoPairs(new[] { true, true, true }, new[] { true, false, false });
        LineagePersistenceResult result = LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2);

        Assert.Equal(2, result.Selected.Count);
        for (int k = 0; k < result.Selected.Count; k++)
            Assert.All(result.Selected[k].Members, p => Assert.Equal(k, result.Assignment.Labels[p]));
    }

    [Fact]
    public void ANonPositiveTemperature_Throws_SinceLifetimesReadLnT()
    {
        var graph = PathGraph();
        double[][] cols = { Column(graph, 0.9, 0.9, 0.9), Column(graph, 0.9, 0.1, 0.9) };
        var stack = DenseTStack.Build(graph, new[] { 0.0, 1.0 }, cols, theta: 0.5);

        Assert.Throws<ArgumentException>(() => LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2));
    }

    [Fact]
    public void LongestLived_HasNoSliceLevel()
    {
        var (graph, stack, cols) = SplitStack();
        Assert.Null(LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2).SliceLevel);
    }

    [Theory]
    [InlineData(0.5)]    // below the coldest level, T = 1
    [InlineData(60.0)]   // above the hottest, T = e⁴ ≈ 54.6
    [InlineData(-1.0)]
    public void Slice_OutsideTheLevelsRead_Throws_RatherThanReadingTheEdge(double tStar)
    {
        var (graph, stack, cols) = SplitStack();
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2 with { SliceTemperature = tStar }));
    }

    [Fact]
    public void Slice_OutsideTheWindow_Throws()
    {
        var (graph, stack, cols) = SplitStack();
        Assert.Throws<ArgumentOutOfRangeException>(() => LineagePersistence.Resolve(graph, stack, cols, q: 20,
            Min2 with { TemperatureWindow = (E(1.5), E(4.5)), SliceTemperature = 1.0 }));
    }

    [Fact]
    public void Slice_ALineageLeftWithNoPoint_IsNotSelected()
    {
        // A chip of 2 on 30 stays erosion, so the main lineage's largest level holds the chip's points too, and
        // the main lineage scores higher. So both are alive at ln T = 2, and the chip claims nothing.
        var (graph, stack, cols) = ChipStack(chip: 2);
        LineagePersistenceResult result = LineagePersistence.Resolve(graph, stack, cols, q: 20,
            new LineagePersistenceConfig { MinClusterSize = 2, SliceTemperature = E(2.0) });

        var only = Assert.Single(result.Selected);
        Assert.Equal(30, only.Members.Length);
        Assert.Equal(1, result.Assignment.Count);
        Assert.Contains(result.AllLineages, l => l.Members.Length == 2);
    }

    [Fact]
    public void ColdCensored_MarksWhatReachesTheColdestLevelRead()
    {
        var (graph, stack, cols) = SplitStack();
        LineagePersistenceResult full = LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2);
        Assert.True(full.AllLineages.Single(l => l.Members.Length == 4).ColdCensored);
        Assert.All(full.AllLineages.Where(l => l.Members.Length == 2), l => Assert.False(l.ColdCensored));

        // A window's cold edge censors just as the grid's does.
        LineagePersistenceResult windowed = LineagePersistence.Resolve(graph, stack, cols, q: 20,
            Min2 with { TemperatureWindow = (E(1.5), E(4.5)) });
        Assert.All(windowed.AllLineages, l => Assert.True(l.ColdCensored));
    }

    [Fact]
    public void ExtendingTheColdEdge_FlipsTheLongestLivedReading_AndTheFlagShowsIt()
    {
        // The dependence on extent, on the design that collapses the toy to its giant: colder levels add
        // lifetime only to the giant, 4.5 e-folds here against the halves' 2.5. Now the giant outranks them and
        // drops them, and it is flagged, while the slice at the halves' level reads the same as before.
        var (graph, stack, cols) = SplitStack(new[] { -3.0, -2.0, -1.0, 0.0, 1.0 }, new[] { 2.0, 3.0, 4.0 });

        LineagePersistenceResult longest = LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2);
        var giant = Assert.Single(longest.Selected);
        Assert.Equal(4, giant.Members.Length);
        Assert.Equal(4.5, giant.Lifetime, precision: 12);
        Assert.True(giant.ColdCensored);

        LineagePersistenceResult slice = LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2 with { SliceTemperature = E(3.0) });
        Assert.Equal(2, slice.Selected.Count);
        Assert.All(slice.Selected, l =>
        {
            Assert.Equal(2, l.Members.Length);
            Assert.False(l.ColdCensored);
        });
    }

    [Fact]
    public void ASingleLevelStack_HasAUnitCell_AndStillResolves()
    {
        var graph = PathGraph();
        double[][] cols = { Column(graph, 0.9, 0.1, 0.9) };
        var stack = DenseTStack.Build(graph, new[] { 2.0 }, cols, theta: 0.5);

        LineagePersistenceResult result = LineagePersistence.Resolve(graph, stack, cols, q: 20, Min2);

        Assert.Equal(2, result.Selected.Count);
        Assert.All(result.Selected, l => Assert.Equal(0.9, l.Score, precision: 12));
    }

    [Fact]
    public void AnEmptyWindow_Abstains()
    {
        var (graph, stack, cols) = SplitStack();
        LineagePersistenceResult result = LineagePersistence.Resolve(graph, stack, cols, q: 20,
            Min2 with { TemperatureWindow = (1000.0, 2000.0), SliceTemperature = 1500.0 });

        Assert.Empty(result.AllLineages);
        Assert.Empty(result.Selected);
        Assert.Null(result.SliceLevel);
        Assert.All(result.Assignment.Labels, l => Assert.Equal(Assignment.Unassigned, l));
    }

    [Theory]
    [InlineData("largest", LineageRepresentative.Largest)]
    [InlineData(" Midpoint ", LineageRepresentative.Midpoint)]
    [InlineData("least-decided", LineageRepresentative.LeastDecided)]
    public void Cli_ParsesTheRepresentative(string token, LineageRepresentative expected)
        => Assert.Equal(expected, SpcCommand.ParseLineageRepresentative(token));

    [Fact]
    public void Cli_RejectsAnUnknownRepresentative()
        => Assert.Throws<ArgumentException>(() => SpcCommand.ParseLineageRepresentative("hottest"));
}
