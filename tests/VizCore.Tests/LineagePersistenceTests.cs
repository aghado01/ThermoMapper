using System.Linq;
using Clustering.Graphical.SPC.Partitions.Hierarchical;
using Clustering.Primitives;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// Track 2 — the lineage-persistence resolver (<see cref="LineagePersistence"/>,
/// wave_clus's insight lifted): lineage tracking across the T-stack, split
/// detection, persistence scoring,
/// and gap selection. Oracles hand-computed (validation independence).
/// </summary>
public sealed class LineagePersistenceTests
{
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
    /// A stack that splits once and the children persist:
    ///   T=1,2: all bonds hot → {0,1,2,3}
    ///   T=3,4,5: the middle bond cold → {0,1} | {2,3} (a clean split)
    /// The two halves outlive the cold giant ⇒ they are the persistent lineages.
    /// </summary>
    private static (CsrGraph Graph, PartitionHierarchy Stack, double[][] Cols) SplitStack()
    {
        var graph = PathGraph();
        var temps = new[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
        var cols = new[]
        {
            Column(graph, 0.9, 0.9, 0.9),  // {0,1,2,3}
            Column(graph, 0.9, 0.9, 0.9),  // {0,1,2,3}
            Column(graph, 0.9, 0.1, 0.9),  // {0,1} {2,3}
            Column(graph, 0.9, 0.1, 0.9),  // {0,1} {2,3}
            Column(graph, 0.9, 0.1, 0.9),  // {0,1} {2,3}
        };
        var stack = DenseTStack.Build(graph, temps, cols, theta: 0.5);
        return (graph, stack, cols);
    }

    [Fact]
    public void SelectsThePersistentHalves_NotTheColdGiant()
    {
        var (graph, stack, cols) = SplitStack();

        LineagePersistenceResult result = LineagePersistence.Resolve(
            graph, stack, cols, q: 20, minClusterSize: 2, splitShare: 0.25);

        // Three lineages tracked: the cold giant (span 1) and the two halves
        // (span 2 each). The halves outscore the giant ⇒ gap selects exactly 2.
        Assert.Equal(2, result.Selected.Count);
        Assert.All(result.Selected, l => Assert.True(l.TSpan >= 2.0 - 1e-9));

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
        LineagePersistenceResult result = LineagePersistence.Resolve(
            graph, stack, cols, q: 20, minClusterSize: 2, splitShare: 0.25);

        // The size-4 giant lineage spans only the two cold levels (T=1,2) — it
        // ends at the split rather than chaining into one of the halves.
        var giant = result.AllLineages.Single(l => l.Members.Length == 4);
        Assert.Equal(1.0, giant.TBirth, precision: 12);
        Assert.Equal(2.0, giant.TDeath, precision: 12);
        Assert.Equal(2, giant.LevelCount);
    }

    [Fact]
    public void TemperatureWindow_BoundsTheAnalyzedLevels()
    {
        var (graph, stack, cols) = SplitStack();

        // Window excludes the cold giant levels (T<3): only the split levels
        // remain, so the two halves are tracked and the giant never appears.
        LineagePersistenceResult result = LineagePersistence.Resolve(
            graph, stack, cols, q: 20, minClusterSize: 2, splitShare: 0.25,
            temperatureWindow: (3.0, 5.0));

        Assert.DoesNotContain(result.AllLineages, l => l.Members.Length == 4);
        Assert.Equal(2, result.Selected.Count);
    }
}
