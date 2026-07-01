using Clustering.Primitives;
using Clustering.Graphical.SPC.Runtime.Core.Solver;
using Graphs.Primitives;
using Xunit;

using Clustering.Graphical.SPC.Runtime.Core;

using Graphs;

namespace Clustering.Graphical.SPC.Tests.Sampler;

/// <summary>
/// The focused PKWang sweep backend: consistency with the single-temperature
/// path, monotone fragmentation over an ascending grid, the geometric grid
/// helper, and traced provenance on the result.
/// </summary>
public sealed class PKWangSweepTests
{
    // Triangle {0,1,2} + pair {3,4,5} joined by a weak bridge (2,3); connected,
    // distinct couplings.
    private static CsrGraph BuildGraph()
    {
        var edges = new[]
        {
            new Edge(0, 1, 10.0),
            new Edge(1, 2, 9.0),
            new Edge(0, 2, 5.0),
            new Edge(3, 4, 8.0),
            new Edge(4, 5, 7.0),
            new Edge(2, 3, 1.0),
        };
        return CsrGraph.FromEdges(edges, nodeCount: 6);
    }

    [Fact]
    public void Sweep_MatchesPerTemperatureCluster()
    {
        CsrGraph g = BuildGraph();
        double[] grid = { 0.7, 4.5, 14.0, 36.0, 72.0 };

        PKWangSweepResult result = PKWangSweep.Run(g, EdgeWeightKind.Coupling, grid, Field.Mean);

        PKWangContext ctx = PKWang.Prepare(g, EdgeWeightKind.Coupling, Field.Mean);
        for (int i = 0; i < grid.Length; i++)
        {
            Assignment expected = PKWang.Cluster(ctx, grid[i]);
            Assert.Equal(expected.Count, result.Partitions[i].Count);
            Assert.Equal(expected.Labels, result.Partitions[i].Labels);
        }
        Assert.Equal(grid, result.Temperatures);
    }

    [Fact]
    public void ClusterCounts_MonotoneOverAscendingGrid()
    {
        CsrGraph g = BuildGraph();
        double[] grid = new[] { 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 100.0 };

        int[] counts = PKWangSweep.Run(g, EdgeWeightKind.Coupling, grid, Field.Local, SymmetrizationRule.Mutual).ClusterCounts();

        for (int i = 1; i < counts.Length; i++)
            Assert.True(counts[i] >= counts[i - 1], $"Cluster count fell at index {i} (T={grid[i]}).");
    }

    [Fact]
    public void Run_RecordsProvenance()
    {
        CsrGraph g = BuildGraph();

        PKWangSweepResult result = PKWangSweep.Run(
            g, EdgeWeightKind.Coupling, new[] { 5.0 }, Field.Local, SymmetrizationRule.Inclusive, theta: 0.4);

        Assert.Equal(Field.Local, result.Field);
        Assert.Equal(SymmetrizationRule.Inclusive, result.Symmetrization);
        Assert.Equal(0.4, result.Theta);
        Assert.Single(result.Partitions);
    }
}
