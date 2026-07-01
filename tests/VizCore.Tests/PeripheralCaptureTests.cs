using System;
using Clustering.Graphical.SPC.Partitions.Strategies;
using Clustering.Primitives;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// P2: peripheral capture — Domany1999 step 2 augmentation.
/// Each node is unioned with its single max-affinity neighbor post-threshold,
/// regardless of θ.
/// </summary>
public sealed class PeripheralCaptureTests
{
    // 4-node graph:  0──1──2  3 (isolated by threshold, but connected to 2 with modest affinity)
    //
    //   edge (0,1) weight=1.0,  affinity=0.90  (above θ=0.5)
    //   edge (1,2) weight=1.0,  affinity=0.80  (above θ=0.5)
    //   edge (2,3) weight=1.0,  affinity=0.20  (below θ=0.5)  ← peripheral capture rescues 3
    //
    // Without capture: {0,1,2}, {3}  (3 is unclassified / singleton)
    // With    capture: {0,1,2,3}     (3's best neighbor is 2 at 0.20, merged in)

    private static (CsrGraph graph, Affinities affinities) BuildGraph()
    {
        var edges = new Edge[]
        {
            new(0, 1, 1.0),
            new(1, 2, 1.0),
            new(2, 3, 1.0),
        };
        var graph = CsrGraph.FromEdges(edges, 4);

        // Affinity array parallel to CSR Targets (symmetric — both directions stored).
        // We need to build the G array so that slot for (i,j) and slot for (j,i) match
        // the edge's affinity. Build from the CSR itself.
        var g = new double[graph.Targets.Length];
        double[] edgeAffinities = { 0.90, 0.80, 0.20 };  // indexed by undirected edge order

        // Assign the affinity value to both CSR slots for each undirected edge.
        // Walk the upper-triangle just like SwendsenWang does (target > source).
        var assigned = new bool[graph.Targets.Length];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
        {
            // Map undirected edge to its affinity by source–target pair.
            double aff = (edge.Source, edge.Target) switch
            {
                (0, 1) => 0.90,
                (1, 2) => 0.80,
                (2, 3) => 0.20,
                _      => 0.0,
            };
            g[edge.Slot] = aff;

            // Also set the reverse slot (j→i) — find it in the j row.
            int reverseSlot = FindSlot(graph, edge.Target, edge.Source);
            g[reverseSlot] = aff;
        }

        return (graph, new Affinities { G = g, Temperature = 0.0 });
    }

    private static int FindSlot(CsrGraph graph, int row, int col)
    {
        for (int slot = graph.RowPointers[row]; slot < graph.RowPointers[row + 1]; slot++)
            if (graph.Targets[slot] == col) return slot;
        throw new InvalidOperationException($"Edge ({row},{col}) not in CSR.");
    }

    [Fact]
    public void WithoutPeripheralCapture_PerimeterNodeIsIsolated()
    {
        var (graph, affinities) = BuildGraph();
        var strategy = new ThresholdBondFrequency { Theta = 0.5, PeripheralCapture = false };

        Assignment result = strategy.Apply(graph, affinities, alignments: null);

        // Nodes 0,1,2 should be in the same cluster; node 3 in its own.
        Assert.Equal(2, result.Count);
        int core = result.Labels[0];
        Assert.Equal(core, result.Labels[1]);
        Assert.Equal(core, result.Labels[2]);
        Assert.NotEqual(core, result.Labels[3]);  // 3 isolated
    }

    [Fact]
    public void WithPeripheralCapture_PerimeterNodeJoinsBestNeighbor()
    {
        var (graph, affinities) = BuildGraph();
        var strategy = new ThresholdBondFrequency { Theta = 0.5, PeripheralCapture = true };

        Assignment result = strategy.Apply(graph, affinities, alignments: null);

        // All 4 nodes should be in the same cluster.
        Assert.Equal(1, result.Count);
        int c = result.Labels[0];
        Assert.Equal(c, result.Labels[1]);
        Assert.Equal(c, result.Labels[2]);
        Assert.Equal(c, result.Labels[3]);
    }

    [Fact]
    public void PeripheralCapture_WithTheta0_SameAsWithout()
    {
        // θ=0 already unions everything above 0; peripheral capture shouldn't break this.
        var (graph, affinities) = BuildGraph();
        var strategyOff = new ThresholdBondFrequency { Theta = 0.0, PeripheralCapture = false };
        var strategyOn  = new ThresholdBondFrequency { Theta = 0.0, PeripheralCapture = true  };

        Assignment off = strategyOff.Apply(graph, affinities, alignments: null);
        Assignment on  = strategyOn .Apply(graph, affinities, alignments: null);

        // Both should merge all 4 nodes (every edge above θ=0).
        Assert.Equal(1, off.Count);
        Assert.Equal(1, on.Count);
    }
}
