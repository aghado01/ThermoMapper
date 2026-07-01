using Clustering.Evaluation.External;
using Clustering.Evaluation.Internal;
using Clustering.Primitives;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// T4: every cluster evaluator scores over the <i>assigned</i> subset only.
/// A point carrying <see cref="Assignment.Unassigned"/> (-1) belongs to no
/// cluster and must be excluded — never densified into a spurious cluster,
/// never used to index a per-cluster array. Each family is checked the same
/// way: scoring an input that contains an unassigned point must (a) not crash
/// and (b) equal scoring the same input with the unassigned points physically
/// removed.
/// </summary>
public sealed class EvaluatorUnassignedTests
{
    // ---- External (predicted vs reference) ----------------------------------

    private static IExternalClusterEvaluator[] Externals() => new IExternalClusterEvaluator[]
    {
        new Purity(),
        new NormalizedMutualInformation(),
        new AdjustedRandIndex(),
        new Homogeneity(),
        new Completeness(),
        new VMeasure(),
    };

    [Fact]
    public void External_DropsUnassignedPredicted_MatchesPhysicalRemoval()
    {
        // Point 4 is unassigned in the prediction; its reference label (2) is a
        // class that exists only on that dropped point.
        int[] predFull = { 0, 0, 1, 1, Assignment.Unassigned };
        int[] refFull  = { 0, 0, 1, 1, 2 };
        int[] predSub  = { 0, 0, 1, 1 };
        int[] refSub   = { 0, 0, 1, 1 };

        foreach (var e in Externals())
        {
            double full = e.Evaluate(predFull, refFull);
            double sub  = e.Evaluate(predSub, refSub);
            Assert.Equal(sub, full, 10);
        }
    }

    [Fact]
    public void External_AllUnassigned_DoesNotThrow()
    {
        int[] pred = { Assignment.Unassigned, Assignment.Unassigned, Assignment.Unassigned };
        int[] reference = { 0, 1, 0 };

        foreach (var e in Externals())
        {
            double score = e.Evaluate(pred, reference);
            Assert.False(double.IsInfinity(score));
        }
    }

    // ---- Internal (point / distance) ----------------------------------------

    private static IInternalClusterEvaluator[] Internals() => new IInternalClusterEvaluator[]
    {
        new SilhouetteEvaluator(),
        new CalinskiHarabaszEvaluator(),
        new DaviesBouldinEvaluator(),
    };

    [Fact]
    public void Internal_DropsUnassignedPoints_MatchesPhysicalRemoval()
    {
        // The unassigned point sits far away — were it not dropped it would
        // wreck every centroid/scatter, so equality proves the exclusion.
        double[][] dataFull =
        {
            new[] { 0.0, 0.0 },
            new[] { 0.0, 1.0 },
            new[] { 10.0, 0.0 },
            new[] { 10.0, 1.0 },
            new[] { 100.0, 100.0 },
        };
        int[] labelsFull = { 0, 0, 1, 1, Assignment.Unassigned };

        double[][] dataSub =
        {
            new[] { 0.0, 0.0 },
            new[] { 0.0, 1.0 },
            new[] { 10.0, 0.0 },
            new[] { 10.0, 1.0 },
        };
        int[] labelsSub = { 0, 0, 1, 1 };

        foreach (var e in Internals())
        {
            double full = e.Evaluate(dataFull, labelsFull);
            double sub  = e.Evaluate(dataSub, labelsSub);
            Assert.Equal(sub, full, 10);
        }
    }

    [Fact]
    public void Internal_AllUnassigned_DoesNotThrow()
    {
        double[][] data =
        {
            new[] { 0.0, 0.0 },
            new[] { 1.0, 1.0 },
            new[] { 2.0, 2.0 },
        };
        int[] labels = { Assignment.Unassigned, Assignment.Unassigned, Assignment.Unassigned };

        foreach (var e in Internals())
        {
            double score = e.Evaluate(data, labels);
            Assert.False(double.IsInfinity(score));
        }
    }

    // ---- Graph (weighted partition) -----------------------------------------

    private static IGraphPartitionEvaluator[] Graphs() => new IGraphPartitionEvaluator[]
    {
        new BondModularity(),
        new BondConductance(),
        new BondCoverage(),
    };

    [Fact]
    public void Graph_SkipsEdgesWithUnassignedEndpoint_MatchesPhysicalRemoval()
    {
        // Node 4 is unassigned. Its two edges (0-4, 3-4) each have one assigned
        // and one unassigned endpoint — assigned-edges-only drops them whole,
        // so the assigned endpoint's degree excludes them too. Placing node 4
        // last keeps nodes 0..3 at the same indices in the reduced graph.
        var full = CsrGraph.FromEdges(new[]
        {
            new Edge(0, 1, 1.0),
            new Edge(1, 2, 0.5),
            new Edge(2, 3, 1.0),
            new Edge(0, 4, 0.3),
            new Edge(3, 4, 0.7),
        }, nodeCount: 5);
        int[] labelsFull = { 0, 0, 1, 1, Assignment.Unassigned };

        var sub = CsrGraph.FromEdges(new[]
        {
            new Edge(0, 1, 1.0),
            new Edge(1, 2, 0.5),
            new Edge(2, 3, 1.0),
        }, nodeCount: 4);
        int[] labelsSub = { 0, 0, 1, 1 };

        foreach (var e in Graphs())
        {
            double scoreFull = e.Evaluate(full, full.Weights, labelsFull, clusterCount: 2);
            double scoreSub  = e.Evaluate(sub, sub.Weights, labelsSub, clusterCount: 2);
            Assert.Equal(scoreSub, scoreFull, 10);
        }
    }

    [Fact]
    public void Graph_AllUnassigned_DoesNotThrow()
    {
        var graph = CsrGraph.FromEdges(new[]
        {
            new Edge(0, 1, 1.0),
            new Edge(1, 2, 0.5),
            new Edge(2, 0, 0.7),
        }, nodeCount: 3);
        int[] labels = { Assignment.Unassigned, Assignment.Unassigned, Assignment.Unassigned };

        foreach (var e in Graphs())
        {
            // clusterCount intentionally non-zero to exercise the edge-skip
            // guard rather than the empty-partition early return.
            double score = e.Evaluate(graph, graph.Weights, labels, clusterCount: 2);
            Assert.Equal(0.0, score, 10);
        }
    }
}
