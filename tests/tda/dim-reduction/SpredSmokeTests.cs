using System;
using Graphs;
using Xunit;

namespace TDA.DimReduction.Tests;

/// <summary>
/// Smoke tests for the SPRED objective + driver on a fixture with a known loop (a unit circle in the
/// xy-plane of R³: β₀=1, β₁=1). A projection onto the circle's plane preserves the H1 loop; a
/// projection onto a plane containing the "flat" axis collapses it to a segment (H1 gone). The
/// objective should score the loop-preserving projection strictly lower — i.e. it does topological
/// work, not just compile.
/// </summary>
public sealed class SpredSmokeTests
{
    private const int N = 120;

    // Unit circle in the xy-plane of R³. Rips over its kNN graph carries one (essential) H1 loop.
    private static double[][] Circle3D(int n)
    {
        var pts = new double[n][];
        for (int i = 0; i < n; i++)
        {
            double t = 2.0 * Math.PI * i / n;
            pts[i] = new[] { Math.Cos(t), Math.Sin(t), 0.0 };
        }
        return pts;
    }

    // Simple, fast, deterministic recipe: kNN + OR-rule, no repair, distance graph. Pure-H1 objective.
    private static PersistenceObjectiveConfig H1Config() => new()
    {
        Graph = new GraphCompilerConfig
        {
            Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 10 },
            Filter = new FilterConfig { Kind = FilterKind.OrRule },
            Repair = new RepairConfig { Kind = RepairKind.NoRepair },
            Projection = new DistanceProjection(),
        },
        Dimensions = [(1, 1.0)],
        MaxDimension = 2,
    };

    // Keeps the loop: project onto the xy-plane (the circle's own plane).
    private static double[][] LoopPreserving() =>
        new[] { new[] { 1.0, 0.0, 0.0 }, new[] { 0.0, 1.0, 0.0 } };

    // Collapses the loop: project onto (z, x) — the circle flattens onto a segment.
    private static double[][] LoopCollapsing() =>
        new[] { new[] { 0.0, 0.0, 1.0 }, new[] { 1.0, 0.0, 0.0 } };

    [Fact]
    public void Objective_RewardsLoopPreservingProjection()
    {
        var objective = new PersistenceObjective(Circle3D(N), H1Config());

        double preserving = objective.Evaluate(LoopPreserving());
        double collapsing = objective.Evaluate(LoopCollapsing());

        Assert.True(preserving < collapsing,
            $"loop-preserving projection ({preserving}) should score below loop-collapsing ({collapsing}).");
    }

    [Fact]
    public void Compute_RunsEndToEnd_ReturnsOrthonormalLoopPreservingProjection()
    {
        double[][] data = Circle3D(N);
        PersistenceObjectiveConfig config = H1Config();

        double[][] proj = Spred.Compute(data, targetDim: 2, config, maxIters: 60, seed: 7);

        Assert.Equal(2, proj.Length);
        for (int i = 0; i < proj.Length; i++)
        {
            double selfDot = 0.0;
            for (int j = 0; j < proj[i].Length; j++) selfDot += proj[i][j] * proj[i][j];
            Assert.InRange(selfDot, 1.0 - 1e-9, 1.0 + 1e-9);

            for (int k = i + 1; k < proj.Length; k++)
            {
                double crossDot = 0.0;
                for (int j = 0; j < proj[i].Length; j++) crossDot += proj[i][j] * proj[k][j];
                Assert.InRange(crossDot, -1e-9, 1e-9);
            }
        }

        // PCA warm-start already lands on the circle's plane, and the anneal only improves; so the
        // result preserves the loop better than a deliberately-collapsed projection.
        var objective = new PersistenceObjective(data, config);
        Assert.True(objective.Evaluate(proj) < objective.Evaluate(LoopCollapsing()));
    }
}
