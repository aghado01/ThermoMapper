using System;
using Graphs;
using Xunit;
using Xunit.Abstractions;

namespace TDA.DimReduction.Tests;

/// <summary>
/// The paper's Figure-1 fixture (Yu &amp; You, 2106.02096 §5): a noisy cylinder S¹×[-2,2] (β₀=1, β₁=1),
/// exercising the **annealer's search** from a loop-losing warm start (PCA's top-2 aligns with the
/// high-variance [-2,2] axis and flattens the loop).
///
/// <para><b>Finding (2026-07-03):</b> the H0+H1 objective preserves *full* persistence, so the naive
/// circle projection (x,y) is <i>not</i> the optimum — projecting the cylinder onto its circle piles
/// every height onto the same θ and destroys the H0 (height) merge structure, scoring it *worse*
/// (~5.8) than even the loop-flattening (h,x) view (~4.5). SPRED instead finds an <i>oblique</i>
/// projection (~1.0) that keeps the loop AND the height spread — better shape preservation than either
/// axis-aligned choice. This test asserts that search wins, not a naive "recover the visible loop".</para>
/// </summary>
public sealed class SpredCylinderTests
{
    private readonly ITestOutputHelper _out;
    public SpredCylinderTests(ITestOutputHelper output) => _out = output;

    private static double[][] Cylinder3D(int n, int seed)
    {
        var rng = new Random(seed);
        var pts = new double[n][];
        for (int i = 0; i < n; i++)
        {
            double t = 2.0 * Math.PI * rng.NextDouble();
            double h = 4.0 * rng.NextDouble() - 2.0;
            pts[i] = new[]
            {
                Math.Cos(t) + 0.03 * (rng.NextDouble() - 0.5),
                Math.Sin(t) + 0.03 * (rng.NextDouble() - 0.5),
                h           + 0.03 * (rng.NextDouble() - 0.5),
            };
        }
        return pts;
    }

    private static PersistenceObjectiveConfig CylinderConfig() => new()
    {
        Graph = new GraphCompilerConfig
        {
            Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 10 },
            Filter = new FilterConfig { Kind = FilterKind.OrRule },
            Repair = new RepairConfig { Kind = RepairKind.NoRepair },
            Projection = new DistanceProjection(),
        },
        Dimensions = [(0, 0.5), (1, 1.0)],   // H0 smooth descent + H1 loop target (paper runs orders 0 and 1)
        MaxDimension = 2,
    };

    [Fact]
    public void Compute_BeatsEveryAxisAlignedProjection_OnCylinder()
    {
        double[][] data = Cylinder3D(n: 70, seed: 3);
        var objective = new PersistenceObjective(data, CylinderConfig());

        double loopLosing = objective.Evaluate(new[] { new[] { 0.0, 0.0, 1.0 }, new[] { 1.0, 0.0, 0.0 } }); // (h,x)
        double loopCircle = objective.Evaluate(new[] { new[] { 1.0, 0.0, 0.0 }, new[] { 0.0, 1.0, 0.0 } }); // (x,y)

        double[][] proj = Spred.Compute(data, targetDim: 2, CylinderConfig(), maxIters: 250, seed: 7);
        double annealed = objective.Evaluate(proj);

        _out.WriteLine($"(h,x)={loopLosing:F4}  (x,y)={loopCircle:F4}  SPRED={annealed:F4}");

        // The anneal finds a projection with strictly better full-persistence fidelity than either
        // axis-aligned view — the search does real work from the loop-losing PCA warm start.
        double bestAxisAligned = Math.Min(loopLosing, loopCircle);
        Assert.True(annealed < 0.75 * bestAxisAligned,
            $"SPRED {annealed:F4} should beat the best axis-aligned projection {bestAxisAligned:F4} " +
            $"(h,x={loopLosing:F4}, x,y={loopCircle:F4}).");
    }
}
