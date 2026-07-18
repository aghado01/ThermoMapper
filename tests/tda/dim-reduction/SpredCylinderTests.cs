using System;
using Graphs;
using Maths.Geometry.DimReduction;
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
        Dimensions = [(0, 0.5), (1, 1.0)],   // paper's multi-order idea, weights deliberately unnormalized: H0 smooth descent + full-strength H1 loop target
        MaxDimension = 2,
        MinPersistence = 0.05,               // prune the ~4n near-diagonal H1 noise loops (P0/P1): huge speedup, denoises
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

    // The forwarding facts live on this fixture because the anneal moves here: on a fixture whose
    // PCA warm start is already optimal, best-so-far tracking returns the warm start under any
    // options, and a dropped annealerOptions parameter would be undetectable.

    /// <summary>Spred.Compute is exactly objective-construction + engine call, options included:
    /// same seed and options must reproduce the direct engine composition bit-for-bit.</summary>
    [Fact]
    public void Compute_AnnealerOptions_ReachTheEngine()
    {
        double[][] data = Cylinder3D(n: 70, seed: 3);
        PersistenceObjectiveConfig config = CylinderConfig();
        var options = new SubspaceAnnealerOptions { IsotropicFraction = 1.0, InitialStep = 0.4 };

        double[][] driver = Spred.Compute(data, targetDim: 2, config, maxIters: 40, seed: 7, options);
        var ph = new PersistenceObjective(data, config);
        double[][] engine = SubspaceAnnealer.Compute(
            data, targetDim: 2, ph.Evaluate, maxIters: 40, seed: 7, options).Projection;
        double[][] defaults = Spred.Compute(data, targetDim: 2, config, maxIters: 40, seed: 7);

        AssertBitIdentical(driver, engine);
        Assert.False(BitIdentical(driver, defaults),
            "Non-default options must alter the seeded trajectory here, or this fact cannot detect a dropped parameter.");
    }

    /// <summary>RunBlock (the diagnostics path shared by every block count) forwards the options to
    /// the engine and reports the engine's own tracked objective.</summary>
    [Fact]
    public void ComputeWithDiagnostics_AnnealerOptions_ReachTheBlockEngine()
    {
        double[][] data = Cylinder3D(n: 70, seed: 3);
        PersistenceObjectiveConfig config = CylinderConfig();
        var options = new SubspaceAnnealerOptions { IsotropicFraction = 1.0, InitialStep = 0.4 };

        DistributedSpredResult run = DistributedSpred.ComputeWithDiagnostics(
            data, targetDim: 2, blockCount: 1, config, maxIters: 40, seed: 7,
            maxDegreeOfParallelism: 1, options);
        var ph = new PersistenceObjective(data, config);
        SubspaceAnnealerResult engine = SubspaceAnnealer.Compute(
            data, targetDim: 2, ph.Evaluate, maxIters: 40, seed: 7, options);

        AssertBitIdentical(run.Projection, engine.Projection);
        Assert.Equal(engine.Objective, run.FullDataObjective);   // single block: local == full == tracked value
    }

    private static void AssertBitIdentical(double[][] a, double[][] b)
    {
        Assert.True(BitIdentical(a, b), "Projections should be bit-identical.");
    }

    private static bool BitIdentical(double[][] a, double[][] b)
    {
        if (a.Length != b.Length) return false;
        for (int i = 0; i < a.Length; i++)
        {
            if (a[i].Length != b[i].Length) return false;
            for (int j = 0; j < a[i].Length; j++)
                if (!a[i][j].Equals(b[i][j])) return false;
        }
        return true;
    }
}
