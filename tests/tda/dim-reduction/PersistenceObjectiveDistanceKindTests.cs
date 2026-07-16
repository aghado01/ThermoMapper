using System;
using Graphs;
using TDA.Ph;
using Xunit;

namespace TDA.DimReduction.Tests;

/// <summary>
/// The DiagramDistance backend selector on the SPRED objective: every backend must preserve the
/// screening ordering (a projection preserving the data's topology scores strictly better than
/// one collapsing it) — the property the ISOLET "H0 matching-cost gate" needs from an
/// approximate distance — and the new config knobs must reject degenerate values at construction.
/// </summary>
public sealed class PersistenceObjectiveDistanceKindTests
{
    [Theory]
    [InlineData(DiagramDistanceKind.Wasserstein)]
    [InlineData(DiagramDistanceKind.SlicedWasserstein)]
    [InlineData(DiagramDistanceKind.SinkhornWasserstein)]
    public void Evaluate_RanksTopologyPreservingProjectionOverCollapsed(DiagramDistanceKind kind)
    {
        var objective = new PersistenceObjective(Circle3D(32), Config(kind));

        double clean = objective.Evaluate(XyPlane());      // keeps the circle
        double collapsed = objective.Evaluate(XzPlane());  // flattens it to a segment

        Assert.True(double.IsFinite(clean));
        Assert.True(double.IsFinite(collapsed));
        Assert.True(clean < collapsed,
            $"{kind}: clean projection scored {clean}, collapsed scored {collapsed}.");
    }

    [Theory]
    [InlineData(DiagramDistanceKind.SlicedWasserstein)]
    [InlineData(DiagramDistanceKind.SinkhornWasserstein)]
    public void Evaluate_ApproximateBackends_Deterministic(DiagramDistanceKind kind)
    {
        var objective = new PersistenceObjective(Circle3D(24), Config(kind));

        Assert.Equal(objective.Evaluate(XyPlane()), objective.Evaluate(XyPlane()), precision: 12);
    }

    [Fact]
    public void Construct_WassersteinOrderBelowOne_Throws()
    {
        // p in (0, 1) previously passed construction and threw deep inside evaluation
        // (DiagramMetrics requires p >= 1) — now rejected up front.
        ArgumentException error = Assert.Throws<ArgumentException>(() =>
            new PersistenceObjective(Circle3D(8), Config(DiagramDistanceKind.Wasserstein) with
            {
                WassersteinOrder = 0.5,
            }));
        Assert.Equal("config", error.ParamName);
    }

    [Fact]
    public void Construct_NonPositiveSlicedDirections_Throws()
    {
        ArgumentException error = Assert.Throws<ArgumentException>(() =>
            new PersistenceObjective(Circle3D(8), Config(DiagramDistanceKind.SlicedWasserstein) with
            {
                SlicedDirections = 0,
            }));
        Assert.Equal("config", error.ParamName);
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(-1.0)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    public void Construct_BadSinkhornEpsilon_Throws(double epsilon)
    {
        ArgumentException error = Assert.Throws<ArgumentException>(() =>
            new PersistenceObjective(Circle3D(8), Config(DiagramDistanceKind.SinkhornWasserstein) with
            {
                SinkhornEpsilon = epsilon,
            }));
        Assert.Equal("config", error.ParamName);
    }

    [Fact]
    public void Construct_NonPositiveSinkhornMaxIters_Throws()
    {
        ArgumentException error = Assert.Throws<ArgumentException>(() =>
            new PersistenceObjective(Circle3D(8), Config(DiagramDistanceKind.SinkhornWasserstein) with
            {
                SinkhornMaxIters = 0,
            }));
        Assert.Equal("config", error.ParamName);
    }

    private static PersistenceObjectiveConfig Config(DiagramDistanceKind kind) => new()
    {
        Graph = new GraphCompilerConfig
        {
            Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 6 },
            Filter = new FilterConfig { Kind = FilterKind.OrRule },
            Repair = new RepairConfig { Kind = RepairKind.NoRepair },
            Projection = new DistanceProjection(),
        },
        Dimensions = [(0, 0.5), (1, 0.5)],
        MaxDimension = 2,
        DiagramDistance = kind,
    };

    private static double[][] Circle3D(int n)
    {
        var pts = new double[n][];
        for (int i = 0; i < n; i++)
        {
            double t = 2.0 * Math.PI * i / n;
            pts[i] = [Math.Cos(t), Math.Sin(t), 0.0];
        }
        return pts;
    }

    private static double[][] XyPlane() =>
    [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ];

    private static double[][] XzPlane() =>
    [
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
}
