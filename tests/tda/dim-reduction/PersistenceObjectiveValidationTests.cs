using System;
using Graphs;
using Xunit;

namespace TDA.DimReduction.Tests;

/// <summary>
/// Rejection facts for <see cref="PersistenceObjective"/> construction — one per validation rule.
/// Every entry point (Spred, DistributedSpred, direct use) constructs through here, so a
/// silently-degenerate recipe (empty diagram, constant objective term, NaN-poisoned weight) fails
/// fast instead of handing the annealer a flat landscape.
/// </summary>
public sealed class PersistenceObjectiveValidationTests
{
    [Fact]
    public void Construct_EmptyDimensions_Throws()
    {
        ArgumentException error = Assert.Throws<ArgumentException>(() =>
            new PersistenceObjective(Circle3D(8), BaseConfig() with { Dimensions = [] }));

        Assert.Equal("config", error.ParamName);
    }

    [Theory]
    [InlineData(-1)]
    [InlineData(2)]
    [InlineData(3)]
    public void Construct_DimensionOutsideRipsRange_Throws(int dim)
    {
        ArgumentException error = Assert.Throws<ArgumentException>(() =>
            new PersistenceObjective(
                Circle3D(8),
                BaseConfig() with { Dimensions = [(dim, 1.0)], MaxDimension = 2 }));

        Assert.Equal("config", error.ParamName);
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(-0.5)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    public void Construct_NonPositiveOrNonFiniteWeight_Throws(double weight)
    {
        ArgumentException error = Assert.Throws<ArgumentException>(() =>
            new PersistenceObjective(
                Circle3D(8),
                BaseConfig() with { Dimensions = [(1, weight)] }));

        Assert.Equal("config", error.ParamName);
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(-2.0)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    public void Construct_NonPositiveOrNonFiniteWassersteinOrder_Throws(double order)
    {
        ArgumentException error = Assert.Throws<ArgumentException>(() =>
            new PersistenceObjective(Circle3D(8), BaseConfig() with { WassersteinOrder = order }));

        Assert.Equal("config", error.ParamName);
    }

    [Theory]
    [InlineData(-0.1)]
    [InlineData(double.NaN)]
    public void Construct_NegativeMinPersistence_Throws(double minPersistence)
    {
        ArgumentException error = Assert.Throws<ArgumentException>(() =>
            new PersistenceObjective(Circle3D(8), BaseConfig() with { MinPersistence = minPersistence }));

        Assert.Equal("config", error.ParamName);
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(-1e6)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    public void Construct_NonPositiveOrNonFinitePathologyPenalty_Throws(double penalty)
    {
        ArgumentException error = Assert.Throws<ArgumentException>(() =>
            new PersistenceObjective(Circle3D(8), BaseConfig() with { PathologyPenalty = penalty }));

        Assert.Equal("config", error.ParamName);
    }

    [Theory]
    [InlineData(0)]
    [InlineData(-1)]
    public void Construct_MaxDimensionBelowOne_Throws(int maxDimension)
    {
        ArgumentException error = Assert.Throws<ArgumentException>(() =>
            new PersistenceObjective(Circle3D(8), BaseConfig() with { MaxDimension = maxDimension }));

        Assert.Equal("config", error.ParamName);
    }

    [Fact]
    public void Construct_RaggedInput_Throws()
    {
        double[][] data =
        [
            [1.0, 2.0, 3.0],
            [4.0, 5.0],
        ];

        ArgumentException error = Assert.Throws<ArgumentException>(() =>
            new PersistenceObjective(data, BaseConfig()));

        Assert.Equal("data", error.ParamName);
    }

    [Fact]
    public void Compute_RaggedInput_Throws()
    {
        double[][] data =
        [
            [1.0, 2.0, 3.0],
            [4.0, 5.0],
        ];

        ArgumentException error = Assert.Throws<ArgumentException>(() =>
            Spred.Compute(data, targetDim: 2, BaseConfig(), maxIters: 0, seed: 1));

        Assert.Equal("data", error.ParamName);
    }

    private static PersistenceObjectiveConfig BaseConfig() => new()
    {
        Graph = new GraphCompilerConfig
        {
            Topology = new TopologyConfig { Kind = TopologyKind.Knn, K = 6 },
            Filter = new FilterConfig { Kind = FilterKind.OrRule },
            Repair = new RepairConfig { Kind = RepairKind.NoRepair },
            Projection = new DistanceProjection(),
        },
        Dimensions = [(1, 1.0)],
        MaxDimension = 2,
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
}
