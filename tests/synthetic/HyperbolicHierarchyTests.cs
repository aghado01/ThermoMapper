using System;
using Synthetic.Manifolds;
using Xunit;

namespace Synthetic.Tests;

public sealed class HyperbolicHierarchyTests
{
    [Fact]
    public void Generate_AllPointsInsideUnitBallAndFinite()
    {
        // The migration's core guarantee: points are real Poincaré-ball points
        // (the old hand-rolled map could drift; the exp map clamps inside).
        var ds = HyperbolicHierarchy.Generate(nPoints: 1500, seed: 3);
        foreach (var p in ds.Features)
        {
            double r2 = 0.0;
            foreach (var c in p)
            {
                Assert.True(double.IsFinite(c));
                r2 += c * c;
            }
            Assert.True(Math.Sqrt(r2) < 1.0, $"point at ball-radius {Math.Sqrt(r2)} escaped the unit ball");
        }
    }

    [Fact]
    public void Generate_PreservesHierarchyLevelsAndCounts()
    {
        var ds = HyperbolicHierarchy.Generate(nPoints: 1500, hierarchyDepth: 3, seed: 3);
        Assert.Equal(1500, ds.Features.Length);
        Assert.NotNull(ds.LabelsByLevel);
        Assert.Equal(3, ds.LabelsByLevel!.Length);
        Assert.True(ds.ClusterCount > 1);
        foreach (var level in ds.LabelsByLevel)
            Assert.Equal(1500, level.Length);
    }

    [Fact]
    public void Generate_IsDeterministicForSameSeed()
    {
        var a = HyperbolicHierarchy.Generate(nPoints: 800, seed: 11);
        var b = HyperbolicHierarchy.Generate(nPoints: 800, seed: 11);
        Assert.Equal(a.Features.Length, b.Features.Length);
        for (int i = 0; i < a.Features.Length; i++)
            for (int d = 0; d < a.Features[i].Length; d++)
                Assert.Equal(a.Features[i][d], b.Features[i][d], 12);
    }
}
