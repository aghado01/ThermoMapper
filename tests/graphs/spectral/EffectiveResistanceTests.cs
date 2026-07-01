#nullable enable
using System;
using System.Collections.Generic;
using Graphs.Primitives;
using Graphs.Spectral;
using Maths.LinAlg;
using Xunit;

namespace Graphs.Spectral.Tests;

public sealed class EffectiveResistanceTests
{
    [Fact]
    public void PathGraph_UnitConductance_ResistanceDistanceIsTwo()
    {
        var g = CsrGraph.FromEdges(new[]
        {
            new Edge(0, 1, 1.0),
            new Edge(1, 2, 1.0),
        }, nodeCount: 3);

        IReadOnlyList<EigenPair> pairs = EffectiveResistance.ComputeEigenpairs(
            g, tailEpsilon: 1e-6, kMax: 8, solverKind: SolverKind.Dense);
        int[] sameComponent = { 0, 0, 0 };

        double r02 = EffectiveResistance.Pair(0, 2, pairs, sameComponent);
        Assert.Equal(2.0, r02, precision: 6);
    }

    [Fact]
    public void CrossComponentPair_ReturnsInfinity()
    {
        var g = CsrGraph.FromEdges(new[]
        {
            new Edge(0, 1, 1.0),
            new Edge(2, 3, 1.0),
        }, nodeCount: 4);

        IReadOnlyList<EigenPair> pairs = EffectiveResistance.ComputeEigenpairs(
            g, kMax: 8, solverKind: SolverKind.Dense);
        int[] components = { 0, 0, 1, 1 };

        double r = EffectiveResistance.Pair(0, 3, pairs, components);
        Assert.True(double.IsPositiveInfinity(r));
    }

    [Fact]
    public void BuildEdgeWeights_SameComponentEdges_AreFiniteAndPositive()
    {
        var g = CsrGraph.FromEdges(new[]
        {
            new Edge(0, 1, 1.0),
            new Edge(1, 2, 1.0),
        }, nodeCount: 3);

        IReadOnlyList<EigenPair> pairs = EffectiveResistance.ComputeEigenpairs(
            g, kMax: 8, solverKind: SolverKind.Dense);
        var weights = EffectiveResistance.BuildEdgeWeights(g, pairs);

        Assert.Equal(2, weights.Count);
        foreach (double w in weights.Values)
        {
            Assert.False(double.IsPositiveInfinity(w));
            Assert.True(w > 0.0);
        }
    }
}
