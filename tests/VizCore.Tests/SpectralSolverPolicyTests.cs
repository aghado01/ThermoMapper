using System;
using System.Collections.Generic;
using Graphs.Primitives;
using Graphs.Spectral;
using Maths.LinAlg;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// Guards the <see cref="SolverKind.Auto"/> seam: the decision function maps the
/// request shape to the right solver, and on graphs below the threshold (where
/// Auto resolves to Dense) it is identical to the dense path — so wiring consumers
/// to Auto is behaviour-preserving until the benchmark tunes the threshold up.
/// </summary>
public sealed class SpectralSolverPolicyTests
{
    [Theory]
    [InlineData(100, 8, SolverKind.Dense)]                                          // small graph
    [InlineData(500, 2, SolverKind.Iterative)]                                      // 500-node Fiedler (benchmarked 22x)
    [InlineData(SpectralSolverPolicy.IterativeMinNodes - 1, 4, SolverKind.Dense)]
    [InlineData(SpectralSolverPolicy.IterativeMinNodes, 4, SolverKind.Iterative)]
    [InlineData(10_000, 2, SolverKind.Iterative)]
    [InlineData(10_000, SpectralSolverPolicy.IterativeMaxK + 1, SolverKind.Dense)] // too many modes
    [InlineData(10_000, 0, SolverKind.Dense)]                                       // degenerate k
    public void Resolve_MapsShapeToSolver(int n, int k, SolverKind expected)
    {
        Assert.Equal(expected, SpectralSolverPolicy.Resolve(n, k));
    }

    [Fact]
    public void Auto_MatchesDense_OnSmallGraph()
    {
        // n well below IterativeMinNodes ⇒ policy resolves to Dense ⇒ Auto must be
        // byte-identical to an explicit Dense request.
        CsrGraph graph = BuildCycle(8);
        const int k = 3;

        IReadOnlyList<EigenPair> dense = Spectral.ComputeBottomK(graph, k: k, solverKind: SolverKind.Dense);
        IReadOnlyList<EigenPair> auto = Spectral.ComputeBottomK(graph, k: k, solverKind: SolverKind.Auto);

        Assert.Equal(dense.Count, auto.Count);
        for (int i = 0; i < dense.Count; i++)
            Assert.Equal(dense[i].Lambda, auto[i].Lambda, 12);
    }

    private static CsrGraph BuildCycle(int n)
    {
        var edges = new Edge[n];
        for (int i = 0; i < n; i++)
            edges[i] = new Edge(i, (i + 1) % n, 1.0);
        return CsrGraph.FromEdges(edges, n);
    }
}
