using System;
using Graphs.Proximity;

namespace Graphs.Pipeline.Generators;

/// <summary>
/// Stage 1 — epsilon-ball generator. Edge (i, j) exists iff
/// <c>dist(i, j) &lt; epsilon</c>. Intrinsically symmetric given a
/// symmetric distance, so the downstream filter (PassThrough is the
/// natural choice) is a no-op.
/// </summary>
/// <remarks>
/// <para>Produces variable-degree graphs: dense regions get more
/// neighbors, sparse regions get fewer. Fragile in high dimensions
/// where distance concentration makes ε-selection brittle —
/// <see cref="KnnGenerator"/> is usually the right choice there.</para>
/// </remarks>
public sealed class EpsilonBallGenerator : ITopologyGenerator
{
    private readonly double _epsilon;

    public EpsilonBallGenerator(double epsilon)
    {
        if (epsilon <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(epsilon), "Epsilon must be positive.");
        _epsilon = epsilon;
    }

    public NeighborSelection Generate(int n, Func<int, int, double> dist)
        => ProximityGraph.SelectEpsilonBall(n, _epsilon, dist);
}
