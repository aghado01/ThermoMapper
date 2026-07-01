using System;
using Graphs.Neighbors;

namespace Graphs.Pipeline.Generators;

/// <summary>
/// Stage 1 — directed K-nearest-neighbor generator. Delegates to
/// <see cref="DirectedKnn.Select"/> and returns candidates before
/// symmetrization.
/// </summary>
public sealed class KnnGenerator : ITopologyGenerator
{
    private readonly int _k;

    public KnnGenerator(int k)
    {
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");
        _k = k;
    }

    public NeighborSelection Generate(int n, Func<int, int, double> dist) =>
        DirectedKnn.Select(n, _k, dist);
}
