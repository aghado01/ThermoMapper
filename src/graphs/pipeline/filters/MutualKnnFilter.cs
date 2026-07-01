using System;
using Graphs.Neighbors;

namespace Graphs.Pipeline.Filters;

/// <summary>
/// Stage 2 — mutual-KNN symmetrizer (AND-rule). Delegates to
/// <see cref="Symmetrization.MutualIntersection"/>.
/// </summary>
public sealed class MutualKnnFilter : ITopologyFilter
{
    private readonly MutualBandwidthSource _bandwidthSource;

    public MutualKnnFilter(MutualBandwidthSource bandwidthSource = MutualBandwidthSource.DirectedKth)
    {
        _bandwidthSource = bandwidthSource;
    }

    public NeighborSelection Filter(NeighborSelection directed, int n, Func<int, int, double> pairDistance) =>
        Symmetrization.MutualIntersection(directed, n, _bandwidthSource);
}
