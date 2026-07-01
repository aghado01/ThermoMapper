using System;
using Graphs.Neighbors;
using Graphs.Pipeline.Filters;

namespace Graphs.Proximity;

public static partial class ProximityGraph
{
    /// <summary>
    /// Standard KNN with OR-symmetrization — thin wrapper over
    /// <see cref="DirectedKnn"/> + <see cref="PassThroughFilter"/>.
    /// </summary>
    public static NeighborSelection SelectKnn(int n, int k, Func<int, int, double> pairDistance)
    {
        NeighborSelection directed = DirectedKnn.Select(n, k, pairDistance);
        return new PassThroughFilter().Filter(directed, n, pairDistance);
    }

    /// <summary>
    /// Mutual KNN (AND-rule) — thin wrapper over
    /// <see cref="DirectedKnn"/> + <see cref="MutualKnnFilter"/>.
    /// </summary>
    public static NeighborSelection SelectMutualKnn(
        int n,
        int k,
        Func<int, int, double> pairDistance,
        MutualBandwidthSource mutualBandwidthSource = MutualBandwidthSource.DirectedKth)
    {
        NeighborSelection directed = DirectedKnn.Select(n, k, pairDistance);
        return new MutualKnnFilter(mutualBandwidthSource).Filter(directed, n, pairDistance);
    }
}
