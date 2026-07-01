using System;
using Graphs.Primitives;

namespace Graphs
{
    public readonly struct NeighborSelection
    {
        public NeighborSelection(
            Neighbor[][] allNeighbors,
            double[] nearestNeighborDistances,
            double[] kthNeighborDistances)
        {
            AllNeighbors = allNeighbors ?? throw new ArgumentNullException(nameof(allNeighbors));
            NearestNeighborDistances = nearestNeighborDistances ?? throw new ArgumentNullException(nameof(nearestNeighborDistances));
            KthNeighborDistances = kthNeighborDistances ?? throw new ArgumentNullException(nameof(kthNeighborDistances));

            int n = allNeighbors.Length;
            if (nearestNeighborDistances.Length != n)
            {
                throw new ArgumentException(
                    "NearestNeighborDistances length must match AllNeighbors row count.",
                    nameof(nearestNeighborDistances));
            }
            if (kthNeighborDistances.Length != n)
            {
                throw new ArgumentException(
                    "KthNeighborDistances length must match AllNeighbors row count.",
                    nameof(kthNeighborDistances));
            }
        }

        public Neighbor[][] AllNeighbors { get; }

        // 1NN distance per node, sampled from the symmetrized neighbor list.
        // Diagnostics that genuinely want the closest-edge view (NeighborhoodScale's
        // median 1NN) read this.
        public double[] NearestNeighborDistances { get; }

        // K-th neighbor distance per node from the directed KNN heap, before
        // mutual pruning or OR-symmetrization. This is the sample bandwidth
        // estimation should use — it matches the "actual kNN search radius"
        // a Gaussian kernel is supposed to scale to. For rules without a true
        // k (EpsilonBall), this carries the per-node furthest-included distance,
        // i.e. the local ball radius.
        public double[] KthNeighborDistances { get; }

        public void Deconstruct(
            out Neighbor[][] allNeighbors,
            out double[] nearestNeighborDistances,
            out double[] kthNeighborDistances)
        {
            allNeighbors = AllNeighbors;
            nearestNeighborDistances = NearestNeighborDistances;
            kthNeighborDistances = KthNeighborDistances;
        }
    }
}
