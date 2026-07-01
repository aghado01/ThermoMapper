using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Graphs;
using Graphs.Primitives;

namespace Graphs.Proximity
{
    public static partial class ProximityGraph
    {
        /// <summary>
        /// Epsilon-ball proximity.
        /// Edge (i,j) exists if d(i,j) &lt; epsilon.
        /// Inherently symmetric given a symmetric distance function.
        /// K parameter is not used — all pairs within radius are included.
        /// Produces a variable-degree graph: dense regions get more edges,
        /// sparse regions get fewer. Can fail in high dimensions where
        /// distance concentration makes epsilon selection fragile.
        ///
        /// Returns per-node neighbor lists and 1-NN distances for auto-delta.
        /// </summary>
        public static NeighborSelection SelectEpsilonBall(
            int n, double epsilon, Func<int, int, double> pairDistance)
        {
            if (epsilon <= 0)
                throw new ArgumentException("Epsilon must be positive for EpsilonBall proximity.");

            // Phase 1: parallel upper-triangle scan — each row i writes only to halfLists[i].
            // No cross-row contention; j's symmetric entry is added in Phase 2.
            var halfLists = new List<Neighbor>[n];
            for (int i = 0; i < n; i++)
                halfLists[i] = new List<Neighbor>();

            Parallel.For(0, n, i =>
            {
                for (int j = i + 1; j < n; j++)
                {
                    double d = pairDistance(i, j);
                    if (d < epsilon)
                        halfLists[i].Add(new Neighbor { Index = j, Distance = d });
                }
            });

            // Phase 2: sequential symmetrization — O(E), no contention.
            var lists = new List<Neighbor>[n];
            for (int i = 0; i < n; i++)
                lists[i] = new List<Neighbor>();

            for (int i = 0; i < n; i++)
            {
                foreach (var nb in halfLists[i])
                {
                    lists[i].Add(nb);
                    lists[nb.Index].Add(new Neighbor { Index = i, Distance = nb.Distance });
                }
            }

            var allNeighbors = new Neighbor[n][];
            var nnDistances = new List<double>(n);
            var kthDistances = new List<double>(n);

            for (int i = 0; i < n; i++)
            {
                lists[i].Sort((a, b) => a.Distance.CompareTo(b.Distance));
                allNeighbors[i] = lists[i].ToArray();

                if (allNeighbors[i].Length > 0)
                {
                    nnDistances.Add(allNeighbors[i][0].Distance);
                    // Per-node furthest-included distance — the analog of the
                    // k-th neighbor distance for an epsilon-ball rule.
                    kthDistances.Add(allNeighbors[i][allNeighbors[i].Length - 1].Distance);
                }
            }

            return new NeighborSelection(allNeighbors, nnDistances.ToArray(), kthDistances.ToArray());
        }
    }
}
