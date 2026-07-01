using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Graphs;
using Graphs.Primitives;

namespace Graphs.Neighbors;

/// <summary>
/// Canonical directed k-nearest-neighbor pass (pre-symmetrization).
/// Used by <see cref="Pipeline.Generators.KnnGenerator"/> on the production
/// path and by the <see cref="Proximity.ProximityGraph"/> KNN convenience
/// wrappers that the test/smoke harnesses build selections with.
/// </summary>
public static class DirectedKnn
{
    public static NeighborSelection Select(int n, int k, Func<int, int, double> dist)
    {
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");
        if (dist is null)
            throw new ArgumentNullException(nameof(dist));

        var directedNeighbors = new Neighbor[n][];
        var kthDistances = new double[n];
        var nnDistances  = new double[n];

        Parallel.For(0, n, i =>
        {
            var heap = new BoundedMinHeap(k);
            for (int j = 0; j < n; j++)
            {
                if (i == j) continue;
                heap.TryAdd(j, dist(i, j));
            }

            Neighbor[] row = heap.GetSorted();

            // Tie-inclusive rank boundary: when candidates tie the K-th
            // distance, a strict top-K keeps an arbitrary subset and the
            // asymmetry breaks mutual-KNN symmetrization at the boundary
            // (quantized data — e.g. Iris' 0.1-step grid — fragments into
            // micro-components). Include EVERY candidate at distance ≤ the
            // K-th, ordered (distance, index) for determinism. No-op when
            // the boundary is tie-free.
            if (row.Length == k && k < n - 1)
            {
                double kth = row[row.Length - 1].Distance;
                var full = new List<Neighbor>(k + 4);
                for (int j = 0; j < n; j++)
                {
                    if (i == j) continue;
                    double d = dist(i, j);
                    if (d <= kth) full.Add(new Neighbor { Index = j, Distance = d });
                }
                if (full.Count > k)
                {
                    full.Sort(static (a, b) =>
                        a.Distance != b.Distance
                            ? a.Distance.CompareTo(b.Distance)
                            : a.Index.CompareTo(b.Index));
                    row = full.ToArray();
                }
            }

            directedNeighbors[i] = row;
            kthDistances[i] = row.Length > 0 ? row[row.Length - 1].Distance : 0.0;
            nnDistances[i]  = row.Length > 0 ? row[0].Distance                : 0.0;
        });

        return new NeighborSelection(directedNeighbors, nnDistances, kthDistances);
    }
}
