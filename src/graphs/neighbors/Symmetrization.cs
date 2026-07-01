using System;
using System.Collections.Generic;
using Graphs;
using Graphs.Primitives;

namespace Graphs.Neighbors;

/// <summary>
/// OR- and mutual-KNN symmetrization over a directed <see cref="NeighborSelection"/>.
/// </summary>
public static class Symmetrization
{
    /// <summary>OR-union: undirected edge (i,j) if i→j or j→i in the directed view.</summary>
    public static NeighborSelection OrUnion(NeighborSelection directed, int n)
    {
        if (directed.AllNeighbors is null)
            throw new ArgumentNullException(nameof(directed));
        if (n < 0)
            throw new ArgumentOutOfRangeException(nameof(n));
        if (directed.AllNeighbors.Length != n)
            throw new ArgumentException("Directed row count must match node count.", nameof(n));

        var pairMap = BuildOrPairMap(directed.AllNeighbors, n);
        return BuildSymmetricSelection(pairMap, n, directed.KthNeighborDistances);
    }

    /// <summary>Mutual intersection: edge (i,j) only when both nominate each other.</summary>
    public static NeighborSelection MutualIntersection(
        NeighborSelection directed,
        int n,
        MutualBandwidthSource bandwidthSource = MutualBandwidthSource.DirectedKth)
    {
        if (directed.AllNeighbors is null)
            throw new ArgumentNullException(nameof(directed));
        if (n < 0)
            throw new ArgumentOutOfRangeException(nameof(n));
        if (directed.AllNeighbors.Length != n)
            throw new ArgumentException("Directed row count must match node count.", nameof(n));

        var directedNeighbors = directed.AllNeighbors;
        var directedSets = new HashSet<int>[n];
        for (int i = 0; i < n; i++)
        {
            var set = new HashSet<int>(directedNeighbors[i].Length);
            foreach (var nb in directedNeighbors[i])
                set.Add(nb.Index);
            directedSets[i] = set;
        }

        var pairMap = new Dictionary<long, double>();
        bool wantsMutualKth = bandwidthSource == MutualBandwidthSource.MutualKth;
        double[]? mutualRowMax = wantsMutualKth ? new double[n] : null;
        bool[]? mutualRowHasNeighbor = wantsMutualKth ? new bool[n] : null;

        for (int i = 0; i < n; i++)
        {
            foreach (var nb in directedNeighbors[i])
            {
                int j = nb.Index;
                if (!directedSets[j].Contains(i)) continue;

                if (mutualRowMax is not null)
                {
                    if (!mutualRowHasNeighbor![i] || nb.Distance > mutualRowMax[i])
                        mutualRowMax[i] = nb.Distance;
                    mutualRowHasNeighbor[i] = true;
                }

                int lo = Math.Min(i, j);
                int hi = Math.Max(i, j);
                long key = (long)lo * n + hi;
                if (!pairMap.TryGetValue(key, out double existing) || nb.Distance < existing)
                    pairMap[key] = nb.Distance;
            }
        }

        var resolvedKth = new double[n];
        if (mutualRowMax is not null)
        {
            for (int i = 0; i < n; i++)
            {
                resolvedKth[i] = mutualRowHasNeighbor![i]
                    ? mutualRowMax[i]
                    : double.PositiveInfinity;
            }
        }
        else
        {
            Array.Copy(directed.KthNeighborDistances, resolvedKth, n);
        }

        return BuildSymmetricSelection(pairMap, n, resolvedKth);
    }

    private static Dictionary<long, double> BuildOrPairMap(Neighbor[][] directedNeighbors, int n)
    {
        var pairMap = new Dictionary<long, double>();
        for (int i = 0; i < n; i++)
        {
            foreach (var nb in directedNeighbors[i])
            {
                int lo = Math.Min(i, nb.Index);
                int hi = Math.Max(i, nb.Index);
                long key = (long)lo * n + hi;
                if (!pairMap.TryGetValue(key, out double existing) || nb.Distance < existing)
                    pairMap[key] = nb.Distance;
            }
        }

        return pairMap;
    }

    private static NeighborSelection BuildSymmetricSelection(
        Dictionary<long, double> pairMap,
        int n,
        double[] kthDistances)
    {
        var lists = new List<Neighbor>[n];
        for (int i = 0; i < n; i++)
            lists[i] = new List<Neighbor>();

        foreach (var kvp in pairMap)
        {
            int lo = (int)(kvp.Key / n);
            int hi = (int)(kvp.Key % n);
            double d = kvp.Value;
            lists[lo].Add(new Neighbor { Index = hi, Distance = d });
            lists[hi].Add(new Neighbor { Index = lo, Distance = d });
        }

        var allNeighbors = new Neighbor[n][];
        var nnDistances  = new double[n];
        for (int i = 0; i < n; i++)
        {
            lists[i].Sort((a, b) => a.Distance.CompareTo(b.Distance));
            allNeighbors[i] = lists[i].ToArray();
            nnDistances[i] = allNeighbors[i].Length > 0
                ? allNeighbors[i][0].Distance
                : double.PositiveInfinity;
        }

        return new NeighborSelection(allNeighbors, nnDistances, kthDistances);
    }
}
