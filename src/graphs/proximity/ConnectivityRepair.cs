using System;
using System.Collections.Generic;
using Graphs.Primitives;
using Graphs.Primitives.Mst;

namespace Graphs.Proximity;

/// <summary>
/// Repairs a (potentially shattered) <see cref="NeighborSelection"/> by
/// injecting the minimal MST bridge edges needed to unify its components
/// into one. The MST-min strategy from the UMAP connectivity literature —
/// only edges that actually bridge disjoint components are added; existing
/// neighborhood structure is preserved.
/// </summary>
/// <remarks>
/// <para>This is a thin shell over <see cref="Boruvka.AddMinimalBridges"/>:
/// the primitive finds the bridges, this class handles the
/// <see cref="NeighborSelection"/> reconstruction. When the 5-stage
/// pipeline lands, this responsibility becomes the
/// <c>MstMinRepair</c> stage implementation.</para>
///
/// <para>The original directed k-th distance sample is passed through
/// unchanged — MST bridge edges are by definition non-local and would
/// skew bandwidth estimation if they leaked into the sample. Downstream
/// consumers of the repaired selection should still see the unpolluted
/// density sample.</para>
/// </remarks>
public static class ConnectivityRepair
{
    /// <summary>
    /// Returns a new <see cref="NeighborSelection"/> with MST-min bridge
    /// edges added so the underlying graph is fully connected. The input
    /// selection is not mutated. <paramref name="pairDistance"/> is
    /// invoked only for nodes in different components.
    /// </summary>
    public static NeighborSelection EnsureConnected(
        NeighborSelection sel, int n, Func<int, int, double> pairDistance)
    {
        var baseNeighbors = sel.AllNeighbors;

        var pairMap = new Dictionary<long, double>();
        var uf = new UnionFind(n);
        for (int i = 0; i < n; i++)
        {
            foreach (var nb in baseNeighbors[i])
            {
                int lo = Math.Min(i, nb.Index);
                int hi = Math.Max(i, nb.Index);
                long key = (long)lo * n + hi;
                pairMap.TryAdd(key, nb.Distance);
                uf.Union(lo, hi);
            }
        }

        // Borůvka's pre-seeded variant emits only the cross-component
        // bridges — exactly the MST-min set.
        var bridges = Boruvka.AddMinimalBridges(n, pairDistance, uf);
        foreach (var bridge in bridges)
        {
            long key = (long)bridge.LoIndex * n + bridge.HiIndex;
            pairMap.TryAdd(key, bridge.Weight);
        }

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
        var nnDistances = new double[n];

        for (int i = 0; i < n; i++)
        {
            lists[i].Sort((a, b) => a.Distance.CompareTo(b.Distance));
            allNeighbors[i] = lists[i].ToArray();
            nnDistances[i] = allNeighbors[i].Length > 0
                ? allNeighbors[i][0].Distance
                : double.PositiveInfinity;
        }

        return new NeighborSelection(allNeighbors, nnDistances, sel.KthNeighborDistances);
    }
}
