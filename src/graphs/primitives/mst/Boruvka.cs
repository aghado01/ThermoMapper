using System;
using System.Collections.Generic;

namespace Graphs.Primitives.Mst;

/// <summary>
/// Borůvka's MST algorithm in its pre-seeded "minimal-bridge" flavor. Each
/// phase performs one O(N²) sweep of cheapest cross-component edges and
/// merges the discovered components. Typical Mutual-kNN inputs (1–3
/// connected components after filtering) finish in 1–2 phases.
/// </summary>
/// <remarks>
/// The full-MST flavor of Borůvka falls out of the same loop by simply
/// seeding the <see cref="UnionFind"/> with every node in its own
/// component. That variant is not exposed yet — <see cref="Prim"/> handles
/// the only current full-MST consumer (HDBSCAN's mutual-reachability
/// path). A sibling <c>Kruskal</c> now joins this namespace for consumers
/// that already materialise a sortable edge list.
/// </remarks>
public static class Boruvka
{
    /// <summary>
    /// One bridge edge between two pre-existing components, as discovered
    /// by Borůvka. <see cref="LoIndex"/> and <see cref="HiIndex"/> are
    /// node indices with <c>LoIndex &lt; HiIndex</c>.
    /// </summary>
    public readonly record struct BridgeEdge(int LoIndex, int HiIndex, double Weight);

    /// <summary>
    /// Finds the minimal set of bridge edges needed to unify the
    /// components encoded by <paramref name="components"/> into a single
    /// connected component. Iterates Borůvka phases until either one
    /// component remains or a phase makes no progress (the graph is
    /// already maximally connectable under <paramref name="pairDistance"/>).
    /// Each phase costs O(N²) <paramref name="pairDistance"/> evaluations;
    /// allocations are bounded by <c>O(numComponents)</c> per phase.
    /// </summary>
    /// <param name="n">Total node count.</param>
    /// <param name="pairDistance">Pairwise distance function — called only
    /// for nodes in different components.</param>
    /// <param name="components">Union-Find pre-seeded with the connected
    /// components of the input graph; mutated in place as bridges are
    /// added.</param>
    /// <returns>The list of bridge edges added across all phases, in
    /// discovery order. Empty when the input was already connected.</returns>
    public static List<BridgeEdge> AddMinimalBridges(
        int                    n,
        Func<int, int, double> pairDistance,
        UnionFind              components)
    {
        var bridges = new List<BridgeEdge>();

        while (true)
        {
            int[] roots = components.GetLabels();
            var rootToIdx = new Dictionary<int, int>();
            int numComp = 0;
            for (int i = 0; i < n; i++)
                if (!rootToIdx.ContainsKey(roots[i]))
                    rootToIdx[roots[i]] = numComp++;
            if (numComp == 1) break;

            // Cheapest outgoing edge per component (no alloc beyond small fixed arrays).
            int[]    cheapSrc  = new int[numComp];
            int[]    cheapDst  = new int[numComp];
            double[] cheapDist = new double[numComp];
            Array.Fill(cheapSrc, -1);
            Array.Fill(cheapDist, double.MaxValue);

            for (int i = 0; i < n; i++)
            {
                int ci = rootToIdx[roots[i]];
                for (int j = i + 1; j < n; j++)
                {
                    int cj = rootToIdx[roots[j]];
                    if (ci == cj) continue;
                    double d = pairDistance(i, j);
                    if (d < cheapDist[ci]) { cheapDist[ci] = d; cheapSrc[ci] = i; cheapDst[ci] = j; }
                    if (d < cheapDist[cj]) { cheapDist[cj] = d; cheapSrc[cj] = j; cheapDst[cj] = i; }
                }
            }

            bool anyAdded = false;
            for (int c = 0; c < numComp; c++)
            {
                if (cheapSrc[c] < 0) continue;
                int a = cheapSrc[c], b = cheapDst[c];
                if (components.Find(a) == components.Find(b)) continue;
                int lo = Math.Min(a, b), hi = Math.Max(a, b);
                bridges.Add(new BridgeEdge(lo, hi, cheapDist[c]));
                components.Union(a, b);
                anyAdded = true;
            }
            if (!anyAdded) break;
        }

        return bridges;
    }
}
