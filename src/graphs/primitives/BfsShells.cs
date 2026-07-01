
using System;
using System.Collections.Generic;

namespace Graphs.Primitives
// ============================================================================
// Graphs.Primitives/BfsShells.cs
// ============================================================================
// Pure-topology BFS primitive. Computes:
//   - per-node hop distance from a seed (with sentinel -1 for unreachable)
//   - per-distance "shells" (each shell = the set of nodes at that hop distance)
//
// This is the canonical BFS primitive for the project, extracted from SPC
// thermo analysis per the SPC maturity plan (B1: BfsShells → src/graphs/).
// No SPC or spin dependency. Consumed by:
//   - SpcThermo.AnalyzeRadialCoherence (when thermo wiring lands)
//   - SpcThermo.BuildLocalHistogram
//   - Graphs.TDA.Mapper.GeodesicDistanceFilter
//   - Future: any radial / shell-decomposition diagnostic
//
// Complexity: O(V + E). Uses level-by-level expansion so shells are emitted
// in distance order. No allocations beyond the two output arrays + the
// current/next level lists.
// ============================================================================

{
    /// <summary>
    /// Result of a single-source BFS shell decomposition over a graph.
    /// </summary>
    public sealed class BfsShellsResult
    {
        /// <summary>Seed node BFS was launched from.</summary>
        public int SeedNode { get; }

        /// <summary>Hop distance per node. <c>Distances[i] == -1</c> means node
        /// <c>i</c> was unreachable from the seed within the (optional) max depth.</summary>
        public int[] Distances { get; }

        /// <summary>Shell decomposition. <c>Shells[d]</c> is the array of node
        /// indices at hop distance <c>d</c> from the seed. <c>Shells[0]</c>
        /// always contains exactly the seed. <c>Shells.Length</c> == MaxDepth + 1.</summary>
        public int[][] Shells { get; }

        /// <summary>Deepest hop distance reached. Equal to <c>Shells.Length - 1</c>.</summary>
        public int MaxDepth { get; }

        /// <summary>Total nodes reached. Equal to the sum of <c>Shells[d].Length</c>.</summary>
        public int ReachableNodeCount { get; }

        internal BfsShellsResult(int seedNode, int[] distances, int[][] shells, int reachableNodeCount)
        {
            SeedNode = seedNode;
            Distances = distances;
            Shells = shells;
            MaxDepth = shells.Length - 1;
            ReachableNodeCount = reachableNodeCount;
        }
    }

    /// <summary>
    /// Single-source BFS over a <see cref="CsrGraph"/>. Returns both per-node
    /// hop distances and the shell decomposition (nodes grouped by distance).
    /// </summary>
    public static class BfsShells
    {
        /// <summary>
        /// Run BFS from <paramref name="seedNode"/> over <paramref name="graph"/>.
        /// Returns hop distances and shells.
        /// </summary>
        /// <param name="graph">Source graph. Edges are traversed unweighted (each
        /// counts as one hop) — for weighted shortest paths use Dijkstra (not yet
        /// in this primitive). The CSR's symmetric edge storage means the BFS
        /// proceeds as on an undirected graph.</param>
        /// <param name="seedNode">Starting node, in [0, NodeCount).</param>
        /// <param name="maxDepth">Optional hop-count cap. If specified, BFS stops
        /// expanding past this depth. <c>null</c> means BFS the entire reachable
        /// set. Useful for local-shell analyses (e.g., k-hop neighborhoods).</param>
        public static BfsShellsResult Compute(CsrGraph graph, int seedNode, int? maxDepth = null)
        {
            int n = graph.NodeCount;
            if (n == 0)
                throw new ArgumentException("Graph is empty.", nameof(graph));
            if (seedNode < 0 || seedNode >= n)
                throw new ArgumentOutOfRangeException(
                    nameof(seedNode),
                    $"seedNode {seedNode} out of range for graph of size {n}.");
            if (maxDepth.HasValue && maxDepth.Value < 0)
                throw new ArgumentOutOfRangeException(nameof(maxDepth), "maxDepth must be >= 0 or null.");

            var distances = new int[n];
            for (int i = 0; i < n; i++) distances[i] = -1;
            distances[seedNode] = 0;

            // Shell 0 is always just the seed.
            var shells = new List<int[]>(capacity: 8) { new[] { seedNode } };
            int reachable = 1;

            // Expand level by level using paired current/next lists.
            // Level d → level d+1: drain every node in `current`, push unvisited
            // neighbors into `next` and mark their distance.
            var current = new List<int>(8) { seedNode };
            var next = new List<int>(16);
            int depth = 0;

            while (current.Count > 0)
            {
                if (maxDepth.HasValue && depth >= maxDepth.Value) break;

                next.Clear();
                foreach (int u in current)
                {
                    int rowStart = graph.RowPointers[u];
                    int rowEnd = graph.RowPointers[u + 1];
                    for (int e = rowStart; e < rowEnd; e++)
                    {
                        int v = graph.Targets[e];
                        if (distances[v] == -1)
                        {
                            distances[v] = depth + 1;
                            next.Add(v);
                        }
                    }
                }

                if (next.Count == 0) break;

                shells.Add(next.ToArray());
                reachable += next.Count;

                // Swap current and next; avoid re-allocating per level.
                (current, next) = (next, current);
                depth++;
            }

            return new BfsShellsResult(seedNode, distances, shells.ToArray(), reachable);
        }

        /// <summary>
        /// Convenience: BFS from seed and return only the per-node hop distances.
        /// Equivalent to <c>Compute(graph, seedNode, maxDepth).Distances</c> but
        /// skips the shell materialization when callers only need distances.
        /// </summary>
        public static int[] ComputeDistances(CsrGraph graph, int seedNode, int? maxDepth = null)
        {
            int n = graph.NodeCount;
            if (n == 0)
                throw new ArgumentException("Graph is empty.", nameof(graph));
            if (seedNode < 0 || seedNode >= n)
                throw new ArgumentOutOfRangeException(
                    nameof(seedNode),
                    $"seedNode {seedNode} out of range for graph of size {n}.");
            if (maxDepth.HasValue && maxDepth.Value < 0)
                throw new ArgumentOutOfRangeException(nameof(maxDepth), "maxDepth must be >= 0 or null.");

            var distances = new int[n];
            for (int i = 0; i < n; i++) distances[i] = -1;
            distances[seedNode] = 0;

            var current = new List<int>(8) { seedNode };
            var next = new List<int>(16);
            int depth = 0;

            while (current.Count > 0)
            {
                if (maxDepth.HasValue && depth >= maxDepth.Value) break;

                next.Clear();
                foreach (int u in current)
                {
                    int rowStart = graph.RowPointers[u];
                    int rowEnd = graph.RowPointers[u + 1];
                    for (int e = rowStart; e < rowEnd; e++)
                    {
                        int v = graph.Targets[e];
                        if (distances[v] == -1)
                        {
                            distances[v] = depth + 1;
                            next.Add(v);
                        }
                    }
                }

                if (next.Count == 0) break;
                (current, next) = (next, current);
                depth++;
            }

            return distances;
        }
    }
}
