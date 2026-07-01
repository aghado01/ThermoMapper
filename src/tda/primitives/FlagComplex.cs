using System;
using System.Collections.Generic;
using Graphs.Primitives;

namespace TDA.Primitives
{
    /// <summary>
    /// Enumerates 2-simplices induced by an undirected graph skeleton.
    /// Returns a flat [T×3] array of sorted vertex triples with v0 &lt; v1 &lt; v2.
    /// </summary>
    public static class FlagComplex
    {
        public static int[] FromEdges(ReadOnlySpan<int> src, ReadOnlySpan<int> dst)
        {
            if (src.Length != dst.Length)
                throw new ArgumentException("Source and destination edge arrays must be the same length.");

            if (src.Length == 0)
                return Array.Empty<int>();

            int vertexCount = GetVertexCount(src, dst);
            var adjacencySets = new HashSet<int>[vertexCount];

            for (int i = 0; i < src.Length; i++)
            {
                int a = src[i];
                int b = dst[i];
                if (a < 0 || b < 0)
                    throw new ArgumentOutOfRangeException(nameof(src), "Vertex indices must be non-negative.");
                if (a == b)
                    continue;

                int u = Math.Min(a, b);
                int v = Math.Max(a, b);

                (adjacencySets[u] ??= new HashSet<int>()).Add(v);
                (adjacencySets[v] ??= new HashSet<int>()).Add(u);
            }

            int[][] adjacency = BuildSortedAdjacency(adjacencySets);
            var triangles = new List<int>();
            AppendTriangles(adjacency, triangles);
            return triangles.ToArray();
        }

        /// <summary>
        /// Counts triangles (3-cycles) in the undirected graph without materializing
        /// the triangle triple list. Delegates to
        /// <see cref="TriangleCount.OnCsrGraph(CsrGraph)"/> — the algorithm
        /// is structural and lives in <c>Graphs.Primitives</c> so graph
        /// diagnostics don't need a TDA-layer dependency.
        /// </summary>
        public static int CountTriangles(CsrGraph graph) => TriangleCount.OnCsrGraph(graph);

        /// <summary>
        /// Enumerates triangle vertex triples (u &lt; v &lt; w) in a symmetric CSR graph.
        /// Uses the same sorted-neighbor merge-walk as <see cref="FromEdges"/>.
        /// </summary>
        public static int[] Triangles(CsrGraph graph)
        {
            int[][] adjacency = BuildSortedAdjacencyFromCsr(graph);
            var triangles = new List<int>();
            AppendTriangles(adjacency, triangles);
            return triangles.ToArray();
        }

        private static int[][] BuildSortedAdjacencyFromCsr(CsrGraph graph)
        {
            var adjacencySets = new HashSet<int>[graph.NodeCount];
            for (int u = 0; u < graph.NodeCount; u++)
            {
                int rowEnd = graph.RowPointers[u + 1];
                for (int e = graph.RowPointers[u]; e < rowEnd; e++)
                {
                    int v = graph.Targets[e];
                    if (v <= u) continue;

                    (adjacencySets[u] ??= new HashSet<int>()).Add(v);
                    (adjacencySets[v] ??= new HashSet<int>()).Add(u);
                }
            }

            return BuildSortedAdjacency(adjacencySets);
        }

        private static int GetVertexCount(ReadOnlySpan<int> src, ReadOnlySpan<int> dst)
        {
            int maxVertex = -1;
            for (int i = 0; i < src.Length; i++)
            {
                maxVertex = Math.Max(maxVertex, src[i]);
                maxVertex = Math.Max(maxVertex, dst[i]);
            }

            return maxVertex + 1;
        }

        private static int[][] BuildSortedAdjacency(HashSet<int>[] adjacencySets)
        {
            var adjacency = new int[adjacencySets.Length][];
            for (int i = 0; i < adjacencySets.Length; i++)
            {
                if (adjacencySets[i] is null || adjacencySets[i].Count == 0)
                {
                    adjacency[i] = Array.Empty<int>();
                    continue;
                }

                int[] neighbors = new int[adjacencySets[i].Count];
                adjacencySets[i].CopyTo(neighbors);
                Array.Sort(neighbors);
                adjacency[i] = neighbors;
            }

            return adjacency;
        }

        private static void AppendTriangles(int[][] adjacency, List<int> triangles)
        {
            for (int u = 0; u < adjacency.Length; u++)
            {
                int[] neighborsOfU = adjacency[u];
                for (int index = 0; index < neighborsOfU.Length; index++)
                {
                    int v = neighborsOfU[index];
                    if (v <= u)
                        continue;

                    AppendSharedHigherNeighbors(u, v, neighborsOfU, adjacency[v], triangles);
                }
            }
        }

        private static void AppendSharedHigherNeighbors(
            int u,
            int v,
            int[] neighborsOfU,
            int[] neighborsOfV,
            List<int> triangles)
        {
            int left = AdvancePast(neighborsOfU, v);
            int right = AdvancePast(neighborsOfV, v);

            while (left < neighborsOfU.Length && right < neighborsOfV.Length)
            {
                int leftValue = neighborsOfU[left];
                int rightValue = neighborsOfV[right];

                if (leftValue == rightValue)
                {
                    triangles.Add(u);
                    triangles.Add(v);
                    triangles.Add(leftValue);
                    left++;
                    right++;
                    continue;
                }

                if (leftValue < rightValue)
                    left++;
                else
                    right++;
            }
        }

        private static int AdvancePast(int[] values, int lowerBoundExclusive)
        {
            int index = 0;
            while (index < values.Length && values[index] <= lowerBoundExclusive)
                index++;
            return index;
        }
    }
}
