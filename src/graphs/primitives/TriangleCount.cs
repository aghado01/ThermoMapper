using System;

namespace Graphs.Primitives
{
    /// <summary>
    /// Pure graph operation: count 3-cycles (triangles) in an undirected
    /// <see cref="CsrGraph"/> via sorted-adjacency neighbor intersection.
    /// Lives in <c>Graphs.Primitives</c> rather than <c>TDA.Ph</c>
    /// because the operation is structural (counts closed triplets in the
    /// adjacency relation) and has no topological content — TDA
    /// consumers depend on it, not the other way around.
    /// </summary>
    /// <remarks>
    /// <para><b>Algorithm.</b> Build sorted adjacency lists per node (CSR
    /// rows are not guaranteed sorted), then for each ordered edge
    /// <c>u &lt; v</c> count shared neighbors <c>w &gt; v</c>. The
    /// per-pair intersection is linear-merge over the two sorted lists.
    /// Each triangle <c>{u, v, w}</c> with <c>u &lt; v &lt; w</c> is
    /// counted exactly once.</para>
    ///
    /// <para><b>Cost.</b> O(Σ_v d(v)²) in the worst case, dominated by
    /// the intersection step on high-degree vertices. Fine for the
    /// proximity graphs SPC and diagnostics actually consume (~50–5000
    /// nodes, bounded k). Not appropriate for dense general graphs;
    /// the <c>Cycles</c> diagnostic gates this behind a node-count cap.</para>
    /// </remarks>
    public static class TriangleCount
    {
        public static int OnCsrGraph(CsrGraph graph)
        {
            if (graph.NodeCount == 0) return 0;

            int[][] adjacency = BuildSortedAdjacency(graph);
            int triangleCount = 0;

            for (int u = 0; u < adjacency.Length; u++)
            {
                int[] neighborsOfU = adjacency[u];
                for (int index = 0; index < neighborsOfU.Length; index++)
                {
                    int v = neighborsOfU[index];
                    if (v <= u)
                        continue;

                    triangleCount += CountSharedHigherNeighbors(v, neighborsOfU, adjacency[v]);
                }
            }
            return triangleCount;
        }

        private static int[][] BuildSortedAdjacency(CsrGraph graph)
        {
            var adjacency = new int[graph.NodeCount][];
            for (int node = 0; node < graph.NodeCount; node++)
            {
                int degree = graph.RowPointers[node + 1] - graph.RowPointers[node];
                if (degree == 0)
                {
                    adjacency[node] = Array.Empty<int>();
                    continue;
                }

                var neighbors = new int[degree];
                Array.Copy(graph.Targets, graph.RowPointers[node], neighbors, 0, degree);
                Array.Sort(neighbors);
                adjacency[node] = neighbors;
            }
            return adjacency;
        }

        private static int CountSharedHigherNeighbors(
            int lowerBoundExclusive, int[] neighborsOfU, int[] neighborsOfV)
        {
            int left  = AdvancePast(neighborsOfU, lowerBoundExclusive);
            int right = AdvancePast(neighborsOfV, lowerBoundExclusive);
            int count = 0;

            while (left < neighborsOfU.Length && right < neighborsOfV.Length)
            {
                int leftValue  = neighborsOfU[left];
                int rightValue = neighborsOfV[right];
                if (leftValue == rightValue)
                {
                    count++;
                    left++;
                    right++;
                    continue;
                }
                if (leftValue < rightValue) left++;
                else                         right++;
            }
            return count;
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
