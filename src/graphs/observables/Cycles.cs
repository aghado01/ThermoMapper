using System;
using System.Collections.Generic;
using Graphs.Primitives;

namespace Graphs.Observables
{
    public readonly record struct CycleReport(
        int CyclomaticComplexity,
        int TriangleCount,
        int Girth,
        double TriangleSaturation,
        double MeanCycleLength,
        int MaxCycleLength);

    public static class Cycles
    {
        public static CycleReport Compute(
            CsrGraph graph,
            ConnectivityReport connectivity,
            int maxNodesForCycleStats = 5000)
        {
            int cyclomaticComplexity = Math.Max(0, connectivity.EdgeCount - connectivity.NodeCount + connectivity.ComponentCount);

            if (cyclomaticComplexity == 0)
            {
                return new CycleReport(
                    CyclomaticComplexity: cyclomaticComplexity,
                    TriangleCount: 0,
                    Girth: -1,
                    TriangleSaturation: 0.0,
                    MeanCycleLength: 0.0,
                    MaxCycleLength: 0);
            }

            if (graph.NodeCount > maxNodesForCycleStats)
            {
                return new CycleReport(
                    CyclomaticComplexity: cyclomaticComplexity,
                    TriangleCount: -1,
                    Girth: -1,
                    TriangleSaturation: double.NaN,
                    MeanCycleLength: double.NaN,
                    MaxCycleLength: -1);
            }

            int triangleCount = TriangleCount.OnCsrGraph(graph);
            double triangleSaturation = triangleCount / (double)Math.Max(1, cyclomaticComplexity);

            HashSet<int> cycleLengths = CollectDistinctCycleLengths(graph);
            if (cycleLengths.Count == 0)
            {
                return new CycleReport(
                    CyclomaticComplexity: cyclomaticComplexity,
                    TriangleCount: triangleCount,
                    Girth: -1,
                    TriangleSaturation: triangleSaturation,
                    MeanCycleLength: 0.0,
                    MaxCycleLength: 0);
            }

            int girth = int.MaxValue;
            int maxCycleLength = 0;
            long totalCycleLength = 0;

            foreach (int cycleLength in cycleLengths)
            {
                if (cycleLength < girth) girth = cycleLength;
                if (cycleLength > maxCycleLength) maxCycleLength = cycleLength;
                totalCycleLength += cycleLength;
            }

            return new CycleReport(
                CyclomaticComplexity: cyclomaticComplexity,
                TriangleCount: triangleCount,
                Girth: girth,
                TriangleSaturation: triangleSaturation,
                MeanCycleLength: totalCycleLength / (double)cycleLengths.Count,
                MaxCycleLength: maxCycleLength);
        }

        public static CycleReport Compute(CsrGraph graph, int maxNodesForCycleStats = 5000)
            => Compute(graph, Connectivity.Validate(graph), maxNodesForCycleStats);

        private static HashSet<int> CollectDistinctCycleLengths(CsrGraph graph)
        {
            var cycleLengths = new HashSet<int>();
            int n = graph.NodeCount;
            var parent = new int[n];
            var depth = new int[n];
            var queue = new int[n];

            for (int root = 0; root < n; root++)
            {
                Array.Fill(parent, -1);
                Array.Fill(depth, -1);

                int head = 0;
                int tail = 0;
                queue[tail++] = root;
                depth[root] = 0;

                while (head < tail)
                {
                    int node = queue[head++];
                    int rowStart = graph.RowPointers[node];
                    int rowEnd = graph.RowPointers[node + 1];

                    for (int edgeIndex = rowStart; edgeIndex < rowEnd; edgeIndex++)
                    {
                        int neighbor = graph.Targets[edgeIndex];
                        if (depth[neighbor] < 0)
                        {
                            depth[neighbor] = depth[node] + 1;
                            parent[neighbor] = node;
                            queue[tail++] = neighbor;
                            continue;
                        }

                        if (parent[node] == neighbor || parent[neighbor] == node || node >= neighbor)
                            continue;

                        int cycleLength = ComputeClosureLength(node, neighbor, parent, depth);
                        if (cycleLength >= 3)
                            cycleLengths.Add(cycleLength);
                    }
                }
            }

            return cycleLengths;
        }

        private static int ComputeClosureLength(int left, int right, int[] parent, int[] depth)
        {
            int cycleLength = 1;
            int leftDepth = depth[left];
            int rightDepth = depth[right];

            while (leftDepth > rightDepth)
            {
                left = parent[left];
                leftDepth--;
                cycleLength++;
            }

            while (rightDepth > leftDepth)
            {
                right = parent[right];
                rightDepth--;
                cycleLength++;
            }

            while (left != right)
            {
                left = parent[left];
                right = parent[right];
                cycleLength += 2;
            }

            return cycleLength;
        }
    }
}
