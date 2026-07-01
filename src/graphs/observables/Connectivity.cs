using System.Collections.Generic;
using Graphs.Primitives;

namespace Graphs.Observables
{
    public readonly record struct ConnectivityReport(
        int ComponentCount,
        int LargestComponent,
        double LargestCoverage,
        int IsolatedNodes,
        int EdgeCount,
        int NodeCount);

    public static class Connectivity
    {
        /// <summary>
        /// Connectivity diagnostics: component count, largest component coverage,
        /// isolated node count. O(E alpha(N)) via union-find.
        /// </summary>
        public static ConnectivityReport Validate(CsrGraph graph)
        {
            int n = graph.NodeCount;
            var uf = new UnionFind(n);

            // Iterate each edge once (source < target direction only)
            for (int i = 0; i < n; i++)
            {
                int rowStart = graph.RowPointers[i];
                int rowEnd = graph.RowPointers[i + 1];
                for (int idx = rowStart; idx < rowEnd; idx++)
                {
                    int j = graph.Targets[idx];
                    if (j > i) uf.Union(i, j);
                }
            }

            int edgeCount = graph.Targets.Length / 2;
            var compSizes = new Dictionary<int, int>(n / 4);
            for (int node = 0; node < n; node++)
            {
                int root = uf.Find(node);
                if (!compSizes.TryGetValue(root, out int size))
                    compSizes[root] = 1;
                else
                    compSizes[root] = size + 1;
            }

            int largest = 0;
            int isolated = 0;
            foreach (int size in compSizes.Values)
            {
                if (size > largest) largest = size;
                if (size == 1) isolated++;
            }

            return new ConnectivityReport(
                ComponentCount: compSizes.Count,
                LargestComponent: largest,
                LargestCoverage: n > 0 ? (double)largest / n : 0.0,
                IsolatedNodes: isolated,
                EdgeCount: edgeCount,
                NodeCount: n);
        }
    }
}
