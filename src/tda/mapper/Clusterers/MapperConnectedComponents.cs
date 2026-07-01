// ============================================================================
// TDA.Mapper — MapperConnectedComponents.cs
// ============================================================================
// IGraphClusterer adapter: connected components of the preimage-induced
// subgraph. The canonical graph MAPPER clusterer (Carrière-Michel-Oudot 2018,
// Hajij-Rosen-Wang 2018) — runs union-find over edges internal to each
// preimage and emits one cluster per connected component.
//
// Pair with the IGraphFilters in Filters/ via Mapper.Build(graph, features,
// filter, cover, clusterer). For thermodynamic graph clustering on the same
// preimage abstraction see MapperSpcClusterer.
// ============================================================================

#nullable enable
using System;
using System.Collections.Generic;
using Graphs.Primitives;
using TDA.Mapper;

namespace TDA.Mapper.Clusterers;

/// <summary>
/// Clusters a preimage by the connected components of the subgraph induced on
/// it. Uses <see cref="UnionFind"/> on the edges internal to the preimage.
///
/// Complexity: O(Σ_u deg(u)) where the sum is over preimage nodes. Sub-linear
/// in the original graph since edges with one endpoint outside the preimage
/// are skipped via O(1) hash-set lookup.
/// </summary>
public sealed class ConnectedComponentsClusterer : IGraphClusterer
{
    public string Name => "Connected components (induced subgraph)";

    public ClusterResult ClusterInduced(CsrGraph graph, IReadOnlyList<int> preimageIndices)
    {
        ArgumentNullException.ThrowIfNull(preimageIndices);
        int k = preimageIndices.Count;
        if (k == 0) return new ClusterResult(Array.Empty<int>(), 0);
        if (k == 1) return new ClusterResult(new[] { 0 }, 1);

        // Map original node id → preimage-local index [0, k).
        // Using a Dictionary is the cleanest way to handle non-contiguous
        // preimage indices; for very large preimages a sparse bool array
        // sized to graph.NodeCount could be faster but is rarely necessary.
        var localOf = new Dictionary<int, int>(k);
        for (int i = 0; i < k; i++) localOf[preimageIndices[i]] = i;

        var uf = new UnionFind(k);

        // Union endpoints of every edge whose other endpoint is also in the preimage.
        for (int i = 0; i < k; i++)
        {
            int u = preimageIndices[i];
            int rowStart = graph.RowPointers[u];
            int rowEnd = graph.RowPointers[u + 1];

            for (int e = rowStart; e < rowEnd; e++)
            {
                int v = graph.Targets[e];
                if (localOf.TryGetValue(v, out int localV))
                {
                    // Avoid double-work: only union when u < v (each undirected edge
                    // appears once in this direction in the symmetric CSR).
                    if (u < v) uf.Union(i, localV);
                }
            }
        }

        // Compress UF roots → sequential cluster labels [0, ccCount).
        var rootToLabel = new Dictionary<int, int>();
        var labels = new int[k];
        int nextLabel = 0;
        for (int i = 0; i < k; i++)
        {
            int root = uf.Find(i);
            if (!rootToLabel.TryGetValue(root, out int label))
            {
                label = nextLabel++;
                rootToLabel[root] = label;
            }
            labels[i] = label;
        }

        return new ClusterResult(labels, nextLabel);
    }
}
