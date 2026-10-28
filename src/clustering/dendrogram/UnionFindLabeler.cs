using System.Collections.Generic;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Dendrograms;

/// <summary>
/// Collapses a <see cref="UnionFind"/>'s sparse root ids into a densely-labeled
/// <see cref="Assignment"/> (labels in <c>[0, Count)</c>, first-seen order).
/// Shared by every union-find partition close: threshold-style cuts and
/// <see cref="Dendrogram"/> cuts.
/// </summary>
public static class UnionFindLabeler
{
    public static Assignment Densify(UnionFind uf, int nodeCount)
    {
        var labels   = new int[nodeCount];
        var labelMap = new Dictionary<int, int>();
        int nextLabel = 0;
        for (int i = 0; i < nodeCount; i++)
        {
            int root = uf.Find(i);
            if (!labelMap.TryGetValue(root, out int dense))
            {
                dense = nextLabel++;
                labelMap[root] = dense;
            }
            labels[i] = dense;
        }
        return new Assignment { Labels = labels, Count = nextLabel };
    }
}
