using System.Collections.Generic;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Strategies;

/// <summary>
/// Internal helper that collapses a <see cref="UnionFind"/>'s sparse
/// root ids into a densely-labeled <see cref="Assignment"/> (labels in
/// <c>[0, Count)</c>). Shared across threshold-style partition
/// strategies that build their cluster structure via union-find.
/// </summary>
internal static class UnionFindLabeler
{
    internal static Assignment Densify(UnionFind uf, int nodeCount)
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
