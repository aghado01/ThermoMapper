using System;
using System.Collections.Generic;
using Graphs.Primitives;

namespace Clustering.Dendrograms;

/// <summary>
/// Co-membership of every node pair, the all-pairs form of the per-edge <see cref="CoMembership"/> currency: the
/// fraction of kept partitions holding i and j in one cluster. It accumulates a sampler's kept draws, and so reaches
/// what an edge column can't: pairs a chain of edges joins with no edge of their own.
/// </summary>
/// <remarks>
/// Holds N(N−1)/2 counts, so about 120 MB at N = 7,800. Filling costs Σ|c|² per kept partition, which a giant cluster
/// makes N²/2.
/// </remarks>
public sealed class DenseCoMembership
{
    private readonly int[] _counts;      // the upper triangle, row by row
    private readonly int[] _rowOffset;   // the slot of (i, j) is _rowOffset[i] + j, for j > i

    private DenseCoMembership(int nodeCount, int partitionCount, int[] counts, int[] rowOffset)
    {
        NodeCount = nodeCount;
        PartitionCount = partitionCount;
        _counts = counts;
        _rowOffset = rowOffset;
    }

    /// <summary>The nodes the partitions label.</summary>
    public int NodeCount { get; }

    /// <summary>The kept partitions accumulated.</summary>
    public int PartitionCount { get; }

    /// <summary>The fraction of kept partitions holding <paramref name="i"/> and <paramref name="j"/> together.</summary>
    public double At(int i, int j)
    {
        if (i == j) return 1.0;
        if (i > j) (i, j) = (j, i);
        return (double)_counts[_rowOffset[i] + j] / PartitionCount;
    }

    /// <summary>
    /// Accumulate kept partitions, given as runs of labels flattened partition by partition with
    /// <paramref name="nodeCount"/> labels each, as a sampler keeps them. Labels lie in [0, nodeCount), and only need
    /// to name blocks within one partition.
    /// </summary>
    public static DenseCoMembership FromPartitions(int nodeCount, IReadOnlyList<int[]> runs)
    {
        ArgumentNullException.ThrowIfNull(runs);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(nodeCount);

        long pairs = (long)nodeCount * (nodeCount - 1) / 2;
        if (pairs > Array.MaxLength)
            throw new ArgumentException($"{nodeCount} nodes have {pairs} pairs, past one array's reach.", nameof(nodeCount));

        var rowOffset = new int[nodeCount];
        long offset = 0;
        for (int i = 0; i < nodeCount; i++)
        {
            rowOffset[i] = (int)(offset - i - 1);
            offset += nodeCount - i - 1;
        }

        var counts = new int[pairs];
        var sizes = new int[nodeCount];
        var starts = new int[nodeCount];
        var members = new int[nodeCount];
        int partitionCount = 0;

        foreach (int[] run in runs)
        {
            ArgumentNullException.ThrowIfNull(run);
            if (run.Length % nodeCount != 0)
                throw new ArgumentException($"A run of {run.Length} labels is not whole partitions of {nodeCount} nodes.", nameof(runs));

            for (int p = 0; p < run.Length; p += nodeCount)
            {
                ReadOnlySpan<int> labels = run.AsSpan(p, nodeCount);

                // A stable counting sort, so each block's members ascend and each pair lands on its (i < j) slot.
                Array.Clear(sizes);
                foreach (int label in labels)
                {
                    if (label < 0 || label >= nodeCount)
                        throw new ArgumentException($"Label {label} lies outside [0, {nodeCount}).", nameof(runs));
                    sizes[label]++;
                }

                int start = 0;
                for (int l = 0; l < nodeCount; l++)
                {
                    starts[l] = start;
                    start += sizes[l];
                }
                for (int i = 0; i < nodeCount; i++)
                    members[starts[labels[i]]++] = i;

                start = 0;
                for (int l = 0; l < nodeCount; l++)
                {
                    int end = start + sizes[l];
                    for (int a = start; a < end - 1; a++)
                    {
                        int row = rowOffset[members[a]];
                        for (int b = a + 1; b < end; b++)
                            counts[row + members[b]]++;
                    }
                    start = end;
                }
                partitionCount++;
            }
        }

        if (partitionCount == 0)
            throw new ArgumentException("At least one kept partition is required.", nameof(runs));
        return new DenseCoMembership(nodeCount, partitionCount, counts, rowOffset);
    }

    /// <summary>The partitions of both accumulated together.</summary>
    public static DenseCoMembership Pool(DenseCoMembership a, DenseCoMembership b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        if (a.NodeCount != b.NodeCount)
            throw new ArgumentException($"Node counts differ: {a.NodeCount} and {b.NodeCount}.");

        var counts = new int[a._counts.Length];
        for (int s = 0; s < counts.Length; s++)
            counts[s] = a._counts[s] + b._counts[s];
        return new DenseCoMembership(a.NodeCount, a.PartitionCount + b.PartitionCount, counts, a._rowOffset);
    }

    /// <summary>
    /// The per-edge column on <paramref name="graph"/>, indexed by CSR slot and meaningful at the <c>j &gt; i</c> slots,
    /// as <see cref="CoMembership.G"/> is.
    /// </summary>
    public double[] EdgeColumn(CsrGraph graph)
    {
        if (graph.NodeCount != NodeCount)
            throw new ArgumentException($"The graph has {graph.NodeCount} nodes, the partitions {NodeCount}.", nameof(graph));

        var column = new double[graph.Targets.Length];
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
            column[edge.Slot] = At(edge.Source, edge.Target);
        return column;
    }
}
