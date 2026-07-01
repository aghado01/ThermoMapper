using System;

namespace TDA.Mapper.Diagnostics;

public readonly record struct MapperNodeReport(
    int NodeCount,
    int MinSize,
    int MaxSize,
    double MedianSize,
    double MeanSize,
    double[] FilterValuePercentiles);

public static class NodeStats
{
    public static MapperNodeReport From(MapperResult result)
    {
        int nodeCount = result.Nodes.Count;
        if (nodeCount == 0)
        {
            return new MapperNodeReport(
                NodeCount: 0,
                MinSize: 0,
                MaxSize: 0,
                MedianSize: 0.0,
                MeanSize: 0.0,
                FilterValuePercentiles: Array.Empty<double>());
        }

        var sizes = new int[nodeCount];
        var filterMeans = new double[nodeCount];
        int minSize = int.MaxValue;
        int maxSize = 0;
        long totalSize = 0;

        for (int i = 0; i < nodeCount; i++)
        {
            var node = result.Nodes[i];
            int size = node.MemberIndices.Length;
            sizes[i] = size;
            filterMeans[i] = node.FilterValueMean;
            totalSize += size;

            if (size < minSize)
                minSize = size;

            if (size > maxSize)
                maxSize = size;
        }

        Array.Sort(sizes);
        Array.Sort(filterMeans);

        return new MapperNodeReport(
            NodeCount: nodeCount,
            MinSize: minSize,
            MaxSize: maxSize,
            MedianSize: ComputeMedian(sizes),
            MeanSize: totalSize / (double)nodeCount,
            FilterValuePercentiles: new[]
            {
                filterMeans[0],
                SamplePercentile(filterMeans, 0.25),
                SamplePercentile(filterMeans, 0.50),
                SamplePercentile(filterMeans, 0.75),
                filterMeans[^1],
            });
    }

    private static double ComputeMedian(int[] sortedSizes)
    {
        int mid = sortedSizes.Length / 2;
        return (sortedSizes.Length & 1) == 0
            ? (sortedSizes[mid - 1] + sortedSizes[mid]) / 2.0
            : sortedSizes[mid];
    }

    private static double SamplePercentile(double[] sortedValues, double percentile)
    {
        int index = (int)Math.Round(percentile * (sortedValues.Length - 1), MidpointRounding.AwayFromZero);
        return sortedValues[index];
    }
}
