using System;
using System.Collections.Generic;
using Synthetic;

namespace Synthetic.Euclidean;

/// <summary>
/// Each cluster lives on its own contiguous range of feature
/// dimensions ("support"), with configurable overlap between adjacent
/// clusters. Feature values are uniform [0.5, 1.5] inside the support
/// and exactly zero outside. Exposes the high-dimensional distance-
/// concentration regime; supportOverlap controls how ambiguous cluster
/// boundaries are.
/// Reference: sparse high-dimensional support benchmark.
/// </summary>
public static class SparseSupports
{
    public static SyntheticDataset Generate(
        int clusterCount = 4,
        int pointsPerCluster = 1500,   // bumped from 50 → 6000 total, plausible high-D scale
        int featureDim = 500,
        int supportSize = 50,
        int supportOverlap = 10,
        int seed = 42)
    {
        if (supportOverlap >= supportSize)
            throw new ArgumentException("supportOverlap must be less than supportSize.");
        int stride = supportSize - supportOverlap;
        int maxStart = (clusterCount - 1) * stride + supportSize;
        if (maxStart > featureDim)
            throw new ArgumentException(
                $"featureDim ({featureDim}) too small for {clusterCount} clusters " +
                $"with supportSize {supportSize} and overlap {supportOverlap}. " +
                $"Need at least {maxStart}.");

        var rng = new Random(seed);
        int n = clusterCount * pointsPerCluster;
        var features = new double[n][];
        var labels = new int[n];

        int idx = 0;
        for (int c = 0; c < clusterCount; c++)
        {
            int supportStart = c * stride;
            for (int p = 0; p < pointsPerCluster; p++)
            {
                var point = new double[featureDim];
                for (int d = supportStart; d < supportStart + supportSize; d++)
                    point[d] = 0.5 + rng.NextDouble();
                features[idx] = point;
                labels[idx] = c;
                idx++;
            }
        }

        return new SyntheticDataset
        {
            Features = features,
            Labels = labels,
            ClusterCount = clusterCount,
            LabelsByLevel = new[] { labels },
            Parameters = new Dictionary<string, object>
            {
                ["clusterCount"] = clusterCount,
                ["pointsPerCluster"] = pointsPerCluster,
                ["featureDim"] = featureDim,
                ["supportSize"] = supportSize,
                ["supportOverlap"] = supportOverlap,
                ["seed"] = seed
            },
            Metadata = new SyntheticDatasetMeta(
                GeneratorName: nameof(SparseSupports),
                GeometryClass: "Euclidean",
                TopologyTag: "sparse",
                HierarchyTag: "none",
                GTNumClusters: clusterCount,
                AmbientDimensionality: featureDim,
                LiteratureReference: "sparse high-dimensional support benchmark")
        };
    }
}


