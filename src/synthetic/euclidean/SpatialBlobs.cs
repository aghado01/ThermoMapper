using System;
using System.Collections.Generic;
using Synthetic;

namespace Synthetic.Euclidean;

/// <summary>
/// Compact Gaussian clusters arranged on a hypersphere so inter-cluster
/// gaps exceed intra-cluster spread. Classic "obvious blobs" case:
/// clear density voids between clusters, no manifold structure,
/// isotropic covariance within each cluster.
/// Reference: classic isotropic Gaussian blobs benchmark.
/// </summary>
public static class SpatialBlobs
{
    public static SyntheticDataset Generate(
        int clusterCount = 4,
        int pointsPerCluster = 1500,   // bumped from 50 → 6000 total, plausible-data scale
        int dimensions = 2,
        double separation = 5.0,
        double spread = 0.5,
        int seed = 42)
    {
        var rng = new Random(seed);
        int n = clusterCount * pointsPerCluster;
        var features = new double[n][];
        var labels = new int[n];
        var centroids = SyntheticData.PlaceCentroidsOnSphere(clusterCount, dimensions, separation, rng);

        int idx = 0;
        for (int c = 0; c < clusterCount; c++)
        {
            for (int p = 0; p < pointsPerCluster; p++)
            {
                var point = new double[dimensions];
                for (int d = 0; d < dimensions; d++)
                    point[d] = centroids[c][d] + spread * SyntheticData.SampleStandardNormal(rng);
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
                ["dimensions"] = dimensions,
                ["separation"] = separation,
                ["spread"] = spread,
                ["seed"] = seed
            },
            Metadata = new SyntheticDatasetMeta(
                GeneratorName: nameof(SpatialBlobs),
                GeometryClass: "Euclidean",
                TopologyTag: "blobs",
                HierarchyTag: "none",
                GTNumClusters: clusterCount,
                AmbientDimensionality: dimensions,
                LiteratureReference: "classic isotropic Gaussian blobs benchmark")
        };
    }
}


