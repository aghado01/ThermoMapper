using System;
using System.Collections.Generic;
using Synthetic;

namespace Synthetic.Euclidean;

/// <summary>
/// Gaussian clusters with anisotropic (elongated) shapes along random
/// rotations. Cluster means are placed close enough that raw
/// feature-space distances overlap substantially between clusters
/// while covariance structure remains distinct. Per-cluster
/// covariance matrices are returned so callers can build a pooled
/// CovarianceInverse (or construct their own) if desired.
/// Reference: anisotropic Gaussian blob benchmark.
/// </summary>
public static class AnisotropicGaussian
{
    public static SyntheticDataset Generate(
        int clusterCount = 3,
        int pointsPerCluster = 2000,   // bumped from 60 → 6000 total, plausible-data scale
        int dimensions = 2,
        double meanSeparation = 2.0,
        double anisotropyRatio = 5.0,
        int seed = 42)
    {
        if (dimensions < 2)
            throw new ArgumentException("AnisotropicGaussian requires dimensions >= 2.");
        var rng = new Random(seed);
        int n = clusterCount * pointsPerCluster;
        var features = new double[n][];
        var labels = new int[n];
        var means = SyntheticData.PlaceCentroidsOnSphere(clusterCount, dimensions, meanSeparation, rng);

        int idx = 0;
        for (int c = 0; c < clusterCount; c++)
        {
            var scale = new double[dimensions];
            scale[0] = anisotropyRatio;
            for (int d = 1; d < dimensions; d++) scale[d] = 1.0;

            var rotation = SyntheticData.RandomRotationMatrix(dimensions, rng);

            for (int p = 0; p < pointsPerCluster; p++)
            {
                var z = new double[dimensions];
                for (int d = 0; d < dimensions; d++)
                    z[d] = scale[d] * SyntheticData.SampleStandardNormal(rng);
                var rotated = SyntheticData.MultiplyMatrixVector(rotation, z);
                var point = new double[dimensions];
                for (int d = 0; d < dimensions; d++)
                    point[d] = means[c][d] + rotated[d];
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
                ["meanSeparation"] = meanSeparation,
                ["anisotropyRatio"] = anisotropyRatio,
                ["seed"] = seed
            },
            Metadata = new SyntheticDatasetMeta(
                GeneratorName: nameof(AnisotropicGaussian),
                GeometryClass: "Euclidean",
                TopologyTag: "anisotropic",
                HierarchyTag: "none",
                GTNumClusters: clusterCount,
                AmbientDimensionality: dimensions,
                LiteratureReference: "anisotropic Gaussian blob benchmark")
        };
    }
}


