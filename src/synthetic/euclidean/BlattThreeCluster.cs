using System;
using System.Collections.Generic;
using Synthetic;

namespace Synthetic.Euclidean;

/// <summary>
/// Three compact Gaussian clusters at fixed corners (-6, 0), (6, 0), (0, 6).
/// Canonical SPC engine-validation dataset.
/// Reference: Blatt, Wiseman, Domany 1996 PRL 76:3251.
/// </summary>
public static class BlattThreeCluster
{
    public static SyntheticDataset Generate(
        int pointsPerCluster = 2000,   // bumped from 100 → 6000 total, plausible-data scale
        double stdDev = 1.0,
        int dimensions = 2,
        int seed = 42)
    {
        if (dimensions < 2)
            throw new ArgumentException("BlattThreeCluster requires dimensions >= 2.");
        var rng = new Random(seed);
        var centers = new double[3][];
        centers[0] = new double[dimensions];
        centers[0][0] = -6.0;
        centers[0][1] = 0.0;
        centers[1] = new double[dimensions];
        centers[1][0] = 6.0;
        centers[1][1] = 0.0;
        centers[2] = new double[dimensions];
        centers[2][0] = 0.0;
        centers[2][1] = 6.0;
        int n = 3 * pointsPerCluster;
        var features = new double[n][];
        var labels = new int[n];

        int idx = 0;
        for (int c = 0; c < 3; c++)
        {
            for (int p = 0; p < pointsPerCluster; p++)
            {
                var point = new double[dimensions];
                point[0] = centers[c][0] + stdDev * SyntheticData.SampleStandardNormal(rng);
                point[1] = centers[c][1] + stdDev * SyntheticData.SampleStandardNormal(rng);
                features[idx] = point;
                labels[idx] = c;
                idx++;
            }
        }

        return new SyntheticDataset
        {
            Features = features,
            Labels = labels,
            ClusterCount = 3,
            LabelsByLevel = new[] { labels },
            Parameters = new Dictionary<string, object>
            {
                ["pointsPerCluster"] = pointsPerCluster,
                ["stdDev"] = stdDev,
                ["dimensions"] = dimensions,
                ["seed"] = seed,
                ["reference"] = "Blatt, Wiseman, Domany 1996 PRL 76:3251"
            },
            Metadata = new SyntheticDatasetMeta(
                GeneratorName: nameof(BlattThreeCluster),
                GeometryClass: "Euclidean",
                TopologyTag: "blobs",
                HierarchyTag: "none",
                GTNumClusters: 3,
                AmbientDimensionality: dimensions,
                LiteratureReference: "Blatt, Wiseman, Domany 1996 PRL 76:3251")
        };
    }
}


