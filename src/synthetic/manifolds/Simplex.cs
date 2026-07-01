using System;
using System.Collections.Generic;
using Synthetic;

namespace Synthetic.Manifolds;

/// <summary>
/// Each data point is a categorical probability vector drawn from a
/// cluster-specific Dirichlet. With disjointSupports = true, each
/// cluster's mass concentrates on a non-overlapping range of
/// categories — producing distribution pairs whose KL / JSD saturate
/// at the metric's maximum. Exposes the stress case for
/// information-geometric metrics and for proximity rules that depend
/// on non-degenerate distance rankings.
/// Reference: Dirichlet-distributed simplex clusters (information-geometric benchmark).
/// </summary>
public static class Simplex
{
    public static SyntheticDataset Generate(
        int clusterCount = 5,
        int pointsPerCluster = 1200,   // bumped from 40 → 6000 total, plausible-data scale
        int categories = 20,
        bool disjointSupports = true,
        double concentration = 20.0,
        int seed = 42)
    {
        var rng = new Random(seed);
        int n = clusterCount * pointsPerCluster;
        var features = new double[n][];
        var labels = new int[n];
        var modes = new double[clusterCount][];

        if (disjointSupports)
        {
            int regionSize = categories / clusterCount;
            if (regionSize < 2)
                throw new ArgumentException(
                    $"categories ({categories}) too small for {clusterCount} disjoint " +
                    $"supports. Need at least 2 * clusterCount categories.");
            for (int c = 0; c < clusterCount; c++)
            {
                var mode = new double[categories];
                int regionStart = c * regionSize;
                for (int k = 0; k < categories; k++)
                    mode[k] = (k >= regionStart && k < regionStart + regionSize)
                        ? 1.0 / regionSize
                        : 1e-6;
                SyntheticData.Normalize(mode);
                modes[c] = mode;
            }
        }
        else
        {
            for (int c = 0; c < clusterCount; c++)
            {
                var mode = new double[categories];
                double peakPos = (double)(c * categories) / clusterCount;
                for (int k = 0; k < categories; k++)
                {
                    double dist = Math.Abs(k - peakPos);
                    mode[k] = Math.Exp(-dist * dist / 4.0);
                }
                SyntheticData.Normalize(mode);
                modes[c] = mode;
            }
        }

        int idx = 0;
        for (int c = 0; c < clusterCount; c++)
        {
            var alpha = new double[categories];
            for (int k = 0; k < categories; k++)
                alpha[k] = concentration * modes[c][k] + 1e-3;
            for (int p = 0; p < pointsPerCluster; p++)
            {
                features[idx] = SyntheticData.SampleDirichlet(alpha, rng);
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
                ["categories"] = categories,
                ["disjointSupports"] = disjointSupports,
                ["concentration"] = concentration,
                ["seed"] = seed
            },
            Metadata = new SyntheticDatasetMeta(
                GeneratorName: nameof(Simplex),
                GeometryClass: "Manifold",
                TopologyTag: "simplex",
                HierarchyTag: "none",
                GTNumClusters: clusterCount,
                AmbientDimensionality: categories,
                LiteratureReference: "Dirichlet-distributed simplex clusters (information-geometric benchmark)",
                FutureMetric: "Wasserstein")
        };
    }
}
