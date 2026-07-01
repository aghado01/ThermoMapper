using System;
using System.Collections.Generic;
using Synthetic;

namespace Synthetic.Euclidean;

/// <summary>
/// Two interleaved half-moons. Non-convex shape; methods that assume
/// convex or spherical clusters will struggle. Matches parameters of
/// sklearn.datasets.make_moons.
/// Reference: scikit-learn make_moons benchmark / classic two moons dataset.
/// </summary>
public static class TwoMoons
{
    public static SyntheticDataset Generate(
        int pointsPerMoon = 3000,   // bumped from 100 → 6000 total, plausible-data scale
        double noise = 0.1,
        int dimensions = 2,
        int seed = 42)
    {
        if (dimensions < 2)
            throw new ArgumentException("TwoMoons requires dimensions >= 2.");
        var rng = new Random(seed);
        int n = 2 * pointsPerMoon;
        var features = new double[n][];
        var labels = new int[n];

        int idx = 0;
        for (int p = 0; p < pointsPerMoon; p++)
        {
            double t = Math.PI * p / (pointsPerMoon - 1);
            var point = new double[dimensions];
            point[0] = Math.Cos(t) + noise * SyntheticData.SampleStandardNormal(rng);
            point[1] = Math.Sin(t) + noise * SyntheticData.SampleStandardNormal(rng);
            features[idx] = point;
            labels[idx] = 0;
            idx++;
        }
        for (int p = 0; p < pointsPerMoon; p++)
        {
            double t = Math.PI * p / (pointsPerMoon - 1);
            var point = new double[dimensions];
            point[0] = 1.0 - Math.Cos(t) + noise * SyntheticData.SampleStandardNormal(rng);
            point[1] = 0.5 - Math.Sin(t) + noise * SyntheticData.SampleStandardNormal(rng);
            features[idx] = point;
            labels[idx] = 1;
            idx++;
        }

        return new SyntheticDataset
        {
            Features = features,
            Labels = labels,
            ClusterCount = 2,
            LabelsByLevel = new[] { labels },
            Parameters = new Dictionary<string, object>
            {
                ["pointsPerMoon"] = pointsPerMoon,
                ["noise"] = noise,
                ["dimensions"] = dimensions,
                ["seed"] = seed
            },
            Metadata = new SyntheticDatasetMeta(
                GeneratorName: nameof(TwoMoons),
                GeometryClass: "Euclidean",
                TopologyTag: "non-convex",
                HierarchyTag: "none",
                GTNumClusters: 2,
                AmbientDimensionality: dimensions,
                LiteratureReference: "scikit-learn make_moons benchmark / classic two moons dataset")
        };
    }
}


