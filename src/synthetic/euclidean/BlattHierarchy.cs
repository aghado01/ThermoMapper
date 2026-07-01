using System;
using System.Collections.Generic;
using Synthetic;

namespace Synthetic.Euclidean;

/// <summary>
/// Coarse &gt; medium &gt; fine Gaussian structure. LabelsByLevel exposes
/// ground truth at each scale so temperature-dependent resolution can
/// be measured against all three simultaneously.
/// Reference: Blatt, Wiseman, Domany 1997 Neural Computation 9:1805.
/// </summary>
public static class BlattHierarchy
{
    public static SyntheticDataset Generate(
        int coarseClusters = 2,
        int mediumPerCoarse = 3,
        int finePerMedium = 4,
        int pointsPerFine = 250,       // bumped from 25 → 6000 total, plausible-data scale
        int dimensions = 2,
        double coarseSeparation = 20.0,
        double mediumSeparation = 4.0,
        double fineSeparation = 0.8,
        double leafSpread = 0.15,
        int seed = 42)
    {
        if (dimensions < 2)
            throw new ArgumentException("BlattHierarchy requires dimensions >= 2.");
        var rng = new Random(seed);
        int mediumClusters = coarseClusters * mediumPerCoarse;
        int fineClusters = mediumClusters * finePerMedium;
        int n = fineClusters * pointsPerFine;

        var features = new double[n][];
        var fineLabels = new int[n];
        var mediumLabels = new int[n];
        var coarseLabels = new int[n];
        var coarseCenters2D = SyntheticData.PlaceCentroidsOnSphere(coarseClusters, 2, coarseSeparation, rng);

        int fineIdx = 0;
        int globalIdx = 0;
        for (int coarse = 0; coarse < coarseClusters; coarse++)
        {
            var mediumCenters = new double[mediumPerCoarse][];
            for (int m = 0; m < mediumPerCoarse; m++)
            {
                mediumCenters[m] = new double[dimensions];
                double angle = 2.0 * Math.PI * m / mediumPerCoarse + rng.NextDouble() * 0.3;
                mediumCenters[m][0] = coarseCenters2D[coarse][0] + mediumSeparation * Math.Cos(angle);
                mediumCenters[m][1] = coarseCenters2D[coarse][1] + mediumSeparation * Math.Sin(angle);
            }

            for (int medium = 0; medium < mediumPerCoarse; medium++)
            {
                int mediumLabel = coarse * mediumPerCoarse + medium;
                for (int fine = 0; fine < finePerMedium; fine++)
                {
                    var fineCenter = new double[dimensions];
                    double angle = 2.0 * Math.PI * fine / finePerMedium + rng.NextDouble() * 0.3;
                    fineCenter[0] = mediumCenters[medium][0] + fineSeparation * Math.Cos(angle);
                    fineCenter[1] = mediumCenters[medium][1] + fineSeparation * Math.Sin(angle);

                    for (int p = 0; p < pointsPerFine; p++)
                    {
                        var point = new double[dimensions];
                        point[0] = fineCenter[0] + leafSpread * SyntheticData.SampleStandardNormal(rng);
                        point[1] = fineCenter[1] + leafSpread * SyntheticData.SampleStandardNormal(rng);
                        features[globalIdx] = point;
                        fineLabels[globalIdx] = fineIdx;
                        mediumLabels[globalIdx] = mediumLabel;
                        coarseLabels[globalIdx] = coarse;
                        globalIdx++;
                    }
                    fineIdx++;
                }
            }
        }

        return new SyntheticDataset
        {
            Features = features,
            Labels = fineLabels,
            ClusterCount = fineClusters,
            LabelsByLevel = new int[][] { coarseLabels, mediumLabels, fineLabels },
            Parameters = new Dictionary<string, object>
            {
                ["coarseClusters"] = coarseClusters,
                ["mediumPerCoarse"] = mediumPerCoarse,
                ["finePerMedium"] = finePerMedium,
                ["pointsPerFine"] = pointsPerFine,
                ["coarseSeparation"] = coarseSeparation,
                ["mediumSeparation"] = mediumSeparation,
                ["fineSeparation"] = fineSeparation,
                ["leafSpread"] = leafSpread,
                ["dimensions"] = dimensions,
                ["seed"] = seed,
                ["reference"] = "Blatt, Wiseman, Domany 1997 Neural Computation 9:1805"
            },
            Metadata = new SyntheticDatasetMeta(
                GeneratorName: nameof(BlattHierarchy),
                GeometryClass: "Euclidean",
                TopologyTag: "hierarchical",
                HierarchyTag: "natural",
                GTNumClusters: fineClusters,
                AmbientDimensionality: dimensions,
                LiteratureReference: "Blatt, Wiseman, Domany 1997 Neural Computation 9:1805")
        };
    }
}


