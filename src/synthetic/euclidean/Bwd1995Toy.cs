using System;
using System.Collections.Generic;
using Synthetic;

namespace Synthetic.Euclidean;

/// <summary>
/// BWD1995 Fig. 1 toy: three dense horizontal stripes (800 points each,
/// uniform) over a sparse uniform background (800 points) — the canonical
/// first SPC data set (Blatt, Wiseman, Domany 1995).
/// </summary>
/// <remarks>
/// <para><b>Geometry provenance.</b> Read off the figure plate
/// (codex-scientiae <c>compendia/clustering/images/BWD1995/imageFile1.png</c>;
/// the docling text-repair is unreliable for figures): background spans
/// <c>[-10,10]×[-3,3]</c>; the three stripes share x ∈ <c>[-7, 7.5]</c> at
/// y ≈ +2 / 0 / −2, height 0.8 each (equal density — the middle stripe only
/// looks taller on the plate because square markers render fatter).</para>
///
/// <para><b>Oracle arithmetic.</b> Each stripe captures background points in
/// proportion to its share of the frame area (~0.097 of 800 ≈ 77), so the
/// expected top-3 cluster sizes are ≈ 877 each — matching the published
/// 900 / 894 / 877 at T_clus = 0.08, θ = 0.5 (plus a 4th cluster of size 2
/// and singletons). The stripes' ~18:1 aspect is part of the point:
/// elongated uniform clusters at ~9:1 density contrast, and equal density
/// keeps all three locked through one shared SP plateau.</para>
/// </remarks>
public static class Bwd1995Toy
{
    public static SyntheticDataset Generate(
        int seed = 42)
    {
        var rng = new Random(seed);
        int nTotal = 3200;
        var features = new double[nTotal][];
        var labels = new int[nTotal];

        int idx = 0;

        // Helper to generate uniform points in a rectangle [xMin, xMax] x [yMin, yMax]
        void AddRectangle(double xMin, double xMax, double yMin, double yMax, int count, int label)
        {
            for (int i = 0; i < count; i++)
            {
                features[idx] = new[]
                {
                    xMin + rng.NextDouble() * (xMax - xMin),
                    yMin + rng.NextDouble() * (yMax - yMin)
                };
                labels[idx] = label;
                idx++;
            }
        }

        // Sparse background over the full frame (label 0).
        AddRectangle(-10, 10, -3, 3, 800, 0);

        // Three dense stripes, top to bottom (labels 1–3). Equal heights:
        // the middle stripe LOOKS taller on the plate but that's marker bias
        // (squares render fatter than crosses/x's); equal density is the
        // design intent — cf. the papers' "all of about the same density",
        // and an unequal-density stripe melts before the others, breaking
        // the SP plateau the method clusters in.
        AddRectangle(-7, 7.5,  1.6,  2.4, 800, 1);
        AddRectangle(-7, 7.5, -0.4,  0.4, 800, 2);
        AddRectangle(-7, 7.5, -2.4, -1.6, 800, 3);

        return new SyntheticDataset
        {
            Features = features,
            Labels = labels,
            ClusterCount = 4, // 3 dense stripes + 1 background
            LabelsByLevel = new[] { labels },
            Parameters = new Dictionary<string, object>
            {
                ["seed"] = seed,
                ["frame"] = "[-10,10]x[-3,3]",
                ["stripe_x"] = "[-7,7.5]",
                ["stripe_y_centers"] = "+2,0,-2",
                ["stripe_heights"] = "0.8,0.8,0.8",
                ["points_per_region"] = 800
            },
            Metadata = new SyntheticDatasetMeta(
                GeneratorName: nameof(Bwd1995Toy),
                GeometryClass: "Euclidean",
                TopologyTag: "stripes",
                HierarchyTag: "none",
                GTNumClusters: 4,
                AmbientDimensionality: 2,
                LiteratureReference: "Blatt, Wiseman, Domany 1995 - Fig 1")
        };
    }
}
