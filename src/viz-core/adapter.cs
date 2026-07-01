using System;
using System.Collections.Generic;
using Synthetic;
using Synthetic.Euclidean;

namespace Viz.Adapters.Synthetic
{
    /// <summary>
    /// Produces a VizDataset from a domain-specific source type.
    /// </summary>
    public interface IVizDatasetAdapter<TSource>
    {
        VizDataset Adapt(TSource source);
    }

    /// <summary>
    /// Adapts a SyntheticDataset into a VizDataset.
    ///
    /// Always produces:
    ///   - PointCloud from Features
    ///   - LabelLayer (GroundTruth) from Labels
    ///
    /// Conditionally produces (from ClusterGeometries when present):
    ///   - GaussianLayer "GT Ellipsoids" for any EllipsoidGeometry entries
    ///   - OverlayGeometryLayer (arc spine samples) for any ArcGeometry entries
    ///     carried as a SpineLayer — spine points are injected as a second
    ///     PointCloud-equivalent via ArcSpineLayer (see below)
    ///
    /// Does NOT compute EdgeLayers — that is the diagnostic harness's job.
    /// The adapter's contract is: given a SyntheticDataset, produce the ground
    /// truth view. Edge computation requires metric + proximity rule choices
    /// that live outside this adapter.
    /// </summary>
    public sealed class SyntheticDatasetAdapter : IVizDatasetAdapter<SyntheticDataset>
    {
        public VizDataset Adapt(SyntheticDataset source)
        {
            if (source is null) throw new ArgumentNullException(nameof(source));
            if (source.Features is null || source.Features.Length == 0)
                throw new ArgumentException("SyntheticDataset has no features.", nameof(source));

            int n = source.Features.Length;
            int d = source.Features[0].Length;

            // Validate all rows are same dimension
            for (int i = 1; i < n; i++)
                if (source.Features[i].Length != d)
                    throw new ArgumentException(
                        $"Feature row {i} has dimension {source.Features[i].Length}, expected {d}.");

            // --- PointCloud ---
            var flat = new double[n * d];
            for (int i = 0; i < n; i++)
                source.Features[i].CopyTo(flat, i * d);

            var cloud = new PointCloud(
                flat, n, d,
                label: BuildCloudLabel(source));

            // --- LabelLayer: ground truth ---
            var gtLayer = new LabelLayer(
                "Ground Truth",
                (int[])source.Labels.Clone(),
                LabelLayerKind.GroundTruth,
                id: "ground-truth");

            // --- ArcSpineLayers: one per ArcGeometry or ManifoldGeometry ---
            // Spine samples are a separate point set, not part of the N-point cloud.
            // We carry them as named SpineLayer objects so the renderer can draw
            // the clean generating curve as an overlay without mixing spine points
            // into the data cloud.
            List<SpineLayer> spineLayers = BuildSpineLayers(source);

            // Note: legacy K=1 "GT Ellipsoids (analytic)" and "Best-Fit Gaussian"
            // layers were removed. GMM overlay generation lives entirely in
            // VizApi.BuildPackage, gated by the user-selected GmmMode.

            string? generatorName = source.Parameters.TryGetValue("generator", out var gn)
                ? gn?.ToString()
                : null;

            return new VizDataset(
                points: cloud,
                labelLayers: new List<LabelLayer> { gtLayer },
                nodeSignalLayers: Array.Empty<NodeSignalLayer>(),
                edgeLayers: Array.Empty<EdgeLayer>(),
                gaussianLayers: Array.Empty<GaussianLayer>(),
                temporalSequences: Array.Empty<TemporalLabelSequence>(),
                spineLayers: spineLayers,
                triangleLayers: Array.Empty<TriangleLayer>(),
                generatorParams: source.Parameters.Count > 0
                    ? source.Parameters
                    : null,
                generatorParamSchema: SchemaCatalog.ForGenerator(generatorName));
        }

        // ---------------------------------------------------------------

        private static string BuildCloudLabel(SyntheticDataset source)
        {
            int n = source.Features.Length;
            int d = source.Features[0].Length;
            string genName = source.Parameters.TryGetValue("generator", out var g)
                ? g?.ToString() ?? "Synthetic"
                : "Synthetic";
            return $"{genName} (N={n}, D={d}, K={source.ClusterCount})";
        }

        private static List<SpineLayer> BuildSpineLayers(SyntheticDataset source)
        {
            var result = new List<SpineLayer>();
            if (source.ClusterGeometries is null) return result;

            for (int i = 0; i < source.ClusterGeometries.Length; i++)
            {
                switch (source.ClusterGeometries[i])
                {
                    case ArcGeometry arc:
                        result.Add(new SpineLayer(
                            name: $"Cluster {i} Arc Spine",
                            clusterIdx: i,
                            kind: SpineLayerKind.Arc,
                            spineSamples: arc.SpineSamples,
                            tangentBases: null,
                            typicalScale: arc.NoiseScale,
                            id: $"cluster-{i}-arc-spine"));
                        break;

                    case ManifoldGeometry mf:
                        result.Add(new SpineLayer(
                            name: $"Cluster {i} Manifold Spine",
                            clusterIdx: i,
                            kind: SpineLayerKind.Manifold,
                            spineSamples: mf.SpineSamples,
                            tangentBases: mf.TangentBases,
                            typicalScale: 0.0,
                            id: $"cluster-{i}-manifold-spine"));
                        break;

                    case MobiusTubeGeometry mob:
                        // LocalFrames[i] = [T, N, B], each double[3] — matches tangentBases M×3×3 layout.
                        result.Add(new SpineLayer(
                            name: $"Cluster {i} Möbius Spine",
                            clusterIdx: i,
                            kind: SpineLayerKind.MobiusTube,
                            spineSamples: mob.SpineSamples,
                            tangentBases: mob.LocalFrames,
                            typicalScale: Math.Sqrt(mob.HalfWidth * mob.HalfThickness),
                            id: $"cluster-{i}-mobius-spine"));
                        break;

                    // EllipsoidGeometry: rendered via GaussianLayer, no spine needed
                    default:
                        break;
                }
            }

            return result;
        }
    }
}
