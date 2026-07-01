using System;
using System.Collections.Generic;
using System.IO;

namespace Viz
{
    /// <summary>
    /// A scene is a named rendering configuration over a VizDataset.
    /// Selects which layers are active, what visual encoding to use, and
    /// what the user-facing title and description are.
    /// </summary>
    public sealed class SceneDescriptor
    {
        public string Title { get; init; } = "";
        public string? Description { get; init; }

        // null = activate all layers of that type
        public IReadOnlyList<string>? ActiveLabelLayers { get; init; }
        public IReadOnlyList<string>? ActiveNodeSignalLayers { get; init; }
        public IReadOnlyList<string>? ActiveEdgeLayers { get; init; }
        public IReadOnlyList<string>? ActiveTriangleLayers { get; init; }
        public IReadOnlyList<string>? ActiveGaussianLayers { get; init; }
        public IReadOnlyList<string>? ActiveSpineLayers { get; init; }
        public IReadOnlyList<string>? ActiveVectorFieldLayers { get; init; }
        public IReadOnlyList<string>? ActiveLineFieldLayers { get; init; }
        public string? ActiveTemporalSequence { get; init; }
        public int TemporalStartIndex { get; init; } = 0;

        public SceneRenderHints Hints { get; init; } = SceneRenderHints.Default;
    }

    public sealed class SceneRenderHints
    {
        public static readonly SceneRenderHints Default = new();

        public bool ShowEdgeWeightsAsOpacity { get; init; } = true;
        public bool ShowGaussianEllipsoids { get; init; } = false;
        public bool ShowSpineOverlays { get; init; } = true;
        public bool ShowTangentBases { get; init; } = false;
        // GmmCluster layers color by cluster id, GmmComponent layers by component id.
        // Keeping both simultaneously exposes the component-to-cluster map visually.
        public bool OverlayComponentAndClusterColoring { get; init; } = false;
        public bool AnnotateSpinColorVsEquilibrium { get; init; } = true;
        // Cross-cluster edges rendered in a distinct warning color when GT labels
        // are present on the EdgeLayer. Within-cluster edges use cluster color.
        public bool HighlightFalseBridges { get; init; } = true;
        public bool ShowVectorField { get; init; } = false;
    }

    /// <summary>
    /// Backend-agnostic intermediate produced by SceneBuilder.
    /// Resolved active layers + render hints. IRenderTarget implementations consume this.
    /// </summary>
    public sealed class ScenePackage
    {
        public string Title { get; }
        public PointCloud Points { get; }
        public IReadOnlyList<LabelLayer> ActiveLabelLayers { get; }
        public IReadOnlyList<NodeSignalLayer> ActiveNodeSignalLayers { get; }
        public IReadOnlyList<EdgeLayer> ActiveEdgeLayers { get; }
        public IReadOnlyList<TriangleLayer> ActiveTriangleLayers { get; }
        public IReadOnlyList<GaussianLayer> ActiveGaussianLayers { get; }
        public IReadOnlyList<SpineLayer> ActiveSpineLayers { get; }
        public IReadOnlyList<VectorFieldLayer> ActiveVectorFieldLayers { get; }
        public IReadOnlyList<LineFieldLayer> ActiveLineFieldLayers { get; }
        public TemporalLabelSequence? ActiveSequence { get; }
        public int TemporalFrameIndex { get; }
        public SceneRenderHints Hints { get; }
        public IReadOnlyDictionary<string, object>? GeneratorParams { get; }
        public GeneratorParamSchema? GeneratorParamSchema { get; }

        public ScenePackage(
            PointCloud points,
            IReadOnlyList<LabelLayer> activeLabelLayers,
            IReadOnlyList<NodeSignalLayer> activeNodeSignalLayers,
            IReadOnlyList<EdgeLayer> activeEdgeLayers,
            IReadOnlyList<TriangleLayer> activeTriangleLayers,
            IReadOnlyList<GaussianLayer> activeGaussianLayers,
            IReadOnlyList<SpineLayer> activeSpineLayers,
            TemporalLabelSequence? activeSequence,
            int temporalFrameIndex,
            SceneRenderHints hints,
            IReadOnlyDictionary<string, object>? generatorParams = null,
            GeneratorParamSchema? generatorParamSchema = null,
            IReadOnlyList<VectorFieldLayer>? activeVectorFieldLayers = null,
            IReadOnlyList<LineFieldLayer>? activeLineFieldLayers = null,
            string title = "")
        {
            Title = title;
            Points = points;
            ActiveLabelLayers = activeLabelLayers;
            ActiveNodeSignalLayers = activeNodeSignalLayers;
            ActiveEdgeLayers = activeEdgeLayers;
            ActiveTriangleLayers = activeTriangleLayers;
            ActiveGaussianLayers = activeGaussianLayers;
            ActiveSpineLayers = activeSpineLayers;
            ActiveSequence = activeSequence;
            TemporalFrameIndex = temporalFrameIndex;
            Hints = hints;
            GeneratorParams = generatorParams;
            GeneratorParamSchema = generatorParamSchema;
            ActiveVectorFieldLayers = activeVectorFieldLayers ?? Array.Empty<VectorFieldLayer>();
            ActiveLineFieldLayers = activeLineFieldLayers ?? Array.Empty<LineFieldLayer>();
        }
    }

    /// <summary>
    /// Resolves a SceneDescriptor against a VizDataset, filtering layers
    /// by name where active lists are provided, passing all through when null.
    /// The resulting ScenePackage is the only thing IRenderTarget sees.
    /// </summary>
    public static class SceneBuilder
    {
        public static ScenePackage Build(VizDataset dataset, SceneDescriptor descriptor)
        {
            if (dataset is null) throw new ArgumentNullException(nameof(dataset));
            if (descriptor is null) throw new ArgumentNullException(nameof(descriptor));

            var labelLayers = Filter(dataset.LabelLayers, descriptor.ActiveLabelLayers);
            var nodeSignalLayers = Filter(dataset.NodeSignalLayers, descriptor.ActiveNodeSignalLayers);
            var edgeLayers = Filter(dataset.EdgeLayers, descriptor.ActiveEdgeLayers);
            var triangleLayers = Filter(dataset.TriangleLayers, descriptor.ActiveTriangleLayers);
            var gaussianLayers = Filter(dataset.GaussianLayers, descriptor.ActiveGaussianLayers);
            var spineLayers = Filter(dataset.SpineLayers, descriptor.ActiveSpineLayers);
            var vectorFieldLayers = Filter(dataset.VectorFieldLayers, descriptor.ActiveVectorFieldLayers);
            var lineFieldLayers = Filter(dataset.LineFieldLayers, descriptor.ActiveLineFieldLayers);

            TemporalLabelSequence? activeSeq = null;
            if (descriptor.ActiveTemporalSequence is not null)
            {
                foreach (var seq in dataset.TemporalSequences)
                    if (seq.Name == descriptor.ActiveTemporalSequence) { activeSeq = seq; break; }
            }

            int frameIndex = activeSeq is not null
                ? Math.Clamp(descriptor.TemporalStartIndex, 0, activeSeq.Frames.Count - 1)
                : 0;

            return new ScenePackage(
                points: dataset.Points,
                activeLabelLayers: labelLayers,
                activeNodeSignalLayers: nodeSignalLayers,
                activeEdgeLayers: edgeLayers,
                activeTriangleLayers: triangleLayers,
                activeGaussianLayers: gaussianLayers,
                activeSpineLayers: spineLayers,
                activeSequence: activeSeq,
                temporalFrameIndex: frameIndex,
                hints: descriptor.Hints,
                generatorParams: dataset.GeneratorParams,
                generatorParamSchema: dataset.GeneratorParamSchema,
                activeVectorFieldLayers: vectorFieldLayers,
                activeLineFieldLayers: lineFieldLayers,
                title: descriptor.Title);
        }

        private static IReadOnlyList<T> Filter<T>(
            IReadOnlyList<T> all,
            IReadOnlyList<string>? activeNames) where T : INamedLayer
        {
            if (activeNames is null) return all;

            var nameSet = new HashSet<string>(activeNames, StringComparer.Ordinal);
            var result = new List<T>(activeNames.Count);
            foreach (var item in all)
            {
                if (nameSet.Contains(item.Name))
                    result.Add(item);
            }
            return result;
        }
    }

    /// <summary>
    /// Renders a ScenePackage to a stream. The format is backend-specific.
    /// Implementations: ThreeJsHtmlRenderTarget, PlotlyJsonRenderTarget, SvgRenderTarget.
    /// </summary>
    public interface IRenderTarget
    {
        void Render(ScenePackage scene, Stream output);
    }
}
