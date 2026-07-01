using System;
using System.Collections.Generic;

namespace Viz;

/// <summary>
/// The complete dataset for a single viz scene. Layers are additive and optional.
/// The scene controller decides which layers to activate and how to render them.
/// </summary>
public sealed class VizDataset
{
    public PointCloud Points { get; }
    public IReadOnlyList<LabelLayer> LabelLayers { get; }
    public IReadOnlyList<NodeSignalLayer> NodeSignalLayers { get; }
    public IReadOnlyList<EdgeLayer> EdgeLayers { get; }
    public IReadOnlyList<TriangleLayer> TriangleLayers { get; }
    public IReadOnlyList<GaussianLayer> GaussianLayers { get; }
    public IReadOnlyList<LineFieldLayer> LineFieldLayers { get; }
    public IReadOnlyList<TemporalLabelSequence> TemporalSequences { get; }
    public IReadOnlyList<SpineLayer> SpineLayers { get; }
    public IReadOnlyList<VectorFieldLayer> VectorFieldLayers { get; }
    public IReadOnlyDictionary<string, object>? GeneratorParams { get; }
    public GeneratorParamSchema? GeneratorParamSchema { get; }

    public VizDataset(
        PointCloud points,
        IReadOnlyList<LabelLayer> labelLayers,
        IReadOnlyList<NodeSignalLayer> nodeSignalLayers,
        IReadOnlyList<EdgeLayer> edgeLayers,
        IReadOnlyList<GaussianLayer> gaussianLayers,
        IReadOnlyList<TemporalLabelSequence> temporalSequences,
        IReadOnlyList<SpineLayer> spineLayers,
        IReadOnlyDictionary<string, object>? generatorParams = null,
        GeneratorParamSchema? generatorParamSchema = null,
        IReadOnlyList<VectorFieldLayer>? vectorFieldLayers = null,
        IReadOnlyList<LineFieldLayer>? lineFieldLayers = null,
        IReadOnlyList<TriangleLayer>? triangleLayers = null)
    {
        Points = points;
        LabelLayers = labelLayers;
        NodeSignalLayers = nodeSignalLayers;
        EdgeLayers = edgeLayers;
        TriangleLayers = triangleLayers ?? Array.Empty<TriangleLayer>();
        GaussianLayers = gaussianLayers;
        TemporalSequences = temporalSequences;
        SpineLayers = spineLayers;
        GeneratorParams = generatorParams;
        GeneratorParamSchema = generatorParamSchema;
        VectorFieldLayers = vectorFieldLayers ?? Array.Empty<VectorFieldLayer>();
        LineFieldLayers = lineFieldLayers ?? Array.Empty<LineFieldLayer>();
    }
}

/// <summary>
/// Ordered snapshots of LabelLayers along a progression axis.
/// Examples: SPC temperature sweep (T-indexed), GMM EM iterations (iter-indexed),
/// recursive GMM depth sweep (depth-indexed).
/// </summary>
public sealed class TemporalLabelSequence
{
    public string Name { get; }
    public TemporalAxis Axis { get; }
    public ReadOnlyMemory<double> IndexValues { get; }
    public IReadOnlyList<LabelLayer> Frames { get; }

    public TemporalLabelSequence(
        string name,
        TemporalAxis axis,
        double[] indexValues,
        IReadOnlyList<LabelLayer> frames)
    {
        Name = name;
        Axis = axis;
        IndexValues = indexValues;
        Frames = frames;
    }
}

public enum TemporalAxis { Temperature, Iteration, Depth, Custom }
