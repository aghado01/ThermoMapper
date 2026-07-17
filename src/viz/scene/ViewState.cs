using System.Collections.Immutable;
using Viz.Contracts;

namespace Viz.Scene;

/// <summary>Durable per-layer state saved independently of scientific artifacts.</summary>
public sealed record LayerViewState(
    LayerId Layer,
    bool Visible,
    ImmutableArray<VisualEncoding> Encodings);

/// <summary>Durable view state for a panel; transient hover/drag state is excluded.</summary>
public sealed record PanelViewState(
    PanelId Panel,
    ArtifactId DisplayCoordinates,
    CameraState Camera,
    ImmutableArray<LayerViewState> Layers,
    ImmutableArray<EntityReference> Selection,
    SemanticId? ScientificAxis,
    double? ScientificAxisCursor);
