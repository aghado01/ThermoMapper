using System.Collections.Immutable;
using Viz.Contracts;

namespace Viz.Scene;

/// <summary>A typed view-only action. Recipe interventions use a different contract.</summary>
public abstract record ViewAction;

public enum SelectionMode
{
    Replace,
    Add,
    Remove,
    Toggle,
}

public sealed record SelectEntities(
    ImmutableArray<EntityReference> Entities,
    SelectionMode Mode) : ViewAction;

public sealed record SetLayerVisibility(LayerId Layer, bool Visible) : ViewAction;

public sealed record SetLayerEncodings(
    LayerId Layer,
    ImmutableArray<VisualEncoding> Encodings) : ViewAction;

public sealed record SetCamera(CameraState Camera) : ViewAction;

public sealed record SetDisplayCoordinates(ArtifactId Coordinates) : ViewAction;

public sealed record SetScientificAxisCursor(
    SemanticId Axis,
    double Value) : ViewAction;

public sealed record FrameSelection : ViewAction;
