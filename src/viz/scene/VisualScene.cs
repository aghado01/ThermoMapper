using System.Collections.Immutable;
using Viz.Contracts;

namespace Viz.Scene;

/// <summary>The bounded visual vocabulary understood by spatial renderers.</summary>
public enum MarkKind
{
    Points,
    Segments,
    Polyline,
    Triangles,
    Vectors,
    LineGlyphs,
    Ellipsoids,
    Text,
}

/// <summary>
/// Connects a visual mark to scientific coordinates, optional topology, and the
/// stable entities returned by picking.
/// </summary>
public sealed record GeometryBinding(
    ArtifactId Coordinates,
    ArtifactId? Topology,
    EntitySetId Entities,
    ArtifactId? InstanceToEntityRelation);

/// <summary>
/// One visual channel. A channel uses either a constant canonical-JSON value or
/// a field artifact, optionally interpreted through a named scale.
/// </summary>
public sealed record VisualEncoding(
    SemanticId Channel,
    string? ConstantJson,
    ArtifactId? Field,
    SemanticId? Scale);

/// <summary>A declarative retained visual layer, independent of Three.js objects.</summary>
public sealed record VisualLayer(
    LayerId Id,
    string Title,
    MarkKind Mark,
    GeometryBinding Geometry,
    ImmutableArray<VisualEncoding> Encodings,
    bool Visible,
    int DrawOrder);

/// <summary>Renderer-ready snapshot for one spatial panel.</summary>
public sealed record SceneSnapshot(
    StudyId Study,
    PanelId Panel,
    CoordinateSpaceId DisplaySpace,
    CameraState Camera,
    ImmutableArray<VisualLayer> Layers,
    ImmutableArray<EntityReference> Selection);
