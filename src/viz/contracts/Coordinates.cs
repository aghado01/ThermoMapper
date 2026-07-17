namespace Viz.Contracts;

/// <summary>A declared geometric coordinate space.</summary>
public sealed record CoordinateSpaceDescriptor(
    CoordinateSpaceId Id,
    SemanticId Geometry,
    int AmbientDimension,
    int? IntrinsicDimension,
    string? Units);

/// <summary>How one coordinate artifact was projected from another.</summary>
public sealed record ProjectionProvenance(
    SemanticId Method,
    ArtifactId SourceCoordinates,
    ProducerProvenance Producer);

/// <summary>Coordinates for an entity set in an explicit geometric space.</summary>
public sealed record CoordinateSetDescriptor(
    ArtifactId Artifact,
    EntitySetId Entities,
    CoordinateSpaceId Space,
    EvidenceRole EvidenceRole,
    ProjectionProvenance? Projection);
