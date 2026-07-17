using System.Collections.Immutable;

namespace Viz.Contracts;

/// <summary>A set of addressable scientific entities.</summary>
public sealed record EntitySetDescriptor(
    EntitySetId Id,
    SemanticId EntityKind,
    string DisplayName,
    long Count);

/// <summary>
/// Storage-neutral description of an artifact payload. Transport interprets the
/// encoding, location, checksum, data type, and shape.
/// </summary>
public sealed record ArtifactPayloadDescriptor(
    SemanticId Encoding,
    SemanticId DataType,
    ImmutableArray<long> Shape,
    string? Location,
    string? Checksum);

/// <summary>One typed, provenance-bearing artifact registered in a study.</summary>
public sealed record ArtifactDescriptor(
    ArtifactId Id,
    SemanticId SemanticType,
    EvidenceRole EvidenceRole,
    ImmutableArray<EntitySetId> Subjects,
    ProducerProvenance Provenance,
    ArtifactPayloadDescriptor Payload);

/// <summary>A relation whose payload maps entities in one set to another.</summary>
public sealed record RelationDescriptor(
    ArtifactId Artifact,
    SemanticId SemanticType,
    EntitySetId Domain,
    EntitySetId Codomain,
    EvidenceRole EvidenceRole,
    ProducerProvenance Provenance);
