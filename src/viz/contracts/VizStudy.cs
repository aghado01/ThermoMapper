using System.Collections.Immutable;

namespace Viz.Contracts;

/// <summary>
/// Durable research document shared by static replay and live execution.
/// Payload bytes are resolved by a transport implementation using artifact
/// descriptors; the study itself remains storage-neutral.
/// </summary>
public sealed record VizStudy(
    ContractVersion ContractVersion,
    StudyId Id,
    string Title,
    ProducerProvenance Provenance,
    ImmutableArray<RunDescriptor> Runs,
    ImmutableArray<EntitySetDescriptor> EntitySets,
    ImmutableArray<CoordinateSpaceDescriptor> CoordinateSpaces,
    ImmutableArray<CoordinateSetDescriptor> CoordinateSets,
    ImmutableArray<ArtifactDescriptor> Artifacts,
    ImmutableArray<RelationDescriptor> Relations,
    ImmutableArray<PanelDescriptor> Panels,
    ImmutableArray<PanelLinkDescriptor> PanelLinks);
