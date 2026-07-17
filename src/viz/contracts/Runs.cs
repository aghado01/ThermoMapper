using System.Collections.Immutable;

namespace Viz.Contracts;

/// <summary>
/// One pinned recipe execution or experimental branch. The recipe itself is a
/// typed artifact owned by orchestration; this descriptor records its outputs.
/// </summary>
public sealed record RunDescriptor(
    RunId Id,
    string Title,
    ArtifactId Recipe,
    ImmutableArray<ArtifactId> Outputs,
    ProducerProvenance Provenance);
