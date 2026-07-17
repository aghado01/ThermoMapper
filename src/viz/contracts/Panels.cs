using System.Collections.Immutable;

namespace Viz.Contracts;

/// <summary>A view intent over a set of study artifacts.</summary>
public sealed record PanelDescriptor(
    PanelId Id,
    SemanticId Kind,
    string Title,
    ImmutableArray<ArtifactId> Artifacts);

/// <summary>A coordination policy between two or more panels.</summary>
public sealed record PanelLinkDescriptor(
    SemanticId Id,
    SemanticId Kind,
    ImmutableArray<PanelId> Panels);
