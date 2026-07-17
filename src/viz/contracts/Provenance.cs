using System.Collections.Immutable;

namespace Viz.Contracts;

/// <summary>The epistemic role an artifact plays in a study.</summary>
public enum EvidenceRole
{
    Observed,
    Oracle,
    Counterfactual,
    DiagnosticDerived,
    PresentationDerived,
}

/// <summary>A configuration value encoded as canonical JSON.</summary>
public sealed record ParameterValue(string Name, string CanonicalJson);

/// <summary>Immutable producer and input provenance attached to an artifact.</summary>
public sealed record ProducerProvenance(
    SemanticId Producer,
    string ProducerVersion,
    ImmutableArray<ArtifactId> Inputs,
    ImmutableArray<ParameterValue> Parameters,
    string Fingerprint,
    ImmutableArray<string> Warnings);
