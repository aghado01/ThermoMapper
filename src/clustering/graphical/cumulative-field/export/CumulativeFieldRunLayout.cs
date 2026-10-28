namespace Clustering.Graphical.CumulativeField.Export;

/// <summary>
/// Declarative schema for a cumulative-field run's artifact bundle — the
/// canonical file and sub-directory names one run emits, as <b>pure data</b> (no
/// path logic, no I/O). Bound to a concrete location by
/// <see cref="CumulativeFieldRunPaths"/> against an <see cref="Archivory.ArtifactScope"/>.
/// The deterministic solver has no checkpoints and no tabular sweep, so the bundle
/// is root metadata plus one CSV sub-directory.
/// </summary>
public sealed record CumulativeFieldRunLayout
{
    /// <summary>Run manifest (config provenance) — root metadata.</summary>
    public string Manifest { get; init; } = "manifest.json";

    /// <summary>Run summary — root metadata.</summary>
    public string Summary { get; init; } = "summary.json";

    /// <summary>Graph-health report — root metadata.</summary>
    public string Health { get; init; } = "graph_health.json";

    /// <summary>Sub-directory for CSV projections; filenames owned by <see cref="CumulativeFieldCsvWriter"/>.</summary>
    public string Csv { get; init; } = "csv";
}
