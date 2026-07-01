namespace Clustering.Graphical.SPC.Export;

/// <summary>
/// Declarative schema for an SPC run's artifact bundle — the canonical file and
/// sub-directory names a single SPC run emits, as <b>pure data</b> (no path logic,
/// no I/O). Bound to a concrete location by <see cref="SpcRunPaths"/> against an
/// <see cref="Archivory.ArtifactScope"/>; the same bundle shape applies whether the
/// run is a standalone root or an arm nested under a larger run.
/// </summary>
/// <remarks>
/// Pure data so it round-trips to JSON (eventually a JSON Schema) once config and
/// code are fully separated. Names are <c>snake_case</c> per the repository artifact
/// convention. Checkpoints are a <b>flat, strategy-agnostic result bag</b> — one
/// directory of per-(T, replica) simulation results, identity carried by the filename
/// (<c>_rep_n</c> vs <c>_final</c>). The scheduler's staging (the former
/// <c>probes/</c>/<c>final/</c> nesting) does not belong on disk; adaptive search
/// lives in the analysis layer reading these results, not in the layout.
/// </remarks>
public sealed record SpcRunLayout
{
    /// <summary>Run manifest (config provenance) — root metadata.</summary>
    public string Manifest { get; init; } = "manifest.json";

    /// <summary>Run summary — root metadata.</summary>
    public string Summary { get; init; } = "summary.json";

    /// <summary>Graph-health report — root metadata.</summary>
    public string Health { get; init; } = "graph_health.json";

    /// <summary>Partition hierarchy (when a hierarchical strategy runs) — root metadata.</summary>
    public string Hierarchy { get; init; } = "partition_hierarchy.json";

    /// <summary>Sub-directory for CSV projections; filenames owned by <see cref="SpcCsvWriter"/>.</summary>
    public string Csv { get; init; } = "csv";

    /// <summary>Sub-directory for tabular (columnar/binary) exports.</summary>
    public string Tabular { get; init; } = "tabular";

    /// <summary>Sub-directory holding the flat bag of SPCX/SPCE simulation results
    /// (identity — replica or final — carried by the filename, not by nesting).</summary>
    public string Checkpoints { get; init; } = "checkpoints";
}
