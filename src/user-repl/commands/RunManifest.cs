using System;
using System.Collections.Generic;
using System.IO;
using System.Text.Json;
using Archivory.Jso;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs;

namespace UserRepl.Commands;

/// <summary>
/// Self-describing artifact written at the root of every CLI run
/// directory. Captures the dataset source, graph construction inputs,
/// sweep/algorithm parameters, and output layout so the
/// <see cref="ExtractCommand"/> (and future tooling) can reconstruct
/// the run end-to-end without re-passing flags.
/// </summary>
/// <remarks>
/// <para><b>Schema versioning.</b> <see cref="SchemaVersion"/> is bumped
/// when the on-disk layout changes in a way that older readers couldn't
/// understand. Readers should warn (not fail) on a newer schema and
/// degrade gracefully when possible.</para>
///
/// <para><b>What's NOT captured.</b> Per-task results
/// (<c>.spcx</c>/<c>.spce</c>) and the post-run CSVs live alongside the
/// manifest but are addressed by convention, not by manifest entries —
/// adding files to the run dir doesn't require manifest churn.</para>
/// </remarks>
public sealed record RunManifest(
    string                    SchemaVersion,
    DateTime                  CreatedUtc,
    string                    Algorithm,        // "spc" | "hdbscan"
    string                    CommandLine,
    DatasetSpec               Dataset,
    GraphSpec?                Graph,            // null for hdbscan
    SweepSpec?                Sweep,            // null for hdbscan
    HdbscanSpec?              Hdbscan,          // null for spc
    OutputSpec                Output,
    RunIdentitySpec?          Identity = null)  // null in pre-Phase-C manifests
{
    public const string CurrentSchemaVersion = "1.0";
    public const string FileName = "manifest.json";

    public static string PathFor(string runDirectory) => Path.Combine(runDirectory, FileName);

    public void WriteTo(string runDirectory)
    {
        string path = PathFor(runDirectory);
        Directory.CreateDirectory(runDirectory);
        File.WriteAllText(path, JsonSerializer.Serialize(this, JsonOptions));
    }

    public static RunManifest ReadFrom(string runDirectory)
    {
        string path = PathFor(runDirectory);
        if (!File.Exists(path))
            throw new FileNotFoundException($"No manifest at {path}.", path);
        string json = File.ReadAllText(path);
        return JsonSerializer.Deserialize<RunManifest>(json, JsonOptions)
            ?? throw new InvalidDataException($"Manifest at {path} deserialized to null.");
    }

    public static readonly JsonSerializerOptions JsonOptions =
        JsonArtifactConventions.Create(UserReplJsonContext.Default,
            indented: true, snakeCase: true, allowNamedFloatingPointLiterals: false);
}

/// <summary>Dataset provenance. Either synthetic-generator info or CSV path is set.</summary>
public sealed record DatasetSpec(
    string                              Source,                 // "synthetic" | "csv"
    string?                             GeneratorName,
    IReadOnlyDictionary<string, object?>? GeneratorParameters,
    int?                                Seed,                   // synthetic-only
    string?                             CsvPath,
    string?                             LabelColumn,
    bool?                               HasHeader,
    string?                             Delimiter);             // serialized as string for readability

/// <summary>Graph construction inputs sufficient to deterministically rebuild the CsrGraph.
/// The L3 <see cref="Projection"/> (Distance / Coupling / Affinity) is embedded directly and
/// carries its own kernel + LMP + bandwidth — no flattened kernel fields, no string
/// discriminator to keep in sync. Topology/filter/repair stay flat (enums/scalars, not
/// unions). Round-trips polymorphically via the <c>kind</c> discriminator on
/// <see cref="IEdgeProjection"/>.</summary>
public sealed record GraphSpec(
    string?          TopologyKind = null,
    string?          FilterKind = null,
    int              K = 0,
    double           Epsilon = 0.0,
    string?          DistanceMetric = null,
    bool             EnsureConnected = false,
    IEdgeProjection? Projection = null);

/// <summary>SPC sweep + cut configuration.</summary>
public sealed record SweepSpec(
    string                              Schedule,               // "FixedGrid" (adaptive scheduling parked)
    string?                             TemperaturesSpec,       // raw --temperatures string for FixedGrid
    int                                 Replicas,
    RunBudget                           SweepBudget,
    RunBudget                           EquilibriumBudget,
    int                                 Q,
    string                              Analyzer,
    string                              PartitionStrategy,
    double                              Theta,
    string?                             TemperaturesResolved = null);

/// <summary>HDBSCAN configuration.</summary>
public sealed record HdbscanSpec(
    int                                 MinPts,
    int?                                MinClusterSize,
    bool                                AllowSingleCluster,
    string                              DistanceMetric);

/// <summary>Run-directory layout (where the artifacts live).</summary>
public sealed record OutputSpec(
    string                              RunDirectory,
    string?                             CheckpointDirectory);

/// <summary>How the run's family folder name was chosen — requested-vs-resolved
/// provenance mirroring the <c>Archivory.RunIdentity</c> resolution.</summary>
public sealed record RunIdentitySpec(
    string                              Family,
    string                              Source,         // "explicit" | "auto:caller={stub}"
    string?                             Requested);
