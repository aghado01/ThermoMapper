using System.Collections.Generic;
using System.Text.Json;
using System.Text.Json.Serialization;
using Archivory.Jso;
using Graphs;
using Graphs.Diagnostics;
using Graphs.Observables;
using Graphs.Distance;
using Clustering.Graphical.SPC;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Partitions.Hierarchical;
using UserRepl.Commands;

namespace UserRepl.Commands;

[JsonSerializable(typeof(GraphHealthReport))]
[JsonSerializable(typeof(RunManifest))]
[JsonSerializable(typeof(GraphConstructionManifest))]
[JsonSerializable(typeof(GraphCompilerConfig))]
[JsonSerializable(typeof(MetricProperties))]
// Registered explicitly (it is also reachable via GraphCompilerConfig) so the
// config fingerprint can serialize the projection polymorphically through the
// interface — deriving the fingerprint from the canonical form rather than a
// hand-rolled kernel switch.
[JsonSerializable(typeof(IEdgeProjection))]
[JsonSerializable(typeof(DiagnosticsLog))]
[JsonSerializable(typeof(DiagnosticMessage))]
[JsonSerializable(typeof(HdbscanPreset))]
[JsonSerializable(typeof(SpcPreset))]
[JsonSerializable(typeof(PartitionHierarchy))]
[JsonSerializable(typeof(SpcRunSummary))]
[JsonSerializable(typeof(SpcSummaryPartition))]
[JsonSerializable(typeof(SpcSummaryRun))]
[JsonSerializable(typeof(HdbscanSummaryPayload))]
[JsonSerializable(typeof(HdbscanRunPayload))]
[JsonSerializable(typeof(HdbscanResultPayload))]
[JsonSerializable(typeof(HdbscanRunPaths))]
[JsonSerializable(typeof(HdbscanDendrogramPayload))]
[JsonSerializable(typeof(HdbscanDendrogramMerge))]
[JsonSerializable(typeof(HdbscanClusterStat))]
internal partial class UserReplJsonContext : JsonSerializerContext
{
}

internal static class UserReplJson
{
    // Option profiles for this assembly's artifacts. Conventions (naming, named-float,
    // indentation) are defined once in Archivory.Jso.JsonArtifactConventions; these are
    // just named flag-combinations bound to UserRepl's source-gen context.
    public static readonly JsonSerializerOptions DefaultOptions =
        JsonArtifactConventions.Create(UserReplJsonContext.Default,
            indented: false, snakeCase: false, allowNamedFloatingPointLiterals: false);

    public static readonly JsonSerializerOptions IndentedOptions =
        JsonArtifactConventions.Create(UserReplJsonContext.Default,
            indented: true, snakeCase: false, allowNamedFloatingPointLiterals: false);

    public static readonly JsonSerializerOptions IndentedSnakeCaseAllowNamedFloatingPointLiteralsOptions =
        JsonArtifactConventions.Create(UserReplJsonContext.Default,
            indented: true, snakeCase: true, allowNamedFloatingPointLiterals: true);

    /// <summary>Shared writer over the canonical (indented + snake_case + named-float)
    /// profile — same conventions as
    /// <see cref="IndentedSnakeCaseAllowNamedFloatingPointLiteralsOptions"/>, plus atomic
    /// file writes. Use for emitting JSON artifacts (summaries, sidecars, manifests).</summary>
    public static readonly JsonArtifactWriter Writer = new(UserReplJsonContext.Default);
}

internal sealed record SpcRunSummary(
    IReadOnlyDictionary<string, object?>? Dataset,
    GraphCompilerConfig Graph,
    string Analyzer,
    string PartitionStrategy,
    int[]? ReferenceLabels,
    SpcSummaryPartition Partition,
    SpcSummaryRun Run);

internal sealed record SpcSummaryPartition(
    int ClusterCount,
    IReadOnlyDictionary<string, double> EvaluatorScores);

internal sealed record SpcSummaryRun(
    string RunDirectory,
    string SweepCsv,
    string PartitionCsv,
    string CriteriaCsv,
    string SessionCsv,
    string ReplicaTracesCsv,
    string? FinalEdgesCsv);

internal sealed record HdbscanSummaryPayload(
    string Algorithm,
    IReadOnlyDictionary<string, object?>? Dataset,
    HdbscanRunPayload Hdbscan,
    int[]? ReferenceLabels,
    HdbscanResultPayload Result,
    HdbscanRunPaths Run);

internal sealed record HdbscanRunPayload(
    int MinPts,
    int? MinClusterSize,
    bool AllowSingleCluster,
    string? Metric);

internal sealed record HdbscanResultPayload(
    int ClusterCount,
    int NoiseCount,
    IReadOnlyDictionary<string, double> EvaluatorScores,
    IReadOnlyList<HdbscanClusterStat> Clusters);

internal sealed record HdbscanRunPaths(
    string RunDirectory,
    string PartitionCsv,
    string DendrogramJson);

internal sealed record HdbscanClusterStat(
    int ClusterId,
    int Size,
    double MeanMembershipProbability);

internal sealed record HdbscanDendrogramPayload(
    int LeafCount,
    string CostAxis,
    HdbscanDendrogramMerge[] Merges);

internal sealed record HdbscanDendrogramMerge(
    int LeftChild,
    int RightChild,
    double Distance,
    int Size,
    double Lambda);
