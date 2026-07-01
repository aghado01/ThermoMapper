using System;
using System.Collections.Generic;
using System.IO;
using Clustering.Primitives;
using Clustering.Graphical.SPC.Partitions;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs.Primitives;
using Archivory.Tabular;

namespace Clustering.Graphical.SPC.Export;

/// <summary>
/// Canonical SPC CSV writer. Owns the per-output filename conventions and
/// the convenience batch shape ("write all standard CSVs to one directory").
/// Each <c>Write*</c> method is intentionally thin: pull a projection from
/// <see cref="SpcTabularProjections"/>, write atomically, return the path.
/// </summary>
/// <remarks>
/// <para>Two-tier access pattern:</para>
/// <list type="bullet">
///   <item>Use these methods when you want the SPC-canonical filename and
///     default CSV shape — the common case.</item>
///   <item>Use <see cref="SpcTabularProjections"/> directly when you need a
///     non-canonical filename, alternate CSV shape, or want to compose the
///     projection with other tabular machinery.</item>
/// </list>
/// <para>Atomic writes are provided by
/// <see cref="TabularProjectionExtensions.WriteToFile"/>; partial writes
/// never leave a half-written file at the canonical path.</para>
/// </remarks>
public static class SpcCsvWriter
{
    // ── Canonical SPC CSV filenames (single source of truth) ─────────────
    public const string SweepFileName             = "spc_sweep.csv";
    public const string PartitionFileName         = "spc_partition.csv";
    public const string CriteriaFileName          = "spc_criteria.csv";
    public const string SessionFileName           = "spc_session.csv";
    public const string PartitionScheduleFileName = "spc_partition_schedule.csv";
    public const string DatasetFileName           = "spc_dataset.csv";
    public const string AnalysisFileName          = "spc_analysis.csv";
    public const string ReplicaTracesFileName     = "spc_replica_traces.csv";
    public const string EquilibriumEdgesFileName        = "spc_equilibrium_edges.csv";

    public static string WriteSweepProfile(SweepProfile profile, string path)
    {
        ArgumentNullException.ThrowIfNull(profile);
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path cannot be null or whitespace.", nameof(path));

        SpcTabularProjections.CreateSweepProfileProjection(profile, "spc_sweep").WriteToFile(path);
        return path;
    }

    public static string WritePartition(Assignment partition, string path, double[][]? features = null, int[]? trueLabels = null)
    {
        ArgumentNullException.ThrowIfNull(partition);
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path cannot be null or whitespace.", nameof(path));

        SpcTabularProjections.CreatePartitionProjection(partition, features, trueLabels, "spc_partition").WriteToFile(path);
        return path;
    }

    public static string WriteCriteria(ProfileCriteria criteria, string path)
    {
        ArgumentNullException.ThrowIfNull(criteria);
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path cannot be null or whitespace.", nameof(path));

        SpcTabularProjections.CreateCriteriaProjection(criteria, "spc_criteria").WriteToFile(path);
        return path;
    }

    public static string WriteSessionSummary(SpcSessionResult result, string path)
    {
        ArgumentNullException.ThrowIfNull(result);
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path cannot be null or whitespace.", nameof(path));

        SpcTabularProjections.CreateSessionSummaryProjection(result, "spc_session").WriteToFile(path);
        return path;
    }

    public static string WritePartitionScheduleRollups(
        IReadOnlyList<SchedulePartitionRollup> rollups,
        string path)
    {
        ArgumentNullException.ThrowIfNull(rollups);
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path cannot be null or whitespace.", nameof(path));

        SpcTabularProjections.CreatePartitionScheduleRollupProjection(rollups, "spc_partition_schedule").WriteToFile(path);
        return path;
    }

    public static string WriteDataset(double[][] features, int[] labels, string path)
    {
        ArgumentNullException.ThrowIfNull(features);
        ArgumentNullException.ThrowIfNull(labels);
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path cannot be null or whitespace.", nameof(path));

        SpcTabularProjections.CreateDatasetProjection(features, labels, "spc_dataset").WriteToFile(path);
        return path;
    }

    /// <summary>
    /// One row per <see cref="SpcRunResult"/> frame — preserves the
    /// (T, replica) trace pre-averaging so variance-band plotting in
    /// Python is straightforward. See
    /// <see cref="SpcTabularProjections.CreateReplicaTracesProjection"/>
    /// for the column schema.
    /// </summary>
    public static string WriteReplicaTraces(IReadOnlyList<SpcRunResult> runs, string path)
    {
        ArgumentNullException.ThrowIfNull(runs);
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path cannot be null or whitespace.", nameof(path));

        SpcTabularProjections.CreateReplicaTracesProjection(runs, "spc_replica_traces").WriteToFile(path);
        return path;
    }

    /// <summary>
    /// One row per undirected edge at the chosen-T equilibrium, carrying
    /// source/target/weight + bond_frequency + spin_agreement. The flat
    /// counterpart to the binary <c>.spce</c> sidecar — preferred for
    /// Python-side edge-weight inspection. Returns the written path; or
    /// <see langword="null"/> when <paramref name="affinities"/> is null
    /// (e.g. a Standard-tier final pass that didn't request edge
    /// observables).
    /// </summary>
    public static string? WriteEquilibriumEdges(
        CsrGraph graph,
        Affinities? affinities,
        Alignments? alignments,
        string path,
        CoMembership? coMembership = null)
    {
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path cannot be null or whitespace.", nameof(path));

        if (affinities is null || affinities.G.Length == 0) return null;

        SpcTabularProjections.CreateEquilibriumEdgesProjection(graph, affinities, alignments, coMembership: coMembership).WriteToFile(path);
        return path;
    }

    /// <summary>
    /// Write the canonical SPC CSV bundle for one session to
    /// <paramref name="runDirectory"/>. Each CSV is written atomically with
    /// its canonical filename (see the <c>*FileName</c> constants). Optional
    /// inputs are only written when supplied.
    /// </summary>
    /// <remarks>
    /// Callers typically create the run directory first via
    /// <see cref="SpcOutputPathHelper.CreateRunDirectory"/>, which applies
    /// the timestamped <c>{yyyyMMdd}_{HHmmss}__{runName}__{guid}</c>
    /// convention so concurrent or sequential runs do not overwrite one
    /// another.
    /// </remarks>
    public static string WriteAllToDirectory(
        SpcSessionResult result,
        string runDirectory,
        double[][]? features = null,
        int[]? trueLabels = null,
        IReadOnlyList<SchedulePartitionRollup>? partitionScheduleRollups = null)
    {
        ArgumentNullException.ThrowIfNull(result);
        if (string.IsNullOrWhiteSpace(runDirectory))
            throw new ArgumentException("Run directory must be provided.", nameof(runDirectory));

        Directory.CreateDirectory(runDirectory);

        WriteSweepProfile(result.Profile, Path.Combine(runDirectory, SweepFileName));
        WritePartition(result.Partition, Path.Combine(runDirectory, PartitionFileName), features, trueLabels);
        WriteCriteria(result.ProfileCriteria, Path.Combine(runDirectory, CriteriaFileName));
        WriteSessionSummary(result, Path.Combine(runDirectory, SessionFileName));

        if (partitionScheduleRollups is not null)
            WritePartitionScheduleRollups(partitionScheduleRollups, Path.Combine(runDirectory, PartitionScheduleFileName));

        if (features is not null && trueLabels is not null)
            WriteDataset(features, trueLabels, Path.Combine(runDirectory, DatasetFileName));

        return runDirectory;
    }
}
