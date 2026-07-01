using System;
using System.Globalization;
using System.IO;
using Clustering.Graphical.SPC.Runtime.Execution;

namespace Clustering.Graphical.SPC.Export;

/// <summary>
/// Directory-layout helpers for SPC runs. Owns the on-disk structure
/// (timestamped run-directory naming, <c>csv/</c> subdirectory,
/// <c>checkpoints/</c> subdirectory) but defers per-file naming to
/// <see cref="SpcCsvWriter"/>, which is the single source of truth for
/// canonical filenames.
/// </summary>
public static class SpcOutputPathHelper
{
    public const string CheckpointDirectoryName = "checkpoints";
    public const string CsvDirectoryName        = "csv";

    /// <summary>
    /// Create a unique run directory under <paramref name="baseDirectory"/>
    /// using the <c>{yyyyMMdd}_{HHmmss}__{runName}__{guid}</c> pattern so
    /// concurrent or sequential runs do not overwrite one another. Set
    /// <paramref name="includeGuid"/> to <see langword="false"/> for a
    /// shorter <c>{yyyyMMdd}_{HHmmss}__{runName}</c> shape when the
    /// timestamp alone is unique enough.
    /// </summary>
    public static string CreateRunDirectory(string baseDirectory, string runName, bool includeGuid = true)
    {
        if (string.IsNullOrWhiteSpace(baseDirectory))
            throw new ArgumentException("Base directory must be provided.", nameof(baseDirectory));

        string normalizedBase = Path.GetFullPath(baseDirectory);
        string timestamp = DateTime.Now.ToString("yyyyMMdd_HHmmss", CultureInfo.InvariantCulture);
        string sanitizedRunName = SanitizePathSegment(runName);

        string directoryName = includeGuid
            ? $"{timestamp}__{sanitizedRunName}__{Guid.NewGuid():N}"
            : $"{timestamp}__{sanitizedRunName}";

        string runDirectory = Path.Combine(normalizedBase, directoryName);
        Directory.CreateDirectory(runDirectory);
        return runDirectory;
    }

    public static string GetCheckpointDirectory(string runDirectory)
        => Path.Combine(runDirectory, CheckpointDirectoryName);

    public static string GetCsvDirectory(string runDirectory)
    {
        string csvDirectory = Path.Combine(runDirectory, CsvDirectoryName);
        Directory.CreateDirectory(csvDirectory);
        return csvDirectory;
    }

    public static string GetCheckpointPath(string runDirectory, double temperature, int replicaIndex)
    {
        string checkpointDirectory = GetCheckpointDirectory(runDirectory);
        Directory.CreateDirectory(checkpointDirectory);
        string fileName = SpcExecutor.GetCheckpointFileName(temperature, replicaIndex);
        return Path.Combine(checkpointDirectory, fileName);
    }

    public static string GetSidecarPath(string checkpointPath)
        => Path.ChangeExtension(checkpointPath, ".spce");

    // ── Per-CSV path helpers — filenames owned by SpcCsvWriter ───────────
    public static string GetSweepCsvPath(string runDirectory)
        => Path.Combine(GetCsvDirectory(runDirectory), SpcCsvWriter.SweepFileName);

    public static string GetPartitionCsvPath(string runDirectory)
        => Path.Combine(GetCsvDirectory(runDirectory), SpcCsvWriter.PartitionFileName);

    public static string GetCriteriaCsvPath(string runDirectory)
        => Path.Combine(GetCsvDirectory(runDirectory), SpcCsvWriter.CriteriaFileName);

    public static string GetSessionCsvPath(string runDirectory)
        => Path.Combine(GetCsvDirectory(runDirectory), SpcCsvWriter.SessionFileName);

    public static string GetPartitionScheduleCsvPath(string runDirectory)
        => Path.Combine(GetCsvDirectory(runDirectory), SpcCsvWriter.PartitionScheduleFileName);

    public static string GetAnalysisCsvPath(string runDirectory)
        => Path.Combine(GetCsvDirectory(runDirectory), SpcCsvWriter.AnalysisFileName);

    public static string GetReplicaTracesCsvPath(string runDirectory)
        => Path.Combine(GetCsvDirectory(runDirectory), SpcCsvWriter.ReplicaTracesFileName);

    public static string GetEquilibriumEdgesCsvPath(string runDirectory)
        => Path.Combine(GetCsvDirectory(runDirectory), SpcCsvWriter.EquilibriumEdgesFileName);

    private static string SanitizePathSegment(string value)
    {
        if (string.IsNullOrWhiteSpace(value))
            return "spc";

        char[] invalid = Path.GetInvalidFileNameChars();
        var builder = new System.Text.StringBuilder(value.Length);
        bool lastWasSeparator = false;

        foreach (char c in value)
        {
            bool isInvalid = Array.IndexOf(invalid, c) >= 0;
            bool isSeparator = isInvalid || char.IsWhiteSpace(c) || c == '.';
            if (isSeparator)
            {
                if (!lastWasSeparator)
                {
                    builder.Append('-');
                    lastWasSeparator = true;
                }
                continue;
            }

            builder.Append(c);
            lastWasSeparator = false;
        }

        string sanitized = builder.ToString().Trim('-');
        return string.IsNullOrWhiteSpace(sanitized) ? "spc" : sanitized;
    }
}
