using System;
using Clustering.Graphical.SPC;
using Clustering.Graphical.SPC.Export;

namespace UserRepl;

public sealed record SpcUserRunResult(
    SpcSessionResult SessionResult,
    SpcUserDataset Dataset,
    string RunDirectory,
    string SweepCsvPath,
    string PartitionCsvPath,
    string CriteriaCsvPath,
    string SessionCsvPath,
    string SummaryJsonPath,
    string ReplicaTracesCsvPath,
    string? EquilibriumEdgesCsvPath = null,
    string? PartitionHierarchyJsonPath = null)
{
    public string WriteTabularExports(string outputDirectory)
    {
        if (string.IsNullOrWhiteSpace(outputDirectory))
            throw new ArgumentException("Output directory must be provided.", nameof(outputDirectory));

        return SpcCsvWriter.WriteAllToDirectory(
            SessionResult,
            outputDirectory,
            features: Dataset.Features,
            trueLabels: Dataset.Labels);
    }
}
