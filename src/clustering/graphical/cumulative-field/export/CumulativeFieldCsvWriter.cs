using System;
using System.Collections.Generic;
using Archivory.Tabular;

namespace Clustering.Graphical.CumulativeField.Export;

/// <summary>
/// The cumulative-field run bundle's CSV writer: owns the bundle's CSV file
/// names and writes the method-native sweep table. The shared currency tables
/// (partition, edge currencies) are written with <c>Clustering.Export</c>'s
/// <c>CurrencyCsvWriter</c> under the names declared here.
/// </summary>
public static class CumulativeFieldCsvWriter
{
    // ── Canonical cumulative-field CSV filenames (single source of truth) ─
    public const string SweepFileName          = "sweep.csv";
    public const string PartitionFileName      = "partition.csv";
    public const string EdgeCurrenciesFileName = "edge_currencies.csv";

    /// <summary>
    /// Write the per-temperature cluster-count sweep, one row per temperature in
    /// caller order: <c>temperature</c>, <c>cluster_count</c>. The draw-free
    /// diagnostic the closed form offers in place of SPC's sweep profile.
    /// </summary>
    public static string WriteSweep(IReadOnlyList<double> temperatures, IReadOnlyList<int> clusterCounts, string path)
    {
        ArgumentNullException.ThrowIfNull(temperatures);
        ArgumentNullException.ThrowIfNull(clusterCounts);
        if (temperatures.Count != clusterCounts.Count)
            throw new ArgumentException(
                $"One cluster count per temperature: {clusterCounts.Count} counts vs {temperatures.Count} temperatures.",
                nameof(clusterCounts));
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path cannot be null or whitespace.", nameof(path));

        TabularProjectionFactory.CreateIndexed("sweep", temperatures.Count)
            .Column("temperature",   i => temperatures[i])
            .Column("cluster_count", i => clusterCounts[i])
            .Build()
            .WriteToFile(path);
        return path;
    }
}
