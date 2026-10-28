using System;
using Archivory.Tabular;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Export;

/// <summary>
/// Writes the shared currency tables of <see cref="CurrencyProjections"/> as CSV,
/// atomically, returning the written path. File names belong to each method's
/// run bundle (SPC's <c>SpcCsvWriter</c>, the cumulative field's
/// <c>CumulativeFieldCsvWriter</c>); this writer owns only the schema.
/// </summary>
public static class CurrencyCsvWriter
{
    /// <summary>Write the resolved partition (see <see cref="CurrencyProjections.CreatePartitionProjection"/>).</summary>
    public static string WritePartition(Assignment partition, string path, double[][]? features = null, int[]? trueLabels = null)
    {
        ArgumentNullException.ThrowIfNull(partition);
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path cannot be null or whitespace.", nameof(path));

        CurrencyProjections.CreatePartitionProjection(partition, features, trueLabels).WriteToFile(path);
        return path;
    }

    /// <summary>
    /// Write the per-edge currencies at one temperature (see
    /// <see cref="CurrencyProjections.CreateEdgeCurrenciesProjection"/>). Returns
    /// <see langword="null"/> without writing when <paramref name="affinities"/>
    /// is absent or empty (e.g. a run that did not request edge currencies).
    /// </summary>
    public static string? WriteEdgeCurrencies(
        CsrGraph graph,
        Affinities? affinities,
        Alignments? alignments,
        string path,
        CoMembership? coMembership = null)
    {
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path cannot be null or whitespace.", nameof(path));

        if (affinities is null || affinities.G.Length == 0) return null;

        CurrencyProjections.CreateEdgeCurrenciesProjection(graph, affinities, alignments, coMembership).WriteToFile(path);
        return path;
    }
}
