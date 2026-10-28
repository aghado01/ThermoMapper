using System;
using System.Collections.Generic;
using System.Linq;
using Archivory.Tabular;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Export;

/// <summary>
/// Tabular projections of the currencies every clustering method shares — the
/// resolved <see cref="Assignment"/> and the per-edge <see cref="Affinities"/>,
/// <see cref="Alignments"/> and <see cref="CoMembership"/> — so each method's run
/// bundle writes them with one schema, named for the currency rather than the
/// method that minted it. Method-specific tables (SPC's sweep profile, the
/// cumulative field's cluster-count sweep) stay with their producers.
/// </summary>
public static class CurrencyProjections
{
    /// <summary>
    /// One row per node: <c>node_index</c>, the feature columns
    /// <c>x0…x{d−1}</c> (when <paramref name="features"/> is given),
    /// <c>true_label</c> (when <paramref name="trueLabels"/> is given), and the
    /// resolved <c>label</c> (<see cref="Assignment.Unassigned"/> for an abstain).
    /// </summary>
    public static TabularProjection CreatePartitionProjection(
        Assignment partition,
        double[][]? features = null,
        int[]? trueLabels = null,
        string tableName = "partition")
    {
        if (partition is null)
            throw new ArgumentNullException(nameof(partition));
        if (features is not null && features.Length != partition.Labels.Length)
            throw new ArgumentException("Feature rows must match partition node count.", nameof(features));
        if (trueLabels is not null && trueLabels.Length != partition.Labels.Length)
            throw new ArgumentException("True label length must match partition node count.", nameof(trueLabels));

        int dimension = features?.Length > 0 ? features[0].Length : 0;
        if (features is not null)
        {
            foreach (var featureRow in features)
            {
                if (featureRow.Length != dimension)
                    throw new ArgumentException("All feature rows must have the same length.", nameof(features));
            }
        }

        var featureCols = features is null
            ? Array.Empty<(string Name, Func<int, object?> Selector)>()
            : Enumerable.Range(0, dimension)
                .Select(d => ($"x{d}", (Func<int, object?>)(i => features[i][d])))
                .ToArray();

        return TabularProjectionFactory.CreateIndexed(tableName, partition.Labels.Length)
            .Column("node_index", i => i)
            .Columns(featureCols)
            .ColumnIf(trueLabels is not null, "true_label", i => trueLabels![i])
            .Column("label", i => partition.Labels[i])
            .Build();
    }

    /// <summary>
    /// One row per undirected edge of <paramref name="graph"/> — the per-edge
    /// currencies at one temperature as a flat table: <c>source</c>,
    /// <c>target</c>, <c>weight</c>, <c>affinity</c>, <c>alignment</c> and
    /// <c>co_membership</c> (NaN when that currency is absent).
    /// </summary>
    /// <remarks>
    /// <para><b>Upper-triangle only.</b> The CSR stores each undirected edge
    /// twice (once per endpoint); this projection emits only <c>i &lt; j</c>, the
    /// slot the currencies are meaningful at.</para>
    /// </remarks>
    public static TabularProjection CreateEdgeCurrenciesProjection(
        CsrGraph graph,
        Affinities affinities,
        Alignments? alignments = null,
        CoMembership? coMembership = null,
        string tableName = "edge_currencies")
    {
        ArgumentNullException.ThrowIfNull(affinities);
        if (affinities.G.Length != graph.Targets.Length)
            throw new ArgumentException(
                $"Affinities.G length ({affinities.G.Length}) does not match " +
                $"CSR slot count ({graph.Targets.Length}).", nameof(affinities));

        // Pre-walk the CSR to build the upper-triangle slot list.
        var rows = new List<(int Slot, int I, int J)>(graph.Targets.Length / 2);
        int n = graph.NodeCount;
        for (int i = 0; i < n; i++)
        {
            int rowEnd = graph.RowPointers[i + 1];
            for (int e = graph.RowPointers[i]; e < rowEnd; e++)
            {
                int j = graph.Targets[e];
                if (j <= i) continue;
                rows.Add((e, i, j));
            }
        }

        return TabularProjectionFactory.CreateIndexed(tableName, rows.Count)
            .Column("source",        k => rows[k].I)
            .Column("target",        k => rows[k].J)
            .Column("weight",        k => graph.Weights[rows[k].Slot])
            .Column("affinity",      k => affinities.G[rows[k].Slot])
            .Column("alignment",     k => alignments   is null ? double.NaN : alignments.G[rows[k].Slot])
            .Column("co_membership", k => coMembership is null ? double.NaN : coMembership.G[rows[k].Slot])
            .Build();
    }
}
