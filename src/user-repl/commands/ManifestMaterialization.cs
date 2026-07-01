using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Graphical.SPC.Runtime.Execution;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;

namespace UserRepl.Commands;

/// <summary>
/// Materialization helpers that hydrate <see cref="RunManifest"/> DTO
/// pieces back into the live objects that the SPC/HDBSCAN code paths
/// consume — datasets, graph configs, etc. Used by
/// <see cref="ExtractCommand"/> (and any future reanalysis tooling)
/// to reconstruct a run from its on-disk manifest without re-passing
/// CLI flags.
/// </summary>
/// <remarks>
/// <para>Kept in a separate file from the DTO declarations so the
/// manifest record types stay declarative — the runtime dependencies
/// (synthetic generator catalog, kernel mixture structs,
/// <see cref="GraphCompilerConfig"/>) only show up at materialization
/// time.</para>
/// </remarks>
public static class ManifestMaterialization
{
    /// <summary>
    /// Reconstructs the dataset from its <see cref="DatasetSpec"/>.
    /// For synthetic datasets, re-invokes the same generator with the
    /// same parameters + seed; for CSV, re-reads the file.
    /// </summary>
    public static SpcUserDataset Materialize(this DatasetSpec spec)
    {
        ArgumentNullException.ThrowIfNull(spec);

        if (string.Equals(spec.Source, "csv", StringComparison.OrdinalIgnoreCase))
        {
            string csvPath = spec.CsvPath
                ?? throw new InvalidOperationException("DatasetSpec.Source='csv' requires CsvPath.");
            return SpcUserSession.FromCsv(
                csvPath,
                spec.LabelColumn,
                spec.HasHeader ?? true,
                ParseDelimiter(spec.Delimiter));
        }

        if (string.Equals(spec.Source, "synthetic", StringComparison.OrdinalIgnoreCase))
        {
            string generatorName = spec.GeneratorName
                ?? throw new InvalidOperationException("DatasetSpec.Source='synthetic' requires GeneratorName.");

            var parameters = new Dictionary<string, object?>(StringComparer.OrdinalIgnoreCase);
            if (spec.Seed.HasValue) parameters["seed"] = spec.Seed.Value;
            if (spec.GeneratorParameters is not null)
                foreach (var kvp in spec.GeneratorParameters) parameters[kvp.Key] = kvp.Value;
            return SpcUserSession.GenerateDataset(generatorName, parameters);
        }

        throw new InvalidOperationException(
            $"Unknown DatasetSpec.Source '{spec.Source}'. Expected 'synthetic' or 'csv'.");
    }

    /// <summary>
    /// Reconstructs a <see cref="GraphCompilerConfig"/> from a
    /// <see cref="GraphSpec"/>, including the distance metric and any
    /// kernel-mixture configuration.
    /// </summary>
    public static GraphCompilerConfig ToConfig(this GraphSpec spec)
    {
        ArgumentNullException.ThrowIfNull(spec);

        return new GraphCompilerConfig
        {
            Topology = new TopologyConfig
            {
                Kind = ResolveTopologyKind(spec),
                K = spec.K,
                Epsilon = spec.Epsilon,
            },
            Filter = new FilterConfig
            {
                Kind = ResolveFilterKind(spec),
            },
            Repair = new RepairConfig
            {
                Kind = spec.EnsureConnected ? RepairKind.MstMin : RepairKind.NoRepair,
            },
            // The projection (kernel + LMP + bandwidth, or distance pass-through) is
            // embedded in the manifest and round-trips polymorphically — no kernel or
            // projection switch to maintain. Older/partial manifests default to a
            // zero-bandwidth Gaussian coupling (bandwidth auto-estimated at build time).
            Projection = spec.Projection ?? new CouplingProjection { Kernel = new Gaussian(0.0) },
        };
    }

    private static TopologyKind ResolveTopologyKind(GraphSpec spec)
        => spec.TopologyKind?.ToLowerInvariant() switch
        {
            "epsilon" => TopologyKind.EpsilonBall,
            "knn" => TopologyKind.Knn,
            _ => TopologyKind.Knn,
        };

    private static FilterKind ResolveFilterKind(GraphSpec spec)
        => spec.FilterKind?.ToLowerInvariant() switch
        {
            "mutualknn" => FilterKind.MutualKnn,
            "or_rule" => FilterKind.OrRule,
            _ => FilterKind.OrRule,
        };

    private static char ParseDelimiter(string? spec)
    {
        if (string.IsNullOrEmpty(spec)) return ',';
        return spec switch
        {
            "tab" or "\\t" => '\t',
            var s when s.Length == 1 => s[0],
            _ => throw new InvalidOperationException($"Invalid delimiter spec '{spec}'."),
        };
    }
}
