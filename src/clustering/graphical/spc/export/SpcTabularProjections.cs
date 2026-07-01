using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Primitives;
using Clustering.Graphical.SPC.Partitions;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Archivory.Tabular;
using Graphs.Primitives;
using Maths.Information;

namespace Clustering.Graphical.SPC.Export;

/// <summary>
/// Domain-specific projections that map SPC types to
/// <see cref="TabularProjection"/> instances. One factory per domain type;
/// column naming and shape are canonical so downstream tooling can rely on
/// stable column names across runs.
/// </summary>
/// <remarks>
/// Consumes the general-purpose
/// <see cref="TabularProjectionFactory.CreateIndexed"/> and
/// <see cref="TabularProjectionFactory.CreateScalar"/> builders for plumbing
/// — column-list assembly, optional-column handling, dictionary flattening
/// — so the methods here express only the SPC-specific decisions about
/// what to surface and what each column is called.
/// </remarks>
public static class SpcTabularProjections
{
    public static TabularProjection CreateSweepProfileProjection(
        SweepProfile profile,
        string tableName = "spc_sweep")
    {
        if (profile is null)
            throw new ArgumentNullException(nameof(profile));

        return TabularProjectionFactory.CreateIndexed(tableName, profile.Count)
            .Column("T",            i => profile.Temperatures[i])
            .Column("Chi",          i => profile.Susceptibility[i])
            .Column("Cv",           i => profile.SpecificHeat[i])
            .Column("LabelEntropy", i => profile.LabelEntropy[i])
            .ColumnIf(profile.BondEntropy != null, "BondEntropy", i => profile.BondEntropy![i])
            .ColumnsFromDictionary(profile.AdditionalChannels, (channel, i) => channel[i])
            .Build();
    }

    public static TabularProjection CreatePartitionProjection(
        Assignment partition,
        double[][]? features = null,
        int[]? trueLabels = null,
        string tableName = "spc_partition")
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
            .Column("spc_label", i => partition.Labels[i])
            .Build();
    }

    public static TabularProjection CreateDatasetProjection(
        double[][] features,
        int[] labels,
        string tableName = "spc_dataset")
    {
        if (features is null)
            throw new ArgumentNullException(nameof(features));
        if (labels is null)
            throw new ArgumentNullException(nameof(labels));
        if (features.Length != labels.Length)
            throw new ArgumentException("Feature rows must match label count.", nameof(labels));

        int dimension = features.Length > 0 ? features[0].Length : 0;
        for (int i = 0; i < features.Length; i++)
        {
            if (features[i] is null)
                throw new ArgumentException($"Feature row {i} cannot be null.", nameof(features));
            if (features[i].Length != dimension)
                throw new ArgumentException("All feature rows must have the same length.", nameof(features));
        }

        var featureCols = Enumerable.Range(0, dimension)
            .Select(d => ($"x{d}", (Func<int, object?>)(i => features[i][d])))
            .ToArray();

        return TabularProjectionFactory.CreateIndexed(tableName, features.Length)
            .Column("Index", i => i)
            .Columns(featureCols)
            .Column("label", i => labels[i])
            .Build();
    }

    public static TabularProjection CreateCriteriaProjection(
        ProfileCriteria criteria,
        string tableName = "spc_criteria")
    {
        if (criteria is null)
            throw new ArgumentNullException(nameof(criteria));

        // Criteria is about T_c estimation diagnostics — partition-quality
        // scores (Purity, modularity, etc.) live in spc_session.csv via
        // SpcSessionResult.EvaluatorScores.
        return TabularProjectionFactory.CreateScalar(tableName)
            .Column("AnchorTemperature",  criteria.AnchorTemperature)
            .Column("AnchorBandLo",       criteria.AnchorBand.Lo)
            .Column("AnchorBandHi",       criteria.AnchorBand.Hi)
            .Column("RefinedTemperature", criteria.RefinedTemperature)
            .Column("CorroborationScore", criteria.CorroborationScore)
            .ColumnsFromDictionary(criteria.Enrichments)
            .Build();
    }

    public static TabularProjection CreateSessionSummaryProjection(
        SpcSessionResult result,
        string tableName = "spc_session")
    {
        if (result is null)
            throw new ArgumentNullException(nameof(result));

        return TabularProjectionFactory.CreateScalar(tableName)
            .Column("AnchorTemperature",  result.ProfileCriteria.AnchorTemperature)
            .Column("AnchorBandLo",       result.ProfileCriteria.AnchorBand.Lo)
            .Column("AnchorBandHi",       result.ProfileCriteria.AnchorBand.Hi)
            .Column("RefinedTemperature", result.ProfileCriteria.RefinedTemperature)
            .Column("CorroborationScore", result.ProfileCriteria.CorroborationScore)
            .Column("ClusterCount",       result.Partition.Count)
            .Column("SubgraphNodes",      result.ScheduleSummary.SubgraphNodes)
            .Column("SubgraphEdges",      result.ScheduleSummary.SubgraphEdges)
            .Column("ChosenTemperature",  result.ScheduleSummary.ChosenTemperature)
            .Column("StabilityScore",     result.ScheduleSummary.StabilityScore)
            .Column("TotalCyclesUsed",    result.ScheduleSummary.TotalCyclesUsed)
            .Column("EarlyStopped",       result.ScheduleSummary.EarlyStopped)
            .Column("ElapsedSeconds",     result.ScheduleSummary.Elapsed.TotalSeconds)
            // Bond entropy at the scheduler's chosen T — a distinct readout from the
            // analyzer's plateau-anchored "BondEntropy" enrichment (at TClus /
            // RefinedTemperature), which ChiPeakSignalAnalyzer always emits into
            // Enrichments below. Namespaced to match the ChosenTemperature column and
            // the GetBondEntropyAtChosenT helper, so the two readouts coexist instead
            // of colliding on a shared "BondEntropy" name.
            .ColumnIf(result.Profile.BondEntropy != null, "BondEntropyAtChosenT", GetBondEntropyAtChosenT(result))
            .ColumnsFromDictionary(result.ProfileCriteria.Enrichments)
            .ColumnsFromDictionary(result.EvaluatorScores)
            .Build();
    }

    private static double GetBondEntropyAtChosenT(SpcSessionResult result)
    {
        var bondEntropy = result.Profile.BondEntropy;
        if (bondEntropy == null || bondEntropy.Count == 0) return double.NaN;

        var temperatures = result.Profile.Temperatures;
        if (temperatures == null || temperatures.Count == 0) return double.NaN;

        // ChosenTemperature anchors the readout; fall back to the coldest grid point
        // when no schedule summary is present (e.g. a profile-only export path).
        double chosenT = result.ScheduleSummary?.ChosenTemperature ?? temperatures[0];

        // Find the nearest grid index, clamped to the bond-entropy series so an
        // off-by-one between the temperature grid and the entropy curve can't index out of range.
        int limit = Math.Min(temperatures.Count, bondEntropy.Count);
        int bestIdx = 0;
        double bestDist = double.MaxValue;
        for (int i = 0; i < limit; i++)
        {
            double dist = Math.Abs(temperatures[i] - chosenT);
            if (dist < bestDist)
            {
                bestDist = dist;
                bestIdx = i;
            }
        }
        return bondEntropy[bestIdx];
    }

    public static TabularProjection CreatePartitionScheduleRollupProjection(
        IReadOnlyList<SchedulePartitionRollup> rollups,
        string tableName = "spc_partition_schedule")
    {
        if (rollups is null)
            throw new ArgumentNullException(nameof(rollups));

        var levelNames = rollups.Count > 0
            ? rollups[0].LevelNames
            : Array.Empty<string>();

        // Validate uniform level structure across rollups before building.
        foreach (var rollup in rollups)
        {
            if (rollup.LevelNames.Count != levelNames.Count)
                throw new ArgumentException("All rollups must have the same level names.", nameof(rollups));
            if (rollup.Purities.Count != levelNames.Count)
                throw new ArgumentException("Rollup purity vector length must match level names.", nameof(rollups));
        }

        var purityCols = Enumerable.Range(0, levelNames.Count)
            .Select(level => ($"purity_{levelNames[level]}", (Func<int, object?>)(i => rollups[i].Purities[level])))
            .ToArray();

        return TabularProjectionFactory.CreateIndexed(tableName, rollups.Count)
            .Column("T",             i => rollups[i].Temperature)
            .Column("replicas",      i => rollups[i].ReplicaCount)
            .Column("clusters",      i => rollups[i].ClusterCount)
            .Column("pooled_cycles", i => rollups[i].PooledCycleCount)
            .Columns(purityCols)
            .Build();
    }

    /// <summary>
    /// One row per <see cref="SpcRunResult"/> frame, preserving the
    /// (T, replica, round) trace before averaging. Designed for
    /// variance-band plotting in Python — the standard
    /// <c>spc_sweep.csv</c> only carries the averaged χ(T) / Cv(T) /
    /// label-entropy curves, so cross-replica variance gets lost. This
    /// projection keeps the raw per-frame signals so the spread can be
    /// reconstructed downstream.
    /// </summary>
    /// <remarks>
    /// <para><b>Columns.</b> <c>temperature</c>, <c>replica_index</c>,
    /// <c>cycle_count</c>, <c>chi_fk</c>, <c>cv</c>, <c>label_entropy</c>,
    /// <c>mean_energy</c>, <c>mean_magnetization</c>,
    /// <c>magnetization_variance</c>. Signals mirror the formulas used
    /// by the <c>Susceptibility</c> / <c>SpecificHeat</c> / <c>LabelEntropy</c>
    /// / <c>Magnetization</c> reductions, computed per-frame rather
    /// than aggregated across the (T, replica) bucket.</para>
    /// </remarks>
    public static TabularProjection CreateReplicaTracesProjection(
        IReadOnlyList<SpcRunResult> runs,
        string tableName = "spc_replica_traces")
    {
        if (runs is null)
            throw new ArgumentNullException(nameof(runs));

        return TabularProjectionFactory.CreateIndexed(tableName, runs.Count)
            .Column("temperature",            i => runs[i].Accumulator.Temperature)
            .Column("replica_index",          i => runs[i].Accumulator.ReplicaIndex)
            .Column("cycle_count",            i => runs[i].Accumulator.DrawCount)
            .Column("chi_fk",                 i => ChiFK(runs[i].Accumulator))
            .Column("cv",                     i => SpecificHeat(runs[i].Accumulator))
            .Column("label_entropy",          i => LabelEntropy(runs[i].Accumulator))
            .Column("mean_energy",            i => MeanEnergy(runs[i].Accumulator))
            .Column("mean_magnetization",     i => MeanMag(runs[i].Accumulator))
            .Column("magnetization_variance", i => MagVariance(runs[i].Accumulator))
            .Build();
    }

    // ── Per-frame signal helpers (delegate to the model reductions) ─────
    // Each computes the same scalar the corresponding model observable in
    // graphs/models/potts/observables/ produces, applied per-frame rather
    // than aggregated across the (T, replica) bucket. The formulas live in
    // one place (the model Reduce); these are thin per-row adapters.

    private static double ChiFK(Accumulator f)
        => Graphs.Models.Potts.Observables.Susceptibility.Reduce(
            f.RunningSumSqClusterSizes, f.DrawCount, f.Spins.Length);

    private static double SpecificHeat(Accumulator f)
        => Graphs.Models.Potts.Observables.SpecificHeat.Reduce(
            f.RunningSumEnergy, f.RunningSumEnergySq, f.DrawCount, f.Temperature);

    private static double LabelEntropy(Accumulator f)
        => Shannon.EntropyNats(f.ClusterSizeHistogram);

    private static double MeanEnergy(Accumulator f)
        => Graphs.Models.Potts.Observables.MeanEnergy.Reduce(f.RunningSumEnergy, f.DrawCount);

    private static double MeanMag(Accumulator f)
        => Graphs.Models.Potts.Observables.Magnetization.Reduce(
            f.RunningSumMag, f.RunningSumMagSq, f.DrawCount).Mean;

    private static double MagVariance(Accumulator f)
        => Graphs.Models.Potts.Observables.Magnetization.Reduce(
            f.RunningSumMag, f.RunningSumMagSq, f.DrawCount).Variance;

    /// <summary>
    /// One row per undirected edge of the input graph, carrying the
    /// CSR weight + the equilibrium bond-frequency and spin-agreement
    /// rates the chosen-T final pass measured. The natural "edge view"
    /// for Python-side inspection — same data the threshold partition
    /// strategies see, surfaced as a flat CSV instead of a
    /// per-checkpoint <c>.spce</c> binary.
    /// </summary>
    /// <remarks>
    /// <para><b>Upper-triangle only.</b> The CSR stores each undirected
    /// edge twice (once per endpoint); this projection emits only
    /// <c>i &lt; j</c> to match the convention threshold partition
    /// strategies use.</para>
    ///
    /// <para><b>Columns.</b> <c>source</c>, <c>target</c>, <c>weight</c>,
    /// <c>bond_frequency</c> (<c>BondFormedCount/CycleCount</c>),
    /// <c>spin_agreement</c> (<c>SpinAgreementCount/CycleCount</c>).</para>
    /// </remarks>
    public static TabularProjection CreateEquilibriumEdgesProjection(
        CsrGraph graph,
        Affinities affinities,
        Alignments? alignments,
        CoMembership? coMembership = null,
        string tableName = "spc_equilibrium_edges")
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
            .Column("source",         k => rows[k].I)
            .Column("target",         k => rows[k].J)
            .Column("weight",         k => graph.Weights[rows[k].Slot])
            .Column("bond_frequency", k => affinities.G[rows[k].Slot])
            .Column("spin_agreement", k => alignments    is null ? double.NaN : alignments.G[rows[k].Slot])
            .Column("co_membership",  k => coMembership  is null ? double.NaN : coMembership.G[rows[k].Slot])
            .Build();
    }
}
