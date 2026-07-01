using System;
using System.Collections.Generic;
using Clustering.Evaluation.External;
using Clustering.Evaluation.Internal;
using Clustering.Primitives;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Partitions;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using Clustering.Graphical.SPC.Partitions.Hierarchical;
using Clustering.Graphical.SPC.Partitions.Strategies;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Profiling.Signals;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC;

/// <summary>
/// Result of one end-to-end SPC run.
/// Carries the partition + sweep diagnostics + chosen-T currencies +
/// a dictionary of evaluator scores keyed by each evaluator's canonical <c>Name</c>.
/// </summary>
/// <remarks>
/// <para><b>Evaluator scores.</b> Both external evaluators (Purity,
/// NMI, ARI, ...) and SPC-aware evaluators (BondModularity,
/// BondCoverage, BondConductance, ...) fold into a single
/// <see cref="EvaluatorScores"/> dictionary keyed by each evaluator's
/// <c>Name</c> property.</para>
///
/// <para><b>Purity convenience.</b> The <see cref="Purity"/> accessor
/// reads the <c>"Purity"</c> key from <see cref="EvaluatorScores"/>
/// for the common case. Returns <see langword="null"/> when no
/// <c>Purity</c> evaluator was supplied.</para>
/// </remarks>
public sealed record SpcSessionResult(
    Assignment Partition,
    SweepProfile Profile,
    ProfileCriteria ProfileCriteria,
    SweepSummary ScheduleSummary,
    CsrGraph Graph,
    Affinities ChosenAffinities,
    Alignments? ChosenAlignments,
    CoMembership? ChosenCoMembership,
    IReadOnlyList<SpcRunResult> SweepRuns,
    IReadOnlyDictionary<string, double> EvaluatorScores,
    PartitionHierarchy? Hierarchy = null)
{
    /// <summary>
    /// Convenience accessor for the <c>"Purity"</c> entry in
    /// <see cref="EvaluatorScores"/>. Returns <see langword="null"/>
    /// when no Purity evaluator was supplied to the session.
    /// </summary>
    public double? Purity =>
        EvaluatorScores.TryGetValue("Purity", out double v) ? v : null;
}

public static class SpcClusteringSession
{
    /// <summary>
    /// End-to-end SPC run on a prebuilt graph.
    /// </summary>
    /// <param name="graph">CSR graph to cluster.</param>
    /// <param name="partitionStrategy">Partition cut to apply to the
    /// chosen-T currencies. Defaults to
    /// <see cref="ThresholdSpinAgreement"/> at θ=0.5.</param>
    /// <param name="analyzer">Signal analyzer that scores the assembled
    /// sweep profile. Defaults to <see cref="ChiPeakSignalAnalyzer"/>.</param>
    /// <param name="sweepStrategy">Sweep-strategy seam (required — the adaptive
    /// default has been parked). Pass a <see cref="FixedGridSweepStrategy"/>
    /// (or any <see cref="ISweepStrategy"/>); see AutoGridFixedSweep for an
    /// auto-bracketed grid.</param>
    /// <param name="hierarchicalStrategy">Optional hierarchical partition strategy.
    /// When supplied, the session produces a supplemental partition hierarchy
    /// alongside the flat cut.</param>
    /// <param name="externalEvaluators">Optional external cluster evaluators.
    /// Evaluators are run only when <paramref name="referenceLabels"/> are provided.</param>
    /// <param name="spcEvaluators">Optional SPC-aware evaluators that score the
    /// partition against the chosen-T bond-frequency field.</param>
    /// <param name="referenceLabels">Optional ground-truth labels for external
    /// evaluator scoring.</param>
    public static SpcSessionResult Run(
        CsrGraph graph,
        IPartitionStrategy? partitionStrategy = null,
        ISignalAnalyzer? analyzer = null,
        ISweepStrategy? sweepStrategy = null,
        IHierarchicalPartitionStrategy? hierarchicalStrategy = null,
        IEnumerable<IExternalClusterEvaluator>? externalEvaluators = null,
        IEnumerable<IGraphPartitionEvaluator>? spcEvaluators = null,
        int[]? referenceLabels = null)
    {
        ISweepStrategy strategy = sweepStrategy ?? throw new ArgumentNullException(
            nameof(sweepStrategy),
            "A sweep strategy is required — the adaptive default has been parked. Pass a " +
            "FixedGridSweepStrategy (see AutoGridFixedSweep for an auto-bracketed grid).");
        SweepResult sweep = strategy.Run(graph, analyzer);

        // A graph too small to sweep (the trivial-result has an empty Affinities)
        // partitions trivially: every node in one cluster.
        var partition = graph.NodeCount < 2
            ? Assignment.FromLabels(new int[graph.NodeCount])
            : (partitionStrategy ?? new ThresholdCoMembership { Theta = 0.5 })
                .Apply(sweep.Graph, sweep.ChosenAffinities, sweep.ChosenAlignments,
                       coMembership: sweep.ChosenCoMembership);

        PartitionHierarchy? hierarchy = graph.NodeCount < 2
            ? null
            : hierarchicalStrategy?.Apply(sweep, graph);

        var scores = new Dictionary<string, double>(StringComparer.Ordinal);
        scores["Coverage"] = partition.Coverage;

        if (externalEvaluators is not null && referenceLabels is not null)
        {
            foreach (var ev in externalEvaluators)
                scores[ev.Name] = ev.Evaluate(partition.Labels, referenceLabels);
        }

        // SPC-aware evaluators score the partition over the chosen-T Affinities field.
        if (spcEvaluators is not null && sweep.ChosenAffinities.G.Length > 0)
        {
            foreach (var ev in spcEvaluators)
                scores[ev.Name] = ev.Evaluate(sweep.Graph, sweep.ChosenAffinities.G, partition.Labels, partition.Count);
        }

        return new SpcSessionResult(
            Partition:           partition,
            Profile:             sweep.Summary.Profile,
            ProfileCriteria:     sweep.ProfileCriteria,
            ScheduleSummary:     sweep.Summary,
            Graph:               sweep.Graph,
            ChosenAffinities:    sweep.ChosenAffinities,
            ChosenAlignments:    sweep.ChosenAlignments,
            ChosenCoMembership:  sweep.ChosenCoMembership,
            SweepRuns:           sweep.SweepRuns,
            EvaluatorScores:     scores,
            Hierarchy:           hierarchy);
    }

    public static SpcSessionResult Run(
        double[][] features,
        IDistanceMetric? metric = null,
        GraphCompilerConfig? graphConfig = null,
        IPartitionStrategy? partitionStrategy = null,
        ISignalAnalyzer? analyzer = null,
        ISweepStrategy? sweepStrategy = null,
        IHierarchicalPartitionStrategy? hierarchicalStrategy = null,
        IEnumerable<IExternalClusterEvaluator>? externalEvaluators = null,
        IEnumerable<IGraphPartitionEvaluator>? spcEvaluators = null,
        int[]? referenceLabels = null)
    {
        GraphCompilerConfig effectiveConfig = graphConfig ?? new GraphCompilerConfig
        {
            Topology = new TopologyConfig
            {
                Kind = TopologyKind.Knn,
                K = 10,
            },
            Projection = new CouplingProjection
            {
                Kernel = new Gaussian(0.5),
                LmpRescale = false,
            },
        };

        CsrGraph graph = SpcGraphBuilder.BuildResult(features, effectiveConfig, metric).Graph;
        return Run(
            graph,
            partitionStrategy,
            analyzer,
            sweepStrategy,
            hierarchicalStrategy,
            externalEvaluators,
            spcEvaluators,
            referenceLabels);
    }
}
