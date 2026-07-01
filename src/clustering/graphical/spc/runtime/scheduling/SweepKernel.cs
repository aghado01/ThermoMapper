using System;
using System.Collections.Generic;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Scheduling;

/// <summary>
/// Stateless per-point sampler helpers shared across sweep strategies.
/// Owns the <c>SwRunSpec</c> assembly + <c>SwRunner.Run</c> call
/// and the post-call <see cref="SpcRunResult"/> wrapping. Strategies
/// (adaptive, fixed-grid, future BARS) call these to keep the per-T
/// sampling contract identical regardless of how the schedule was
/// produced.
/// </summary>
public static class SweepKernel
{
    /// <summary>
    /// Run a single equilibrium probe at temperature <paramref name="T"/>
    /// collecting scalar moments only — the cheap measurement every
    /// χ(T)-style sweep uses.
    /// </summary>
    public static SpcRunResult RunProbePoint(
        CsrGraph  graph,
        double    T,
        int       q,
        RunBudget budget,
        int?      seed)
    {
        var runResult = SwRunner.Run(new SwRunSpec
        {
            Graph        = graph,
            Temperature  = T,
            Q            = q,
            Accumulation = AccumulationSpec.None,
            Seed         = seed,
            Budget       = budget,
            ReplicaIndex = 0,
        });

        return new SpcRunResult
        {
            Graph       = graph,
            Accumulator = runResult.Accumulator,
        };
    }

    /// <summary>
    /// Run the equilibrium pass at the chosen temperature with the
    /// caller-specified <paramref name="accumulation"/>. Used by
    /// <see cref="Hierarchical.BlattPartitionStrategy"/> for its per-phase
    /// equilibria; typically called with
    /// <see cref="AccumulationSpec.Currencies"/> so the cut policy has
    /// access to both currency precursors.
    /// </summary>
    public static SpcRunResult RunEquilibrium(
        CsrGraph         graph,
        double           T,
        int              q,
        RunBudget        budget,
        AccumulationSpec accumulation,
        int?             seed)
    {
        var runResult = SwRunner.Run(new SwRunSpec
        {
            Graph        = graph,
            Temperature  = T,
            Q            = q,
            Accumulation = accumulation,
            Seed         = seed,
            Budget       = budget,
            ReplicaIndex = 0,
        });

        return new SpcRunResult
        {
            Graph       = graph,
            Accumulator = runResult.Accumulator,
        };
    }

    /// <summary>
    /// Canonical empty-result shape for graphs with fewer than 2 nodes —
    /// every node is its own cluster downstream. Strategies short-circuit
    /// to this so their main pipelines never see <c>N&lt;2</c>.
    /// </summary>
    public static SweepResult BuildTrivialResult(
        CsrGraph graph,
        int      q,
        TimeSpan elapsed)
    {
        return new SweepResult
        {
            Summary = new SweepSummary
            {
                SubgraphNodes = graph.NodeCount,
                SubgraphEdges = 0,
                Elapsed       = elapsed,
                Profile       = SweepProfile.Empty,
            },
            SweepRuns = Array.Empty<SpcRunResult>(),
            ProfileCriteria = new ProfileCriteria(
                AnchorTemperature:   0.0,
                AnchorBand:          (0.0, 0.0),
                RefinedTemperature:  0.0,
                CorroborationScore:  0.0,
                Enrichments:         new Dictionary<string, double>()),
            Graph = graph,
            ChosenAffinities = new Affinities
            {
                Temperature = 0.0,
                G           = Array.Empty<double>(),
            },
            ChosenAlignments = null,
        };
    }
}
