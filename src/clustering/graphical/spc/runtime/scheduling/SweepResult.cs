using System;
using System.Collections.Generic;
using Clustering.Graphical.SPC.Profiling;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Scheduling;

/// <summary>
/// Summary from one <see cref="ISweepStrategy.Run"/> call. Intended
/// for both runtime introspection (debug output, profiling) and analytic
/// post-processing.
/// </summary>
/// <remarks>
/// Sampler-agnostic and strategy-agnostic — shared by every
/// <see cref="ISweepStrategy"/> implementation. (Previously co-located with
/// the adaptive scheduler config; rehoused here when the adaptive strategy
/// was parked so the fixed-grid path and the partition/session layers no
/// longer depend on an adaptive-named file.)
/// </remarks>
public sealed class SweepSummary
{
    public required int      SubgraphNodes     { get; init; }
    public required int      SubgraphEdges     { get; init; }
    public          double   ChosenTemperature { get; init; }
    public          double   StabilityScore    { get; init; }
    public          int      TotalCyclesUsed   { get; init; }
    public          bool     EarlyStopped      { get; init; }
    public          TimeSpan Elapsed           { get; init; }

    /// <summary>
    /// Per-stage wall-clock breakdown. Keys are strategy-defined stage
    /// names (e.g. <c>"sweep"</c>, <c>"equilibrium"</c> for fixed-grid).
    /// The sum should approximate <see cref="Elapsed"/> minus negligible
    /// bookkeeping overhead. Empty for trivial graphs or strategies that
    /// don't track per-stage timing.
    /// </summary>
    public IReadOnlyDictionary<string, TimeSpan> StageTimings { get; init; }
        = new Dictionary<string, TimeSpan>();

    /// <summary>
    /// The full χ(T) profile assembled from the sweep. Sorted ascending in
    /// T. Useful for plotting the susceptibility curve.
    /// </summary>
    public SweepProfile Profile { get; init; } = SweepProfile.Empty;
}

/// <summary>
/// Result of one sweep-strategy run. Carries the probe traces that
/// assembled the χ(T) profile and the model-agnostic currencies minted
/// at the chosen temperature — the narrow waist between the sweep and any
/// downstream partition or evaluation step.
/// </summary>
/// <remarks>
/// Sampler-agnostic and strategy-agnostic — any <see cref="ISweepStrategy"/>
/// (SW, PKWang, future estimators) populates the same shape.
/// <see cref="ChosenAlignments"/> is <see langword="null"/> when the sweep
/// did not collect spin-agreement counts (i.e.
/// <see cref="Core.Sampler.AccumulationSpec.Alignments"/> was false).
/// </remarks>
public sealed class SweepResult
{
    public required SweepSummary Summary { get; init; }
    public required IReadOnlyList<SpcRunResult> SweepRuns { get; init; }
    public required ProfileCriteria ProfileCriteria { get; init; }
    public required CsrGraph Graph { get; init; }

    /// <summary>Bond-survival affinities minted at the chosen temperature — always present.</summary>
    public required Affinities ChosenAffinities { get; init; }

    /// <summary>Spin-agreement alignments minted at the chosen temperature — null unless
    /// <see cref="Core.Sampler.AccumulationSpec.Alignments"/> was set.</summary>
    public Alignments? ChosenAlignments { get; init; }

    /// <summary>Co-membership frequencies minted at the chosen temperature — null unless
    /// <see cref="Core.Sampler.AccumulationSpec.CoMembership"/> was set.</summary>
    public CoMembership? ChosenCoMembership { get; init; }

    /// <summary>Convenience accessor: temperature at which the currencies were minted.</summary>
    public double ChosenTemperature => Summary.ChosenTemperature;
}
