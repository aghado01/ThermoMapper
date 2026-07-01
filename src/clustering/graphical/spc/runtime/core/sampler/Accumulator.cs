using System;
using Graphs.Models.Potts;

namespace Clustering.Graphical.SPC.Runtime.Core.Sampler;

/// <summary>
/// The Swendsen–Wang sampler's <b>accumulator</b>: the running field of sufficient statistics it
/// folds each <c>Draw</c> into, plus the resume state (spins + RNG words). The model observables
/// (<see cref="Graphs.Models.Potts.Observables.Susceptibility"/> &amp; siblings) reduce the scalar
/// moments to χ/Cv/M/⟨E⟩; the optional per-edge counts reduce to the
/// <see cref="Graphs.Primitives.Affinities"/> / <see cref="Graphs.Primitives.Alignments"/>
/// currencies. A solver (PKWang) has no accumulator — no stream to fold.
/// </summary>
/// <remarks>
/// Replaces the old <c>PottsModelCoreObservables</c> + <c>PottsModelEdgeObservables</c> pair — the
/// <c>Core</c>/<c>Edge</c> split conflated <i>form-degree</i> (scalar moments vs the per-edge 1-form)
/// with <i>capability-gating</i> (always-on vs optional). Here the degree lives in the reductions
/// (<see cref="IReduction"/>) and the gate is just the nullable per-edge arrays.
/// </remarks>
public sealed record Accumulator
{
    // ── Identifying (validated on Restore) ───────────────────────────────
    public required double Temperature  { get; init; }
    public required int    Q            { get; init; }
    public          int    ReplicaIndex { get; init; }

    // ── Draw position (was CycleCount) ───────────────────────────────────
    public required int DrawCount { get; init; }

    // ── Configuration + resume state ─────────────────────────────────────
    public required int[]  Spins                { get; init; }
    public required int[]  ClusterSizeHistogram { get; init; }
    public required ulong  RngState0            { get; init; }
    public required ulong  RngState1            { get; init; }
    public required ulong  RngState2            { get; init; }
    public required ulong  RngState3            { get; init; }

    // ── Scalar sufficient-statistics (0-form → χ/Cv/M/⟨E⟩) ───────────────
    public required double RunningSumSqClusterSizes     { get; init; }
    public required double RunningSumSqClusterSizesExcl { get; init; }
    public required double RunningSumEnergy             { get; init; }
    public required double RunningSumEnergySq           { get; init; }
    public required double RunningSumMag                { get; init; }
    public required double RunningSumMagSq              { get; init; }

    // ── Per-edge sufficient-statistics (1-form → Affinities / Alignments / CoMembership; capability-gated, null when off) ──
    public int[]? BondFormedCount    { get; init; }
    public int[]? SpinAgreementCount { get; init; }
    public int[]? CoMembershipCount  { get; init; }

    // ── Per-node sufficient-statistics (un-reduced 0-form → MeanClusterSize / GiantParticipation
    //    landscapes via SwLandscapes; capability-gated, null when off) ──
    // The same per-node quantities the sweep collapses to global χ (Σ|c|²) and the M order
    // parameter, kept un-collapsed as per-node fields a downstream resolution step can ascend.
    public double[]? SumClusterSizePerNode    { get; init; }
    public double[]? SumInGiantClusterPerNode { get; init; }

    // ── Provenance ───────────────────────────────────────────────────────
    public DateTime TimestampUtc { get; init; } = DateTime.UtcNow;
}
