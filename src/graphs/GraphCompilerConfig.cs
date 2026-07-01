using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Pipeline;

namespace Graphs;

// =========================================================================
// GraphCompilerConfig — the declarative DTO that drives GraphCompiler.Build
// =========================================================================
// One record per pipeline stage. Strategy enums carry Auto as default at
// every stage EXCEPT Topology (which is a mental-model choice the engine
// can't guess — Knn vs EpsilonBall). Build() evaluates the whole config
// holistically before executing, so order-of-method-calls and
// X-interrupts-Y questions go away (LMP renders global bandwidth dead,
// for example — the engine sees the whole board).
//
// This file declares the config types. Build() consumption lives in
// GraphCompiler.cs.
// =========================================================================

/// <summary>
/// Top-level configuration for <see cref="GraphCompiler.Build"/>.
/// Deeply immutable; a 1:1 snapshot ends up in the
/// <c>GraphConstructionManifest</c> for reproducibility.
/// </summary>
public sealed record GraphCompilerConfig
{
    /// <summary><b>Required.</b> No sensible engine default for the
    /// fundamental "what does local mean" choice.</summary>
    public required TopologyConfig Topology { get; init; }

    public FilterConfig             Filter      { get; init; } = new();
    public RepairConfig             Repair      { get; init; } = new();
    public RefinementConfig         Refinement  { get; init; } = new();

    /// <summary><b>Required.</b> L3 projection — declares the output weight
    /// kind the caller wants. <see cref="DistanceProjection"/> stops at the
    /// conditioned distance graph; <see cref="CouplingProjection"/> applies a
    /// kernel (SPC / Potts). The engine emits exactly what is requested and
    /// echoes <see cref="IEdgeProjection.Kind"/> onto
    /// <see cref="GraphBuildResult.WeightKind"/>.</summary>
    public required IEdgeProjection Projection { get; init; }

    public PathologyInterruptConfig Interrupts  { get; init; } = new();
}

// ── Stage 1: Topology ───────────────────────────────────────────────────

public sealed record TopologyConfig
{
    public required TopologyKind Kind { get; init; }

    /// <summary>Honored when <see cref="Kind"/> is
    /// <see cref="TopologyKind.Knn"/>. Defaults to 10 when omitted.</summary>
    public int? K { get; init; }

    /// <summary>Honored when <see cref="Kind"/> is
    /// <see cref="TopologyKind.EpsilonBall"/>. Required on that path.</summary>
    public double? Epsilon { get; init; }
}

public enum TopologyKind { Knn, EpsilonBall }

// ── Stage 2: Filter ─────────────────────────────────────────────────────

public sealed record FilterConfig
{
    public FilterKind Kind { get; init; } = FilterKind.Auto;

    /// <summary>Honored when <see cref="Kind"/> resolves to
    /// <see cref="FilterKind.MutualKnn"/> (either explicit or
    /// auto-picked). Null falls back to
    /// <see cref="MutualBandwidthSource.DirectedKth"/>.</summary>
    public MutualBandwidthSource? MutualBandwidthSource { get; init; }
}

public enum FilterKind
{
    /// <summary>Diagnostic-driven: hubness skewness &gt; 3.0 →
    /// MutualKnn; else OrRule. The auto-pick decision is recorded in
    /// the DiagnosticsLog.</summary>
    Auto,
    /// <summary>OR-rule: edge (i,j) exists if either side nominated
    /// the other. Denser; hub-prone in high D.</summary>
    OrRule,
    /// <summary>AND-rule: edge (i,j) only if both nominated. Suppresses
    /// hubs; shatters in high D — usually paired with MstMin repair.</summary>
    MutualKnn,
}

// ── Stage 3: Repair ─────────────────────────────────────────────────────

public sealed record RepairConfig
{
    public RepairKind Kind { get; init; } = RepairKind.Auto;
}

public enum RepairKind
{
    /// <summary>Diagnostic-driven: largest-component coverage &lt; 0.95 →
    /// MstMin; else NoRepair.</summary>
    Auto,
    NoRepair,
    /// <summary>Minimal Borůvka bridges across disconnected components.</summary>
    MstMin,
    // MstAll deferred — no current use case.
}

// ── Stage 4: Refinement ─────────────────────────────────────────────────

public sealed record RefinementConfig
{
    public RefinementKind Kind { get; init; } = RefinementKind.Auto;

    /// <summary>
    /// Optional upper-bound on geodesic path distances used by the
    /// <see cref="RefinementKind.PathNeighbor"/> refiner. When specified,
    /// refinement stops exploring paths longer than this distance.
    /// </summary>
    public double? MaxDistance { get; init; }
}

public enum RefinementKind
{
    /// <summary>Resolves to <see cref="PathNeighbor"/> — bounded SSSP with
    /// edge-local early exit via <c>PathNeighborRefiner</c>. Heavier than
    /// <see cref="Euclidean"/> pass-through; sweep configs should treat this
    /// as the engine default.</summary>
    Auto,
    /// <summary>Explicit pass-through: edge distances are not refined over
    /// the repaired topology.</summary>
    Euclidean,
    /// <summary>Geodesic path-neighbor refinement over the Stage 3 topology.</summary>
    PathNeighbor,
}

// ── Stage 5: Projection (L3) ─────────────────────────────────────────────
// The former ScalingConfig + BandwidthConfig are gone. Coupling — kernel,
// LMP, and bandwidth-strategy override — now lives on CouplingProjection
// (see GraphProjection.cs); DistanceProjection stops at the L2 distance graph.

// ── Pathology interrupts ────────────────────────────────────────────────

/// <summary>
/// Fatal-threshold settings. Any non-null field activates a fail-fast
/// check at the relevant diagnostic gate inside <see cref="GraphCompiler.Build"/>.
/// Crossing the threshold raises
/// <see cref="GraphPathologyException"/>.
/// </summary>
/// <remarks>
/// Engine-default heuristics still emit <c>Warning</c>-level messages
/// into the <c>DiagnosticsLog</c> when something looks suspicious —
/// these interrupt thresholds escalate specific signals from warning
/// to fatal.
/// </remarks>
public sealed record PathologyInterruptConfig
{
    /// <summary>Maximum allowed in-degree skewness on the directed
    /// KNN graph (Stage 1 output). Cross → Fatal.</summary>
    public double? MaxHubnessSkewness { get; init; }

    /// <summary>Minimum allowed coverage of the largest connected
    /// component after Stage 3. Below → Fatal.</summary>
    public double? MinLargestCoverage { get; init; }

    /// <summary>Minimum median edge weight on the final CsrGraph
    /// (Stage 5 output). Below → bandwidth collapsed → Fatal.</summary>
    public double? MinMedianWeight { get; init; }

    /// <summary>Maximum median edge weight on the final CsrGraph
    /// (Stage 5 output). Above → bandwidth too tight (saturation)
    /// → Fatal.</summary>
    public double? MaxMedianWeight { get; init; }

    /// <summary>Maximum allowed fraction of edges whose coupling weight
    /// rounds to near-zero (&lt; 1e-8). Caught via Stage 5's speculative
    /// pre-pass before the final weight array is allocated. Above →
    /// "delta-collapse" pathology → Fatal. Recommended setting for
    /// MST-repaired graphs: 0.05 (5%).</summary>
    public double? MaxNearZeroEdgeRatio { get; init; }
}
