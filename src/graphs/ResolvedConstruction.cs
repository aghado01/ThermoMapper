using Graphs.Distance;

namespace Graphs;

/// <summary>
/// What a build actually did, beside the config that requested it: what each <c>Auto</c> resolved to,
/// the statistic behind it, and the units and connectivity of the graph returned. The units are the
/// construction's own declaration, so a consumer reads them here and never from a compiler option.
/// Statistics are null when not finite.
/// </summary>
public sealed record ResolvedConstruction(
    FilterKind Filter,
    double? InDegreeSkewness,
    RepairKind Repair,
    double? CoverageBeforeRepair,
    int RepairEdges,
    RefinementKind Refinement,
    int Components,
    double? LargestCoverage,
    CouplingResolution? Coupling);

/// <summary>The coupling stage's resolution. Null on a distance build.</summary>
/// <param name="Kernel">The kernel descriptor's name (Gaussian, Cauchy, Mixture, ...).</param>
/// <param name="Bandwidth">The bandwidth used, in distance units; the Gaussian component for a mixture.</param>
/// <param name="Lmp">Whether LMP ran (an explicit setting, or what the auto trigger decided).</param>
/// <param name="LmpNearZeroRatio">The statistic the auto trigger read; null when LMP was set explicitly.</param>
/// <param name="CutEdges">Edges the coupling floor dropped.</param>
/// <param name="RestoredEdges">Edges under the floor kept all the same, the strongest across each split dropping them would have opened.</param>
/// <param name="GaugeFactor">The factor applied last to every coupling; 1 for no gauge.</param>
/// <param name="CoupledComponents">Components over the edges of positive coupling.</param>
/// <param name="RepairEdgesDropped">Repair edges absent from the returned graph.</param>
/// <param name="RepairEdgesZeroed">Repair edges present with zero coupling.</param>
/// <param name="HasBwdUnits">The couplings are <c>J = K-hat⁻¹·exp(−d²/2a²)</c> with <c>a</c> the mean pair
/// distance and no LMP, the kernel for which BWD's q-only T_ps(q) is derived.</param>
public sealed record CouplingResolution(
    string Kernel,
    BandwidthStrategy BandwidthStrategy,
    double? Bandwidth,
    CouplingFidelity Fidelity,
    SphericalIntrinsicMode SphericalMode,
    bool Lmp,
    double? LmpNearZeroRatio,
    LmpRankRule LmpRank,
    double CouplingFloor,
    int CutEdges,
    int RestoredEdges,
    GaugeKind Gauge,
    double GaugeFactor,
    ScaleEdgeSet ScaleEdges,
    int? ScaleRank,
    int CoupledComponents,
    int ZeroCouplingEdges,
    int RepairEdgesDropped,
    int RepairEdgesZeroed,
    bool HasBwdUnits);
