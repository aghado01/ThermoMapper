using System.Collections.Generic;
using Graphs.Coupling;
using Graphs.Primitives;

namespace Graphs.Pipeline;

/// <summary>
/// Stage 5 of the graph-construction pipeline: convert refined
/// per-edge distances into kernel coupling weights, producing the
/// final <see cref="CsrGraph"/>.
/// </summary>
/// <remarks>
/// <para>Two implementation families:</para>
/// <list type="bullet">
///   <item><c>GlobalBandwidthScaler</c> — Gaussian/Cauchy/Laplacian/
///     Linear (or mixture) kernels with a MAD-estimated global
///     bandwidth drawn from the selection's k-th neighbor distance
///     sample (the kNN search radius).</item>
///   <item><c>LocalMutualProximityScaler</c> — bypasses global
///     bandwidth entirely; weights are derived from local
///     neighborhood probability counts (SIMD popcount-friendly).
///     Right when Stage 3 injected long bridges that a global kernel
///     would collapse to near-zero.</item>
/// </list>
/// </remarks>
public interface IEdgeScaler
{
    /// <summary>
    /// Convert distances to coupling weights. Returns the assembled
    /// graph plus the resolved bandwidth metadata so the build result
    /// can capture provenance without re-deriving it.
    /// </summary>
    ScalerResult Scale(NeighborSelection refined, int n);
}

/// <summary>
/// Output of an <see cref="IEdgeScaler.Scale"/> call. Exactly one of
/// <see cref="SingleBandwidth"/> and <see cref="MixtureBandwidth"/> is
/// non-null — single-kernel scalers populate the former, mixture
/// scalers populate the latter. Local Mutual Proximity scalers
/// populate neither (no global bandwidth concept). <see cref="CutEdges"/> counts the
/// undirected edges the kernel produced and the coupling floor dropped;
/// <see cref="GaugeFactor"/> is the global factor applied last (1 for none); <see cref="Lmp"/> is set by the LMP scalers;
/// <see cref="RestoredEdges"/> are the edges under the floor kept across each split it would have opened (null for none).
/// </summary>
public readonly record struct ScalerResult(
    CsrGraph          Graph,
    double?           SingleBandwidth,
    MixtureBandwidth? MixtureBandwidth,
    int               CutEdges = 0,
    double            GaugeFactor = 1.0,
    LmpOutcome?       Lmp = null,
    IReadOnlySet<(int Lo, int Hi)>? RestoredEdges = null);

/// <summary>Whether Local Mutual Proximity ran, and for the auto trigger the near-zero ratio it read.</summary>
public readonly record struct LmpOutcome(bool Applied, double? NearZeroRatio);
