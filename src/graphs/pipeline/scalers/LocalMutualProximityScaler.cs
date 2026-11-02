using System;
using System.Collections.Generic;
using Graphs;
using Graphs.Diagnostics;
using Graphs.Observables;
using Graphs.Primitives;
using Graphs.Proximity;

namespace Graphs.Pipeline.Scalers;

/// <summary>
/// Stage 5 — Local Mutual Proximity scaler. Computes coupling weights
/// from local neighborhood probability counts (SIMD popcount-friendly)
/// rather than a global MAD bandwidth. The compiler protects the repair
/// stage's MST bridges, which keep their kernel coupling; any other edge
/// that is a node's weakest gets an exact 0 under the strict rank rule.
/// </summary>
/// <remarks>
/// <para><b>Inner scaler.</b> LMP operates over an existing CsrGraph's
/// weight distribution — it needs *some* initial edge weights to count
/// against. The constructor takes an <c>innerScaler</c>
/// (default: Gaussian global with auto-estimated bandwidth, preserving
/// today's behavior). The inner scaler produces the initial weighted
/// graph; LMP re-weights it in-place using
/// <see cref="LocalMutualProximity.ApplyLocalScaling"/>.</para>
///
/// <para>The resulting <see cref="ScalerResult.SingleBandwidth"/> is
/// inherited from the inner scaler — it's still useful for provenance
/// even though LMP's final couplings are not derived from it.</para>
/// </remarks>
public sealed class LocalMutualProximityScaler : IEdgeScaler
{
    private readonly IEdgeScaler _innerScaler;
    private readonly bool        _weightsAreCouplings;
    private readonly Func<NeighborSelection, int, IReadOnlySet<(int Lo, int Hi)>?>? _protectedProvider;
    private readonly LmpRankRule _rankRule;

    public LocalMutualProximityScaler(
        IEdgeScaler? innerScaler         = null,
        bool         weightsAreCouplings = true,
        Func<NeighborSelection, int, IReadOnlySet<(int Lo, int Hi)>?>? protectedEdgeProvider = null,
        LmpRankRule  rankRule            = LmpRankRule.Strict)
    {
        _rankRule = rankRule;
        // Default inner: Gaussian global scaler with auto bandwidth. The
        // compiler composes LMP after the global scaler in
        // GraphCompiler.BuildScaler; this default just keeps the scaler
        // usable standalone.
        _innerScaler = innerScaler ?? new GlobalBandwidthScaler(kernel: KernelType.Gaussian);
        _weightsAreCouplings = weightsAreCouplings;
        _protectedProvider = protectedEdgeProvider;
    }

    public ScalerResult Scale(NeighborSelection refined, int n)
    {
        ScalerResult inner = _innerScaler.Scale(refined, n);
        IReadOnlySet<(int Lo, int Hi)>? protectedEdges = _protectedProvider?.Invoke(refined, n);
        CsrGraph rescaled = LocalMutualProximity.ApplyLocalScaling(
            inner.Graph,
            weightsAreCouplings: _weightsAreCouplings,
            protectedEdges: protectedEdges,
            rankRule: _rankRule);

        return inner with { Graph = rescaled, Lmp = new LmpOutcome(true, null) };
    }
}
