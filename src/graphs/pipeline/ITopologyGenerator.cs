using System;

namespace Graphs.Pipeline;

/// <summary>
/// Stage 1 of the graph-construction pipeline: emit the raw, **directed**
/// (pre-symmetrization) candidate edges. The output is a
/// <see cref="NeighborSelection"/> whose <see cref="NeighborSelection.AllNeighbors"/>
/// rows are the per-node top-k candidates from the generator's own
/// perspective — no OR/AND symmetrization has been applied.
/// </summary>
/// <remarks>
/// <para><b>Why directed.</b> Splitting the directed pass from the
/// symmetrization step gives Stage 2 (<see cref="ITopologyFilter"/>)
/// something meaningful to do — the OR rule (PassThrough) versus the
/// AND rule (MutualKnn) becomes a per-filter choice rather than a
/// monolithic generator choice. This matches the
/// Dalmia &amp; Sia / UMAP-style decomposition where the prior on
/// connectivity is a swappable component.</para>
///
/// <para>For generators that have no meaningful "directed phase"
/// (e.g. EpsilonBall, which is intrinsically symmetric given a
/// symmetric distance), the implementation should populate the
/// symmetric view in <see cref="NeighborSelection.AllNeighbors"/> and
/// let the filter stage no-op via PassThrough.</para>
/// </remarks>
public interface ITopologyGenerator
{
    /// <summary>
    /// Generate the per-node candidate edges. <paramref name="dist"/> is
    /// the pairwise distance function; <paramref name="n"/> is the node
    /// count.
    /// </summary>
    NeighborSelection Generate(int n, Func<int, int, double> dist);
}
