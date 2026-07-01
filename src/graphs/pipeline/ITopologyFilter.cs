using System;

namespace Graphs.Pipeline;

/// <summary>
/// Stage 2 of the graph-construction pipeline: apply a topology prior
/// (OR vs AND) to the directed selection from Stage 1
/// (<see cref="ITopologyGenerator"/>) and emit a symmetric
/// <see cref="NeighborSelection"/>.
/// </summary>
/// <remarks>
/// <para><b>The OR/AND prior.</b> The filter encodes whether an edge
/// (i, j) exists when only one side selected the other (PassThrough,
/// the OR rule — denser graph, hub-prone in high dimensions) or only
/// when both sides selected each other (MutualKnn, the AND rule —
/// hub-suppressing but shattering-prone). The choice changes the
/// downstream graph's degree distribution and connectivity character
/// fundamentally; see the design doc for the academic context.</para>
///
/// <para>For generators that are intrinsically symmetric (EpsilonBall),
/// PassThrough is the natural choice — it preserves the selection
/// unchanged.</para>
/// </remarks>
public interface ITopologyFilter
{
    /// <summary>
    /// Apply the filter rule to the directed selection. <paramref name="n"/>
    /// is passed so filters can use index-derived hashing if needed;
    /// <paramref name="pairDistance"/> is supplied for filters that
    /// re-evaluate distances during symmetrization (most implementations
    /// will not need it).
    /// </summary>
    NeighborSelection Filter(NeighborSelection directed, int n, Func<int, int, double> pairDistance);
}
