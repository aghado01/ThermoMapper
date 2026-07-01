using System;

namespace Graphs.Pipeline;

/// <summary>
/// Stage 3 of the graph-construction pipeline: restore global
/// reachability to a (potentially shattered) symmetric
/// <see cref="NeighborSelection"/> from Stage 2
/// (<see cref="ITopologyFilter"/>).
/// </summary>
/// <remarks>
/// <para><b>When repair runs.</b> The MutualKnn filter in high
/// dimensions reliably shatters the graph into disconnected components.
/// Without repair, spectral methods, BFS-based diagnostics, and
/// downstream clustering all degenerate. The repair stage decides how
/// to bridge the components (MST-min vs MST-all) — or whether to leave
/// the graph shattered for callers that genuinely want per-component
/// behavior.</para>
///
/// <para>Default for the auto-pipeline:
/// <c>MstMinRepair</c> when <see cref="Graphs.Diagnostics.Connectivity"/>
/// reports more than one component, else <c>NoRepair</c>.</para>
/// </remarks>
public interface ITopologyRepair
{
    /// <summary>
    /// Repair connectivity if needed. <paramref name="pairDistance"/> is
    /// required because the repair may need to evaluate cross-component
    /// distances that the generator never sampled.
    /// </summary>
    NeighborSelection Repair(NeighborSelection input, int n, Func<int, int, double> pairDistance);
}
