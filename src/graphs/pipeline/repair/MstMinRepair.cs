using System;
using Graphs.Proximity;

namespace Graphs.Pipeline.Repair;

/// <summary>
/// Stage 3 — MST-min connectivity repair. Adds only the minimal set of
/// MST bridge edges needed to unify the input's disconnected components
/// into one. Thin wrapper over
/// <see cref="ConnectivityRepair.EnsureConnected"/>, which uses the
/// pre-seeded Borůvka primitive
/// (<see cref="Graphs.Primitives.Mst.Boruvka.AddMinimalBridges"/>) under
/// the hood.
/// </summary>
/// <remarks>
/// <para>The "MST-min" name reflects the academic distinction with
/// MST-all (inject the entire spanning tree of the global graph).
/// MST-min preserves the local topology produced by Stage 2 and only
/// adds bridges where they're strictly required for global
/// reachability; MST-all over-saturates with long-distance edges that
/// distort downstream weight distributions.</para>
/// </remarks>
public sealed class MstMinRepair : ITopologyRepair
{
    public NeighborSelection Repair(NeighborSelection input, int n, Func<int, int, double> pairDistance)
        => ConnectivityRepair.EnsureConnected(input, n, pairDistance);
}
