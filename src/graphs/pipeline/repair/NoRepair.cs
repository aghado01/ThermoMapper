using System;

namespace Graphs.Pipeline.Repair;

/// <summary>
/// Stage 3 — null repair. Pass-through; the caller has explicitly opted
/// out of connectivity repair (e.g. they want per-component clustering
/// behavior or know the input is already connected).
/// </summary>
public sealed class NoRepair : ITopologyRepair
{
    public NeighborSelection Repair(NeighborSelection input, int n, Func<int, int, double> pairDistance)
        => input;
}
