namespace Graphs.Pipeline;

/// <summary>
/// What the compiler hands a scaler beyond the selection it scales: the repair edges
/// (kept past the floor), the coupling floor, and the neighbour sample the bandwidth reads.
/// </summary>
/// <param name="RepairEdges">Edges the repair stage added; null for none.</param>
/// <param name="CouplingFloor">A kernel coupling at or under this is dropped, except on a repair edge or the strongest across a split dropping it would open; 0 keeps every edge.</param>
/// <param name="ScaleSelection">The edges the mean-edge-distance bandwidth averages; null = the selection being scaled.</param>
/// <param name="DirectedSelection">The directed k-NN rows, read when <paramref name="ScaleRank"/> is set.</param>
/// <param name="ScaleRank">Read the bandwidth from each node's first r directed neighbours; null = as configured by the strategy.</param>
public sealed record ScalingContext(
    RepairEdgeSet? RepairEdges = null,
    double CouplingFloor = 1e-10,
    NeighborSelection? ScaleSelection = null,
    NeighborSelection? DirectedSelection = null,
    int? ScaleRank = null)
{
    /// <summary>The smallest coupling a repair edge is kept at when its kernel underflows below it.</summary>
    public const double RepairEdgeMinimum = 1e-300;
}
