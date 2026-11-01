namespace Graphs.Pipeline;

/// <summary>
/// What the compiler hands a scaler beyond the selection it scales: the repair edges
/// (kept past the floor) and the coupling floor itself.
/// </summary>
/// <param name="RepairEdges">Edges the repair stage added; null for none.</param>
/// <param name="CouplingFloor">A kernel coupling at or under this is dropped, except on a repair edge; 0 keeps every edge.</param>
public sealed record ScalingContext(
    RepairEdgeSet? RepairEdges = null,
    double CouplingFloor = 1e-10)
{
    /// <summary>The smallest coupling a repair edge is kept at when its kernel underflows below it.</summary>
    public const double RepairEdgeMinimum = 1e-300;
}
