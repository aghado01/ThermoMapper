namespace Graphs.Pipeline;

/// <summary>
/// Stage 4 of the graph-construction pipeline: re-evaluate distances
/// over the (possibly repaired) topology before they get fed into a
/// kernel. <see cref="Graphs.Pipeline.Refinement.PathNeighborRefiner"/>
/// recomputes geodesic distances via parallel bounded SSSP, which matters
/// when Stage 3 injected long MST bridges (raw ambient distance across a
/// manifold void is a poor proxy for local proximity).
/// </summary>
/// <remarks>
/// <para><see cref="RefinementKind.Auto"/> in <see cref="GraphCompilerConfig"/>
/// resolves to <see cref="Graphs.Pipeline.Refinement.PathNeighborRefiner"/>
/// and the explicit Euclidean pass-through refiner is retired.</para>
/// </remarks>
public interface IMetricRefiner
{
    /// <summary>
    /// Re-weight the edges in <paramref name="input"/> according to the
    /// refiner's distance interpretation, returning a new
    /// <see cref="NeighborSelection"/> with updated per-edge distances.
    /// The topology (edge set) is preserved.
    /// </summary>
    NeighborSelection Refine(NeighborSelection input, int n);
}
