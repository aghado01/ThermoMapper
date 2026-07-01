using Graphs.Diagnostics;
using Graphs.Observables;

namespace TDA.Mapper.Diagnostics;

public readonly record struct NerveTopologyReport(
    int NerveNodeCount,
    int NerveEdgeCount,
    int ConnectedComponents,
    int LargestComponentSize,
    int LoopCount,
    bool IsTreeLike,
    int TriangleCount,
    int Girth);

public static class NerveTopology
{
    public static NerveTopologyReport From(MapperResult result)
    {
        var connectivity = Connectivity.Validate(result.Nerve);
        var cycles = Cycles.Compute(result.Nerve, connectivity);

        return new NerveTopologyReport(
            NerveNodeCount: connectivity.NodeCount,
            NerveEdgeCount: connectivity.EdgeCount,
            ConnectedComponents: connectivity.ComponentCount,
            LargestComponentSize: connectivity.LargestComponent,
            LoopCount: cycles.CyclomaticComplexity,
            IsTreeLike: cycles.CyclomaticComplexity == 0,
            TriangleCount: cycles.TriangleCount,
            Girth: cycles.Girth);
    }
}
