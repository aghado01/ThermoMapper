#nullable enable
using System.Collections.Generic;
using Graphs.Primitives;
using Maths.Topology;

namespace TDA.Ph;

/// <summary>
/// H1 loop edges from involuted persistence on a distance-weighted skeleton —
/// the protect-set callers inject into GraphCompiler.Build (as its ProtectedEdgeSource) so LMP demotion spares load-bearing loops.
/// </summary>
public static class H1CycleEdges
{
    /// <summary>
    /// Union of all H1 representative cycle edges on <paramref name="distanceGraph"/>
    /// (normalized <c>(Lo, Hi)</c> vertex pairs).
    /// </summary>
    public static HashSet<(int Lo, int Hi)> FromDistanceGraph(CsrGraph distanceGraph)
    {
        if (distanceGraph.NodeCount == 0)
            return new HashSet<(int, int)>();

        SimplicialFiltration filtration = RipsFiltration.GraphRips(
            distanceGraph,
            FiltrationWeights.RawDistance,
            maxDimension: 2,
            label: "load-bearing");

        Barcode barcode = PersistentInvolutedHomology.Compute(filtration, representatives: true);

        var set = new HashSet<(int, int)>();
        foreach (var (_, edges) in BarCycleEdges.H1Loops(barcode, filtration))
        {
            foreach (TDA.Ph.UndirectedEdge e in edges)
                set.Add((e.Lo, e.Hi));
        }

        return set;
    }
}
