using System;
using System.Collections.Generic;
using Graphs.Primitives;
using Graphs.Primitives.Mst;

namespace Clustering.Dendrograms;

/// <summary>
/// The thermal dendrogram producer: merge heights are <b>decoupling
/// temperatures</b> read off per-edge affinity curves G_e(T). An edge couples
/// at T_e = the hottest grid temperature with G_e(T) ≥ θ; single-linkage over
/// descending T_e (the negation trick orders the build; the raw temperatures
/// are re-stamped as heights) yields the thermal merge tree — hot singletons
/// at the leaves, the cold root last. <c>CostAxis = "temperature"</c>, so the
/// axis-alignment law holds against thermal landscapes directly. Heights
/// DESCEND in build order (<see cref="LandscapeWalk"/> detects the
/// orientation); <see cref="Dendrogram.CutAt"/> assumes ascending heights and
/// does not apply — cut thermally with <see cref="Dendrogram.CutToK"/>.
/// </summary>
/// <remarks>
/// Producer-agnostic over the currency: feed PKWang's exact closed-form G(T)
/// columns or SW's sampled co-membership/affinity columns. Sampled Ĝ(T) is
/// only noisily monotone — the BARS monotonizer slots in UPSTREAM (replace
/// the raw columns with monotonized posterior curves; same signature).
/// Heights are grid-quantized: T_e resolves to grid points — refine the grid
/// (or BARS-interpolate) for finer merge structure.
/// </remarks>
public static class ThermalDendrogram
{
    /// <summary>
    /// Builds the thermal merge FOREST from grid-major per-edge field columns
    /// (one column per ascending grid temperature, indexed by CSR slot — the
    /// same convention as <c>Landscape.ValuesByGridPoint</c>). Edges whose
    /// field never reaches θ within the grid do not couple — their endpoints
    /// are <b>thermal outliers in the observed window</b> (e.g. a sparse
    /// background whose ordering temperature lies below the grid floor), and
    /// the result is a forest: multiple roots, with never-coupled leaves
    /// belonging to no merge at all. The walk handles forests natively and
    /// such leaves resolve to <c>Assignment.Unassigned</c> — the honest
    /// abstain, not an error.
    /// </summary>
    public static Dendrogram FromEdgeCurves(
        CsrGraph graph,
        IReadOnlyList<double> temperatures,
        IReadOnlyList<double[]> edgeFieldByGridPoint,
        double theta = 0.5)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(temperatures);
        ArgumentNullException.ThrowIfNull(edgeFieldByGridPoint);
        if (graph.NodeCount < 2)
            throw new ArgumentOutOfRangeException(nameof(graph), "A thermal dendrogram needs at least 2 nodes.");
        if (temperatures.Count == 0)
            throw new ArgumentException("At least one grid temperature is required.", nameof(temperatures));
        if (edgeFieldByGridPoint.Count != temperatures.Count)
            throw new ArgumentException(
                $"One edge-field column per grid temperature: {edgeFieldByGridPoint.Count} columns vs {temperatures.Count} temperatures.",
                nameof(edgeFieldByGridPoint));
        for (int t = 1; t < temperatures.Count; t++)
            if (temperatures[t] <= temperatures[t - 1])
                throw new ArgumentException("Grid temperatures must be strictly ascending.", nameof(temperatures));
        int slots = graph.Targets.Length;
        for (int t = 0; t < edgeFieldByGridPoint.Count; t++)
            if (edgeFieldByGridPoint[t].Length != slots)
                throw new ArgumentException(
                    $"Edge-field column {t} length ({edgeFieldByGridPoint[t].Length}) does not match CSR slot count ({slots}).",
                    nameof(edgeFieldByGridPoint));

        int n = graph.NodeCount;
        var coupled = new List<MstEdge>();
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
        {
            // T_e = the hottest grid temperature where the bond still holds.
            for (int t = temperatures.Count - 1; t >= 0; t--)
            {
                if (edgeFieldByGridPoint[t][edge.Slot] >= theta)
                {
                    coupled.Add(new MstEdge(edge.Source, edge.Target, -temperatures[t])); // negate: hottest first
                    break;
                }
            }
        }

        MstEdge[] sorted = coupled.ToArray();
        Array.Sort(sorted); // ascending negated weight = descending decoupling temperature

        var mst = new MstEdge[n - 1];
        int treeEdges = Kruskal.BuildMinimumSpanningTree(sorted, n, mst);

        // Forest-tolerant single-linkage build (the shared DendrogramBuilder
        // requires a spanning tree; thermal forests are first-class here).
        // Same id convention: internal ids n .. n+treeEdges-1 in build order.
        var uf = new UnionFind(2 * n - 1);
        uf.Reset();
        var merges = new DendrogramNode[treeEdges];
        int nextId = n;
        for (int i = 0; i < treeEdges; i++)
        {
            MstEdge e = mst[i];
            int ra = uf.Find(e.U);
            int rb = uf.Find(e.V);
            int sizeA = uf.Size(ra);
            int sizeB = uf.Size(rb);
            uf.Union(ra, rb);
            uf.Reroot(ra, nextId, sizeA + sizeB);
            // The negation ordered the build; heights stay on the physical axis.
            merges[i] = new DendrogramNode(ra, rb, -e.Weight, sizeA + sizeB);
            nextId++;
        }

        return new Dendrogram(merges, n, CostAxis: "temperature");
    }
}
