// ============================================================================
// TDA.Mapper — Mapper.cs
// ============================================================================
// Core MAPPER (Singh, Mémoli, Carlsson 2007) interfaces, DTOs, and orchestrator.
//
// Pipeline (Mapper.Build):
//   1. Apply filter (lens function) to every data point (or graph node)
//   2. Cover the filter range with overlapping bins (intervals or product cubes)
//   3. For each cover bin, locally cluster the preimage
//   4. Build nerve graph: nodes are (bin, local-cluster); edges connect nodes
//      whose member sets overlap (weighted by overlap count)
//
// The nerve is materialized as a CsrGraph from the Graphs.Primitives layer
// so it plugs into the same diagnostic / viz infrastructure as proximity graphs.
//
// Three parallel input paths:
//   IFilter / ICover            — 1-D point-cloud lens (the common case)
//   IMultiFilter / IMultiCover  — multi-D point-cloud lens (product covers)
//   IGraphFilter / IGraphClust. — graph-input MAPPER (Carrière-Michel-Oudot 2018,
//                                 Hajij-Rosen-Wang 2018); clusters preimages via
//                                 connected components of the induced subgraph
//
// All three Build overloads share the inverted-index nerve construction via an
// internal helper. The per-bin clustering step is parameterized via a delegate
// so the graph and point-cloud paths share the same outer scaffold.
// ============================================================================

#nullable enable
using System;
using System.Collections.Generic;
using Graphs.Primitives;
using TDA.Mapper.Cover;

namespace TDA.Mapper;

// ── Filter interfaces ───────────────────────────────────────────────────────

/// <summary>
/// 1-D lens function on point-cloud data: maps each data point to a scalar
/// filter value. For multi-D point-cloud lenses see <see cref="IMultiFilter"/>;
/// for graph-input lenses see <see cref="IGraphFilter"/>.
/// </summary>
public interface IFilter
{
    double[] Apply(double[][] data);
    string Name { get; }
}

/// <summary>
/// Multi-D lens function on point-cloud data: maps each data point to a
/// Dimension-length filter vector. Used for product-cover MAPPER.
/// </summary>
public interface IMultiFilter
{
    int Dimension { get; }
    double[][] Apply(double[][] data);
    string Name { get; }
}

/// <summary>
/// Graph-input lens function: maps each graph node to a scalar filter value.
/// May use graph structure (degree, Fiedler vector, geodesic distance from
/// a seed), node features, or both. Features are optional — purely graph-
/// topological filters pass null.
/// </summary>
public interface IGraphFilter
{
    double[] Apply(CsrGraph graph, double[][]? features = null);
    string Name { get; }
}

// ── Clusterer interfaces ────────────────────────────────────────────────────

/// <summary>
/// Point-cloud clusterer applied to each cover preimage.
/// </summary>
public interface IClusterer
{
    ClusterResult Cluster(double[][] subset);
    string Name { get; }
}

/// <summary>
/// Graph-input clusterer: clusters a preimage by the connected components
/// (or other graph-topological criterion) of the subgraph induced on the
/// preimage. This is the defining difference of graph MAPPER — two nodes
/// are clustered together iff they are connected through other preimage
/// nodes, not merely close in some metric.
/// </summary>
public interface IGraphClusterer
{
    ClusterResult ClusterInduced(CsrGraph graph, IReadOnlyList<int> preimageIndices);
    string Name { get; }
}

public sealed record ClusterResult(int[] Labels, int K);

// ── DTOs ────────────────────────────────────────────────────────────────────

public readonly record struct MapperNode(
    int BinId,
    int LocalClusterId,
    int[] MemberIndices,
    double FilterValueMean,
    double FilterValueMin,
    double FilterValueMax);

public sealed class MapperResult
{
    public required IReadOnlyList<MapperNode> Nodes { get; init; }
    public required CsrGraph Nerve { get; init; }
    public required string FilterName { get; init; }
    public required string CoverName { get; init; }
    public required string ClustererName { get; init; }
    public required int EmptyBinCount { get; init; }
}

// ── Orchestrator ────────────────────────────────────────────────────────────

public static class Mapper
{
    // ── 1-D point-cloud MAPPER ──────────────────────────────────────────────

    public static MapperResult Build(
        double[][] data,
        IFilter filter,
        ICover cover,
        IClusterer clusterer)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(filter);
        ArgumentNullException.ThrowIfNull(cover);
        ArgumentNullException.ThrowIfNull(clusterer);

        if (data.Length == 0)
            return EmptyResult(filter.Name, cover.Name, clusterer.Name, "Empty input data");

        double[] filterValues = filter.Apply(data);
        var coverResult = cover.Generate(filterValues);

        // Per-bin clustering: extract feature-vector subset, call point-cloud clusterer.
        return BuildFromBins(
            bins: coverResult.Bins,
            scalarFilter: filterValues,
            pointCount: data.Length,
            clusterBin: bin =>
                            {
                                var subset = new double[bin.PointIndices.Count][];
                                for (int i = 0; i < bin.PointIndices.Count; i++)
                                    subset[i] = data[bin.PointIndices[i]];
                                return clusterer.Cluster(subset);
                            },
            filterName: filter.Name,
            coverName: cover.Name,
            clustererName: clusterer.Name);
    }

    // ── Multi-D point-cloud MAPPER ──────────────────────────────────────────

    public static MapperResult Build(
        double[][] data,
        IMultiFilter filter,
        IMultiCover cover,
        IClusterer clusterer)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(filter);
        ArgumentNullException.ThrowIfNull(cover);
        ArgumentNullException.ThrowIfNull(clusterer);

        if (filter.Dimension != cover.Dimension)
            throw new ArgumentException(
                $"Filter dimension ({filter.Dimension}) does not match cover dimension ({cover.Dimension}).");

        if (data.Length == 0)
            return EmptyResult(filter.Name, cover.Name, clusterer.Name, "Empty input data");

        double[][] filterValues = filter.Apply(data);
        var coverResult = cover.Generate(filterValues);

        // For multi-D filters, use the first dimension's filter values for per-node
        // diagnostics. Downstream viz can recompute per-dim statistics if needed.
        var firstDimValues = new double[data.Length];
        for (int i = 0; i < data.Length; i++)
            firstDimValues[i] = filterValues[i].Length > 0 ? filterValues[i][0] : 0.0;

        return BuildFromBins(
            bins: coverResult.Bins,
            scalarFilter: firstDimValues,
            pointCount: data.Length,
            clusterBin: bin =>
                            {
                                var subset = new double[bin.PointIndices.Count][];
                                for (int i = 0; i < bin.PointIndices.Count; i++)
                                    subset[i] = data[bin.PointIndices[i]];
                                return clusterer.Cluster(subset);
                            },
            filterName: filter.Name,
            coverName: cover.Name,
            clustererName: clusterer.Name);
    }

    // ── Graph-input MAPPER ──────────────────────────────────────────────────

    /// <summary>
    /// Run MAPPER with a graph-input lens and graph-induced-subgraph clusterer.
    ///
    /// This is the variant that consumes an existing <see cref="CsrGraph"/> (e.g.,
    /// the proximity graph SPC will run on) and produces a nerve summarizing its
    /// topology. The clustering step uses connected components of the subgraph
    /// induced on each preimage, so two nodes share a cluster iff they're
    /// connected through other preimage nodes — a stronger statement than
    /// "close in metric coordinates."
    ///
    /// <paramref name="features"/> may be null for purely graph-topological
    /// filters (degree, Fiedler vector, geodesic distance from a seed). Pass
    /// features when using <c>PointFilterAdapter</c> to lift a feature-space
    /// filter (e.g., <c>HyperbolicFilters.PoincareRadial</c>) onto graph nodes.
    /// </summary>
    public static MapperResult Build(
        CsrGraph graph,
        double[][]? features,
        IGraphFilter filter,
        ICover cover,
        IGraphClusterer clusterer)
    {
        ArgumentNullException.ThrowIfNull(filter);
        ArgumentNullException.ThrowIfNull(cover);
        ArgumentNullException.ThrowIfNull(clusterer);

        int n = graph.NodeCount;
        if (n == 0)
            return EmptyResult(filter.Name, cover.Name, clusterer.Name, "Empty input graph");

        double[] filterValues = filter.Apply(graph, features);
        var coverResult = cover.Generate(filterValues);

        // Per-bin clustering: connected components of the subgraph induced on
        // the preimage. No feature extraction step — graph topology drives the cut.
        return BuildFromBins(
            bins: coverResult.Bins,
            scalarFilter: filterValues,
            pointCount: n,
            clusterBin: bin => clusterer.ClusterInduced(graph, bin.PointIndices),
            filterName: filter.Name,
            coverName: cover.Name,
            clustererName: clusterer.Name);
    }

    // ── Shared build path ───────────────────────────────────────────────────

    /// <summary>
    /// Given cover bins and a per-bin clustering delegate, build the MAPPER
    /// nodes and nerve graph. The clustering step is fully abstracted via
    /// <paramref name="clusterBin"/>, so this method is identical for
    /// point-cloud and graph variants.
    /// </summary>
    private static MapperResult BuildFromBins(
        IReadOnlyList<CoverBin> bins,
        double[] scalarFilter,
        int pointCount,
        Func<CoverBin, ClusterResult> clusterBin,
        string filterName,
        string coverName,
        string clustererName)
    {
        var nodes = new List<MapperNode>();
        int emptyBins = 0;

        foreach (var bin in bins)
        {
            if (bin.PointIndices.Count == 0)
            {
                emptyBins++;
                continue;
            }

            var cluster = clusterBin(bin);

            var clusterMembers = new Dictionary<int, List<int>>(cluster.K);
            for (int i = 0; i < bin.PointIndices.Count; i++)
            {
                int label = cluster.Labels[i];
                if (!clusterMembers.TryGetValue(label, out var list))
                {
                    list = new List<int>();
                    clusterMembers[label] = list;
                }
                list.Add(bin.PointIndices[i]);
            }

            foreach (var (localClusterId, members) in clusterMembers)
            {
                if (members.Count == 0) continue;

                double sum = 0.0, min = double.PositiveInfinity, max = double.NegativeInfinity;
                foreach (int idx in members)
                {
                    double v = scalarFilter[idx];
                    sum += v;
                    if (v < min) min = v;
                    if (v > max) max = v;
                }

                nodes.Add(new MapperNode(
                    BinId: bin.BinId,
                    LocalClusterId: localClusterId,
                    MemberIndices: members.ToArray(),
                    FilterValueMean: sum / members.Count,
                    FilterValueMin: min,
                    FilterValueMax: max));
            }
        }

        var nerve = BuildNerveInverted(nodes, pointCount);

        return new MapperResult
        {
            Nodes = nodes,
            Nerve = nerve,
            FilterName = filterName,
            CoverName = coverName,
            ClustererName = clustererName,
            EmptyBinCount = emptyBins,
        };
    }

    /// <summary>
    /// Inverted-index nerve construction. Complexity: O(Σ_p k_p²) where k_p is
    /// the number of nodes containing point p (typically small, bounded by
    /// cover overlap × clusters per bin).
    /// </summary>
    private static CsrGraph BuildNerveInverted(IReadOnlyList<MapperNode> nodes, int pointCount)
    {
        int m = nodes.Count;
        if (m == 0) return CsrGraph.FromEdges(Array.Empty<Edge>(), 0);

        var pointToNodes = new List<int>[pointCount];
        for (int p = 0; p < pointCount; p++) pointToNodes[p] = new List<int>(2);

        for (int nodeId = 0; nodeId < m; nodeId++)
            foreach (int p in nodes[nodeId].MemberIndices)
                pointToNodes[p].Add(nodeId);

        var nodeOverlaps = new Dictionary<int, int>[m];
        for (int u = 0; u < m; u++) nodeOverlaps[u] = new Dictionary<int, int>();

        for (int p = 0; p < pointCount; p++)
        {
            var pn = pointToNodes[p];
            int count = pn.Count;
            if (count < 2) continue;

            for (int i = 0; i < count; i++)
            {
                int u = pn[i];
                for (int j = i + 1; j < count; j++)
                {
                    int v = pn[j];
                    int lo = u < v ? u : v;
                    int hi = u < v ? v : u;
                    nodeOverlaps[lo].TryGetValue(hi, out int cur);
                    nodeOverlaps[lo][hi] = cur + 1;
                }
            }
        }

        int totalEdges = 0;
        for (int u = 0; u < m; u++) totalEdges += nodeOverlaps[u].Count;

        var edges = new Edge[totalEdges];
        int eIdx = 0;
        for (int u = 0; u < m; u++)
            foreach (var kvp in nodeOverlaps[u])
                edges[eIdx++] = new Edge(u, kvp.Key, kvp.Value);

        return CsrGraph.FromEdges(edges, m);
    }

    private static MapperResult EmptyResult(
        string filterName, string coverName, string clustererName, string warning) =>
        new()
        {
            Nodes = Array.Empty<MapperNode>(),
            Nerve = CsrGraph.FromEdges(Array.Empty<Edge>(), 0),
            FilterName = filterName,
            CoverName = coverName,
            ClustererName = clustererName,
            EmptyBinCount = 0,
        };
}
