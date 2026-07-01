// ============================================================================
// TDA.Mapper — GraphFilters.cs
// ============================================================================
// Graph-input lens functions (IGraphFilter). Two families:
//
//   Graph-topological filters  — use only graph structure (Degree,
//                                WeightedDegree, GeodesicDistance).
//                                Pass features = null at the Mapper.Build call.
//
//   Lifted feature filters     — wrap a point-cloud IFilter to operate on
//                                graph nodes (PointFilterAdapter). Pass
//                                features as the node feature vectors.
//
// FiedlerVectorFilter lives in its own file because it carries a non-trivial
// dependency chain (Maths.LinAlg + Graphs.Diagnostics), mirroring how
// Pca1DFilter is split out from DefaultFilters.cs.
//
// References:
//   Carrière, Michel, Oudot (2018) "Statistical analysis and parameter
//     selection for Mapper" — graph MAPPER as a variant
//   Hajij, Rosen, Wang (2018) "Mapper on graphs for network visualization"
// ============================================================================

#nullable enable
using System;
using Graphs.Primitives;
using TDA.Mapper;

namespace TDA.Mapper.Filters;

// ── Factory ─────────────────────────────────────────────────────────────────

/// <summary>Factory for graph-input filter implementations.</summary>
public static class GraphFilters
{
    /// <summary>Unweighted node degree: number of neighbors per node.</summary>
    public static IGraphFilter Degree => new DegreeFilter();

    /// <summary>Weighted node degree: sum of edge weights (J values) per node.
    /// More informative than unweighted degree on kernel-weighted graphs.</summary>
    public static IGraphFilter WeightedDegree => new WeightedDegreeFilter();

    /// <summary>BFS hop distance from <paramref name="seedNode"/> to every other
    /// node. Unreachable nodes get +Infinity. The natural graph-native analog
    /// of "hierarchy depth from root" — for HyperbolicHierarchy with the root
    /// at the graph node closest to origin, this is a strong topology probe.
    /// Pass <paramref name="maxDepth"/> to cap the BFS (e.g., for k-hop local
    /// neighborhood filters).</summary>
    public static IGraphFilter GeodesicDistance(int seedNode, int? maxDepth = null)
        => new GeodesicDistanceFilter(seedNode, maxDepth);

    /// <summary>Fiedler vector (eigenvector of the second-smallest eigenvalue
    /// of the weighted graph Laplacian). Captures the dominant connectivity
    /// gradient — the "principal axis" of the graph. Requires a connected
    /// graph (use <c>ConnectivityRepair.EnsureConnected</c> for MST-bridged repair).
    /// Uses dense eigendecomposition via <c>Maths.LinAlg.DenseEigen.DecomposeSymmetric</c>;
    /// suitable for graphs up to ~few thousand nodes.</summary>
    public static IGraphFilter FiedlerVector => new FiedlerVectorFilter();

    /// <summary>Lift a point-cloud <see cref="IFilter"/> to operate on graph
    /// nodes by applying it to the node feature vectors. The graph topology
    /// is unused; only the lifted filter's per-point logic runs.
    ///
    /// Use this to compose graph MAPPER with existing feature-space filters
    /// like <c>HyperbolicFilters.PoincareRadial</c>. Pass non-null
    /// <c>features</c> to <c>Mapper.Build(graph, features, filter, cover, clusterer)</c>.</summary>
    public static IGraphFilter FromPointFilter(IFilter pointFilter) => new PointFilterAdapter(pointFilter);
}

// ── DegreeFilter / WeightedDegreeFilter ─────────────────────────────────────

internal sealed class DegreeFilter : IGraphFilter
{
    public string Name => "Node degree";

    public double[] Apply(CsrGraph graph, double[][]? features = null)
    {
        int n = graph.NodeCount;
        var degrees = new double[n];
        for (int i = 0; i < n; i++)
            degrees[i] = graph.RowPointers[i + 1] - graph.RowPointers[i];
        return degrees;
    }
}

internal sealed class WeightedDegreeFilter : IGraphFilter
{
    public string Name => "Weighted node degree (Σ_j J_ij)";

    public double[] Apply(CsrGraph graph, double[][]? features = null)
    {
        int n = graph.NodeCount;
        var weighted = new double[n];
        for (int i = 0; i < n; i++)
        {
            int rowStart = graph.RowPointers[i];
            int rowEnd = graph.RowPointers[i + 1];
            double sum = 0.0;
            for (int e = rowStart; e < rowEnd; e++) sum += graph.Weights[e];
            weighted[i] = sum;
        }
        return weighted;
    }
}

// ── GeodesicDistanceFilter ──────────────────────────────────────────────────

internal sealed class GeodesicDistanceFilter : IGraphFilter
{
    private readonly int _seedNode;
    private readonly int? _maxDepth;

    public string Name => _maxDepth.HasValue
        ? $"Geodesic distance from node {_seedNode} (BFS hops, max depth {_maxDepth.Value})"
        : $"Geodesic distance from node {_seedNode} (BFS hops)";

    public GeodesicDistanceFilter(int seedNode, int? maxDepth = null)
    {
        if (seedNode < 0) throw new ArgumentOutOfRangeException(nameof(seedNode), "seedNode must be >= 0");
        if (maxDepth.HasValue && maxDepth.Value < 0)
            throw new ArgumentOutOfRangeException(nameof(maxDepth), "maxDepth must be >= 0 or null");
        _seedNode = seedNode;
        _maxDepth = maxDepth;
    }

    public double[] Apply(CsrGraph graph, double[][]? features = null)
    {
        // Delegate to the canonical BFS primitive. ComputeDistances returns int
        // hop counts with -1 sentinel for unreachable; convert to double with
        // +Infinity for unreachable so cover filtering excludes them automatically.
        int[] hops = BfsShells.ComputeDistances(graph, _seedNode, _maxDepth);

        int n = hops.Length;
        var distances = new double[n];
        for (int i = 0; i < n; i++)
            distances[i] = hops[i] == -1 ? double.PositiveInfinity : (double)hops[i];
        return distances;
    }
}

// ── PointFilterAdapter ──────────────────────────────────────────────────────

/// <summary>
/// Lifts a point-cloud <see cref="IFilter"/> to operate on graph nodes by
/// applying the underlying filter to the node feature vectors. Graph topology
/// is unused — this is the bridge for "I want PoincareRadial on the graph,
/// please."
/// </summary>
internal sealed class PointFilterAdapter : IGraphFilter
{
    private readonly IFilter _pointFilter;

    public string Name => $"Point→Graph: {_pointFilter.Name}";

    public PointFilterAdapter(IFilter pointFilter)
    {
        ArgumentNullException.ThrowIfNull(pointFilter);
        _pointFilter = pointFilter;
    }

    public double[] Apply(CsrGraph graph, double[][]? features = null)
    {
        if (features is null)
            throw new ArgumentException(
                $"PointFilterAdapter (wrapping {_pointFilter.Name}) requires non-null features — " +
                "the wrapped point-cloud filter operates on node feature vectors.",
                nameof(features));

        if (features.Length != graph.NodeCount)
            throw new ArgumentException(
                $"Features length ({features.Length}) does not match graph node count ({graph.NodeCount}).",
                nameof(features));

        return _pointFilter.Apply(features);
    }
}
