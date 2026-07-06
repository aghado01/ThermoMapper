#nullable enable
using System;
using System.Collections.Generic;
using Graphs.Primitives;

namespace TDA.Ph;

/// <summary>
/// Full Vietoris-Rips builder for small point clouds: construct the complete distance graph
/// (optionally thresholded) and delegate to the graph-restricted Rips materializer.
/// </summary>
public static class FullRips
{
    /// <summary>
    /// Build a threshold-bounded full Vietoris-Rips filtration from a Euclidean point cloud.
    /// Edges with distance greater than <paramref name="threshold"/> are omitted; the default
    /// includes every pair. The current Rips materializer emits the 2-skeleton, so
    /// <paramref name="maxDimension"/> above 2 is accepted but still capped by that path.
    /// </summary>
    public static SimplicialFiltration Build(
        double[][] points,
        int maxDimension = 2,
        double threshold = double.PositiveInfinity,
        string label = "FullRips")
    {
        ArgumentNullException.ThrowIfNull(points);
        if (maxDimension < 0)
            throw new ArgumentOutOfRangeException(nameof(maxDimension));
        if (double.IsNaN(threshold) || threshold < 0.0)
            throw new ArgumentOutOfRangeException(nameof(threshold), "Full Rips threshold must be non-negative.");

        int n = points.Length;
        if (n == 0)
        {
            var emptyGraph = new CsrGraph
            {
                Targets = Array.Empty<int>(),
                Weights = Array.Empty<double>(),
                RowPointers = new[] { 0 },
                NodeCount = 0,
            };
            return RipsFiltration.RipsFromGraph(emptyGraph, FiltrationWeights.RawDistance, maxDimension, label);
        }

        int dimension = points[0]?.Length ?? throw new ArgumentException("Point rows must not be null.", nameof(points));
        var edges = new List<Edge>(n * (n - 1) / 2);
        for (int i = 0; i < n; i++)
        {
            double[] pi = RequirePoint(points, i, dimension);
            for (int j = i + 1; j < n; j++)
            {
                double dist = EuclideanDistance(pi, RequirePoint(points, j, dimension));
                if (dist <= threshold)
                    edges.Add(new Edge(i, j, dist));
            }
        }

        CsrGraph graph = CsrGraph.FromEdges(edges.ToArray(), n);
        return RipsFiltration.RipsFromGraph(graph, FiltrationWeights.RawDistance, maxDimension, label);
    }

    private static double[] RequirePoint(double[][] points, int index, int dimension)
    {
        double[]? point = points[index];
        if (point is null)
            throw new ArgumentException("Point rows must not be null.", nameof(points));
        if (point.Length != dimension)
            throw new ArgumentException("Point cloud rows must have a consistent dimension.", nameof(points));
        return point;
    }

    private static double EuclideanDistance(double[] a, double[] b)
    {
        double sum = 0.0;
        for (int k = 0; k < a.Length; k++)
        {
            double d = a[k] - b[k];
            sum += d * d;
        }
        return Math.Sqrt(sum);
    }
}
