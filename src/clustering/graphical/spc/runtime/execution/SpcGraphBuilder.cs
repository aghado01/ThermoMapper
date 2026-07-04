using System;
using Graphs;
using Graphs.Coupling;
using Graphs.Distance;
using Graphs.Primitives;
using Graphs.Proximity;

namespace Clustering.Graphical.SPC.Runtime.Execution;

public static class SpcGraphBuilder
{
    /// <summary>
    /// Build the SPC graph for a feature matrix, returning the rich
    /// <see cref="GraphBuildResult"/> that surfaces the directed
    /// <see cref="NeighborSelection"/>, resolved bandwidth, and pre-repair
    /// graph alongside the final <see cref="CsrGraph"/>. Callers that only
    /// want the graph should pick <c>.Graph</c> off the result; downstream
    /// diagnostics (<see cref="Graphs.Diagnostics.GraphHealth"/>) consume
    /// the full shape without re-running KNN.
    /// </summary>
    public static GraphBuildResult BuildResult(
        double[][] features,
        GraphCompilerConfig config,
        IDistanceMetric? metric = null,
        ProtectedEdgeSource? protectedEdges = null)
    {
        if (features is null)
            throw new ArgumentNullException(nameof(features));
        if (config is null)
            throw new ArgumentNullException(nameof(config));

        if (features.Length == 0)
            throw new ArgumentException("Features must contain at least one observation.", nameof(features));

        int dimension = features[0]?.Length ?? throw new ArgumentException(
            "Feature vectors must not be null.", nameof(features));
        for (int i = 1; i < features.Length; i++)
        {
            if (features[i] is null)
                throw new ArgumentException($"Feature vector at index {i} is null.", nameof(features));
            if (features[i].Length != dimension)
                throw new ArgumentException(
                    $"Feature vector at index {i} has dimension {features[i].Length}, expected {dimension}.",
                    nameof(features));
        }

        return GraphCompiler.Build(config, features.Length, GraphMetric.FromFeatures(features, metric), protectedEdges);
    }
}
