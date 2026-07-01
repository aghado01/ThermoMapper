using System;
using System.Collections.Generic;

namespace Viz;

/// <summary>
/// A named sparse graph over the N points: k-NN, mutual k-NN, MST edges, etc.
/// Used to render connectivity and highlight false bridges.
/// </summary>
public sealed class EdgeLayer : INamedLayer
{
    public string Id { get; }
    public string Name { get; }
    public ReadOnlyMemory<int> Src { get; }
    public ReadOnlyMemory<int> Dst { get; }
    public ReadOnlyMemory<double> Weight { get; }
    public string? Metric { get; }
    public ProximitySpec? Proximity { get; }
    public ReadOnlyMemory<int>? EdgeClusterSrc { get; }
    public ReadOnlyMemory<int>? EdgeClusterDst { get; }

    public EdgeLayer(
        string name,
        int[] src,
        int[] dst,
        double[] weight,
        string? metric = null,
        ProximitySpec? proximity = null,
        int[]? edgeClusterSrc = null,
        int[]? edgeClusterDst = null,
        string? id = null)
    {
        Id = LayerIdentity.Resolve(id, name);
        Name = name;
        Src = src;
        Dst = dst;
        Weight = weight;
        Metric = metric;
        Proximity = proximity;
        EdgeClusterSrc = edgeClusterSrc;
        EdgeClusterDst = edgeClusterDst;
    }
}

public sealed class TriangleLayer : INamedLayer
{
    public string Id { get; }
    public string Name { get; }
    public TriangleSource Source { get; }
    public string SourceEdgeLayerId { get; }
    public ReadOnlyMemory<int> Vertices { get; }

    public TriangleLayer(
        string name,
        TriangleSource source,
        string sourceEdgeLayerId,
        int[] vertices,
        string? id = null)
    {
        if (vertices is null) throw new ArgumentNullException(nameof(vertices));
        if (vertices.Length % 3 != 0)
            throw new ArgumentException("Triangle vertices must be a flat array of length 3*T.", nameof(vertices));

        Id = LayerIdentity.Resolve(id, name);
        Name = name;
        Source = source;
        SourceEdgeLayerId = sourceEdgeLayerId;
        Vertices = vertices;
    }

    public static TriangleLayer FromFlagComplex(
        EdgeLayer edges,
        int[] vertices,
        string? name = null,
        string? id = null)
    {
        if (edges is null) throw new ArgumentNullException(nameof(edges));

        return new TriangleLayer(
            name ?? $"{edges.Name} Triangles",
            TriangleSource.FlagComplex,
            edges.Id,
            vertices,
            id);
    }
}

public enum TriangleSource
{
    FlagComplex,
    Filtered,
    Custom,
}

public sealed class GaussianLayer : INamedLayer
{
    public string Id { get; }
    public string Name { get; }
    public ReadOnlyMemory<double> Means { get; }
    public ReadOnlyMemory<double> Covariances { get; }
    public ReadOnlyMemory<double> Weights { get; }
    public int K { get; }
    public int D { get; }
    public ReadOnlyMemory<int>? ComponentToClusterMap { get; }

    public GaussianLayer(
        string name,
        double[] means,
        double[] covariances,
        double[] weights,
        int k,
        int d,
        int[]? componentToClusterMap = null,
        string? id = null)
    {
        Id = LayerIdentity.Resolve(id, name);
        Name = name;
        Means = means;
        Covariances = covariances;
        Weights = weights;
        K = k;
        D = d;
        ComponentToClusterMap = componentToClusterMap;
    }
}

public sealed class SpineLayer : INamedLayer
{
    public string Id { get; }
    public string Name { get; }
    public int ClusterIdx { get; }
    public SpineLayerKind Kind { get; }
    public double[][] SpineSamples { get; }
    public double[][][]? TangentBases { get; }
    public double TypicalScale { get; }

    public SpineLayer(
        string name,
        int clusterIdx,
        SpineLayerKind kind,
        double[][] spineSamples,
        double[][][]? tangentBases,
        double typicalScale = 0.0,
        string? id = null)
    {
        Id = LayerIdentity.Resolve(id, name);
        Name = name;
        ClusterIdx = clusterIdx;
        Kind = kind;
        SpineSamples = spineSamples;
        TangentBases = tangentBases;
        TypicalScale = typicalScale;
    }
}

public enum SpineLayerKind { Arc, Manifold, MobiusTube }

public sealed class VectorFieldLayer : INamedLayer
{
    public string Id { get; }
    public string Name { get; }
    public ReadOnlyMemory<double> Vectors { get; }
    public int N { get; }
    public int D { get; }

    public VectorFieldLayer(string name, double[] vectors, int n, int d, string? id = null)
    {
        Id = LayerIdentity.Resolve(id, name);
        Name = name;
        Vectors = vectors;
        N = n;
        D = d;
    }
}
