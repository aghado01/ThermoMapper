using System;
using System.Collections.Generic;

namespace Viz;

/// <summary>
/// The invariant core: N points in D dimensions, optionally named,
/// with no opinion about where they came from.
/// </summary>
public sealed class PointCloud
{
    public ReadOnlyMemory<double> Features { get; }  // N x D, row-major
    public int N { get; }
    public int D { get; }
    public string? Label { get; }

    public PointCloud(double[] features, int n, int d, string? label = null)
    {
        Features = features;
        N = n;
        D = d;
        Label = label;
    }
}

/// <summary>
/// Shared contract for all named layer types.
/// Enables type-safe filtering in SceneBuilder without dynamic dispatch.
/// </summary>
public interface INamedLayer
{
    string Id { get; }
    string Name { get; }
}

internal static class LayerIdentity
{
    public static string Resolve(string? id, string name) =>
        string.IsNullOrWhiteSpace(id) ? name : id;
}

/// <summary>
/// A named assignment of each of the N points to an integer label.
/// "True labels", "SPC @ T=0.1", "GMM component", "GMM cluster (via map)" are
/// all just LabelLayers with different names and sources.
/// Negative values are legal: -1 = unassigned / noise.
/// </summary>
public sealed class LabelLayer : INamedLayer
{
    public string Id { get; }
    public string Name { get; }
    public ReadOnlyMemory<int> Labels { get; }
    public LabelLayerKind Kind { get; }

    public LabelLayer(string name, int[] labels, LabelLayerKind kind, string? id = null)
    {
        Id = LayerIdentity.Resolve(id, name);
        Name = name;
        Labels = labels;
        Kind = kind;
    }
}

public enum LabelLayerKind
{
    GroundTruth,
    SpinColor,
    EquilibriumCluster,
    GmmComponent,
    GmmCluster,
    Custom,
}

/// <summary>
/// A named per-node scalar function on the graph 0-skeleton: filter functions,
/// eigenfunctions, susceptibility, Mahalanobis distance, coherence, responsibility[k],
/// log-likelihood, etc.
/// </summary>
public sealed class NodeSignalLayer : INamedLayer
{
    public string Id { get; }
    public string Name { get; }
    public ReadOnlyMemory<double> Values { get; }
    public ScalarSource Source { get; }

    public NodeSignalLayer(string name, double[] values, ScalarSource source, string? id = null)
    {
        Id = LayerIdentity.Resolve(id, name);
        Name = name;
        Values = values;
        Source = source;
    }
}

public enum ScalarSource
{
    FilterFunction,
    Eigenfunction,
    Susceptibility,
    CoherenceScore,
    MahalanobisDistance,
    PercolationArrival,
    Responsibility,
    LogLikelihood,
    Custom,
}

/// <summary>
/// A named per-node unoriented direction field on the graph 0-skeleton.
/// d and -d are mathematically equivalent: this is a line in tangent space, not a vector.
/// </summary>
public sealed class LineFieldLayer : INamedLayer
{
    public string Id { get; }
    public string Name { get; }
    public ReadOnlyMemory<double> Directions { get; }
    public int N { get; }
    public int D { get; }
    public LineFieldSource Source { get; }

    public LineFieldLayer(string name, double[] directions, int n, int d, LineFieldSource source, string? id = null)
    {
        Id = LayerIdentity.Resolve(id, name);
        Name = name;
        Directions = directions;
        N = n;
        D = d;
        Source = source;
    }
}

public enum LineFieldSource
{
    SpectralGradient,
    LocalPca,
    Custom,
}
