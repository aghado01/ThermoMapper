// ============================================================================
// TDA.Mapper — default_filters.cs
// ============================================================================
// Standard MAPPER filters (geometry-agnostic).
//
// Hyperbolic-specialized filters live in HyperbolicFilters.cs.
// PCA-based filters live in pca_filter.cs (they depend on the LinearAlgebra primitive).
// ============================================================================

#nullable enable
using System;
using TDA.Mapper;

namespace TDA.Mapper.Filters;

/// <summary>Factory for standard (geometry-agnostic) filter implementations.</summary>
public static class Filters
{
    /// <summary>Use the first coordinate of each point as the filter value.
    /// Useful when data is already projected to a meaningful 1-D axis upstream.</summary>
    public static IFilter Identity => new IdentityFilter();

    /// <summary>Euclidean L2 norm — distance from origin in flat coordinates.
    /// For data lying in B^n, prefer <see cref="HyperbolicFilters.PoincareRadial"/>
    /// — Euclidean norm and Poincaré radial distance are different functions.</summary>
    public static IFilter EuclideanNorm => new EuclideanNormFilter();

    /// <summary>First principal component projection. Computes PCA once over
    /// the full dataset, then projects each point onto PC1.</summary>
    public static IFilter Pca1D(double[][] data) => new Pca1DFilter(data);
}

internal sealed class IdentityFilter : IFilter
{
    public string Name => "Identity (first coordinate)";

    public double[] Apply(double[][] data)
    {
        ArgumentNullException.ThrowIfNull(data);
        var values = new double[data.Length];
        for (int i = 0; i < data.Length; i++)
            values[i] = data[i].Length > 0 ? data[i][0] : 0.0;
        return values;
    }
}

internal sealed class EuclideanNormFilter : IFilter
{
    public string Name => "Euclidean norm";

    public double[] Apply(double[][] data)
    {
        ArgumentNullException.ThrowIfNull(data);
        var values = new double[data.Length];
        for (int i = 0; i < data.Length; i++)
        {
            double sumSq = 0.0;
            for (int d = 0; d < data[i].Length; d++)
                sumSq += data[i][d] * data[i][d];
            values[i] = Math.Sqrt(sumSq);
        }
        return values;
    }
}
