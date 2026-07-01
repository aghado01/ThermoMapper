// ============================================================================
// TDA.Mapper — pca_filter.cs
// ============================================================================
// PCA-based filter implementation. Wires to the LinearAlgebra.Pca primitive in
// src/dev/linalg/PCA.cs (which delegates eigendecomposition to LinearAlgebra.Eigen).
// Separated from default_filters.cs because it carries a non-trivial dependency.
// ============================================================================

#nullable enable
using System;
using Maths.LinAlg;
using TDA.Mapper;

namespace TDA.Mapper.Filters;

/// <summary>
/// Projects each point onto the first principal component of the full dataset.
/// PCA is computed once at construction; Apply is O(N·d) per call.
/// </summary>
internal sealed class Pca1DFilter : IFilter
{
    private readonly double[] _direction;
    private readonly double[] _mean;

    public string Name => "PCA-1D (first principal component)";

    public Pca1DFilter(double[][] data)
    {
        ArgumentNullException.ThrowIfNull(data);
        if (data.Length < 2)
            throw new ArgumentException("PCA filter requires at least 2 data points.", nameof(data));

        var pca = Pca.Compute(data, numComponents: 1);
        _direction = pca.Components[0];     // first principal component
        _mean = pca.Mean;
    }

    public double[] Apply(double[][] data)
    {
        ArgumentNullException.ThrowIfNull(data);
        var values = new double[data.Length];

        for (int i = 0; i < data.Length; i++)
        {
            double proj = 0.0;
            int dim = Math.Min(data[i].Length, _direction.Length);
            for (int d = 0; d < dim; d++)
                proj += (data[i][d] - _mean[d]) * _direction[d];
            values[i] = proj;
        }
        return values;
    }
}
