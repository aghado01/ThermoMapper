// ============================================================================
// TDA.Mapper — HyperbolicFilters.cs
// ============================================================================
// Hyperbolic-aware MAPPER filters for data living in B^n (open Poincaré ball).
//
// These filters use geodesic distances rather than Euclidean. The geometric
// difference matters for hierarchical hyperbolic data: Euclidean norm on B^n
// coordinates is *not* a proxy for hierarchy depth, but Poincaré radial distance
// from origin *is* — by construction in the HyperbolicHierarchy generator,
// hierarchy depth maps onto radial position.
//
// PoincareRadial is the headline filter for the HyperbolicHierarchy diagnostic:
// a tree-structured dataset under this filter should produce a tree-like nerve.
// ============================================================================

#nullable enable
using System;
using TDA.Mapper;

namespace TDA.Mapper.Filters;

/// <summary>Factory for hyperbolic-aware filter implementations.</summary>
public static class HyperbolicFilters
{
    /// <summary>
    /// Poincaré geodesic distance from origin in B^n: <c>d(0, x) = 2·arctanh(||x||)</c>.
    /// The most informative MAPPER filter for HyperbolicHierarchy and similar
    /// tree-structured hyperbolic data — hierarchy depth maps onto this radial
    /// distance by construction.
    ///
    /// Points on or past the unit-ball boundary (||x|| ≥ 1 − ε) emit +Infinity;
    /// they are excluded from all cover bins by <see cref="UniformCover"/>.
    /// </summary>
    public static IFilter PoincareRadial => new PoincareRadialFilter();

    /// <summary>
    /// Local density estimate via Poincaré-distance k-NN: density(x) = 1 / mean k-NN distance.
    /// Useful when hierarchical structure is irregular (different branches have
    /// different densities). For tree-structured data, PoincareRadial is usually
    /// a stronger lens.
    /// </summary>
    public static IFilter PoincareLocalDensity(int k = 10) => new PoincareLocalDensityFilter(k);
}

internal sealed class PoincareRadialFilter : IFilter
{
    public string Name => "Poincaré radial (2·arctanh(||x||))";

    public double[] Apply(double[][] data)
    {
        ArgumentNullException.ThrowIfNull(data);
        const double boundaryEps = 1e-12;

        var values = new double[data.Length];
        for (int i = 0; i < data.Length; i++)
        {
            double normSq = 0.0;
            for (int d = 0; d < data[i].Length; d++)
                normSq += data[i][d] * data[i][d];
            double norm = Math.Sqrt(normSq);

            // Poincaré geodesic distance from origin: 2·arctanh(||x||).
            // Defined only for ||x|| < 1; boundary points emit +∞ (excluded by cover).
            values[i] = norm < 1.0 - boundaryEps
                ? 2.0 * Math.Atanh(norm)
                : double.PositiveInfinity;
        }
        return values;
    }
}

internal sealed class PoincareLocalDensityFilter : IFilter
{
    private readonly int _k;

    public string Name => $"Poincaré local density (k={_k})";

    public PoincareLocalDensityFilter(int k)
    {
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k), "k must be >= 1");
        _k = k;
    }

    public double[] Apply(double[][] data)
    {
        ArgumentNullException.ThrowIfNull(data);
        int n = data.Length;
        if (n == 0) return Array.Empty<double>();

        var values = new double[n];
        var distances = new double[n];

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                distances[j] = i == j ? double.PositiveInfinity : PoincareDistance(data[i], data[j]);

            // Partial sort: smallest k via in-place selection.
            // For sketch simplicity, full sort — Apply is one-shot, not hot-path.
            Array.Sort(distances);

            int count = Math.Min(_k, n - 1);
            double sumNn = 0.0;
            int finite = 0;
            for (int kk = 0; kk < count; kk++)
            {
                if (double.IsFinite(distances[kk]))
                {
                    sumNn += distances[kk];
                    finite++;
                }
            }

            if (finite == 0)
                values[i] = 0.0;          // isolated point in finite-distance sense
            else
            {
                double meanNn = sumNn / finite;
                values[i] = meanNn > 0 ? 1.0 / meanNn : double.PositiveInfinity;
            }
        }
        return values;
    }

    /// <summary>
    /// Poincaré distance on B^n: <c>d(x, y) = arcosh(1 + 2·||x − y||² / ((1 − ||x||²)(1 − ||y||²)))</c>.
    /// Returns +∞ if either point is at or past the boundary.
    /// </summary>
    private static double PoincareDistance(double[] x, double[] y)
    {
        double normXSq = 0, normYSq = 0, diffSq = 0;
        int dim = Math.Min(x.Length, y.Length);
        for (int d = 0; d < dim; d++)
        {
            normXSq += x[d] * x[d];
            normYSq += y[d] * y[d];
            double diff = x[d] - y[d];
            diffSq += diff * diff;
        }
        double denom = (1.0 - normXSq) * (1.0 - normYSq);
        if (denom <= 1e-12) return double.PositiveInfinity;

        double arg = 1.0 + 2.0 * diffSq / denom;
        return Math.Acosh(arg);
    }
}
