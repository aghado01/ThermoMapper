using System;
using Graphs.Distance;
using Graphs.Distance.Geodesic;
using Maths.Geometry;

namespace Graphs;

/// <summary>
/// Data-substrate binding for graph construction: a dataset-bound distance function
/// paired with the <see cref="MetricProperties"/> describing the geometry it realizes.
/// The two faces of one metric travel together so they cannot drift. Passed to
/// <see cref="GraphCompiler.Build(GraphCompilerConfig, int, GraphMetric)"/> alongside the
/// declarative <see cref="GraphCompilerConfig"/>.
/// </summary>
/// <remarks>
/// <para><see cref="Properties"/> is null for metric-less construction (raw distance,
/// test presets, sweeps); strategy then resolves from the projection's
/// <c>BandwidthOverride</c> or falls to the historical MAD default.</para>
/// </remarks>
public sealed record GraphMetric(
    Func<int, int, double> Distance,
    MetricProperties? Properties = null,
    int? AmbientDimension = null,
    double[][]? Features = null,
    IRiemannianManifold? Manifold = null)
{
    /// <summary>
    /// Bind an <see cref="IDistanceMetric"/> to a concrete dataset, producing the
    /// index-level closure plus the carried properties. Null metric keeps the
    /// historical inline Euclidean fallback with null provenance.
    /// </summary>
    public static GraphMetric FromFeatures(double[][] features, IDistanceMetric? metric = null)
    {
        if (features is null) throw new ArgumentNullException(nameof(features));

        int dimension = features.Length > 0 ? features[0].Length : 0;
        metric?.Properties.ValidateDimension(dimension);

        Func<int, int, double> distance = metric is not null
            ? (i, j) => metric.Distance(features[i], features[j])
            : (i, j) => Euclidean(features[i], features[j], dimension);

        IRiemannianManifold? manifold = metric switch
        {
            PoincareMetric => new PoincareBallManifold(dimension),
            SphericalGeodesicMetric => new SphericalManifold(dimension),
            FisherRaoSimplexMetric => new SphericalManifold(dimension),
            _ => null,
        };

        return new GraphMetric(distance, metric?.Properties, dimension, features, manifold);
    }

    private static double Euclidean(ReadOnlySpan<double> xi, ReadOnlySpan<double> xj, int dim)
    {
        double sum = 0.0;
        for (int d = 0; d < dim; d++)
        {
            double diff = xi[d] - xj[d];
            sum += diff * diff;
        }

        return Math.Sqrt(sum);
    }
}
