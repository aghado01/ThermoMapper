using System;
using Graphs.Distance;
using Graphs.Distance.Geodesic;
using Maths.Distance.Geodesic;
using Maths.Geometry;
using Xunit;

namespace VizCore.Tests;

public sealed class DistanceGeodesicTests
{
    [Fact]
    public void SphericalGeodesic_SameVectors_ReturnsZero()
    {
        double[] a = { 1.0, 2.0, 3.0 };
        double result = SphericalGeodesic.Distance(a, a);
        Assert.Equal(0.0, result, 10);
    }

    [Fact]
    public void SphericalGeodesic_OrthogonalVectors_ReturnsPiOver2()
    {
        double[] a = { 1.0, 0.0, 0.0 };
        double[] b = { 0.0, 1.0, 0.0 };
        double result = SphericalGeodesic.Distance(a, b);
        Assert.Equal(Math.PI / 2.0, result, 10);
    }

    [Fact]
    public void Wasserstein1_PmfVectors_CorrectDistance()
    {
        double[] p = { 0.2, 0.5, 0.3 };
        double[] q = { 0.1, 0.3, 0.6 };

        // CDF_p = [0.2, 0.7, 1.0]
        // CDF_q = [0.1, 0.4, 1.0]
        // sum |CDF_p - CDF_q| = 0.1 + 0.3 + 0.0 = 0.4
        double expected = 0.4;
        double result = Wasserstein1.Distance(p, q);

        Assert.Equal(expected, result, 10);
    }

    [Fact]
    public void Wasserstein1Metric_ImplementsDistanceMetricContract()
    {
        var metric = new Wasserstein1Metric();
        double[] pArray = { 0.1, 0.9 };
        double[] qArray = { 0.5, 0.5 };

        double result = metric.Distance(pArray.AsSpan(), qArray.AsSpan());
        Assert.Equal(0.4, result, 10);
    }

    [Fact]
    public void PoincareMetric_ImplementsDistanceMetricContract()
    {
        var metric = new PoincareMetric();
        double[] a = { 0.1, 0.2, 0.3 };
        double[] b = { 0.2, 0.3, 0.4 };

        double result = metric.Distance(a.AsSpan(), b.AsSpan());
        Assert.Equal(Poincare.Distance(a, b), result, 10);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(8)]
    public void PoincareBallManifold_LogExp_RoundTrips_AcrossDimensions(int dimension)
    {
        var manifold = new PoincareBallManifold(dimension);
        double[] p = CreatePointFixture(dimension, 0.10, -0.15, 0.20, 0.05);
        double[] q = CreatePointFixture(dimension, 0.20, -0.05, 0.10, -0.03);
        var tangent = new double[dimension];
        var recovered = new double[dimension];

        manifold.LogMap(p, q, tangent);
        manifold.ExpMap(p, tangent, recovered);

        for (int i = 0; i < q.Length; i++)
            Assert.Equal(q[i], recovered[i], 9);

        Assert.Equal(manifold.Distance(p, q), manifold.Norm(p, tangent), 9);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(8)]
    public void PoincareBallManifold_LogMap_RiemannianNormMatchesDistance_NotAmbientNorm(int dimension)
    {
        var manifold = new PoincareBallManifold(dimension);
        double[] p = CreatePointFixture(dimension, 0.22, -0.18, 0.11, 0.04);
        double[] q = CreatePointFixture(dimension, -0.05, 0.08, -0.02, 0.03);
        var tangent = new double[dimension];

        manifold.LogMap(p, q, tangent);

        double distance = manifold.Distance(p, q);
        double riemannianNorm = manifold.Norm(p, tangent);
        double ambientNorm = EuclideanNorm(tangent);
        double conformalFactor = 2.0 / (1.0 - SquaredNorm(p));

        Assert.Equal(distance, riemannianNorm, 9);
        Assert.Equal(distance, conformalFactor * ambientNorm, 9);
        Assert.True(
            Math.Abs(distance - ambientNorm) > 1e-3,
            $"Expected the ambient tangent norm to differ from geodesic distance away from the origin. distance={distance}, ambient={ambientNorm}");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(8)]
    public void PoincareBallManifold_DirectedRadii_AgreeWithGeodesicDistance(int dimension)
    {
        var manifold = new PoincareBallManifold(dimension);
        double[] p = CreatePointFixture(dimension, 0.16, -0.12, 0.09, 0.02);
        double[] q = CreatePointFixture(dimension, -0.09, 0.14, -0.04, -0.01);
        var forward = new double[dimension];
        var backward = new double[dimension];

        manifold.LogMap(p, q, forward);
        manifold.LogMap(q, p, backward);

        double distance = manifold.Distance(p, q);
        double forwardRadius = manifold.Norm(p, forward);
        double backwardRadius = manifold.Norm(q, backward);

        Assert.Equal(distance, forwardRadius, 9);
        Assert.Equal(distance, backwardRadius, 9);
        Assert.Equal(forwardRadius, backwardRadius, 9);
    }

    [Fact]
    public void FisherRaoSimplexMetric_ImplementsDistanceMetricContract()
    {
        var metric = new FisherRaoSimplexMetric();
        double[] pArray = { 0.2, 0.5, 0.3 };
        double[] qArray = { 0.1, 0.3, 0.6 };

        double result = metric.Distance(pArray.AsSpan(), qArray.AsSpan());
        Assert.Equal(FisherRaoSimplex.Distance(pArray, qArray), result, 10);
    }

    [Fact]
    public void FisherRaoHalfPlaneMetric_ImplementsDistanceMetricContract()
    {
        var metric = new FisherRaoHalfPlaneMetric();
        double[] p = { 0.0, Math.Log(1.0) };
        double[] q = { 1.0, Math.Log(2.0) };

        double result = metric.Distance(p.AsSpan(), q.AsSpan());
        Assert.Equal(FisherRaoHalfPlane.Distance(p, q), result, 10);
    }

    private static double[] CreatePointFixture(int dimension, params double[] seed)
    {
        var point = new double[dimension];
        for (int i = 0; i < dimension; i++)
            point[i] = i < seed.Length ? seed[i] : 0.01 * (i + 1);

        return point;
    }

    private static double SquaredNorm(double[] vector)
    {
        double sum = 0.0;
        for (int i = 0; i < vector.Length; i++)
            sum += vector[i] * vector[i];

        return sum;
    }

    private static double EuclideanNorm(double[] vector)
        => Math.Sqrt(SquaredNorm(vector));
}
