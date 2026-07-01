using System;
using System.Runtime.CompilerServices;
using Graphs.Distance;
using Maths.Distance.Geodesic;

namespace Graphs.Distance.Geodesic;

public readonly struct SphericalGeodesicMetric : IDistanceMetric
{
    private static readonly MetricProperties _props = new(
        IsBounded:           true,
        MaxValue:            Math.PI,
        RequiresProbability: false,
        RequiresUnitNorm:    true,
        FixedDimension:      null,
        Geometry:            SpaceGeometry.Spherical,
        BandwidthStrategy:   BandwidthStrategy.QuantileNormalized,
        Name:                "SphericalGeodesic");

    public MetricProperties Properties => _props;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
        => SphericalGeodesic.Distance(a, b);
}

public readonly struct PoincareMetric : IDistanceMetric
{
    private static readonly MetricProperties _props = new(
        IsBounded:           false,
        MaxValue:            0.0,
        RequiresProbability: false,
        RequiresUnitNorm:    false,
        FixedDimension:      null,
        Geometry:            SpaceGeometry.Hyperbolic,
        BandwidthStrategy:   BandwidthStrategy.LogScaleHyperbolic,
        Name:                "Poincare");

    public MetricProperties Properties => _props;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
        => Poincare.Distance(a, b);
}

public readonly struct FisherRaoSimplexMetric : IDistanceMetric
{
    private static readonly MetricProperties _props = new(
        IsBounded:           true,
        MaxValue:            Math.PI,
        RequiresProbability: true,
        RequiresUnitNorm:    false,
        FixedDimension:      null,
        Geometry:            SpaceGeometry.Spherical,
        BandwidthStrategy:   BandwidthStrategy.QuantileNormalized,
        Name:                "FisherRaoSimplex");

    public MetricProperties Properties => _props;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
        => FisherRaoSimplex.Distance(a, b);
}

public readonly struct FisherRaoHalfPlaneMetric : IDistanceMetric
{
    private static readonly MetricProperties _props = new(
        IsBounded:           false,
        MaxValue:            0.0,
        RequiresProbability: false,
        RequiresUnitNorm:    false,
        FixedDimension:      2,
        Geometry:            SpaceGeometry.Hyperbolic,
        BandwidthStrategy:   BandwidthStrategy.LogScaleHyperbolic,
        Name:                "FisherRaoHalfPlane");

    public MetricProperties Properties => _props;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
        => FisherRaoHalfPlane.Distance(a, b);
}

public readonly struct Wasserstein1Metric : IDistanceMetric
{
    private static readonly MetricProperties _props = new(
        IsBounded:           false,
        MaxValue:            0.0,
        RequiresProbability: false,
        RequiresUnitNorm:    false,
        FixedDimension:      null,
        Geometry:            SpaceGeometry.Euclidean,
        BandwidthStrategy:   BandwidthStrategy.MadConsistencyFactor,
        Name:                "Wasserstein1");

    public MetricProperties Properties => _props;

    public double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
        => Wasserstein1.Distance(a, b);
}
