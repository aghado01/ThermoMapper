using System;
using System.Runtime.CompilerServices;
using Graphs.Distance;
using Maths.Distance.Euclidean;

namespace Graphs.Distance.Euclidean;

public readonly struct ManhattanMetric : IDistanceMetric
{
    private static readonly MetricProperties _props = new(
        IsBounded:           false,
        MaxValue:            0.0,
        RequiresProbability: false,
        RequiresUnitNorm:    false,
        FixedDimension:      null,
        Geometry:            SpaceGeometry.Euclidean,
        BandwidthStrategy:   BandwidthStrategy.MadConsistencyFactor,
        Name:                "Manhattan");

    public MetricProperties Properties => _props;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
        => Minkowski.Distance(a, b, 1.0);
}

public readonly struct EuclideanMetric : IDistanceMetric
{
    private static readonly MetricProperties _props = new(
        IsBounded:           false,
        MaxValue:            0.0,
        RequiresProbability: false,
        RequiresUnitNorm:    false,
        FixedDimension:      null,
        Geometry:            SpaceGeometry.Euclidean,
        BandwidthStrategy:   BandwidthStrategy.MadConsistencyFactor,
        Name:                "Euclidean");

    public MetricProperties Properties => _props;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
        => Minkowski.Distance(a, b, 2.0);
}

public readonly struct MinkowskiMetric : IDistanceMetric
{
    public double P { get; }

    public MinkowskiMetric(double p)
    {
        if (p > 2.0)
            Console.Error.WriteLine(
                $"[warn] Minkowski p={p} > 2: clamped to 2. Lᵖ stays a valid metric for p≥1, but " +
                "distance concentration grows with p and erodes the cluster contrast SPC needs; prefer p∈[1,2].");
        P = Minkowski.ClampExponent(p);
    }

    public MetricProperties Properties => new(
        IsBounded:           false,
        MaxValue:            0.0,
        RequiresProbability: false,
        RequiresUnitNorm:    false,
        FixedDimension:      null,
        Geometry:            SpaceGeometry.Euclidean,
        BandwidthStrategy:   BandwidthStrategy.MadConsistencyFactor,
        Name:                $"Minkowski(p={P})");

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
        => Minkowski.Distance(a, b, P);
}

/// <summary>
/// Symbol-Hamming over <see cref="ReadOnlySpan{T}"/> of double coordinates.
/// </summary>
public readonly struct HammingMetric : IDistanceMetric
{
    private static readonly MetricProperties _props = new(
        IsBounded:           false,
        MaxValue:            0.0,
        RequiresProbability: false,
        RequiresUnitNorm:    false,
        FixedDimension:      null,
        Geometry:            SpaceGeometry.Euclidean,
        BandwidthStrategy:   BandwidthStrategy.MadConsistencyFactor,
        Name:                "Hamming");

    public MetricProperties Properties => _props;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public double Distance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        if (a.Length != b.Length)
            throw new ArgumentException("Vectors must be of the same length.");

        int diffs = 0;
        for (int i = 0; i < a.Length; i++)
            if (a[i] != b[i]) diffs++;
        return diffs;
    }
}
