using System;

namespace Maths.Regression.Spline;

/// <summary>
/// The carrier the free-knot sampler walks: a knot configuration <c>(k, ξ)</c> — the interior knot locations
/// of a spline on the unit interval, a Field over [0,1]. Degree is a fixed model hyperparameter held by the
/// <see cref="SplineBasis"/>, not part of the per-step state. The reversible-jump moves never mutate a config;
/// each produces a fresh one, so the backing array is treated as immutable by convention.
/// </summary>
public sealed class KnotConfig
{
    /// <summary>Interior knots, strictly increasing, in the open interval (0,1). Length is k.</summary>
    public double[] InteriorKnots { get; }

    /// <summary>Number of interior knots, k.</summary>
    public int Count => InteriorKnots.Length;

    public KnotConfig(double[] interiorKnots)
    {
        ArgumentNullException.ThrowIfNull(interiorKnots);
        InteriorKnots = interiorKnots;
    }
}
