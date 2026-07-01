using System;

namespace Graphs.Distance;

/// <summary>
/// How bandwidth should be estimated for a given metric — consumed by
/// the graph-construction layer's bandwidth dispatch (Pass 3). The
/// enum lives here in Pass 2 so every metric struct can declare its
/// preference up front; Pass 3 wires the dispatch in
/// <see cref="T:Graphs.BandwidthEstimation"/>.
/// </summary>
public enum BandwidthStrategy
{
    /// <summary>
    /// MAD of nearest-neighbor distances rescaled by the kernel's
    /// MAD-to-σ consistency factor (the historical default). Correct
    /// for Euclidean / Manhattan / Minkowski distances in moderate
    /// dimensions where NN distances are approximately distributed as
    /// half-normal / exponential / Cauchy draws.
    /// </summary>
    MadConsistencyFactor,

    /// <summary>
    /// Normalize NN distances to <c>[0, 1]</c> by dividing by a high
    /// quantile (e.g. 95th percentile), then apply the standard MAD
    /// dispatch in normalized space. Right for bounded metrics
    /// (JSD, FisherRao/Simplex, Cosine) where MAD-on-raw-distances is
    /// scale-mismatched against the consistency factors.
    /// </summary>
    QuantileNormalized,

    /// <summary>
    /// MAD computed on <c>log(d)</c> of NN distances, then
    /// back-transformed. Right for hyperbolic metrics (Poincaré)
    /// where the NN distance distribution is roughly log-normal and
    /// raw MAD is dominated by boundary-region outliers.
    /// </summary>
    LogScaleHyperbolic,

    /// <summary>
    /// Mean of the sample. Used for exact replication of BWD/DMY papers
    /// where bandwidth is defined as the mean of the K-th neighbor distances.
    /// </summary>
    MeanEdgeDistance,
}

/// <summary>
/// The metric's manifold curvature class. Selects the bandwidth
/// consistency factor at graph-build time (Pass 3). Distinct from
/// <see cref="MetricProperties.IsBounded"/>, which is a range property,
/// not a geometry.
/// </summary>
public enum SpaceGeometry
{
    /// <summary>
    /// Flat (κ = 0). Linear MAD→σ consistency factors apply directly
    /// (the historical Euclidean default).
    /// </summary>
    Euclidean,

    /// <summary>
    /// Negatively curved (κ &lt; 0): Poincaré ball, Fisher-Rao half-plane.
    /// Hyperbolic-ball volume grows exponentially with radius, so NN
    /// distances are roughly log-normal — bandwidth is estimated in log
    /// space and the linear MAD→σ factor does not apply.
    /// </summary>
    Hyperbolic,

    /// <summary>
    /// Positively curved (κ &gt; 0): spherical-geodesic (cosine),
    /// Fisher-Rao simplex. Bounded, compact; bandwidth uses quantile
    /// normalization, factor handled in the normalized linear range.
    /// </summary>
    Spherical,
}

/// <summary>
/// Companion descriptor for <see cref="IDistanceMetric"/>. Encodes
/// the metric's domain constraints, range, and the right way to
/// estimate kernel bandwidth from its NN distance sample. Consulted
/// once at graph-build time — never on the hot path.
/// </summary>
/// <param name="IsBounded">True when distance values are bounded
/// above (cosine ≤ π, JSD ≤ 1, etc.). When true, <see cref="MaxValue"/>
/// is the upper bound; when false, <see cref="MaxValue"/> is
/// meaningless.</param>
/// <param name="MaxValue">Upper bound of <c>Distance(a, b)</c> when
/// <see cref="IsBounded"/> is true. Ignored otherwise.</param>
/// <param name="RequiresProbability">True when inputs must be valid
/// probability mass functions (non-negative, sum to 1). The
/// FisherRao-simplex and JSD metrics require this.</param>
/// <param name="RequiresUnitNorm">True when inputs are expected to be
/// L2-normalized. The cosine spherical-geodesic family requires this.</param>
/// <param name="FixedDimension">When non-null, the metric only accepts
/// vectors of exactly this length (e.g. FisherRao/half-plane is 2D
/// only, Wasserstein-1 is 1D only).</param>
/// <param name="Geometry">The metric's manifold curvature class. Used to
/// resolve geometry-aware bandwidth consistency factors at graph-build time.</param>
/// <param name="BandwidthStrategy">The recommended bandwidth strategy
/// (see <see cref="Graphs.Distance.BandwidthStrategy"/>). Pass 3 wires
/// the dispatch.</param>
/// <param name="Name">Stable human-readable identifier for logging /
/// validation messages.</param>
public readonly record struct MetricProperties(
    bool              IsBounded,
    double            MaxValue,
    bool              RequiresProbability,
    bool              RequiresUnitNorm,
    int?              FixedDimension,
    SpaceGeometry     Geometry,
    BandwidthStrategy BandwidthStrategy,
    string            Name)
{
    /// <summary>
    /// Fallback descriptor: unbounded, unconstrained, MAD-consistency-factor
    /// bandwidth, name <c>"unknown"</c>. Use only as a placeholder; every
    /// real metric should declare its own properties.
    /// </summary>
    public static MetricProperties Default => new(
        IsBounded:           false,
        MaxValue:            0.0,
        RequiresProbability: false,
        RequiresUnitNorm:    false,
        FixedDimension:      null,
        Geometry:            SpaceGeometry.Euclidean,
        BandwidthStrategy:   BandwidthStrategy.MadConsistencyFactor,
        Name:                "unknown");

    /// <summary>
    /// Validates that an input batch matches the metric's contract.
    /// Currently checks <see cref="FixedDimension"/>; richer checks
    /// (RequiresProbability, RequiresUnitNorm) land with the bandwidth
    /// dispatch work in Pass 3. Throws <see cref="ArgumentException"/>
    /// on the first violation with a message identifying the metric
    /// and the offending dimension.
    /// </summary>
    /// <param name="dimension">Common dimension of every input vector.</param>
    public void ValidateDimension(int dimension)
    {
        if (FixedDimension is int d && dimension != d)
            throw new ArgumentException(
                $"Metric '{Name}' requires inputs of dimension {d}; got {dimension}.");
    }
}
